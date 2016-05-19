#!/usr/bin/env python2
#
# This file is part of the dune-hdd project:
#   https://github.com/pymor/dune-hdd
# Copyright Holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

from __future__ import division, print_function

from functools import partial
from tempfile import NamedTemporaryFile
from scipy.sparse import coo_matrix, bmat

import numpy as np

from pymor.core.defaults import set_defaults
set_defaults({'pymor.core.cache.default_regions.disk_max_size': 107374182400})

from pymor.algorithms.basisextension import pod_basis_extension
from pymor.algorithms.greedy import greedy
from pymor.algorithms.timestepping import ImplicitEulerTimeStepper
from pymor.core.logger import getLogger
from pymor.discretizations.basic import InstationaryDiscretization
from pymor.grids.oned import OnedGrid
from pymor.operators.constructions import LincombOperator
from pymor.operators.numpy import NumpyMatrixOperator
from pymor.parameters.spaces import CubicParameterSpace
from pymor.playground.algorithms.blockbasisextension import pod_block_basis_extension
from pymor.playground.operators.block import BlockOperator
from pymor.playground.reductors import GenericBlockRBReconstructor
from pymor.reductors.basic import reduce_generic_rb
from pymor.vectorarrays.block import BlockVectorArray
from pymor.vectorarrays.list import ListVectorArray
from pymor.vectorarrays.numpy import NumpyVectorArray
import pymor.core.logger

from dune.pymor.la.container import make_listvectorarray

from simdb.run import new_dataset, add_values, add_logfile

logfile = NamedTemporaryFile(delete=False).name
pymor.core.logger.FILENAME = logfile

logger = getLogger('.morepas3.main')
logger.setLevel('INFO')


class InstationaryDuneVisualizer(object):

    def __init__(self, disc, prefix):
        self.disc = disc
        self.prefix = prefix

    def visualize(self, U, *args, **kwargs):
        import numpy as np
        dune_disc = self.disc._impl
        assert isinstance(U, ListVectorArray)
        filename = kwargs['filename'] if 'filename' in kwargs else self.prefix
        size = len(U)
        pad = len(str(size))
        for ss in np.arange(size):
            dune_disc.visualize(U._list[ss]._impl,
                                filename + '_' + str(ss).zfill(pad),
                                'solution',
                                False) # do not add dirichlet shift

config= {'subdomains': '[1 1 1]',
         'end_time' : 0.01,
         'nt' : 10,
         'num_training_samples' : 1,
         'max_rb_size' : 500,
         'target_error': 1e-1,
         'initial_data': 'exp(-((x[0]-0.5)*(x[0]-0.5)+(x[1]-0.5)*(x[1]-0.5))/0.01)'}
new_dataset('morepas3', **config)

logger.info('initializing DUNE module ...')
from generic_multiscale import dune_module, examples, wrapper
Example = examples[2]['aluconformgrid']['fem']['istl']

# all but the 'type' will be discarded, so no point in setting more details here
solver_options = Example.solver_options('bicgstab.amg.ilu0')

def create_example(ExampleType, num_grid_elements):
    logger_cfg = ExampleType.logger_options()
    logger_cfg.set('info', 0, True)
    logger_cfg.set('info_color', 'blue', True)

    grid_cfg = ExampleType.grid_options('grid.multiscale.provider.cube')
    grid_cfg.set('lower_left',     '[0 0]', True)
    grid_cfg.set('upper_right',    '[1 1]', True)
    grid_cfg.set('num_elements',   num_grid_elements, True)
    grid_cfg.set('num_partitions', '[2 1]', True)

    boundary_cfg = ExampleType.boundary_options('stuff.grid.boundaryinfo.alldirichlet')

    problem_cfg = ExampleType.problem_options('hdd.linearelliptic.problem.thermalblock')
    problem_cfg.set('diffusion_factor.num_elements', '[2 1]', True)

    return ExampleType(logger_cfg, grid_cfg, boundary_cfg, problem_cfg)


def create_discretization(example, config, name, nt):
    stat_blocked_disc = wrapper[example.discretization()]
    stat_nonblocked_disc = stat_blocked_disc.as_nonblocked()

    return (stat_blocked_disc, InstationaryDiscretization(T=config['end_time'],
                                                          initial_data=make_listvectorarray(wrapper[example.project(config['initial_data'])]),
                                                          operator=stat_nonblocked_disc.operator,
                                                          rhs=stat_nonblocked_disc.rhs,
                                                          mass=stat_nonblocked_disc.products['l2'],
                                                          time_stepper=ImplicitEulerTimeStepper(nt, invert_options=solver_options.get_str('type')),
                                                          products=stat_nonblocked_disc.products,
                                                          operators=stat_nonblocked_disc.operators,
                                                          functionals=stat_nonblocked_disc.functionals,
                                                          vector_operators=stat_nonblocked_disc.vector_operators,
                                                          visualizer=InstationaryDuneVisualizer(stat_nonblocked_disc,
                                                                                                name + '.solution'),
                                                          parameter_space=CubicParameterSpace(stat_nonblocked_disc.parameter_type, 0.1, 10),
                                                          cache_region='disk',
                                                          name=name))


# create reference discretization
example_ref = create_example(Example, '[32 32]')
disc_ref_stat, disc_ref = create_discretization(example_ref, config, 'reference discretization', config['nt']*2)
logger.info('  grid has {} subdomains'.format(disc_ref_stat.num_subdomains))
logger.info('  parameter type is {}'.format(disc_ref.parameter_type))
logger.info('  reference discretization has {} DoFs'.format(disc_ref.solution_space.dim))

# create discretization
example = create_example(Example, '[8 8]')
disc_stat, disc = create_discretization(example, config, 'detailed discretization', config['nt'])
logger.info('  discretization has {} DoFs'.format(disc.solution_space.dim))

# logger.info('visualizing grid and data functions ...')
# example_ref.visualize('example')
# disc_ref_stat._impl.visualize(example_ref.project(config['initial_data']), 'initial_data.reference', 'initial_data')
# example.visualize('example.reference')
# disc_stat._impl.visualize(example.project(config['initial_data']), 'initial_data', 'initial_data')

# # compute sample trajectories
# mu = (0.1, 1.0)
# U = disc.solve(mu)
# disc.visualize(U)
# U_ref = disc_ref.solve(mu)
# disc_ref.visualize(U_ref)

# create product and norm
mu_fixed = disc.operator.parse_parameter((1, 1))
mu_fixed = wrapper.dune_parameter(mu_fixed)
elliptic_product = disc_ref.products['elliptic']
elliptic_product = wrapper[elliptic_product._impl.freeze_parameter(mu_fixed)]


def space_norm_squared(U):
    return elliptic_product.apply2(U, U)


def norm(U, order=2):
    '''
      L^2-in-time, X-in-space
    '''
    T = config['end_time']
    nt = len(U)
    time_grid = OnedGrid(domain=(0., T), num_intervals=nt-1)
    assert len(U) == time_grid.size(1)
    qq = time_grid.quadrature_points(0, order=order)
    integral = 0.
    for entity in np.arange(time_grid.size(0)):
        # get quadrature
        qq_e = qq[entity] # points
        ww = time_grid.reference_element.quadrature(order)[1] # weights
        ie = time_grid.integration_elements(0)[entity] # integration element
        # create shape function evaluations
        a = time_grid.centers(1)[entity]
        b = time_grid.centers(1)[entity + 1]
        SF = np.array((1./(a - b)*qq_e[..., 0] - b/(a - b),
                       1./(b - a)*qq_e[..., 0] - a/(b - a)))
        U_a = U._list[entity]
        U_b = U._list[entity + 1]
        values = np.zeros(len(qq_e))
        for ii in np.arange(len(qq_e)):
            # compute U(t)
            U_t = U_a.copy()
            U_t.scal(SF[0][ii])
            U_t.axpy(SF[1][ii], U_b)
            # compute the X-norm of U(t)
            values[ii] = space_norm_squared(make_listvectorarray(U_t))
        integral += np.dot(values, ww)*ie
    return np.sqrt(integral)


def project(example_ref, nt_ref, U):
    T = config['end_time']
    time_grid_ref = OnedGrid(domain=(0., T), num_intervals=nt_ref)
    time_grid = OnedGrid(domain=(0., T), num_intervals=len(U)-1)
    for t_n in time_grid_ref.centers(1):
        coarse_entity = (time_grid.centers(1) <= t_n).nonzero()[0][-1]
        a = time_grid.centers(1)[coarse_entity]
        b = time_grid.centers(1)[coarse_entity + 1]
        SF = np.array((1./(a - b)*t_n - b/(a - b),
                       1./(b - a)*t_n - a/(b - a)))
        import ipdb; ipdb.set_trace()


U_proj = project(example_ref, len(U_ref) - 1, U)

norm(U_ref)

a = b

class Reconstructor(object):

    def __init__(self, disc, RB):
        self._disc = disc
        self._RB = RB

    def reconstruct(self, U):
        if U.dim == 0:
            return self._disc.globalize_vectors(self._disc.solution_space.zeros(1))
        else:
            assert U.dim == np.sum(len(RB) for RB in self._RB)
            local_lens = [len(RB) for RB in self._RB]
            local_starts = np.insert(np.cumsum(local_lens), 0, 0)[:-1]
            localized_U = [NumpyVectorArray(U._array[:, local_starts[ii]:(local_starts[ii] + local_lens[ii])])
                           for ii in np.arange(len(self._RB))]
            U = GenericBlockRBReconstructor(self._RB).reconstruct(BlockVectorArray(localized_U))
            return self._disc.globalize_vectors(U)

    def restricted_to_subbasis(self, dim):
        if not isinstance(dim, tuple):
            dim = len(self.RB)*[dim]
        assert all([dd <= len(rb) for dd, rb in izip(dim, self.RB)])
        return Reconstructor(self._disc, [rb.copy(ind=range(dd)) for rb, dd in izip(self.RB, dim)])


def reductor(discretization, RB, vector_product=None, disable_caching=True, extends=None):
    if RB is None:
        RB = [stat_blocked_disc.local_operator(ss).source.empty() for ss in np.arange(stat_blocked_disc.num_subdomains)]
    rd, rc, reduction_data = reduce_generic_rb(stat_blocked_disc, RB, vector_product, disable_caching, extends)

    def unblock_op(op):
        assert op._blocks[0][0] is not None
        if isinstance(op._blocks[0][0], LincombOperator):
            coefficients = op._blocks[0][0].coefficients
            operators = [None for kk in np.arange(len(op._blocks[0][0].operators))]
            for kk in np.arange(len(op._blocks[0][0].operators)):
                ops = [[op._blocks[ii][jj].operators[kk]
                        if op._blocks[ii][jj] is not None else None
                        for jj in np.arange(op.num_source_blocks)]
                       for ii in np.arange(op.num_range_blocks)]
                operators[kk] = unblock_op(BlockOperator(ops))
            return LincombOperator(operators=operators, coefficients=coefficients)
        else:
            assert all(all([isinstance(block, NumpyMatrixOperator) if block is not None else True
                           for block in row])
                       for row in op._blocks)
            if op.source.dim == 0 and op.range.dim == 0:
                return NumpyMatrixOperator(np.zeros((0, 0)))
            elif op.source.dim == 1:
                mat = np.concatenate([op._blocks[ii][0]._matrix
                                      for ii in np.arange(op.num_range_blocks)],
                                     axis=1)
            elif op.range.dim == 1:
                mat = np.concatenate([op._blocks[0][jj]._matrix
                                      for jj in np.arange(op.num_source_blocks)],
                                     axis=1)
            else:
                mat = bmat([[coo_matrix(op._blocks[ii][jj]._matrix)
                             if op._blocks[ii][jj] is not None else coo_matrix((op._range_dims[ii], op._source_dims[jj]))
                             for jj in np.arange(op.num_source_blocks)]
                            for ii in np.arange(op.num_range_blocks)]).toarray()
            return NumpyMatrixOperator(mat)

    reduced_op = unblock_op(rd.operator)
    reduced_rhs = unblock_op(rd.rhs)
    return (InstationaryDiscretization(T=config['end_time'],
                                       initial_data=reduced_op.source.zeros(1),
                                       operator=reduced_op,
                                       rhs=unblock_op(rd.rhs),
                                       mass=unblock_op(rd.products['l2']),
                                       time_stepper=ImplicitEulerTimeStepper(config['nt']),
                                       products={kk: unblock_op(rd.products[kk]) for kk in rd.products.keys()},
                                       operators={kk: unblock_op(rd.operators[kk])
                                                  for kk in rd.operators.keys() if kk != 'operator'},
                                       functionals={kk: unblock_op(rd.functionals[kk])
                                                    for kk in rd.functionals.keys() if kk != 'rhs'},
                                       vector_operators={kk: unblock_op(rd.vector_operators[kk])
                                                         for kk in rd.vector_operators.keys()},
                                       parameter_space=rd.parameter_space,
                                       cache_region='disk',
                                       name='reduced non-blocked discretization'),
            Reconstructor(stat_blocked_disc, RB),
            reduction_data)


def extension(basis, U):
    if not isinstance(U, BlockVectorArray):
        U = BlockVectorArray([stat_blocked_disc.localize_vector(U, ii)
                              for ii in np.arange(stat_blocked_disc.num_subdomains)])
    return pod_block_basis_extension(basis,
                                     U,
                                     count=1,
                                     product=[stat_blocked_disc.local_product(ss, 'h1')
                                              for ss in np.arange(stat_blocked_disc.num_subdomains)])


greedy_data = greedy(nonblocked_disc,
                     reductor,
                     nonblocked_disc.parameter_space.sample_uniformly(config['num_training_samples']),
                     use_estimator=False,
                     error_norm=norm,
                     extension_algorithm=extension,
                     max_extensions=config['max_rb_size'],
                     target_error=config['target_error'])
rd, rc = greedy_data['reduced_discretization'], greedy_data['reconstructor']
RB = greedy_data['basis']

logger.info(' ')

add_values(time=greedy_data['time'],
           max_err_mus=greedy_data['max_err_mus'],
           extensions=greedy_data['extensions'],
           max_errs=greedy_data['max_errs'],
           basis_sizes=[len(local_RB) for local_RB in RB])
add_logfile(logfile)

