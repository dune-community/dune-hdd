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

logger = getLogger('.battery.main')
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
         'end_time' : 0.001,
         'nt' : 25,
         'num_training_samples' : 50,
         'max_rb_size' : 500,
         'target_error': 1e-10}
new_dataset('parabolic_LRBMS_ENUMATH2015', **config)

logger.info('initializing DUNE module ...')
from generic import dune_module, examples, wrapper
Example = examples[3]['aluconformgrid']['fem']['istl']

logger_cfg = Example.logger_options()
logger_cfg.set('info', 99, True)
logger_cfg.set('info_color', 'blue', True)

# multibat-pymor/problem-data/46x20x20_h4e-6m
battery_geometry = {'lower_left':   '[0      0     0]',
                    'upper_right':  '[0.0184 0.008 0.008]',
                    'num_elements': '[46     20    20]',
                    'separator':    '[0.0084 0.01; 0 0.008; 0 0.008]',
                    'filename': 'geometry__46x20x20_h4e-6m'}
# multibat-pymor/problem-data/ellisoid_5_5_16.8
# battery_geometry = {'lower_left':   '[0      0     0]',
#                     'upper_right':  '[0.0246 0.006 0.006]',
#                     'num_elements': '[246    60    60]',
#                     'separator':    '[]',
#                     'filename': 'geometry__ellisoid_5_5_16.8'}

grid_cfg = Example.grid_options('grid.multiscale.provider.cube')
grid_cfg.set('lower_left', battery_geometry['lower_left'], True)
grid_cfg.set('upper_right', battery_geometry['upper_right'], True)
grid_cfg.set('num_elements', battery_geometry['num_elements'], True)
grid_cfg.set('num_partitions', config['subdomains'], True)

boundary_cfg = dune_module.Dune.Stuff.Common.Configuration()
boundary_cfg.set('type', 'stuff.grid.boundaryinfo.normalbased')
boundary_cfg.set('default', 'neumann')
boundary_cfg.set('dirichlet.0', '[-1 0 0]')
boundary_cfg.set('dirichlet.1', '[1 0 0]')

dirichlet_value = '0'
problem_cfg = Example.problem_options('hdd.linearelliptic.problem.battery')
problem_cfg.set('diffusion_factor.lower_left', battery_geometry['lower_left'], True)
problem_cfg.set('diffusion_factor.upper_right', battery_geometry['upper_right'], True)
problem_cfg.set('diffusion_factor.num_elements', battery_geometry['num_elements'], True)
problem_cfg.set('diffusion_factor.separator', battery_geometry['separator'], True)
problem_cfg.set('diffusion_factor.filename', battery_geometry['filename'], True)
problem_cfg.set('dirichlet.expression', dirichlet_value, True)
problem_cfg.set('neumann.expression', '0', True)
problem_cfg.set('force.value', '1000', True)

# all but the 'type' will be discarded, so no point in setting more details here
solver_options = Example.solver_options('bicgstab.amg.ilu0')


example = Example(logger_cfg, grid_cfg, boundary_cfg, problem_cfg)
stat_blocked_disc = wrapper[example.discretization()]
stat_nonblocked_disc = stat_blocked_disc.as_nonblocked()
logger.info('  grid has {} subdomains'.format(stat_blocked_disc.num_subdomains))
logger.info('  discretization has {} DoFs'.format(stat_nonblocked_disc.solution_space.dim))
logger.info('  parameter type is {}'.format(stat_nonblocked_disc.parameter_type))
# logger.info('visualizing grid and data functions ...')
# example.visualize(problem_cfg.get_str('type'))

nonblocked_disc = InstationaryDiscretization(T=config['end_time'],
                                             initial_data=make_listvectorarray(wrapper[stat_nonblocked_disc._impl.create_vector()]),
                                             operator=stat_nonblocked_disc.operator,
                                             rhs=stat_nonblocked_disc.rhs,
                                             mass=stat_nonblocked_disc.products['l2'],
                                             time_stepper=ImplicitEulerTimeStepper(config['nt'], invert_options=solver_options.get_str('type')),
                                             products=stat_nonblocked_disc.products,
                                             operators=stat_nonblocked_disc.operators,
                                             functionals=stat_nonblocked_disc.functionals,
                                             vector_operators=stat_nonblocked_disc.vector_operators,
                                             visualizer=InstationaryDuneVisualizer(stat_nonblocked_disc,
                                                                                   problem_cfg.get_str('type') + '.solution'),
                                             parameter_space=CubicParameterSpace(stat_nonblocked_disc.parameter_type, 0.1, 10),
                                             cache_region='disk',
                                             name='detailed non-blocked discretization')

# mu = 0.6
# U = nonblocked_disc.solve(mu)
# logger.info('visualizing trajectory ...')
# nonblocked_disc.visualize(U)



def norm(U):
    return np.max(nonblocked_disc.h1_norm(U))


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

