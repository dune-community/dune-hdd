#!/usr/bin/env python2
#
# This file is part of the dune-hdd project:
#   https://github.com/pymor/dune-hdd
# Copyright Holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

from __future__ import division, print_function

import numpy as np

from functools import partial

from pymor.algorithms.timestepping import ImplicitEulerTimeStepper
from pymor.core.logger import getLogger
from pymor.discretizations.basic import InstationaryDiscretization
from pymor.grids.oned import OnedGrid
from pymor.parameters.spaces import CubicParameterSpace
from pymor.vectorarrays.list import ListVectorArray

from dune.pymor.la.container import make_listvectorarray

from generic_multiscale import dune_module, examples, wrapper


logger = getLogger('.morepas3.prepare')
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


def bochner_norm(T, space_norm2, U, mu=None, order=2):
    '''
      L^2-in-time, X-in-space
    '''
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
            values[ii] = space_norm2(make_listvectorarray(U_t), mu)
        integral += np.dot(values, ww)*ie
    return np.sqrt(integral)


def discretize(num_elements, num_refinements, T, nt, initial_data, parameter_range, name='detailed discretization'):
    Example = examples[2]['aluconformgrid']['fem']['istl']

    # all but the 'type' will be discarded, so no point in setting more details here
    solver_options = Example.solver_options('bicgstab.amg.ilu0')

    logger_cfg = Example.logger_options()
    logger_cfg.set('info', -1, True)
    logger_cfg.set('info_color', 'blue', True)

    grid_cfg = Example.grid_options('grid.multiscale.provider.cube')
    grid_cfg.set('lower_left',     '[0 0]', True)
    grid_cfg.set('upper_right',    '[1 1]', True)
    grid_cfg.set('num_elements',   num_elements, True)
    grid_cfg.set('num_partitions', '[2 1]', True)
    grid_cfg.set('num_refinements', num_refinements, True)

    boundary_cfg = Example.boundary_options('stuff.grid.boundaryinfo.alldirichlet')

    problem_cfg = Example.problem_options('hdd.linearelliptic.problem.thermalblock')
    problem_cfg.set('diffusion_factor.num_elements', '[1 1]', True)

    example = Example(logger_cfg, grid_cfg, boundary_cfg, problem_cfg)

    elliptic_LRBMS_disc = wrapper[example.discretization()]
    elliptic_disc = elliptic_LRBMS_disc.as_nonblocked()

    def prolong(coarse_disc, coarse_U):
        time_grid_ref = OnedGrid(domain=(0., T), num_intervals=nt)
        time_grid = OnedGrid(domain=(0., T), num_intervals=(len(coarse_U) - 1))
        U_fine = [None for ii in time_grid_ref.centers(1)]
        for n in np.arange(len(time_grid_ref.centers(1))):
            t_n = time_grid_ref.centers(1)[n]
            coarse_entity = min((time_grid.centers(1) <= t_n).nonzero()[0][-1],
                                time_grid.size(0) - 1)
            a = time_grid.centers(1)[coarse_entity]
            b = time_grid.centers(1)[coarse_entity + 1]
            SF = np.array((1./(a - b)*t_n - b/(a - b),
                           1./(b - a)*t_n - a/(b - a)))
            U_t = coarse_U.copy(ind=coarse_entity)
            U_t.scal(SF[0][0])
            U_t.axpy(SF[1][0], coarse_U, x_ind=(coarse_entity + 1))
            U_fine[n] = wrapper[example.prolong(coarse_disc._impl, U_t._list[0]._impl)]
        return make_listvectorarray(U_fine)

    if isinstance(initial_data, str):
        initial_data = make_listvectorarray(wrapper[example.project(initial_data)])
    else:
        coarse_disc = initial_data[0]
        initial_data = initial_data[1]
        assert len(initial_data) == 1
        initial_data = example.prolong(coarse_disc._impl, initial_data._list[0]._impl)
        initial_data = make_listvectorarray(wrapper[initial_data])

    parabolic_disc = InstationaryDiscretization(T=T,
                                                initial_data=initial_data,
                                                operator=elliptic_disc.operator,
                                                rhs=elliptic_disc.rhs,
                                                mass=elliptic_disc.products['l2'],
                                                time_stepper=ImplicitEulerTimeStepper(nt, invert_options=solver_options.get_str('type')),
                                                products=elliptic_disc.products,
                                                operators=elliptic_disc.operators,
                                                functionals=elliptic_disc.functionals,
                                                vector_operators=elliptic_disc.vector_operators,
                                                visualizer=InstationaryDuneVisualizer(elliptic_disc, 'dune_discretization.solution'),
                                                parameter_space=CubicParameterSpace(elliptic_disc.parameter_type,
                                                                                    parameter_range[0],
                                                                                    parameter_range[1]),
                                                cache_region='disk',
                                                name=name)

    return {'example': example,
            'initial_data': initial_data,
            'wrapper': wrapper,
            'elliptic_LRBMS_disc': elliptic_LRBMS_disc,
            'elliptic_disc': elliptic_disc,
            'parabolic_disc': parabolic_disc,
            'prolongator': prolong}


def prepare(cfg):

    detailed_data = discretize(cfg['dune_num_elements'], cfg['dune_num_refinements'], cfg['end_time'], cfg['nt'],
                               cfg['initial_data'], (cfg['mu_min'], cfg['mu_max']))
    wrapper, parabolic_disc = detailed_data['wrapper'], detailed_data['parabolic_disc']

    logger.info('creating products and norms ...')

    for tp in ('mu_bar', 'mu_hat', 'mu_tilde', 'mu_min', 'mu_max'):
        detailed_data[tp] = parabolic_disc.parse_parameter(cfg[tp])

    space_products = {}
    for kk, prod in parabolic_disc.products.items():
        space_products[kk] = prod
        if prod.parametric:
            for tp in 'mu_bar', 'mu_hat', 'mu_tilde':
                mu = wrapper.dune_parameter(detailed_data[tp])
                space_products['{}_{}'.format(kk, tp)] = wrapper[prod._impl.freeze_parameter(mu)]

    def create_norm2(prod):
        def norm2(U, mu=None):
            return prod.apply2(U, U, mu=mu)[0][0]
        return norm2

    space_norms2 = {kk: create_norm2(prod)
                    for kk, prod in space_products.items()}

    def create_bochner_norm(space_norm2):
        return partial(bochner_norm, cfg['end_time'], space_norm2, order=cfg['integration_order_time'])

    bochner_norms = {kk: create_bochner_norm(space_norm2)
                     for kk, space_norm2 in space_norms2.items()}

    detailed_data['space_products'] = space_products
    detailed_data['space_norms2'] = space_norms2
    detailed_data['bochner_norms'] = bochner_norms

    return detailed_data


