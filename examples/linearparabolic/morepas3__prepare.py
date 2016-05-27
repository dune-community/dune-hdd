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


def discretize(num_elements, num_partitions, T, nt, initial_data, parameter_range, name='detailed discretization'):
    Example = examples[2]['aluconformgrid']['fem']['istl']

    # all but the 'type' will be discarded, so no point in setting more details here
    solver_options = Example.solver_options('bicgstab.amg.ilu0')

    logger_cfg = Example.logger_options()
    logger_cfg.set('info', -1, True)
    logger_cfg.set('info_color', 'blue', True)

    grid_cfg = Example.grid_options('grid.multiscale.provider.cube')
    grid_cfg.set('lower_left',     '[0 0]', True)
    grid_cfg.set('upper_right',    '[5 1]', True)
    grid_cfg.set('num_elements',   num_elements, True)
    grid_cfg.set('num_partitions', num_partitions, True)

    boundary_cfg = Example.boundary_options('stuff.grid.boundaryinfo.alldirichlet')

    problem_cfg = Example.problem_options('hdd.linearelliptic.problem.OS2015.spe10model1')
    problem_cfg.set('parametric_channel', 'true', True)
    problem_cfg.set('channel_boundary_layer', '0', True)
    problem_cfg.set('filename', 'perm_case1.dat', True)
    problem_cfg.set('lower_left', '[0 0]', True)
    problem_cfg.set('upper_right', '[5 1]', True)
    problem_cfg.set('num_elements', '[100 20]', True)
    problem_cfg.set('forces.0.domain', '[0.95 1.10; 0.30 0.45]', True)
    problem_cfg.set('forces.0.value', '2000', True)
    problem_cfg.set('forces.1.domain', '[3.00 3.15; 0.75 0.90]', True)
    problem_cfg.set('forces.1.value', '-1000', True)
    problem_cfg.set('forces.2.domain', '[4.25 4.40; 0.25 0.40]', True)
    problem_cfg.set('forces.2.value', '-1000', True)
    problem_cfg.set('channel.0.value', '-1.07763239495', True)
    problem_cfg.set('channel.1.value', '-1.07699512772', True)
    problem_cfg.set('channel.2.value', '-1.07356156439', True)
    problem_cfg.set('channel.3.value', '-1.06602281736', True)
    problem_cfg.set('channel.4.value', '-1.06503683743', True)
    problem_cfg.set('channel.5.value', '-1.07974870426', True)
    problem_cfg.set('channel.6.value', '-1.05665895923', True)
    problem_cfg.set('channel.7.value', '-1.08310334837', True)
    problem_cfg.set('channel.8.value', '-1.05865484973', True)
    problem_cfg.set('channel.9.value', '-1.05871039535', True)
    problem_cfg.set('channel.10.value', '-1.08136695901', True)
    problem_cfg.set('channel.11.value', '-1.08490172721', True)
    problem_cfg.set('channel.12.value', '-1.06641120758', True)
    problem_cfg.set('channel.13.value', '-1.06812773298', True)
    problem_cfg.set('channel.14.value', '-1.07695652049', True)
    problem_cfg.set('channel.15.value', '-1.08630079205', True)
    problem_cfg.set('channel.16.value', '-1.08273722112', True)
    problem_cfg.set('channel.17.value', '-1.07500402155', True)
    problem_cfg.set('channel.18.value', '-1.08607142562', True)
    problem_cfg.set('channel.19.value', '-1.07268761799', True)
    problem_cfg.set('channel.20.value', '-1.08537037362', True)
    problem_cfg.set('channel.21.value', '-1.08466927273', True)
    problem_cfg.set('channel.22.value', '-1.08444661815', True)
    problem_cfg.set('channel.23.value', '-1.08957037967', True)
    problem_cfg.set('channel.24.value', '-1.08047394052', True)
    problem_cfg.set('channel.25.value', '-1.08221229083', True)
    problem_cfg.set('channel.26.value', '-1.08568599863', True)
    problem_cfg.set('channel.27.value', '-1.08428347872', True)
    problem_cfg.set('channel.28.value', '-1.09104098734', True)
    problem_cfg.set('channel.29.value', '-1.09492700673', True)
    problem_cfg.set('channel.30.value', '-1.09760440537', True)
    problem_cfg.set('channel.31.value', '-1.09644989453', True)
    problem_cfg.set('channel.32.value', '-1.09441681025', True)
    problem_cfg.set('channel.33.value', '-1.09533290654', True)
    problem_cfg.set('channel.34.value', '-1.1001430808', True)
    problem_cfg.set('channel.35.value', '-1.10065627621', True)
    problem_cfg.set('channel.36.value', '-1.10125877186', True)
    problem_cfg.set('channel.37.value', '-1.10057485893', True)
    problem_cfg.set('channel.38.value', '-1.10002261906', True)
    problem_cfg.set('channel.39.value', '-1.10219154209', True)
    problem_cfg.set('channel.40.value', '-1.09994463801', True)
    problem_cfg.set('channel.41.value', '-1.10265630533', True)
    problem_cfg.set('channel.42.value', '-1.10448566526', True)
    problem_cfg.set('channel.43.value', '-1.10735820121', True)
    problem_cfg.set('channel.44.value', '-1.1070022367', True)
    problem_cfg.set('channel.45.value', '-1.10777650387', True)
    problem_cfg.set('channel.46.value', '-1.10892785562', True)
    problem_cfg.set('channel.0.domain', '[1.7 1.75; 0.5 0.55]', True)
    problem_cfg.set('channel.1.domain', '[1.75 1.8; 0.5 0.55]', True)
    problem_cfg.set('channel.2.domain', '[1.8 1.85; 0.5 0.55]', True)
    problem_cfg.set('channel.3.domain', '[1.85 1.9; 0.5 0.55]', True)
    problem_cfg.set('channel.4.domain', '[1.9 1.95; 0.5 0.55]', True)
    problem_cfg.set('channel.5.domain', '[1.95 2.0; 0.5 0.55]', True)
    problem_cfg.set('channel.6.domain', '[2.0 2.05; 0.5 0.55]', True)
    problem_cfg.set('channel.7.domain', '[2.05 2.1; 0.5 0.55]', True)
    problem_cfg.set('channel.8.domain', '[2.1 2.15; 0.5 0.55]', True)
    problem_cfg.set('channel.9.domain', '[2.15 2.2; 0.5 0.55]', True)
    problem_cfg.set('channel.10.domain', '[2.2 2.25; 0.5 0.55]', True)
    problem_cfg.set('channel.11.domain', '[2.25 2.3; 0.5 0.55]', True)
    problem_cfg.set('channel.12.domain', '[2.3 2.35; 0.5 0.55]', True)
    problem_cfg.set('channel.13.domain', '[2.35 2.4; 0.5 0.55]', True)
    problem_cfg.set('channel.14.domain', '[2.4 2.45; 0.5 0.55]', True)
    problem_cfg.set('channel.15.domain', '[2.45 2.5; 0.5 0.55]', True)
    problem_cfg.set('channel.16.domain', '[2.5 2.55; 0.5 0.55]', True)
    problem_cfg.set('channel.17.domain', '[2.55 2.6; 0.5 0.55]', True)
    problem_cfg.set('channel.18.domain', '[2.6 2.65; 0.5 0.55]', True)
    problem_cfg.set('channel.19.domain', '[2.65 2.7; 0.5 0.55]', True)
    problem_cfg.set('channel.20.domain', '[2.7 2.75; 0.5 0.55]', True)
    problem_cfg.set('channel.21.domain', '[2.75 2.8; 0.5 0.55]', True)
    problem_cfg.set('channel.22.domain', '[2.8 2.85; 0.5 0.55]', True)
    problem_cfg.set('channel.23.domain', '[2.85 2.9; 0.5 0.55]', True)
    problem_cfg.set('channel.24.domain', '[2.9 2.95; 0.5 0.55]', True)
    problem_cfg.set('channel.25.domain', '[2.95 3.0; 0.5 0.55]', True)
    problem_cfg.set('channel.26.domain', '[3.0 3.05; 0.5 0.55]', True)
    problem_cfg.set('channel.27.domain', '[3.05 3.1; 0.5 0.55]', True)
    problem_cfg.set('channel.28.domain', '[3.1 3.15; 0.5 0.55]', True)
    problem_cfg.set('channel.29.domain', '[3.15 3.2; 0.5 0.55]', True)
    problem_cfg.set('channel.30.domain', '[3.2 3.25; 0.5 0.55]', True)
    problem_cfg.set('channel.31.domain', '[3.25 3.3; 0.5 0.55]', True)
    problem_cfg.set('channel.32.domain', '[3.3 3.35; 0.5 0.55]', True)
    problem_cfg.set('channel.33.domain', '[3.35 3.4; 0.5 0.55]', True)
    problem_cfg.set('channel.34.domain', '[3.4 3.45; 0.5 0.55]', True)
    problem_cfg.set('channel.35.domain', '[3.45 3.5; 0.5 0.55]', True)
    problem_cfg.set('channel.36.domain', '[3.5 3.55; 0.5 0.55]', True)
    problem_cfg.set('channel.37.domain', '[3.55 3.6; 0.5 0.55]', True)
    problem_cfg.set('channel.38.domain', '[3.6 3.65; 0.5 0.55]', True)
    problem_cfg.set('channel.39.domain', '[3.65 3.7; 0.5 0.55]', True)
    problem_cfg.set('channel.40.domain', '[3.7 3.75; 0.5 0.55]', True)
    problem_cfg.set('channel.41.domain', '[3.75 3.8; 0.5 0.55]', True)
    problem_cfg.set('channel.42.domain', '[3.8 3.85; 0.5 0.55]', True)
    problem_cfg.set('channel.43.domain', '[3.85 3.9; 0.5 0.55]', True)
    problem_cfg.set('channel.44.domain', '[3.9 3.95; 0.5 0.55]', True)
    problem_cfg.set('channel.45.domain', '[3.95 4.0; 0.5 0.55]', True)
    problem_cfg.set('channel.46.domain', '[4.0 4.05; 0.5 0.55]', True)
    problem_cfg.set('channel.47.value', '-1.10372589211', True)
    problem_cfg.set('channel.48.value', '-1.1020889988', True)
    problem_cfg.set('channel.49.value', '-1.09806955069', True)
    problem_cfg.set('channel.50.value', '-1.10000902421', True)
    problem_cfg.set('channel.51.value', '-1.08797468724', True)
    problem_cfg.set('channel.52.value', '-1.08827472176', True)
    problem_cfg.set('channel.53.value', '-1.08692237109', True)
    problem_cfg.set('channel.54.value', '-1.07893190093', True)
    problem_cfg.set('channel.55.value', '-1.08748373853', True)
    problem_cfg.set('channel.56.value', '-1.07445197324', True)
    problem_cfg.set('channel.57.value', '-1.08246613163', True)
    problem_cfg.set('channel.58.value', '-1.06726790504', True)
    problem_cfg.set('channel.59.value', '-1.07891217847', True)
    problem_cfg.set('channel.60.value', '-1.07260827126', True)
    problem_cfg.set('channel.61.value', '-1.07094062748', True)
    problem_cfg.set('channel.62.value', '-1.0692399429', True)
    problem_cfg.set('channel.63.value', '-1.00099885701', True)
    problem_cfg.set('channel.64.value', '-1.00109544002', True)
    problem_cfg.set('channel.65.value', '-0.966491003242', True)
    problem_cfg.set('channel.66.value', '-0.802284684014', True)
    problem_cfg.set('channel.67.value', '-0.980790923021', True)
    problem_cfg.set('channel.68.value', '-0.614478271687', True)
    problem_cfg.set('channel.69.value', '-0.288129858959', True)
    problem_cfg.set('channel.70.value', '-0.929509396842', True)
    problem_cfg.set('channel.71.value', '-0.992376505995', True)
    problem_cfg.set('channel.72.value', '-0.968162494855', True)
    problem_cfg.set('channel.73.value', '-0.397316938901', True)
    problem_cfg.set('channel.74.value', '-0.970934956609', True)
    problem_cfg.set('channel.75.value', '-0.784344730096', True)
    problem_cfg.set('channel.76.value', '-0.539725422323', True)
    problem_cfg.set('channel.77.value', '-0.915632282372', True)
    problem_cfg.set('channel.78.value', '-0.275089177273', True)
    problem_cfg.set('channel.79.value', '-0.949684959286', True)
    problem_cfg.set('channel.80.value', '-0.936132529794', True)
    problem_cfg.set('channel.47.domain', '[2.6 2.65; 0.45 0.50]', True)
    problem_cfg.set('channel.48.domain', '[2.65 2.7; 0.45 0.50]', True)
    problem_cfg.set('channel.49.domain', '[2.7 2.75; 0.45 0.50]', True)
    problem_cfg.set('channel.50.domain', '[2.75 2.8; 0.45 0.50]', True)
    problem_cfg.set('channel.51.domain', '[2.8 2.85; 0.45 0.50]', True)
    problem_cfg.set('channel.52.domain', '[2.85 2.9; 0.45 0.50]', True)
    problem_cfg.set('channel.53.domain', '[2.9 2.95; 0.45 0.50]', True)
    problem_cfg.set('channel.54.domain', '[2.95 3.0; 0.45 0.50]', True)
    problem_cfg.set('channel.55.domain', '[3.0 3.05; 0.45 0.50]', True)
    problem_cfg.set('channel.56.domain', '[3.05 3.1; 0.45 0.50]', True)
    problem_cfg.set('channel.57.domain', '[3.1 3.15; 0.45 0.50]', True)
    problem_cfg.set('channel.58.domain', '[3.15 3.2; 0.45 0.50]', True)
    problem_cfg.set('channel.59.domain', '[3.2 3.25; 0.45 0.50]', True)
    problem_cfg.set('channel.60.domain', '[3.25 3.3; 0.45 0.50]', True)
    problem_cfg.set('channel.61.domain', '[3.3 3.35; 0.45 0.50]', True)
    problem_cfg.set('channel.62.domain', '[3.35 3.4; 0.45 0.50]', True)
    problem_cfg.set('channel.63.domain', '[3.4 3.45; 0.45 0.50]', True)
    problem_cfg.set('channel.64.domain', '[3.45 3.5; 0.45 0.50]', True)
    problem_cfg.set('channel.65.domain', '[3.5 3.55; 0.45 0.50]', True)
    problem_cfg.set('channel.66.domain', '[3.55 3.6; 0.45 0.50]', True)
    problem_cfg.set('channel.67.domain', '[3.6 3.65; 0.45 0.50]', True)
    problem_cfg.set('channel.68.domain', '[3.65 3.7; 0.45 0.50]', True)
    problem_cfg.set('channel.69.domain', '[3.7 3.75; 0.45 0.50]', True)
    problem_cfg.set('channel.70.domain', '[3.75 3.8; 0.45 0.50]', True)
    problem_cfg.set('channel.71.domain', '[3.8 3.85; 0.45 0.50]', True)
    problem_cfg.set('channel.72.domain', '[3.85 3.9; 0.45 0.50]', True)
    problem_cfg.set('channel.73.domain', '[3.9 3.95; 0.45 0.50]', True)
    problem_cfg.set('channel.74.domain', '[3.95 4.0; 0.45 0.50]', True)
    problem_cfg.set('channel.75.domain', '[4.0 4.05; 0.45 0.50]', True)
    problem_cfg.set('channel.76.domain', '[4.05 4.1; 0.45 0.50]', True)
    problem_cfg.set('channel.77.domain', '[4.1 4.15; 0.45 0.50]', True)
    problem_cfg.set('channel.78.domain', '[4.15 4.2; 0.45 0.50]', True)
    problem_cfg.set('channel.79.domain', '[4.2 4.25; 0.45 0.50]', True)
    problem_cfg.set('channel.80.domain', '[4.25 4.3; 0.45 0.50]', True)
    problem_cfg.set('channel.81.value', '-1.10923642795', True)
    problem_cfg.set('channel.82.value', '-1.10685618623', True)
    problem_cfg.set('channel.83.value', '-1.1057800376', True)
    problem_cfg.set('channel.84.value', '-1.10187723629', True)
    problem_cfg.set('channel.85.value', '-1.10351710464', True)
    problem_cfg.set('channel.86.value', '-1.10037551137', True)
    problem_cfg.set('channel.87.value', '-1.09724407076', True)
    problem_cfg.set('channel.88.value', '-1.09604600208', True)
    problem_cfg.set('channel.89.value', '-1.09354469656', True)
    problem_cfg.set('channel.90.value', '-1.08934455354', True)
    problem_cfg.set('channel.91.value', '-1.08155476586', True)
    problem_cfg.set('channel.92.value', '-1.07815397899', True)
    problem_cfg.set('channel.93.value', '-1.09174062023', True)
    problem_cfg.set('channel.94.value', '-1.07433616068', True)
    problem_cfg.set('channel.95.value', '-1.08030587701', True)
    problem_cfg.set('channel.81.domain', '[1.95 2.0; 0.40 0.45]', True)
    problem_cfg.set('channel.82.domain', '[2.0 2.05; 0.40 0.45]', True)
    problem_cfg.set('channel.83.domain', '[2.05 2.1; 0.40 0.45]', True)
    problem_cfg.set('channel.84.domain', '[2.1 2.15; 0.40 0.45]', True)
    problem_cfg.set('channel.85.domain', '[2.15 2.2; 0.40 0.45]', True)
    problem_cfg.set('channel.86.domain', '[2.2 2.25; 0.40 0.45]', True)
    problem_cfg.set('channel.87.domain', '[2.25 2.3; 0.40 0.45]', True)
    problem_cfg.set('channel.88.domain', '[2.3 2.35; 0.40 0.45]', True)
    problem_cfg.set('channel.89.domain', '[2.35 2.4; 0.40 0.45]', True)
    problem_cfg.set('channel.90.domain', '[2.4 2.45; 0.40 0.45]', True)
    problem_cfg.set('channel.91.domain', '[2.45 2.5; 0.40 0.45]', True)
    problem_cfg.set('channel.92.domain', '[2.5 2.55; 0.40 0.45]', True)
    problem_cfg.set('channel.93.domain', '[2.55 2.6; 0.40 0.45]', True)
    problem_cfg.set('channel.94.domain', '[2.6 2.65; 0.40 0.45]', True)
    problem_cfg.set('channel.95.domain', '[2.65 2.7; 0.40 0.45]', True)
    problem_cfg.set('channel.96.value', '-1.00032869407', True)
    problem_cfg.set('channel.97.value', '-1.01175908905', True)
    problem_cfg.set('channel.98.value', '-1.04954395793', True)
    problem_cfg.set('channel.99.value', '-1.017967697', True)
    problem_cfg.set('channel.100.value', '-1.04647184091', True)
    problem_cfg.set('channel.101.value', '-1.01911894831', True)
    problem_cfg.set('channel.102.value', '-1.00699340158', True)
    problem_cfg.set('channel.103.value', '-0.995492960025', True)
    problem_cfg.set('channel.104.value', '-1.0373059007', True)
    problem_cfg.set('channel.96.domain', '[2.25 2.3; 0.35 0.40]', True)
    problem_cfg.set('channel.97.domain', '[2.3 2.35; 0.35 0.40]', True)
    problem_cfg.set('channel.98.domain', '[2.35 2.4; 0.35 0.40]', True)
    problem_cfg.set('channel.99.domain', '[2.4 2.45; 0.35 0.40]', True)
    problem_cfg.set('channel.100.domain', '[2.45 2.5; 0.35 0.40]', True)
    problem_cfg.set('channel.101.domain', '[2.5 2.55; 0.35 0.40]', True)
    problem_cfg.set('channel.102.domain', '[2.55 2.6; 0.35 0.40]', True)
    problem_cfg.set('channel.103.domain', '[2.6 2.65; 0.35 0.40]', True)
    problem_cfg.set('channel.104.domain', '[2.65 2.7; 0.35 0.4]', True)

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
        # initial_data = elliptic_disc.operator.apply_inverse(initial_data, mu=(1, 1))
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
                                                name='{} ({} DoFs)'.format(name, elliptic_disc.solution_space.dim))

    return {'example': example,
            'initial_data': initial_data,
            'wrapper': wrapper,
            'elliptic_LRBMS_disc': elliptic_LRBMS_disc,
            'elliptic_disc': elliptic_disc,
            'parabolic_disc': parabolic_disc,
            'prolongator': prolong}


def prepare(cfg):

    detailed_data = discretize(cfg['dune_num_elements'], cfg['dune_num_partitions'], cfg['end_time'], cfg['nt'],
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


