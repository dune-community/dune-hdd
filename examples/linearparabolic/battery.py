#!/usr/bin/env python2
#
# This file is part of the dune-hdd project:
#   https://github.com/pymor/dune-hdd
# Copyright Holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

from __future__ import division, print_function

from pymor.algorithms.timestepping import ImplicitEulerTimeStepper
from pymor.core.logger import getLogger
from pymor.discretizations.basic import InstationaryDiscretization
from pymor.parameters.spaces import CubicParameterSpace
from pymor.vectorarrays.list import ListVectorArray

from dune.pymor.la.container import make_listvectorarray

logger = getLogger('.battery.main')
logger.setLevel('INFO')


class InstationaryDuneVisualizer(object):

    def __init__(self, disc, prefix):
        self.disc = disc
        self.prefix = prefix

    def visualize(self, U, *args, **kwargs):
        dune_disc = self.disc._impl
        import numpy as np
        assert isinstance(U, ListVectorArray)
        filename = kwargs['filename'] if 'filename' in kwargs else self.prefix
        size = len(U)
        pad = len(str(size))
        for ss in np.arange(size):
            dune_U = U._list[ss]._impl
            dune_disc.visualize(U._list[ss]._impl,
                                filename + '_' + str(ss).zfill(pad),
                                'solution',
                                False) # do not add dirichlet shift


logger.info('initializing DUNE module:')
from generic import dune_module, examples, wrapper
Example = examples[3]['aluconformgrid']['fem']['istl']

logger_cfg = Example.logger_options()
logger_cfg.set('info', 99, True)
logger_cfg.set('info_color', 'blue', True)

# multibat-pymor/problem-data/46x20x20_h4e-6m
battery_geometry = {'lower_left':   '[0      0     0]',
                    'upper_right':  '[0.0184 0.008 0.008]',
                    'num_elements': '[46     20    20]'}

grid_cfg = Example.grid_options('stuff.grid.provider.cube')
grid_cfg.set('lower_left', battery_geometry['lower_left'], True)
grid_cfg.set('upper_right', battery_geometry['upper_right'], True)
grid_cfg.set('num_elements', battery_geometry['num_elements'], True)

# boundary_cfg = Example.boundary_options('stuff.grid.boundaryinfo.alldirichlet')
boundary_cfg = dune_module.Dune.Stuff.Common.Configuration()
boundary_cfg.set('type', 'stuff.grid.boundaryinfo.normalbased')
boundary_cfg.set('default', 'neumann')
boundary_cfg.set('dirichlet.0', '[-1 0 0]')
boundary_cfg.set('dirichlet.1', '[1 0 0]')

dirichlet_value = '298'
problem_cfg = Example.problem_options('hdd.linearelliptic.problem.battery')
problem_cfg.set('dirichlet.expression,', dirichlet_value, True)
problem_cfg.set('neumann.expression', '0', True)
problem_cfg.set('force.expression', '0', True)

# all but the 'type' will be discarded, so no point in setting more details here
solver_options = Example.solver_options('bicgstab.amg.ilu0')

example = Example(logger_cfg, grid_cfg, boundary_cfg, problem_cfg)
stationary_discretization = wrapper[example.discretization()]
stationary_discretization.unlock()
stationary_discretization.solver_options = solver_options
stationary_discretization.lock()

logger.info('  discretization has {} DoFs'.format(stationary_discretization.solution_space.dim))
logger.info('  parameter type is {}'.format(stationary_discretization.parameter_type))
logger.info('visualizing grid and data functions ...')
example.visualize(problem_cfg.get_str('type'))

logger.info('projecting initial values ...')
dirichlet_shift = make_listvectorarray(stationary_discretization.vector_operators['dirichlet'])
initial_values = make_listvectorarray(wrapper[example.project(dirichlet_value)])
initial_values -= dirichlet_shift

discretization = InstationaryDiscretization(T=0.0000001,
                                            initial_data=initial_values,
                                            operator=stationary_discretization.operator,
                                            rhs=stationary_discretization.rhs,
                                            mass=stationary_discretization.products['l2_0'],
                                            time_stepper=ImplicitEulerTimeStepper(100, invert_options=solver_options.get_str('type')),
                                            products=stationary_discretization.products,
                                            operators=stationary_discretization.operators,
                                            functionals=stationary_discretization.functionals,
                                            vector_operators=stationary_discretization.vector_operators,
                                            visualizer=InstationaryDuneVisualizer(stationary_discretization,
                                                                                  problem_cfg.get_str('type') + '.solution'),
                                            parameter_space=CubicParameterSpace(stationary_discretization.parameter_type,
                                                                                [0.5, 0.5, 0.5, 0.5, 0.1, 0.1],
                                                                                [2.0, 2.0, 2.0, 2.0, 10., 1.0]))

# mu = {'ANODE': 1.04, 'CATHODE': 1.58, 'CC_ANODE': 238, 'CC_CATHODE': 398, 'ELECTROLYTE': 0.6, 'SEPARATOR': 0.3344}
mu = {'ANODE': 1.04, 'CATHODE': 1.58, 'CC_ANODE': 1.04, 'CC_CATHODE': 1.58, 'ELECTROLYTE': 0.6, 'SEPARATOR': 0.3344}
U = discretization.solve(mu)
U += dirichlet_shift
logger.info('visualizing trajectory ...')
discretization.visualize(U)

logger.info(' ')

