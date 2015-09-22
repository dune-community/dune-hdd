#!/usr/bin/env python2
#
# This file is part of the dune-hdd project:
#   https://github.com/pymor/dune-hdd
# Copyright Holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

from __future__ import division, print_function

import numpy as np

from pymor.algorithms.timestepping import ImplicitEulerTimeStepper
from pymor.core.logger import getLogger
from pymor.discretizations.basic import InstationaryDiscretization
from pymor.operators.constructions import induced_norm
from pymor.operators.constructions import LincombOperator
from pymor.parameters.spaces import CubicParameterSpace
from pymor.vectorarrays.list import ListVectorArray

from dune.pymor.core import wrap_module
from dune.pymor.la.container import make_listvectorarray

import linearparabolicthermalblockexample as dune_module


def prepare(cfg):
    logger = getLogger('.thermalblock.prepare')
    logger.setLevel('INFO')

    logger.info('Initializing DUNE module ({}):'.format(cfg['dune_example']))
    Example = dune_module.__dict__[cfg['dune_example']]
    example = Example(num_blocks=cfg['dune_num_blocks'],
                      num_grid_elements=cfg['dune_num_grid_elements'],
                      products=cfg['dune_products'],
                      info_log_levels=cfg['dune_log_info_level'],
                      debug_log_levels=cfg['dune_log_debug_level'],
                      enable_warnings=cfg['dune_log_enable_warnings'],
                      enable_colors=True,
                      info_color='blue')
    _, wrapper = wrap_module(dune_module)
    stationary_discretization = wrapper[example.discretization()]
    logger.info('  discretization has {} DoFs'.format(stationary_discretization.solution_space.dim))
    logger.info('  parameter type is {}'.format(stationary_discretization.parameter_type))

    logger.info('Visualizing grid and data functions ...')
    example.visualize(cfg['dune_example'])

    logger.info('Projecting initial values ...')
    initial_values = make_listvectorarray(wrapper[example.project(cfg['initial_values_expr'])])

    class InstationaryDuneVisualizer(object):

        def __init__(self, dune_disc, prefix):
            self.dune_disc = dune_disc
            self.prefix = prefix

        def visualize(self, U, *args, **kwargs):
            assert isinstance(U, ListVectorArray)
            filename = kwargs['filename'] if 'filename' in kwargs else self.prefix
            size = len(U)
            pad = len(str(size))
            for ss in np.arange(size):
                self.dune_disc.visualize(U._list[ss]._impl,
                                         filename + '_' + str(ss).zfill(pad),
                                         'solution')

    discretization = InstationaryDiscretization(T=cfg['end_time'],
                                                initial_data=initial_values,
                                                operator=stationary_discretization.operator,
                                                rhs=stationary_discretization.rhs,
                                                mass=stationary_discretization.products['l2'],
                                                time_stepper=ImplicitEulerTimeStepper(100),
                                                products=stationary_discretization.products,
                                                operators=stationary_discretization.operators,
                                                functionals=stationary_discretization.functionals,
                                                vector_operators=stationary_discretization.vector_operators,
                                                visualizer=InstationaryDuneVisualizer(stationary_discretization._impl,
                                                                                      cfg['dune_example'] + '__solution'),
                                                parameter_space=CubicParameterSpace(stationary_discretization.parameter_type,
                                                                                    cfg['parameter_range'][0],
                                                                                    cfg['parameter_range'][1]))

    return {'example': example,
            'discretization': discretization,
            'wrapper': wrapper}

