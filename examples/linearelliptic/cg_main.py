#! /usr/bin/env python
# This file is part of the dune-hdd project:
#   http://users.dune-project.org/projects/dune-hdd
# Copyright holders: Felix Albrecht
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

from __future__ import print_function

import os

try:
    import linearellipticexamplecg as example_module

    # create example and write settings file
    example = example_module.LinearellipticExampleCG__DuneSGrid__lt___2__2___gt___1()
    settingsFilename = example.static_id() + '.settings'
    if not os.path.exists(settingsFilename):
        example.write_settings_file()
    # initialize everything (grid, problem, discretization)
    example.initialize([os.getcwd()])
    discretization = example.discretization()
    # inspect the discretization
    available_operators = list(discretization.available_operators())
    print('discretization provides:')
    for operator_id in available_operators:
        operator = discretization.get_operator(operator_id)
        print('  operator \'{}\' with parameter type: {}'.format(operator_id,
                                                                 operator.parameter_type().report()))
    available_functionals = list(discretization.available_functionals())
    for functional_id in available_functionals:
        functional = discretization.get_functional(functional_id)
        print('  functional \'{}\' with parameter type: {}'.format(functional_id,
                                                                   functional.parameter_type().report()))
    available_solver_contexts = list(discretization.solver_options())
    for solver_context in available_solver_contexts:
        solver_options = discretization.solver_options(context)
        print('  solver context \'{}\': \'{}\''.format(solver_context, solver_options))
    # solve
    mu = example_module.Dune.Pymor.Parameter('diffusion', [0.1, 1.0, 0.1, 1.0])
    print('')
    print('solving for parameter mu = {}... '.format(mu.report()), end='')
    vector = discretization.create_vector()
    discretization.solve(vector, mu)
    print('done')
    filename = 'solution'
    print('writing to {}.vtu... '.format(filename), end='')
    discretization.visualize(vector, filename, 'CG solution')
    print('done')

except ImportError:
    print('There was an error importing \'linearellipticexamplecg\'!')
    print('Did you build the example successfully,')
    print('i.e. did \'make linearellipticexamplecg\' succeed?')
    raise
