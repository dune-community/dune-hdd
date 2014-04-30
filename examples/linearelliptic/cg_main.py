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
    example = example_module.LinearellipticExampleCG__DuneSGrid__lt___2__2___gt__()
    settingsFilename = example.static_id() + '.cfg'
    if not os.path.exists(settingsFilename):
        example.write_config_file()
    # initialize everything (grid, problem, discretization)
    example.initialize([os.getcwd()])
    discretization = example.discretization()
    if discretization.parametric():
        # solve
        mu = example_module.Dune.Pymor.Parameter('diffusion_factor', [0.1, 1.0, 0.1, 1.0])
        print('')
        print('solving for parameter mu = {}... '.format(mu.report()), end='')
        vector = discretization.create_vector()
        discretization.solve(vector, mu)
        print('done')
        filename = 'solution'
        print('writing to {}.vtu... '.format(filename), end='')
        discretization.visualize(vector, filename, 'CG solution')
        print('done')
    else:
        print('solving ... ', end='')
        vector = discretization.create_vector()
        discretization.solve(vector)
        print('done')
        filename = 'solution'
        print('writing to {}.vtu... '.format(filename), end='')
        discretization.visualize(vector, filename, 'CG solution')
        print('done')

except ImportError:
    print('There was an error importing \'linearellipticexamplecg\'!')
    print('Did you build the example successfully, i.e. did')
    print('  make linearellipticexamplecg')
    print('succeed?')
    raise
