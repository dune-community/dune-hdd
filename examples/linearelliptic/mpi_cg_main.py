#! /usr/bin/env python
# This file is part of the dune-hdd project:
#   http://users.dune-project.org/projects/dune-hdd
# Copyright holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

from __future__ import print_function

import os

from mpi4py import MPI
comm = MPI.COMM_WORLD 
import mpi_cg_bindings 

# create example and write settings file
example = mpi_cg_bindings.MpiCGExample()
settingsFilename = 'foo'
if not os.path.exists(settingsFilename):
    example.write_config_file()

discretization = example.discretization()
if discretization.parametric():
    # solve
    mu = mpi_cg_bindings.Dune.Pymor.Parameter('mu', [0.1])
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

