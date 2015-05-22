#! /usr/bin/env python
# This file is part of the dune-hdd project:
#   http://users.dune-project.org/projects/dune-hdd
# Copyright holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

from __future__ import print_function, division

import os
import sys

from mpi4py import MPI
comm = MPI.COMM_WORLD
from mpi_cg_bindings import MpiCGExample, Dune


try:
    fn = sys.argv[1]
    del sys.argv[1]
    cfg = Dune.Stuff.Common.Configuration(sys.argv, fn)
except Exception as e:
    print('Could not load parameter file, defaulting')
    cfg = Dune.Stuff.Common.Configuration('grids.refinements', 2)

example = MpiCGExample(cfg.get_int("grids.refinements"))
discretization = example.discretization()
if discretization.parametric():
    # solve
    mu = Dune.Pymor.Parameter('mu', [0.1])
    print('')
    print('solving for parameter mu = {}... '.format(mu.report()), end='')
    vector = discretization.create_vector()
    discretization.solve(cfg.sub('solver'), vector, mu)
    print('done')
    filename = 'solution'
    print('writing to {}.vtu... '.format(filename), end='')
    discretization.visualize(vector, filename, 'CG_solution')
    print('done')
else:
    print('solving ... ', end='')
    vector = discretization.create_vector()
    discretization.solve(cfg.sub('solver'), vector)
    print('done')
    filename = 'solution'
    print('writing to {}.vtu... '.format(filename), end='')
    discretization.visualize(vector, filename, 'CG_solution')
    print('done')

