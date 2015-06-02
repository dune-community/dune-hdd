#! /usr/bin/env python
# This file is part of the dune-hdd project:
#   http://users.dune-project.org/projects/dune-hdd
# Copyright holders: Rene Milk
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

from __future__ import print_function, division

import os
import sys

import mpi_cg_main
import pymor.tools.mpi

print(os.getcwd())
fn = '/home/r_milk01/projekte/uni/dune/build/dune-pymor-paper/gcc-release/dune-hdd/examples/linearelliptic/mpi_cg.ini'
pymor.tools.mpi.call(mpi_cg_main.run_cg, fn, [sys.argv[0]])
