#! /usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the dune-hdd project:
#   http://users.dune-project.org/projects/dune-hdd
# Copyright holders: Rene Milk
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

from __future__ import print_function, division

import os
import sys

import mpi_cg_main
from mpi_cg_bindings import MpiCGExample, Dune

from pymor.parameters.base import Parameter
import pymor.tools.mpi as pmpi
from pymor.discretizations.mpi import mpi_wrap_discretization
from pymor.vectorarrays.mpi import MPIVectorArrayAutoComm, MPIVectorArrayNoComm

print(os.getcwd())
try:
    config_file = sys.argv[2]
except IndexError as ix:
    config_file = '/home/r_milk01/projekte/uni/dune/build/dune-pymor-paper/gcc-release/dune-hdd/examples/linearelliptic/mpi_cg.ini'

example_id = pmpi.call(pmpi.function_call_manage, mpi_cg_main.init_example, config_file, [sys.argv[0]])
disc_id = pmpi.call(pmpi.function_call_manage, mpi_cg_main.discretize, example_id, config_file)
d = mpi_wrap_discretization(disc_id, use_with=False, array_type=MPIVectorArrayNoComm)

'''Bei array_type sollten wir es mit MPIVectorArrayNoComm versuchen. Außerdem würde ich es erstmal 'use_with=False'
wählen. Dies gibt dir eine MPIDiscretization, die für 'solve' MPI-verteilt das 'solve' der
gewrappten Diskretisierungen aufruft. Mit 'use_with=True' sollte man eine gewöhnliche Diskretisierung erhalten,
in der nur die Operatoren ausgetauscht sind und in der dann eine generisch implementierte 'solve'-Methode läuft.
Schließlich solltest Du noch 'with_apply2=True' angeben, was den 'apply2'-Call ebenfalls als MPI-Call implementiert.
(Dafür müssen wir in dune.pymor.operators allerdings noch apply2 in der Wrapper-Klasse implementieren,
was aber einfach sein sollte ..)
'''

mu = Parameter({'diffusion': [0.15227525, 0.87955853, 0.24041678, 0.24039507, 1, 1, 1, 1]})
# mu = Parameter({'diffusion': [1,1,1,1, 1, 1, 1, 1]})
u = d.solve(mu)
# print(u.l2_norm())
d.visualize(u, file_name='solution.vtu', delete=False)
