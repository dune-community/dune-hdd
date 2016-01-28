#!/usr/bin/env python2
#
# This file is part of the dune-hdd project:
#   https://github.com/pymor/dune-hdd
# Copyright Holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

from __future__ import division, print_function

from dune.pymor.core import wrap_module

import linearellipticexamplespe10 as dune_module
_, wrapper = wrap_module(dune_module)

examples = {}

for grid in ('Yasp', 'AluConform'):
    for space in ('fem', 'pdelab'):
        for la in ('istl', 'eigen'):
            name = 'Spe10Example_3d{}Grid_{}_{}'.format(grid, space, la)
            lowergrid = grid.lower() + 'grid'
            if name in dune_module.__dict__:
                if lowergrid not in examples:
                    examples[lowergrid] = {}
                if space not in examples[lowergrid]:
                    examples[lowergrid][space] = {}
                examples[lowergrid][space][la] = dune_module.__dict__[name]

