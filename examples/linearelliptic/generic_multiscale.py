#!/usr/bin/env python2
#
# This file is part of the dune-hdd project:
#   https://github.com/pymor/dune-hdd
# Copyright Holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

from __future__ import division, print_function

from dune.pymor.core import wrap_module

import linearellipticexamplegeneric as dune_module
_, wrapper = wrap_module(dune_module)

examples = {}

for dd in (1, 2, 3):
    for grid in ('Yasp', 'AluConform', 'Sp', 'S'):
        for space in ('fem', 'pdelab'):
            for la in ('istl', 'eigen'):
                name = 'GenericLinearellipticMultiscaleExample_{}d{}Grid_{}_{}'.format(dd, grid, space, la)
                lowergrid = grid.lower() + 'grid'
                if name in dune_module.__dict__:
                    if dd not in examples:
                        examples[dd] = {}
                    if lowergrid not in examples[dd]:
                        examples[dd][lowergrid] = {}
                    if space not in examples[dd][lowergrid]:
                        examples[dd][lowergrid][space] = {}
                    examples[dd][lowergrid][space][la] = dune_module.__dict__[name]

