# This file is part of the dune-hdd project:
#   http://users.dune-project.org/projects/dune-hdd
# Copyright holders: Felix Albrecht
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

add_library(linearelliptic_discretization_blockswipdg_fem_eigensparse_2daluconform
            OBJECT EXCLUDE_FROM_ALL
            block-swipdg_fem_eigensparse_2daluconform.cc)

add_library(linearelliptic_discretization_blockswipdg_estimator_fem_eigensparse_2daluconform
            OBJECT EXCLUDE_FROM_ALL
            block-swipdg-estimator_fem_eigensparse_2daluconform.cc)
