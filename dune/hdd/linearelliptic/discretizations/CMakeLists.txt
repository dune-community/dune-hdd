# This file is part of the dune-hdd project:
#   http://users.dune-project.org/projects/dune-hdd
# Copyright holders: Felix Albrecht
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

add_library(linearelliptic_discretization_swipdg_leaf_fem_eigensparse_2daluconform
            OBJECT EXCLUDE_FROM_ALL
            swipdg_leaf_fem_eigensparse_2daluconform.cc)
add_library(linearelliptic_discretization_swipdg_level_fem_eigensparse_2daluconform
            OBJECT EXCLUDE_FROM_ALL
            swipdg_level_fem_eigensparse_2daluconform.cc)
            
add_library(linearelliptic_discretization_swipdg__estimator_level_fem_eigensparse_2daluconform
            OBJECT EXCLUDE_FROM_ALL
            swipdg_estimator_level_fem_eigensparse_2daluconform.cc)
