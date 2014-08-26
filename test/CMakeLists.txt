# This file is part of the dune-hdd project:
#   http://users.dune-project.org/projects/dune-hdd
# Copyright holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

add_library(test_linearelliptic_swipg_esv2007_fem_eigensparse_2daluconform
            OBJECT
            EXCLUDE_FROM_ALL
            linearelliptic-swipdg_ESV2007_fem_eigensparse_2daluconform.cc)

add_library(test_linearelliptic_blockswipg_esv2007_fem_eigensparse_2daluconform
            OBJECT
            EXCLUDE_FROM_ALL
            linearelliptic-block-swipdg_ESV2007_fem_eigensparse_2daluconform.cc)
add_library(test_linearelliptic_blockswipg_os2014_fem_eigensparse_2daluconform
            OBJECT
            EXCLUDE_FROM_ALL
            linearelliptic-block-swipdg_OS2014_fem_eigensparse_2daluconform.cc)

add_library(test_linearelliptic_swipg
            STATIC
            EXCLUDE_FROM_ALL
            $<TARGET_OBJECTS:test_linearelliptic_blockswipg_os2014_fem_eigensparse_2daluconform>
            $<TARGET_OBJECTS:test_linearelliptic_blockswipg_esv2007_fem_eigensparse_2daluconform>
            $<TARGET_OBJECTS:test_linearelliptic_swipg_esv2007_fem_eigensparse_2daluconform>
            $<TARGET_OBJECTS:linearelliptic_problems_common_2daluconform>
            $<TARGET_OBJECTS:linearelliptic_problem_esv2007_2daluconform>
            $<TARGET_OBJECTS:linearelliptic_testcase_esv2007_2daluconform>
            $<TARGET_OBJECTS:linearelliptic_testcase_esv2007_multiscale_2daluconform>
            $<TARGET_OBJECTS:linearelliptic_testcase_os2014_multiscale_2daluconform>
            $<TARGET_OBJECTS:linearelliptic_discretization_swipdg_level_fem_eigensparse_2daluconform>
            $<TARGET_OBJECTS:linearelliptic_discretization_blockswipdg_fem_eigensparse_2daluconform>
            $<TARGET_OBJECTS:linearelliptic_discretization_swipdg_estimator_level_fem_eigensparse_2daluconform>
            $<TARGET_OBJECTS:linearelliptic_discretization_blockswipdg_estimator_fem_eigensparse_2daluconform>)
add_dune_alugrid_flags(test_linearelliptic_swipg)

enable_testing()

DEPENDENCYCHECK(${hdd})

BEGIN_TESTCASES()
END_TESTCASES()

if (NOT CMAKE_WITH_AUTOTOOLS)
  target_link_libraries(test_OS2014 test_linearelliptic_swipg)
endif (NOT CMAKE_WITH_AUTOTOOLS)
