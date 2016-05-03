# This file is part of the dune-hdd project:
#   http://users.dune-project.org/projects/dune-hdd
# Copyright holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

# test infrastructure
enable_testing()
DEPENDENCYCHECK()
BEGIN_TESTCASES()
END_TESTCASES()

# expectation targets
add_library(obj_test_linearelliptic-block-swipdg-expectations_os2015_2daluconform
            OBJECT
            EXCLUDE_FROM_ALL
            linearelliptic-block-swipdg-expectations_os2015_2daluconform.cxx)
add_library(obj_test_linearelliptic-block-swipdg-expectations_os2014_2daluconform
            OBJECT
            EXCLUDE_FROM_ALL
            linearelliptic-block-swipdg-expectations_os2014_2daluconform.cxx)
add_library(obj_test_linearelliptic-swipdg-expectations_os2014_2daluconform
            OBJECT
            EXCLUDE_FROM_ALL
            linearelliptic-swipdg-expectations_os2014_2daluconform.cxx)
add_library(obj_test_linearelliptic-swipdg-expectations_os2014_2dsgrid
            OBJECT
            EXCLUDE_FROM_ALL
            linearelliptic-swipdg-expectations_os2014_2dsgrid.cxx)

# test target libs
add_library(lib_test_OS2015_SISC__6_1__academic_example__estimator_study
            STATIC
            EXCLUDE_FROM_ALL
            $<TARGET_OBJECTS:obj_test_linearelliptic-block-swipdg-expectations_os2015_2daluconform>)
add_dune_alugrid_flags(lib_test_OS2015_SISC__6_1__academic_example__estimator_study)

add_library(lib_test_OS2015_SISC__6_1__multiscale_example__estimator_study
            STATIC
            EXCLUDE_FROM_ALL
            $<TARGET_OBJECTS:obj_test_linearelliptic-block-swipdg-expectations_os2015_2daluconform>)
add_dune_alugrid_flags(lib_test_OS2015_SISC__6_1__multiscale_example__estimator_study)
add_library(lib_test_OS2014_FVCA7__estimator_study
            STATIC
            EXCLUDE_FROM_ALL
            $<TARGET_OBJECTS:obj_test_linearelliptic-swipdg-expectations_os2014_2daluconform>
            $<TARGET_OBJECTS:obj_test_linearelliptic-block-swipdg-expectations_os2014_2daluconform>)
add_dune_alugrid_flags(lib_test_OS2014_FVCA7__estimator_study)

# test targets
target_link_libraries(test_OS2015_SISC__6_1__academic_example__estimator_study
                      lib_test_OS2015_SISC__6_1__academic_example__estimator_study)
target_link_libraries(test_OS2015_SISC__6_1__multiscale_example__estimator_study
                      lib_test_OS2015_SISC__6_1__multiscale_example__estimator_study)
target_link_libraries(test_OS2014_FVCA7__estimator_study
                      lib_test_OS2014_FVCA7__estimator_study)

# clang analizer targets (should be guarded)
#add_analyze(OS2014.cc)
#target_link_libraries(analyze_OS2014.cc lib_test_OS2014)

# link spe10 data file if present
if (NOT ${SPE10MODEL1DATA} STREQUAL "SPE10MODEL1DATA-NOTFOUND")
  set (SPE10MODEL1DATA_TARGET_FILENAME "${CMAKE_CURRENT_BINARY_DIR}/perm_case1.dat")
  if (NOT EXISTS ${SPE10MODEL1DATA_TARGET_FILENAME})
    execute_process(COMMAND ln -s "${SPE10MODEL1DATA}" "${SPE10MODEL1DATA_TARGET_FILENAME}")
  endif(NOT EXISTS ${SPE10MODEL1DATA_TARGET_FILENAME})
endif (NOT ${SPE10MODEL1DATA} STREQUAL "SPE10MODEL1DATA-NOTFOUND")
