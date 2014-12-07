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
add_library(obj_test_linearelliptic-block-swipdg-expectations_os2014_2daluconform
            OBJECT
            EXCLUDE_FROM_ALL
            linearelliptic-block-swipdg-expectations_os2014_2daluconform.cxx)
add_library(obj_test_linearelliptic-block-swipdg-expectations_esv2007_2daluconform
            OBJECT
            EXCLUDE_FROM_ALL
            linearelliptic-block-swipdg-expectations_esv2007_2daluconform.cxx)
add_library(obj_test_linearelliptic-block-swipdg-expectations_spe10_2daluconform
            OBJECT
            EXCLUDE_FROM_ALL
            linearelliptic-block-swipdg-expectations_spe10_2daluconform.cxx)

add_library(obj_test_linearelliptic-swipdg-expectations_esv2007_2daluconform
            OBJECT
            EXCLUDE_FROM_ALL
            linearelliptic-swipdg-expectations_esv2007_2daluconform.cxx)
add_library(obj_test_linearelliptic-swipdg-expectations_esv2007_2dsgrid
            OBJECT
            EXCLUDE_FROM_ALL
            linearelliptic-swipdg-expectations_esv2007_2dsgrid.cxx)
add_library(obj_test_linearelliptic-swipdg-expectations_spe10_2daluconform
            OBJECT
            EXCLUDE_FROM_ALL
            linearelliptic-swipdg-expectations_spe10_2daluconform.cxx)
add_library(obj_test_linearelliptic-swipdg-expectations_spe10_2dsgrid
            OBJECT
            EXCLUDE_FROM_ALL
            linearelliptic-swipdg-expectations_spe10_2dsgrid.cxx)

# test target libs
add_library(obj_test_os2014_nonparametric_convergence_study
            OBJECT
            EXCLUDE_FROM_ALL
            OS2014_nonparametric_convergence_study.cxx)
add_library(lib_test_OS2014_nonparametric_convergence_study
            STATIC
            EXCLUDE_FROM_ALL
            $<TARGET_OBJECTS:obj_test_os2014_nonparametric_convergence_study>
            $<TARGET_OBJECTS:obj_test_linearelliptic-swipdg-expectations_esv2007_2daluconform>
            $<TARGET_OBJECTS:obj_test_linearelliptic-block-swipdg-expectations_os2014_2daluconform>
            $<TARGET_OBJECTS:obj_test_linearelliptic-block-swipdg-expectations_esv2007_2daluconform>)
add_dune_alugrid_flags(lib_test_OS2014_nonparametric_convergence_study)

add_library(lib_test_OS2014_parametric_convergence_study
            STATIC
            EXCLUDE_FROM_ALL
            $<TARGET_OBJECTS:obj_test_linearelliptic-block-swipdg-expectations_os2014_2daluconform>)
add_dune_alugrid_flags(lib_test_OS2014_parametric_convergence_study)

add_library(lib_test_OS2014_parametric_multiscale_convergence_study
            STATIC
            EXCLUDE_FROM_ALL
            $<TARGET_OBJECTS:obj_test_linearelliptic-block-swipdg-expectations_spe10_2daluconform>)
add_dune_alugrid_flags(lib_test_OS2014_parametric_multiscale_convergence_study)

add_library(lib_test_linearelliptic-swipdg
            STATIC
            EXCLUDE_FROM_ALL
            $<TARGET_OBJECTS:obj_test_linearelliptic-swipdg-expectations_esv2007_2daluconform>
            $<TARGET_OBJECTS:obj_test_linearelliptic-swipdg-expectations_esv2007_2dsgrid>
            $<TARGET_OBJECTS:obj_test_linearelliptic-swipdg-expectations_spe10_2daluconform>
            $<TARGET_OBJECTS:obj_test_linearelliptic-swipdg-expectations_spe10_2dsgrid>)
add_dune_alugrid_flags(lib_test_linearelliptic-swipdg)

add_library(lib_test_linearelliptic-block-swipdg
            STATIC
            EXCLUDE_FROM_ALL
            $<TARGET_OBJECTS:obj_test_linearelliptic-block-swipdg-expectations_esv2007_2daluconform>)

add_library(lib_test_OS2014-FVCA7-poster
            STATIC
            EXCLUDE_FROM_ALL
            $<TARGET_OBJECTS:obj_test_linearelliptic-block-swipdg-expectations_esv2007_2daluconform>)

# test targets
target_link_libraries(test_OS2014_parametric_convergence_study lib_test_OS2014_parametric_convergence_study)
target_link_libraries(test_OS2014_parametric_multiscale_convergence_study lib_test_OS2014_parametric_multiscale_convergence_study)
target_link_libraries(test_OS2014_nonparametric_convergence_study lib_test_OS2014_nonparametric_convergence_study)
target_link_libraries(test_linearelliptic-block-swipdg lib_test_linearelliptic-block-swipdg)
target_link_libraries(test_OS2014-FVCA7-poster lib_test_linearelliptic-block-swipdg)
target_link_libraries(test_linearelliptic-swipdg lib_test_linearelliptic-swipdg)

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
