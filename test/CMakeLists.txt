# This file is part of the dune-hdd project:
#   http://users.dune-project.org/projects/dune-hdd
# Copyright holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

add_library(obj_test_linearelliptic-block-swipdg-expectations_os2014_2daluconform
            OBJECT
            EXCLUDE_FROM_ALL
            linearelliptic-block-swipdg-expectations_os2014_2daluconform.cxx)
add_library(obj_test_linearelliptic-block-swipdg-expectations_esv2007_2daluconform
            OBJECT
            EXCLUDE_FROM_ALL
            linearelliptic-block-swipdg-expectations_esv2007_2daluconform.cxx)
add_library(obj_test_os2014
            OBJECT
            EXCLUDE_FROM_ALL
            OS2014.cxx)
add_library(lib_test_OS2014
            STATIC
            EXCLUDE_FROM_ALL
            $<TARGET_OBJECTS:obj_test_os2014>
            $<TARGET_OBJECTS:obj_test_linearelliptic-block-swipdg-expectations_os2014_2daluconform>
            $<TARGET_OBJECTS:obj_test_linearelliptic-block-swipdg-expectations_esv2007_2daluconform>)
add_dune_alugrid_flags(lib_test_OS2014)

add_library(lib_test_linearelliptic-block-swipdg
            STATIC
            EXCLUDE_FROM_ALL
            $<TARGET_OBJECTS:obj_test_linearelliptic-block-swipdg-expectations_esv2007_2daluconform>)

add_library(lib_test_OS2014-FVCA7-poster
            STATIC
            EXCLUDE_FROM_ALL
            $<TARGET_OBJECTS:obj_test_linearelliptic-block-swipdg-expectations_esv2007_2daluconform>)

enable_testing()

DEPENDENCYCHECK()

BEGIN_TESTCASES()
END_TESTCASES()

target_link_libraries(test_OS2014 lib_test_OS2014)
target_link_libraries(test_linearelliptic-block-swipdg lib_test_linearelliptic-block-swipdg)
target_link_libraries(test_OS2014-FVCA7-poster lib_test_linearelliptic-block-swipdg)

add_analyze(OS2014.cc)
target_link_libraries(analyze_OS2014.cc lib_test_OS2014)

# link spe10 data file if present
if (NOT ${SPE10MODEL1DATA} STREQUAL "SPE10MODEL1DATA-NOTFOUND")
  set (SPE10MODEL1DATA_TARGET_FILENAME "${CMAKE_CURRENT_BINARY_DIR}/perm_case1.dat")
  if (NOT EXISTS ${SPE10MODEL1DATA_TARGET_FILENAME})
    execute_process(COMMAND ln -s "${SPE10MODEL1DATA}" "${SPE10MODEL1DATA_TARGET_FILENAME}")
  endif(NOT EXISTS ${SPE10MODEL1DATA_TARGET_FILENAME})
endif (NOT ${SPE10MODEL1DATA} STREQUAL "SPE10MODEL1DATA-NOTFOUND")
