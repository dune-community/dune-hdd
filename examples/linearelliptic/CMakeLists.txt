# This file is part of the dune-pymor project:
#   https://github.com/pyMor/dune-pymor
# Copyright Holders: Felix Albrecht, Stephan Rave
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

add_library(dunehdd-linearelliptic-cg
            STATIC
            "cg.cc")

add_library(dunehdd-linearelliptic-swipdg
            STATIC
            "swipdg.cc")

add_library(dunehdd-linearelliptic-thermalblock
            STATIC
            "thermalblock.cc")


ADD_EXECUTABLE(example_linearelliptic_block_swipdg "block-swipdg_main.cc" ${COMMON_HEADER})
TARGET_LINK_LIBRARIES(example_linearelliptic_block_swipdg
                      ${COMMON_LIBS} ${Boost_FILESYSTEM_LIBRARY})

ADD_EXECUTABLE(example_linearelliptic_cg "cg_main.cc" ${COMMON_HEADER})
TARGET_LINK_LIBRARIES(example_linearelliptic_cg
                      dunehdd-linearelliptic-cg
                      ${COMMON_LIBS} ${Boost_FILESYSTEM_LIBRARY})

ADD_EXECUTABLE(example_linearelliptic_swipdg "swipdg_main.cc" ${COMMON_HEADER})
TARGET_LINK_LIBRARIES(example_linearelliptic_swipdg
                      dunehdd-linearelliptic-swipdg
                      ${COMMON_LIBS} ${Boost_FILESYSTEM_LIBRARY})

ADD_EXECUTABLE(example_linearelliptic_thermalblock "thermalblock_main.cc" ${COMMON_HEADER})
TARGET_LINK_LIBRARIES(example_linearelliptic_thermalblock
                      dunehdd-linearelliptic-thermalblock
                      ${COMMON_LIBS} ${Boost_FILESYSTEM_LIBRARY})

add_python_bindings(linearellipticexamplecg cg_bindings_generator.py cg.hh dunehdd-linearelliptic-cg)
add_python_bindings(linearellipticexamplethermalblock thermalblock_bindings_generator.py thermalblock.hh dunehdd-linearelliptic-thermalblock)

if(NOT CMAKE_WITH_AUTOTOOLS)
  add_dune_alugrid_flags(example_linearelliptic_block_swipdg)
  add_dune_alugrid_flags(example_linearelliptic_cg)
  add_dune_alugrid_flags(example_linearelliptic_swipdg)
  add_dune_alugrid_flags(example_linearelliptic_thermalblock)
endif(NOT CMAKE_WITH_AUTOTOOLS)


# link spe10 data file if present
if (NOT ${SPE10MODEL1DATA} STREQUAL "SPE10MODEL1DATA-NOTFOUND")
  set (SPE10MODEL1DATA_TARGET_FILENAME "${CMAKE_CURRENT_BINARY_DIR}/perm_case1.dat")
  if (NOT EXISTS ${SPE10MODEL1DATA_TARGET_FILENAME})
    execute_process(COMMAND ln -s "${SPE10MODEL1DATA}" "${SPE10MODEL1DATA_TARGET_FILENAME}")
  endif(NOT EXISTS ${SPE10MODEL1DATA_TARGET_FILENAME})
endif (NOT ${SPE10MODEL1DATA} STREQUAL "SPE10MODEL1DATA-NOTFOUND")
