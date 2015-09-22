# This file is part of the dune-hdd project:
#   https://github.com/pymor/dune-hdd
# Copyright Holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ADD_EXECUTABLE(example_linearelliptic_cg "cg_main.cc" ${COMMON_HEADER})
#TARGET_LINK_LIBRARIES(example_linearelliptic_cg
#                      dunehdd-linearelliptic-cg
#                      ${COMMON_LIBS} ${Boost_FILESYSTEM_LIBRARY})

add_python_bindings(linearparabolicthermalblockexample
                    thermalblock_bindings_generator.py
                    thermalblock_bindings_generator.hh
                    "")

foreach(target
        linearparabolicthermalblockexample
        )
  add_dune_parmetis_flags(${target})
  add_dune_superlu_flags(${target})
  add_dune_umfpack_flags(${target})
  add_dune_tbb_flags(${target})
  add_dune_alugrid_flags(${target})
  add_dune_mpi_flags(${target})
endforeach(target)

# link spe10 data file if present
if (NOT ${SPE10MODEL1DATA} STREQUAL "SPE10MODEL1DATA-NOTFOUND")
  set (SPE10MODEL1DATA_TARGET_FILENAME "${CMAKE_CURRENT_BINARY_DIR}/perm_case1.dat")
  if (NOT EXISTS ${SPE10MODEL1DATA_TARGET_FILENAME})
    execute_process(COMMAND ln -s "${SPE10MODEL1DATA}" "${SPE10MODEL1DATA_TARGET_FILENAME}")
  endif(NOT EXISTS ${SPE10MODEL1DATA_TARGET_FILENAME})
endif (NOT ${SPE10MODEL1DATA} STREQUAL "SPE10MODEL1DATA-NOTFOUND")
if (NOT ${SPE10MODEL2DATA} STREQUAL "SPE10MODEL2DATA-NOTFOUND")
  set (SPE10MODEL2DATA_TARGET_FILENAME "${CMAKE_CURRENT_BINARY_DIR}/perm_case2a.dat")
  if (NOT EXISTS ${SPE10MODEL2DATA_TARGET_FILENAME})
    execute_process(COMMAND ln -s "${SPE10MODEL2DATA}" "${SPE10MODEL2DATA_TARGET_FILENAME}")
  endif(NOT EXISTS ${SPE10MODEL2DATA_TARGET_FILENAME})
endif (NOT ${SPE10MODEL2DATA} STREQUAL "SPE10MODEL2DATA-NOTFOUND")

foreach(file
        "thermalblock_main.py"
        "thermalblock_prepare.py")
  if (NOT EXISTS ${file})
    execute_process(COMMAND ln -s "${CMAKE_CURRENT_SOURCE_DIR}/${file}" "${CMAKE_CURRENT_BINARY_DIR}/${file}")
  endif (NOT EXISTS ${file})
endforeach(file)
