# This file is part of the dune-hdd project:
#   https://github.com/pymor/dune-hdd
# Copyright Holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

add_library(dunehdd-linearelliptic-cg
            STATIC
            "cg.cc")

#add_library(dunehdd-linearelliptic-swipdg
#            STATIC
#            "swipdg.cc")


#ADD_EXECUTABLE(example_linearelliptic_block_swipdg "block-swipdg_main.cc" ${COMMON_HEADER})
#TARGET_LINK_LIBRARIES(example_linearelliptic_block_swipdg
#                      ${COMMON_LIBS} ${Boost_FILESYSTEM_LIBRARY})

ADD_EXECUTABLE(example_linearelliptic_cg "cg_main.cc" ${COMMON_HEADER})
TARGET_LINK_LIBRARIES(example_linearelliptic_cg
                      dunehdd-linearelliptic-cg
                      ${COMMON_LIBS} ${Boost_FILESYSTEM_LIBRARY})

ADD_EXECUTABLE(mpi_cg_main "mpi_cg_main.cc" ${COMMON_HEADER})
TARGET_LINK_LIBRARIES(mpi_cg_main
                      dunehdd-linearelliptic-cg 
                      ${COMMON_LIBS} )


#ADD_EXECUTABLE(example_linearelliptic_swipdg "swipdg_main.cc" ${COMMON_HEADER})
#TARGET_LINK_LIBRARIES(example_linearelliptic_swipdg
#                      dunehdd-linearelliptic-swipdg
#                      ${COMMON_LIBS} ${Boost_FILESYSTEM_LIBRARY})

ADD_EXECUTABLE(example_linearelliptic_thermalblock "thermalblock_main.cc" ${COMMON_HEADER})
TARGET_LINK_LIBRARIES(example_linearelliptic_thermalblock
                     ${COMMON_LIBS} ${Boost_FILESYSTEM_LIBRARY})

#ADD_EXECUTABLE(example_linearelliptic_OS2014_spe10_visualization
#               "OS2014_spe10_visualization.cc" ${COMMON_HEADER})
#TARGET_LINK_LIBRARIES(example_linearelliptic_OS2014_spe10_visualization
#                      ${COMMON_LIBS} ${Boost_FILESYSTEM_LIBRARY})

add_python_bindings(linearellipticexamplecg
                    cg_bindings_generator.py
                    cg.hh
                    dunehdd-linearelliptic-cg)
add_python_bindings(linearellipticexamplethermalblock
                    thermalblock_bindings_generator.py
                    thermalblock_bindings_generator.hh
                    "")
add_python_bindings(mpi_cg_bindings
                    mpi_cg_bindings_generator.py
                    mpi_cg_example.hh
                    dunehdd-linearelliptic-cg)
add_python_bindings(linearellipticexampleOS2014
                    OS2014_bindings_generator.py
                    OS2014_bindings_generator.hh
                    "")

#        example_linearelliptic_block_swipdg
#        example_linearelliptic_OS2014_spe10_visualization
#        example_linearelliptic_swipdg
#        example_linearelliptic_thermalblock
#        linearellipticexampleOS2014)
foreach(target
        mpi_cg_main
        mpi_cg_bindings
        example_linearelliptic_cg
        linearellipticexamplecg
        linearellipticexamplethermalblock
        linearellipticexampleOS2014
        )
  add_dune_parmetis_flags(${target})
  add_dune_superlu_flags(${target})
  add_dune_umfpack_flags(${target})
  add_dune_tbb_flags(${target})
  add_dune_alugrid_flags(${target})
  add_dune_mpi_flags(${target})
endforeach(_target)

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
foreach(_file "cg_main.py" "OS2014_academic.py" "OS2014_estimators.py" "OS2014_offline.py" "OS2014_online.py" "OS2014_prepare.py" "OS2014_spe10.py" "thermalblock_main.py")
  if (NOT EXISTS ${_file})
    execute_process(COMMAND ln -s "${CMAKE_CURRENT_SOURCE_DIR}/${_file}" "${CMAKE_CURRENT_BINARY_DIR}//${_file}")
  endif (NOT EXISTS ${_file})
endforeach(_file)
