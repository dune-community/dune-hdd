# This file is part of the dune-hdd project:
#   https://github.com/pymor/dune-hdd
# Copyright Holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

add_executable(example_linearelliptic_kuehlergrill "kuehlergrill.cc" ${COMMON_HEADER})
target_link_libraries(example_linearelliptic_kuehlergrill ${COMMON_LIBS})

add_executable(example_linearelliptic_OS2014_FVCA7__visualization "OS2014_FVCA7__visualization.cc" ${COMMON_HEADER})
target_link_libraries(example_linearelliptic_OS2014_FVCA7__visualization ${COMMON_LIBS})

add_executable(example_linearelliptic_OS2015_SISC__6_1__multiscale_example__visualizations
               "OS2015_SISC__6_1__multiscale_example__visualizations.cc" ${COMMON_HEADER})
target_link_libraries(example_linearelliptic_OS2015_SISC__6_1__multiscale_example__visualizations
                      ${COMMON_LIBS} ${Boost_FILESYSTEM_LIBRARY})

add_python_bindings(linearellipticexamplegenericmultiscale
                    generic_multiscale_bindings_generator.py
                    generic_multiscale_bindings_generator.hh "")
add_python_bindings(linearellipticexampleMRS201643
                    MRS2016__4_3_bindings_generator.py
                    MRS2016__4_3_bindings_generator.hh "")
add_python_bindings(linearellipticexampleMRS201651
                    MRS2016__5_1_bindings_generator.py
                    MRS2016__5_1_bindings_generator.hh "")
add_python_bindings(linearellipticexampleOS2015SISC62
                    OS2015_SISC__6_2__bindings_generator.py
                    OS2015_SISC__6_2__bindings_generator.hh "")

foreach(target
        example_linearelliptic_kuehlergrill
        example_linearelliptic_OS2014_FVCA7__visualization
        example_linearelliptic_OS2015_SISC__6_1__multiscale_example__visualizations
        linearellipticexamplegenericmultiscale
        linearellipticexampleMRS201651
        linearellipticexampleMRS201643
        linearellipticexampleOS2015SISC62)
  add_dune_alugrid_flags(${target})
  add_dune_mpi_flags(${target})
  add_dune_parmetis_flags(${target})
  add_dune_superlu_flags(${target})
  add_dune_tbb_flags(${target})
  add_dune_umfpack_flags(${target})
endforeach(target)

foreach(file
        "generic_multiscale.py"
        "OS2015_SISC__6_2__academic_example.py"
        "OS2015_SISC__6_2__multiscale_example.py")
  if (NOT EXISTS ${file})
    execute_process(COMMAND ln -s "${CMAKE_CURRENT_SOURCE_DIR}/${file}" "${CMAKE_CURRENT_BINARY_DIR}/${file}")
  endif (NOT EXISTS ${file})
endforeach(file)

# link spe10 data file if present
if (NOT ${SPE10MODEL1DATA} STREQUAL "SPE10MODEL1DATA-NOTFOUND")
  set (SPE10MODEL1DATA_TARGET_FILENAME "${CMAKE_CURRENT_BINARY_DIR}/perm_case1.dat")
  if (NOT EXISTS ${SPE10MODEL1DATA_TARGET_FILENAME})
    execute_process(COMMAND ln -s "${SPE10MODEL1DATA}" "${SPE10MODEL1DATA_TARGET_FILENAME}")
  endif(NOT EXISTS ${SPE10MODEL1DATA_TARGET_FILENAME})
endif (NOT ${SPE10MODEL1DATA} STREQUAL "SPE10MODEL1DATA-NOTFOUND")
if (NOT ${SPE10MODEL2DATA} STREQUAL "SPE10MODEL2DATA-NOTFOUND")
  set (SPE10MODEL2DATA_TARGET_FILENAME "${CMAKE_CURRENT_BINARY_DIR}/spe_perm.dat")
  if (NOT EXISTS ${SPE10MODEL2DATA_TARGET_FILENAME})
    execute_process(COMMAND ln -s "${SPE10MODEL2DATA}" "${SPE10MODEL2DATA_TARGET_FILENAME}")
  endif(NOT EXISTS ${SPE10MODEL2DATA_TARGET_FILENAME})
endif (NOT ${SPE10MODEL2DATA} STREQUAL "SPE10MODEL2DATA-NOTFOUND")

