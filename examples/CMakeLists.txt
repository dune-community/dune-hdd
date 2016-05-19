# This file is part of the dune-hdd project:
#   https://github.com/pymor/dune-hdd
# Copyright Holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

add_subdirectory(linearelliptic)
add_subdirectory(linearparabolic)

add_executable(example_empty "empty.cc" ${COMMON_HEADER})
target_link_libraries(example_empty ${COMMON_LIBS} dunehdd)

foreach(target
        example_empty)
  add_dune_parmetis_flags(${target})
  add_dune_superlu_flags(${target})
  add_dune_umfpack_flags(${target})
  add_dune_tbb_flags(${target})
  add_dune_alugrid_flags(${target})
  add_dune_mpi_flags(${target})
endforeach(target)
