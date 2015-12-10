# This file is part of the dune-hdd project:
#   https://github.com/pymor/dune-hdd
# Copyright Holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

add_subdirectory(linearelliptic)
add_subdirectory(linearparabolic)

ADD_EXECUTABLE(example_empty "empty.cc" ${COMMON_HEADER})
TARGET_LINK_LIBRARIES(example_empty ${COMMON_LIBS})
add_dune_alugrid_flags(example_empty)
