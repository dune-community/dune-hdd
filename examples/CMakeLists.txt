# This file is part of the dune-pymor project:
#   https://github.com/pyMor/dune-pymor
# Copyright Holders: Felix Albrecht, Stephan Rave
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

add_subdirectory(linearelliptic)

ADD_EXECUTABLE(example_empty "empty.cc" ${COMMON_HEADER})
TARGET_LINK_LIBRARIES(example_empty ${COMMON_LIBS})
add_dune_alugrid_flags(example_empty)
