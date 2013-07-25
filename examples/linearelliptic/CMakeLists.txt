# This file is part of the dune-pymor project:
#   https://github.com/pyMor/dune-pymor
# Copyright Holders: Felix Albrecht, Stephan Rave
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

add_library( dunehdd-example-linearelliptic-cg-with-dune-gdt STATIC "cg-with-dune-gdt.cc" )

ADD_EXECUTABLE( example_linearelliptic_cg_with_dune_gdt "cg-with-dune-gdt_main.cc" ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES( example_linearelliptic_cg_with_dune_gdt dunehdd-example-linearelliptic-cg-with-dune-gdt ${COMMON_LIBS} )
