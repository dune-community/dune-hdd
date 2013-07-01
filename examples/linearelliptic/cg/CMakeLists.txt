# This file is part of the dune-hdd project:
#   http://users.dune-project.org/projects/dune-hdd
# Copyright holders: Felix Albrecht
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

ADD_EXECUTABLE( examples_linearelliptic_cg_gdt "gdt.cc" ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES( examples_linearelliptic_cg_gdt ${COMMON_LIBS} )

