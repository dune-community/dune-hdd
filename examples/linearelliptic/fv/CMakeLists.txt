# This file is part of the dune-hdd project:
#   http://users.dune-project.org/projects/dune-hdd
# Copyright holders: Felix Albrecht
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

ADD_EXECUTABLE( examples_linearelliptic_fv_pdelab "pdelab.cc" ${COMMON_HEADERS} )
TARGET_LINK_LIBRARIES( examples_linearelliptic_fv_pdelab ${COMMON_LIBS} )

