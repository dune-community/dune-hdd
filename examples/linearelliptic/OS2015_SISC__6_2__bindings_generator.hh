// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_EXAMPLES_LINEARELLIPTIC_OS2015_SISC_6_2_BINDINGS_GENERATOR_HH
#define DUNE_HDD_EXAMPLES_LINEARELLIPTIC_OS2015_SISC_6_2_BINDINGS_GENERATOR_HH

#ifdef ENABLE_MPI
# undef ENABLE_MPI
#endif
#define ENABLE_MPI 0

#ifdef ENABLE_PARMETIS
# undef ENABLE_PARMETIS
#endif
#define ENABLE_PARMETIS 0

#include "config.h"

#include "OS2015_SISC__6_2.hh"

#include <dune/stuff/common/disable_warnings.hh>
# if HAVE_ALUGRID
#   include <dune/grid/alugrid.hh>
# endif

#endif // DUNE_HDD_EXAMPLES_LINEARELLIPTIC_OS2015_SISC_6_2_BINDINGS_GENERATOR_HH
