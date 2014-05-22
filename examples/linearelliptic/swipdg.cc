// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#ifdef HAVE_DUNE_GRID_MULTISCALE
#undef HAVE_DUNE_GRID_MULTISCALE
#endif

#include "swipdg.hh"

#if HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
# define ENABLE_ALUGRID 1
# include <dune/grid/alugrid.hh>

template class LinearellipticExampleSWIPDG< Dune::ALUConformGrid< 2, 2 > >;

#endif // HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
