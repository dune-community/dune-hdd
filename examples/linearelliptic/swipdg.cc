// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include "swipdg.hh"

#   include <dune/stuff/common/disable_warnings.hh>
#if HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
#   include <dune/stuff/common/reenable_warnings.hh>
# define ENABLE_ALUGRID 1
# include <dune/grid/alugrid.hh>
#endif // HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
#include <dune/stuff/common/disable_warnings.hh>
# include <dune/grid/sgrid.hh>
#include <dune/stuff/common/reenable_warnings.hh>


template class LinearellipticExampleSWIPDG< Dune::SGrid< 1, 1 > >;

#if HAVE_ALUGRID

template class LinearellipticExampleSWIPDG< Dune::ALUConformGrid< 2, 2 > >;

#endif // HAVE_ALUGRID
