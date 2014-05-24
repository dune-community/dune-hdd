// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#   include <dune/stuff/common/disable_warnings.hh>
#if HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
#   include <dune/stuff/common/reenable_warnings.hh>
# define ENABLE_ALUGRID 1
# include <dune/grid/alugrid.hh>
#endif // HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
#include <dune/stuff/common/disable_warnings.hh>
# include <dune/grid/sgrid.hh>
#include <dune/stuff/common/reenable_warnings.hh>

# include "cg.hh"
#include <dune/stuff/common/reenable_warnings.hh> // <- here for the python bindings!


template class LinearellipticExampleCG< Dune::SGrid< 1, 1 > >;
template class LinearellipticExampleCG< Dune::SGrid< 2, 2 > >;
template class LinearellipticExampleCG< Dune::SGrid< 3, 3 > >;

#   include <dune/stuff/common/disable_warnings.hh>
#if HAVE_ALUGRID
#   include <dune/stuff/common/reenable_warnings.hh>

template class LinearellipticExampleCG< Dune::ALUConformGrid< 2, 2 > >;

#endif // HAVE_ALUGRID
