// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#if HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
# define ENABLE_ALUGRID 1
# include <dune/grid/alugrid.hh>
#endif

# include "thermalblock.hh"
#include <dune/stuff/common/reenable_warnings.hh> // <- here for the python bindings!

#if HAVE_ALUGRID


template class ThermalblockExample< Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming > >;


#endif // HAVE_ALUGRID
