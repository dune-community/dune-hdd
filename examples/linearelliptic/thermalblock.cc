// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include "thermalblock.hh"


// sgrid
#include <dune/stuff/common/disable_warnings.hh>
# include <dune/grid/sgrid.hh>
#include <dune/stuff/common/reenable_warnings.hh>

template class ThermalblockExample< Dune::SGrid< 1, 1 > >;
template class ThermalblockExample< Dune::SGrid< 2, 2 > >;
template class ThermalblockExample< Dune::SGrid< 3, 3 > >;

// yaspgrid
#include <dune/grid/yaspgrid.hh>

template class ThermalblockExample< Dune::YaspGrid< 1 > >;
template class ThermalblockExample< Dune::YaspGrid< 2 > >;
template class ThermalblockExample< Dune::YaspGrid< 3 > >;
