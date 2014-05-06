// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

# include "cg.hh"
#include <dune/stuff/common/reenable_warnings.hh> // <- here for the python bindings!

// sgrid
#include <dune/stuff/common/disable_warnings.hh>
# include <dune/grid/sgrid.hh>
#include <dune/stuff/common/reenable_warnings.hh>

template class LinearellipticExampleCG< Dune::SGrid< 1, 1 > >;
template class LinearellipticExampleCG< Dune::SGrid< 2, 2 > >;
template class LinearellipticExampleCG< Dune::SGrid< 3, 3 > >;

// yaspgrid
#include <dune/grid/yaspgrid.hh>

template class LinearellipticExampleCG< Dune::YaspGrid< 1 > >;
template class LinearellipticExampleCG< Dune::YaspGrid< 2 > >;
template class LinearellipticExampleCG< Dune::YaspGrid< 3 > >;
