// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include "OS2014.hh"

#if HAVE_ALUGRID && HAVE_DUNE_GRID_MULTISCALE

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace TestCases {


template class OS2014Multiscale< ALUGrid< 2, 2, simplex, conforming > >;


} // namespace TestCases
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // HAVE_ALUGRID && HAVE_DUNE_GRID_MULTISCALE
