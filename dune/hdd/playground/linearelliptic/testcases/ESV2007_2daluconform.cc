// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include "ESV2007.hh"

#if HAVE_ALUGRID

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace TestCases {


template class ESV2007< ALUGrid< 2, 2, simplex, conforming > >;


} // namespace TestCases
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // HAVE_ALUGRID
