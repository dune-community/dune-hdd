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
namespace Problems {


template class ESV2007< typename ALUGrid< 2, 2, simplex, conforming >::template Codim< 0 >::Entity,
                        double, 2, double, 1 >;


} // namespace Problems
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // HAVE_ALUGRID
