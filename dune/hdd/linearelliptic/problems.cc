// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#if HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
# define ENABLE_ALUGRID 1
# include <dune/grid/alugrid.hh>
#endif // HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H

#include "problems.hh"

namespace Dune {
namespace HDD {
namespace internal {


void lib_exists() {}


} // namespace internal
namespace LinearElliptic {

#if HAVE_ALUGRID


template class ProblemProvider< typename ALUGrid< 2, 2, Dune::simplex, Dune::conforming >::template Codim< 0 >::Entity,
                                double, 2, double, 1 >;


#endif // HAVE_ALUGRID

} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune
