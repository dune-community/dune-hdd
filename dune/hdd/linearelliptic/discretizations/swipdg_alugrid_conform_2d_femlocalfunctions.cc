// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#if HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
# define ENABLE_ALUGRID 1
# include <dune/grid/alugrid.hh>
#endif // HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H

#include "swipdg.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Discretizations {

#if HAVE_ALUGRID && HAVE_DUNE_FEM_LOCALFUNCTIONS


template class SWIPDG< ALUConformGrid< 2, 2 >,
                   Stuff::Grid::ChooseLayer::leaf,
                   double, 1, 1, GDT::ChooseSpaceBackend::fem_localfunction,
                   Stuff::LA::ChooseBackend::common_dense >;

template class SWIPDG< ALUConformGrid< 2, 2 >,
                   Stuff::Grid::ChooseLayer::level,
                   double, 1, 1, GDT::ChooseSpaceBackend::fem_localfunction,
                   Stuff::LA::ChooseBackend::common_dense >;

# if HAVE_DUNE_ISTL

template class SWIPDG< ALUConformGrid< 2, 2 >,
                   Stuff::Grid::ChooseLayer::leaf,
                   double, 1, 1, GDT::ChooseSpaceBackend::fem_localfunction,
                   Stuff::LA::ChooseBackend::istl_sparse >;

template class SWIPDG< ALUConformGrid< 2, 2 >,
                   Stuff::Grid::ChooseLayer::level,
                   double, 1, 1, GDT::ChooseSpaceBackend::fem_localfunction,
                   Stuff::LA::ChooseBackend::istl_sparse >;

# endif // HAVE_DUNE_ISTL
# if HAVE_EIGEN

template class SWIPDG< ALUConformGrid< 2, 2 >,
                   Stuff::Grid::ChooseLayer::leaf,
                   double, 1, 1, GDT::ChooseSpaceBackend::fem_localfunction,
                   Stuff::LA::ChooseBackend::eigen_dense >;

template class SWIPDG< ALUConformGrid< 2, 2 >,
                   Stuff::Grid::ChooseLayer::level,
                   double, 1, 1, GDT::ChooseSpaceBackend::fem_localfunction,
                   Stuff::LA::ChooseBackend::eigen_dense >;

template class SWIPDG< ALUConformGrid< 2, 2 >,
                   Stuff::Grid::ChooseLayer::leaf,
                   double, 1, 1, GDT::ChooseSpaceBackend::fem_localfunction,
                   Stuff::LA::ChooseBackend::eigen_sparse >;

template class SWIPDG< ALUConformGrid< 2, 2 >,
                   Stuff::Grid::ChooseLayer::level,
                   double, 1, 1, GDT::ChooseSpaceBackend::fem_localfunction,
                   Stuff::LA::ChooseBackend::eigen_sparse >;


# endif // HAVE_EIGEN
#endif // HAVE_ALUGRID && HAVE_DUNE_FEM_LOCALFUNCTIONS

} // namespace Discretizations
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune
