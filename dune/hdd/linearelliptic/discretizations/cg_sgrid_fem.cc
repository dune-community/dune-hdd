// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include <dune/stuff/common/disable_warnings.hh>
# include <dune/grid/sgrid.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#include "cg.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Discretizations {

#if HAVE_DUNE_FEM


template class CG< Dune::SGrid< 1, 1 >, Stuff::Grid::ChooseLayer::leaf, double, 1, 1,
                   GDT::ChooseSpaceBackend::fem, Stuff::LA::ChooseBackend::common_dense >;
template class CG< Dune::SGrid< 2, 2 >, Stuff::Grid::ChooseLayer::leaf, double, 1, 1,
                   GDT::ChooseSpaceBackend::fem, Stuff::LA::ChooseBackend::common_dense >;
template class CG< Dune::SGrid< 3, 3 >, Stuff::Grid::ChooseLayer::leaf, double, 1, 1,
                   GDT::ChooseSpaceBackend::fem, Stuff::LA::ChooseBackend::common_dense >;

template class CG< Dune::SGrid< 1, 1 >, Stuff::Grid::ChooseLayer::level, double, 1, 1,
                   GDT::ChooseSpaceBackend::fem, Stuff::LA::ChooseBackend::common_dense >;
template class CG< Dune::SGrid< 2, 2 >, Stuff::Grid::ChooseLayer::level, double, 1, 1,
                   GDT::ChooseSpaceBackend::fem, Stuff::LA::ChooseBackend::common_dense >;
template class CG< Dune::SGrid< 3, 3 >, Stuff::Grid::ChooseLayer::level, double, 1, 1,
                   GDT::ChooseSpaceBackend::fem, Stuff::LA::ChooseBackend::common_dense >;

# if HAVE_DUNE_ISTL

template class CG< Dune::SGrid< 1, 1 >, Stuff::Grid::ChooseLayer::leaf, double, 1, 1,
                   GDT::ChooseSpaceBackend::fem, Stuff::LA::ChooseBackend::istl_sparse >;
template class CG< Dune::SGrid< 2, 2 >, Stuff::Grid::ChooseLayer::leaf, double, 1, 1,
                   GDT::ChooseSpaceBackend::fem, Stuff::LA::ChooseBackend::istl_sparse >;
template class CG< Dune::SGrid< 3, 3 >, Stuff::Grid::ChooseLayer::leaf, double, 1, 1,
                   GDT::ChooseSpaceBackend::fem, Stuff::LA::ChooseBackend::istl_sparse >;

template class CG< Dune::SGrid< 1, 1 >, Stuff::Grid::ChooseLayer::level, double, 1, 1,
                   GDT::ChooseSpaceBackend::fem, Stuff::LA::ChooseBackend::istl_sparse >;
template class CG< Dune::SGrid< 2, 2 >, Stuff::Grid::ChooseLayer::level, double, 1, 1,
                   GDT::ChooseSpaceBackend::fem, Stuff::LA::ChooseBackend::istl_sparse >;
template class CG< Dune::SGrid< 3, 3 >, Stuff::Grid::ChooseLayer::level, double, 1, 1,
                   GDT::ChooseSpaceBackend::fem, Stuff::LA::ChooseBackend::istl_sparse >;

# endif // HAVE_DUNE_ISTL
# if HAVE_EIGEN

template class CG< Dune::SGrid< 1, 1 >, Stuff::Grid::ChooseLayer::leaf, double, 1, 1,
                   GDT::ChooseSpaceBackend::fem, Stuff::LA::ChooseBackend::eigen_dense >;
template class CG< Dune::SGrid< 2, 2 >, Stuff::Grid::ChooseLayer::leaf, double, 1, 1,
                   GDT::ChooseSpaceBackend::fem, Stuff::LA::ChooseBackend::eigen_dense >;
template class CG< Dune::SGrid< 3, 3 >, Stuff::Grid::ChooseLayer::leaf, double, 1, 1,
                   GDT::ChooseSpaceBackend::fem, Stuff::LA::ChooseBackend::eigen_dense >;

template class CG< Dune::SGrid< 1, 1 >, Stuff::Grid::ChooseLayer::level, double, 1, 1,
                   GDT::ChooseSpaceBackend::fem, Stuff::LA::ChooseBackend::eigen_dense >;
template class CG< Dune::SGrid< 2, 2 >, Stuff::Grid::ChooseLayer::level, double, 1, 1,
                   GDT::ChooseSpaceBackend::fem, Stuff::LA::ChooseBackend::eigen_dense >;
template class CG< Dune::SGrid< 3, 3 >, Stuff::Grid::ChooseLayer::level, double, 1, 1,
                   GDT::ChooseSpaceBackend::fem, Stuff::LA::ChooseBackend::eigen_dense >;

template class CG< Dune::SGrid< 1, 1 >, Stuff::Grid::ChooseLayer::leaf, double, 1, 1,
                   GDT::ChooseSpaceBackend::fem, Stuff::LA::ChooseBackend::eigen_sparse >;
template class CG< Dune::SGrid< 2, 2 >, Stuff::Grid::ChooseLayer::leaf, double, 1, 1,
                   GDT::ChooseSpaceBackend::fem, Stuff::LA::ChooseBackend::eigen_sparse >;
template class CG< Dune::SGrid< 3, 3 >, Stuff::Grid::ChooseLayer::leaf, double, 1, 1,
                   GDT::ChooseSpaceBackend::fem, Stuff::LA::ChooseBackend::eigen_sparse >;

template class CG< Dune::SGrid< 1, 1 >, Stuff::Grid::ChooseLayer::level, double, 1, 1,
                   GDT::ChooseSpaceBackend::fem, Stuff::LA::ChooseBackend::eigen_sparse >;
template class CG< Dune::SGrid< 2, 2 >, Stuff::Grid::ChooseLayer::level, double, 1, 1,
                   GDT::ChooseSpaceBackend::fem, Stuff::LA::ChooseBackend::eigen_sparse >;
template class CG< Dune::SGrid< 3, 3 >, Stuff::Grid::ChooseLayer::level, double, 1, 1,
                   GDT::ChooseSpaceBackend::fem, Stuff::LA::ChooseBackend::eigen_sparse >;


# endif // HAVE_EIGEN
#endif // HAVE_DUNE_FEM

} // namespace Discretizations
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune
