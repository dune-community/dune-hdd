// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include "swipdg.hh"

#if HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_EIGEN

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Discretizations {


template class SWIPDG< ALUGrid< 2, 2, simplex, conforming >,
                       Stuff::Grid::ChooseLayer::leaf,
                       double, 1, 1,
                       GDT::ChooseSpaceBackend::fem,
                       Stuff::LA::ChooseBackend::eigen_sparse >;


} // namespace Discretizations
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_EIGEN

