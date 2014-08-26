// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include "block-swipdg.hh"

#if HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_EIGEN && HAVE_DUNE_GRID_MULTISCALE

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Discretizations {


template class BlockSWIPDGEstimator<
    GDT::Spaces::Block< GDT::Spaces::DiscontinuousLagrange::FemBased<
        grid::Part::Local::IndexBased::Const< grid::Part::Leaf::Const< ALUGrid< 2, 2, simplex, conforming > > >,
        1, double, 1, 1 > >,
    Dune::Stuff::LA::EigenDenseVector< double >,
    ProblemInterface< typename ALUGrid< 2, 2, simplex, conforming >::template Codim< 0 >::Entity, double, 2, double, 1 >,
    ALUGrid< 2, 2, simplex, conforming > >;


} // namespace Discretizations
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_EIGEN && HAVE_DUNE_GRID_MULTISCALE
