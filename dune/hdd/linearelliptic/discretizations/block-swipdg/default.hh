// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_BLOCK_SWIPDG_HH
#define DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_BLOCK_SWIPDG_HH

#include <dune/grid/multiscale/provider.hh>

#include "base.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Discretizations {
namespace BlockSwipdg {


// forward, needed in the Traits
template <class G, class R = double, int r = 1, int p = 1, Stuff::LA::ChooseBackend la = Stuff::LA::default_sparse_backend>
class DefaultMultiscaleGrid;


namespace internal {


template <class G, class R, int r, int p, Stuff::LA::ChooseBackend la>
class DefaultMultiscaleGridTraits
{
public:
  typedef DefaultMultiscaleGrid<G, R, r, p, la> derived_type;
};


} // namespace internal


/**
 * \attention The given problem is replaced by a Problems::ZeroBoundary.
 * \attention The given boundary info config is replaced by a Stuff::Grid::BoundaryInfos::AllDirichlet.
 * \attention The boundary info for the local oversampled discretizations is hardwired to dirichlet zero atm!
 */
template <class G, class R, int r, int p, Stuff::LA::ChooseBackend la>
class DefaultMultiscaleGrid
  : public internal::Base<internal::DefaultMultiscaleGridTraits<G, R, r, p, la>, G, R, r, p, la>

{
  typedef internal::Base<internal::DefaultMultiscaleGridTraits<G, R, r, p, la>, G, R, r, p, la> BaseType;
public:
  using typename BaseType::GridProviderType;
  using typename BaseType::ProblemType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".default-multiscalegrid";
  }

  DefaultMultiscaleGrid(const GridProviderType& grid_provider,
                        const Stuff::Common::Configuration& bound_inf_cfg,
                        const ProblemType& prob,
                        const std::vector< std::string >& only_these_products = {})
    : BaseType(grid_provider, bound_inf_cfg, prob, only_these_products)
  {}
}; // DefaultMultiscaleGrid


#if HAVE_ALUGRID && HAVE_DUNE_FEM
# if HAVE_DUNE_ISTL

extern template class DefaultMultiscaleGrid< ALUGrid< 2, 2, simplex, conforming, No_Comm >,
                                             double,
                                             1,
                                             1,
                                             Stuff::LA::ChooseBackend::istl_sparse >;

#   if HAVE_MPI

extern template class DefaultMultiscaleGrid< ALUGrid< 2, 2, simplex, conforming, MPI_Comm >,
                                             double,
                                             1,
                                             1,
                                             Stuff::LA::ChooseBackend::istl_sparse >;

#   endif // HAVE_MPI
# endif // HAVE_DUNE_ISTL
# if HAVE_EIGEN

extern template class DefaultMultiscaleGrid< ALUGrid< 2, 2, simplex, conforming, No_Comm >,
                                             double,
                                             1,
                                             1,
                                             Stuff::LA::ChooseBackend::eigen_sparse >;

#   if HAVE_MPI

extern template class DefaultMultiscaleGrid< ALUGrid< 2, 2, simplex, conforming, MPI_Comm >,
                                             double,
                                             1,
                                             1,
                                             Stuff::LA::ChooseBackend::eigen_sparse >;

#   endif // HAVE_MPI
# endif // HAVE_EIGEN
#endif // HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_ISTL


template< Stuff::LA::ChooseBackend la, class G, class R = double, int r = 1, int p = 1 >
DefaultMultiscaleGrid< G, R, r, p, la > make_default(const grid::Multiscale::ProviderInterface< G >& grid_provider,
                                                     const DSC::Configuration& boundary_info,
                                                     const ProblemInterface< typename G::template Codim< 0 >::Entity,
                                                                             typename G::ctype, G::dimension,
                                                                             R, r >& problem,
                                                     const std::vector< std::string >& only_these_products = {})
{
  return DefaultMultiscaleGrid< G, R, r, p, la >(grid_provider, boundary_info, problem, only_these_products);
}


} // namespace BlockSwipdg
} // namespace Discretizations
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_BLOCK_SWIPDG_HH
