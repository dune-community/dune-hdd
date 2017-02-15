// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_SWIPDG_HH
#define DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_SWIPDG_HH

#include <memory>
#include <vector>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/timer.hh>
#include <dune/common/static_assert.hh>

#include <dune/stuff/common/disable_warnings.hh>
# if HAVE_ALUGRID
#   include <dune/grid/alugrid.hh>
# endif
#include <dune/stuff/common/reenable_warnings.hh>

#if HAVE_DUNE_GRID_MULTISCALE
# include <dune/grid/multiscale/provider/interface.hh>
#endif

#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/timedlogging.hh>
#include <dune/stuff/grid/layers.hh>
#include <dune/stuff/grid/provider.hh>
#include <dune/stuff/la/container.hh>
#include <dune/stuff/la/solver.hh>
#include <dune/stuff/grid/walker/functors.hh>

#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/assembler/tmp-storage.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/functionals/l2.hh>
#include <dune/gdt/operators/oswaldinterpolation.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/playground/functionals/swipdg.hh>
#include <dune/gdt/playground/operators/elliptic-swipdg.hh>
#include <dune/gdt/playground/operators/fluxreconstruction.hh>
#include <dune/gdt/playground/products/swipdgpenalty.hh>
#include <dune/gdt/spaces/fv/default.hh>
#include <dune/gdt/spaces/rt/pdelab.hh>
#include <dune/gdt/spaces/dg.hh>

#include "base.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Discretizations {


// forward, for friendlyness
template< class GridImp, class RangeFieldImp, int rangeDim, int polynomialOrder, Stuff::LA::ChooseBackend la_backend >
class BlockSWIPDG;


// forward, needed in the Traits
template< class GridImp, Stuff::Grid::ChooseLayer layer, class RangeFieldImp, int rangeDim, int polynomialOrder = 1,
#if HAVE_DUNE_FEM
          GDT::ChooseSpaceBackend space_backend = GDT::ChooseSpaceBackend::fem,
#else
# error No suitable space backend available!
#endif
          Stuff::LA::ChooseBackend la_backend  = Stuff::LA::default_sparse_backend >
class SWIPDG;


namespace internal {


template< class GridImp, Stuff::Grid::ChooseLayer layer, class RangeFieldImp, int rangeDim, int polynomialOrder,
          GDT::ChooseSpaceBackend space_backend,
          Stuff::LA::ChooseBackend la_backend >
class SWIPDGTraits
  : public internal::ContainerBasedDefaultTraits< typename Stuff::LA::Container< RangeFieldImp, la_backend >::MatrixType,
                                                  typename Stuff::LA::Container< RangeFieldImp, la_backend >::VectorType>
{
public:
  typedef SWIPDG< GridImp, layer, RangeFieldImp, rangeDim, polynomialOrder, space_backend, la_backend > derived_type;
  typedef GridImp GridType;
  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRange = rangeDim;
  static const unsigned int polOrder = polynomialOrder;

private:
  typedef GDT::Spaces::DGProvider< GridType, layer, space_backend, polOrder, RangeFieldType, dimRange > SpaceProvider;

  friend class SWIPDG< GridImp, layer, RangeFieldImp, rangeDim, polynomialOrder, space_backend, la_backend >;

public:
  typedef typename SpaceProvider::Type TestSpaceType;
  typedef TestSpaceType AnsatzSpaceType;
  typedef typename TestSpaceType::GridViewType GridViewType;
}; // class SWIPDGTraits


} // namespace internal


template< class GridImp, Stuff::Grid::ChooseLayer layer, class RangeFieldImp, int rangeDim, int polynomialOrder,
          GDT::ChooseSpaceBackend space_backend,
          Stuff::LA::ChooseBackend la_backend >
class SWIPDG
  : public ContainerBasedDefault< internal::SWIPDGTraits< GridImp, layer, RangeFieldImp, rangeDim, polynomialOrder,
                                                          space_backend, la_backend > >
{
  typedef ContainerBasedDefault< internal::SWIPDGTraits< GridImp, layer, RangeFieldImp, rangeDim, polynomialOrder,
                                                         space_backend, la_backend > > BaseType;
  typedef SWIPDG< GridImp, layer, RangeFieldImp, rangeDim, polynomialOrder, space_backend, la_backend > ThisType;
public:
  typedef internal::SWIPDGTraits< GridImp, layer, RangeFieldImp, rangeDim, polynomialOrder, space_backend, la_backend >
      Traits;
  typedef GridImp GridType;
  using typename BaseType::ProblemType;
  using typename BaseType::BoundaryInfoType;
  using typename BaseType::TestSpaceType;
  using typename BaseType::AnsatzSpaceType;
  using typename BaseType::MatrixType;
  using typename BaseType::VectorType;
  using typename BaseType::GridViewType;
  using typename BaseType::RangeFieldType;

  typedef typename TestSpaceType::PatternType PatternType;

  static const unsigned int dimDomain = BaseType::dimDomain;
  static const unsigned int dimRange = BaseType::dimRange;

private:
  typedef typename Traits::SpaceProvider SpaceProvider;

  typedef Stuff::Grid::ProviderInterface< GridType >      GridProviderType;
#if HAVE_DUNE_GRID_MULTISCALE
  typedef grid::Multiscale::ProviderInterface< GridType > MsGridProviderType;
#endif
  using typename BaseType::AffinelyDecomposedMatrixType;
  using typename BaseType::AffinelyDecomposedVectorType;

  typedef typename ProblemType::DiffusionFactorType::NonparametricType DiffusionFactorType;
  typedef typename ProblemType::DiffusionTensorType::NonparametricType DiffusionTensorType;
  typedef GDT::Operators::EllipticSWIPDG< DiffusionFactorType, MatrixType
                                        , TestSpaceType, AnsatzSpaceType
                                        , GridViewType, DiffusionTensorType > EllipticOperatorType;

public:
  static std::string static_id();

  template< Stuff::Grid::ChooseLayer lr = layer >
  SWIPDG(GridProviderType& grid_provider,
         const Stuff::Common::Configuration& bound_inf_cfg,
         const ProblemType& prob,
         const int level_or_subdomain = 0,
         const std::vector< std::string >& only_these_products = {}
#if HAVE_DUNE_GRID_MULTISCALE
       , typename std::enable_if<    (lr != Stuff::Grid::ChooseLayer::local)
                                  && (lr != Stuff::Grid::ChooseLayer::local_oversampled), void >::type* /*disable_for_local_grid_parts*/ = nullptr
#endif
        );

#if HAVE_DUNE_GRID_MULTISCALE

  SWIPDG(const MsGridProviderType& grid_provider,
         const Stuff::Common::Configuration& bound_inf_cfg,
         const ProblemType& prob,
         const int level_or_subdomain = 0,
         const std::vector< std::string >& only_these_products = {});

#endif // HAVE_DUNE_GRID_MULTISCALE

  void init(const bool prune = false);

private:
  friend class BlockSWIPDG< GridImp, RangeFieldImp, rangeDim, polynomialOrder, la_backend >;

  const RangeFieldType beta_;
  using BaseType::pattern_;
  const std::vector< std::string > only_these_products_;
}; // class SWIPDG

} // namespace Discretizations
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#include "swipdg.hxx"
#endif // DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_SWIPDG_HH
