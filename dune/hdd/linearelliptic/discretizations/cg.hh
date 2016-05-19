// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_CG_HH
#define DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_CG_HH

#include <memory>
#include <vector>

#include <dune/common/timer.hh>

#if HAVE_ALUGRID
# include <dune/grid/alugrid.hh>
#endif

#if HAVE_DUNE_GRID_MULTISCALE
# include <dune/grid/multiscale/provider/interface.hh>
#endif

#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/timedlogging.hh>
#include <dune/stuff/grid/layers.hh>
#include <dune/stuff/grid/provider.hh>
#include <dune/stuff/la/container.hh>
#include <dune/stuff/la/solver.hh>

#include <dune/gdt/spaces/cg.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/playground/operators/elliptic-cg.hh>
#include <dune/gdt/functionals/l2.hh>

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
#elif HAVE_DUNE_PDELAB
          GDT::ChooseSpaceBackend space_backend = GDT::ChooseSpaceBackend::pdelab,
#else
# error No suitable space backend available!
#endif
          Stuff::LA::ChooseBackend la_backend = Stuff::LA::default_sparse_backend >
class CG;


namespace internal {


template< class GridImp, Stuff::Grid::ChooseLayer layer, class RangeFieldImp, int rangeDim, int polynomialOrder,
          GDT::ChooseSpaceBackend space_backend,
          Stuff::LA::ChooseBackend la_backend >
class CGTraits
  : public internal::ContainerBasedDefaultTraits< typename Stuff::LA::Container< RangeFieldImp, la_backend >::MatrixType,
                                                  typename Stuff::LA::Container< RangeFieldImp, la_backend >::VectorType >
{
  typedef internal::ContainerBasedDefaultTraits< typename Stuff::LA::Container< RangeFieldImp, la_backend >::MatrixType,
  typename Stuff::LA::Container< RangeFieldImp, la_backend >::VectorType > BaseType;
public:
  typedef CG< GridImp, layer, RangeFieldImp, rangeDim, polynomialOrder, space_backend, la_backend > derived_type;
  typedef GridImp GridType;
  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRange = rangeDim;
  static const unsigned int polOrder = polynomialOrder;

private:
  typedef GDT::Spaces::CGProvider< GridType, layer, space_backend, polOrder, RangeFieldType, dimRange > SpaceProvider;

  friend class CG< GridImp, layer, RangeFieldImp, rangeDim, polynomialOrder, space_backend, la_backend >;

public:
  typedef typename SpaceProvider::Type TestSpaceType;
  typedef TestSpaceType AnsatzSpaceType;
  typedef typename TestSpaceType::GridViewType GridViewType;
}; // class CGTraits


} // namespace internal


template< class GridImp, Stuff::Grid::ChooseLayer layer, class RangeFieldImp, int rangeDim, int polynomialOrder,
          GDT::ChooseSpaceBackend space_backend,
          Stuff::LA::ChooseBackend la_backend >
class CG
  : public ContainerBasedDefault< internal::CGTraits< GridImp, layer, RangeFieldImp, rangeDim, polynomialOrder, space_backend, la_backend > >
{
  typedef ContainerBasedDefault< internal::CGTraits< GridImp, layer, RangeFieldImp, rangeDim, polynomialOrder, space_backend, la_backend > > BaseType;
public:
  typedef internal::CGTraits< GridImp, layer, RangeFieldImp, rangeDim, polynomialOrder, space_backend, la_backend > Traits;
  typedef GridImp GridType;
  using typename BaseType::ProblemType;
  using typename BaseType::BoundaryInfoType;
  using typename BaseType::TestSpaceType;
  using typename BaseType::AnsatzSpaceType;
  using typename BaseType::MatrixType;
  using typename BaseType::VectorType;
  using typename BaseType::GridViewType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::PatternType;

protected:
  typedef typename Traits::SpaceProvider SpaceProvider;

  typedef Stuff::Grid::ProviderInterface< GridType >      GridProviderType;
#if HAVE_DUNE_GRID_MULTISCALE
  typedef grid::Multiscale::ProviderInterface< GridType > MsGridProviderType;
#endif
  using typename BaseType::AffinelyDecomposedMatrixType;
  using typename BaseType::AffinelyDecomposedVectorType;

  typedef typename ProblemType::DiffusionFactorType::NonparametricType DiffusionFactorType;
  typedef typename ProblemType::DiffusionTensorType::NonparametricType DiffusionTensorType;
  typedef GDT::Operators::EllipticCG< DiffusionFactorType, MatrixType, TestSpaceType,
                                      AnsatzSpaceType, GridViewType, DiffusionTensorType > EllipticOperatorType;

public:
  static std::string static_id();

  CG(GridProviderType& grid_provider,
     Stuff::Common::Configuration bound_inf_cfg,
     const ProblemType& prob,
     const int level = 0,
     const std::vector< std::string >& only_these_products = {});

#if HAVE_DUNE_GRID_MULTISCALE

  CG(MsGridProviderType& grid_provider,
     const Stuff::Common::Configuration& bound_inf_cfg,
     const ProblemType& prob,
     const int level_or_subdomain = 0,
     const std::vector< std::string >& only_these_products = {});

#endif // HAVE_DUNE_GRID_MULTISCALE

  void init(const bool prune = false);

private:
  friend class BlockSWIPDG< GridImp, RangeFieldImp, rangeDim, polynomialOrder, la_backend >;

  template< class SetConstraints, class ClearConstraints >
  std::shared_ptr< AffinelyDecomposedMatrixType >
  make_zero_dirichlet_product(const std::shared_ptr< AffinelyDecomposedMatrixType >& prod,
                              const SetConstraints& set,
                              const ClearConstraints& clear)
  {
    auto ret = std::make_shared<AffinelyDecomposedMatrixType>(prod->copy());
    for (ssize_t qq = 0; qq < ret->num_components(); ++qq) {
      // clear rows
      clear.apply(*ret->component(qq));
      // clear columns
      for (const auto& column : clear.dirichlet_DoFs())
        ret->component(qq)->unit_col(column);
    }
    if (!ret->has_affine_part()) {
      // we can assume that there is at least one component
      assert(ret->num_components() > 0);
      ret->register_affine_part(new MatrixType(ret->component(0)->copy()));
      *ret->affine_part() *= 0;
    }
    // set rows
    set.apply(*ret->affine_part());
    // set columns
    for (const auto& column : set.dirichlet_DoFs())
      ret->affine_part()->unit_col(column);
    return ret;
  } // ... make_zero_dirichlet_product(...)

  using BaseType::pattern_;
  const std::vector< std::string > only_these_products_;
}; // class CG


//#if HAVE_ALUGRID

//extern template class CG< ALUGrid< 2, 2, simplex, conforming >,
//                          Stuff::Grid::ChooseLayer::leaf,
//                          double,
//                          1,
//                          1,
//                          GDT::ChooseSpaceBackend::fem,
//                          Stuff::LA::ChooseBackend::istl_sparse >;

//#endif // HAVE_ALUGRID


} // namespace Discretizations
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_CG_HH
