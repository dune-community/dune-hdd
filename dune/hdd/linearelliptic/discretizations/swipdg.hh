// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_SWIPDG_HH
#define DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_SWIPDG_HH

#include <memory>
#include <vector>

#include <dune/common/timer.hh>
#include <dune/common/static_assert.hh>

#if HAVE_ALUGRID
# include <dune/grid/alugrid.hh>
#endif
#include <dune/stuff/common/disable_warnings.hh>
# include <dune/grid/sgrid.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#if HAVE_DUNE_GRID_MULTISCALE
# include <dune/grid/multiscale/provider/interface.hh>
#endif

#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/grid/layers.hh>
#include <dune/stuff/grid/provider.hh>
#include <dune/stuff/la/container.hh>
#include <dune/stuff/la/solver.hh>
#include <dune/stuff/playground/functions/ESV2007.hh>

#include <dune/gdt/spaces/discontinuouslagrange.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/playground/assembler/functors.hh>
#include <dune/gdt/playground/operators/elliptic-swipdg.hh>
#include <dune/gdt/functionals/l2.hh>
#include <dune/gdt/playground/functionals/swipdg.hh>
#include <dune/gdt/products/l2.hh>
#include <dune/gdt/products/h1.hh>
#include <dune/gdt/playground/products/elliptic.hh>
#include <dune/gdt/operators/oswaldinterpolation.hh>
#include <dune/gdt/playground/spaces/finitevolume/default.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/playground/spaces/raviartthomas/pdelab.hh>
#include <dune/gdt/playground/operators/fluxreconstruction.hh>
#include <dune/gdt/playground/products/ESV2007.hh>
#include <dune/gdt/assembler/tmp-storage.hh>

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
          GDT::ChooseSpaceBackend space_backend = GDT::ChooseSpaceBackend::fem,
          Stuff::LA::ChooseBackend la_backend = Stuff::LA::ChooseBackend::istl_sparse >
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
  typedef GDT::Spaces::DiscontinuousLagrangeProvider< GridType, layer, space_backend,
                                                      polOrder, RangeFieldType, dimRange > SpaceProvider;

  friend class SWIPDG< GridImp, layer, RangeFieldImp, rangeDim, polynomialOrder, space_backend, la_backend >;

public:
  typedef typename SpaceProvider::Type TestSpaceType;
  typedef TestSpaceType AnsatzSpaceType;
  typedef typename TestSpaceType::GridViewType GridViewType;
}; // class SWIPDGTraits


namespace Estimators {


static const size_t over_integrate = 2;


class LocalNonconformityESV2007Base
{
public:
  static std::string id() { return "eta_NC_ESV2007"; }
};


template< class DiscType, class GridType, int dimRange >
class LocalNonconformityESV2007
  : public LocalNonconformityESV2007Base
{
public:
  static const bool available = false;
};

//#if HAVE_ALUGRID

/**
 *  \brief computes the local nonconformity estimator as defined in ESV2007
 */
template< class DiscType >
class LocalNonconformityESV2007< DiscType, ALUGrid< 2, 2, simplex, conforming >, 1 >
  : public LocalNonconformityESV2007Base
  , public GDT::Functor::Codim0< typename DiscType::GridViewType >
{
  typedef LocalNonconformityESV2007< DiscType, ALUGrid< 2, 2, simplex, conforming >, 1 > ThisType;
  typedef GDT::Functor::Codim0< typename DiscType::GridViewType > FunctorBaseType;
public:
  static const bool available = true;

  typedef typename FunctorBaseType::GridViewType GridViewType;
  typedef typename FunctorBaseType::EntityType   EntityType;

  typedef typename DiscType::RangeFieldType RangeFieldType;
  typedef typename DiscType::VectorType     VectorType;

private:
  typedef typename DiscType::AnsatzSpaceType AnsatzSpaceType;
  typedef GDT::ConstDiscreteFunction< AnsatzSpaceType, VectorType > ConstDiscreteFunctionType;
  typedef GDT::DiscreteFunction< AnsatzSpaceType, VectorType > DiscreteFunctionType;
  typedef typename ConstDiscreteFunctionType::DifferenceType DifferenceType;

  typedef typename DiscType::ProblemType ProblemType;
  typedef typename ProblemType::DiffusionFactorType::NonparametricType DiffusionFactorType;
  typedef typename ProblemType::DiffusionTensorType::NonparametricType DiffusionTensorType;

  typedef GDT::LocalOperator::Codim0Integral< GDT::LocalEvaluation::Elliptic< DiffusionFactorType,
                                                                              DiffusionTensorType > > LocalOperatorType;
  typedef GDT::TmpStorageProvider::Matrices< RangeFieldType > TmpStorageProviderType;

  static const DiscType& assert_disc(const DiscType& disc)
  {
    if (disc.parametric())
      DUNE_THROW(NotImplemented, "Not implemented yet for parametric discretizations!");
    assert(disc.problem().diffusion_factor()->has_affine_part());
    assert(disc.problem().diffusion_tensor()->has_affine_part());
    return disc;
  } // ... assert_disc(...)

public:
  static RangeFieldType estimate(const DiscType& disc, const VectorType& vector)
  {
    ThisType estimator(disc, vector);
    GDT::GridWalker< GridViewType > grid_walker(*disc.grid_view());
    grid_walker.add(estimator);
    grid_walker.walk();
    return std::sqrt(estimator.result_);
  } // ... estimate(...)

  LocalNonconformityESV2007(const DiscType& disc, const VectorType& vector)
    : disc_(assert_disc(disc))
    , vector_(vector)
    , discrete_solution_(*disc_.ansatz_space(), vector_)
    , oswald_interpolation_(*disc.ansatz_space())
    , difference_(Stuff::Common::make_unique< DifferenceType >(discrete_solution_ - oswald_interpolation_))
    , local_operator_(Estimators::over_integrate,
                      *disc.problem().diffusion_factor()->affine_part(),
                      *disc.problem().diffusion_tensor()->affine_part())
    , tmp_local_matrices_({1, local_operator_.numTmpObjectsRequired()}, 1, 1)
    , result_(0.0)
  {}

  virtual void prepare()
  {
    const GDT::Operators::OswaldInterpolation< GridViewType > oswald_interpolation_operator(*disc_.grid_view());
    oswald_interpolation_operator.apply(discrete_solution_, oswald_interpolation_);
    result_ = 0.0;
  } // ... prepare(...)

  RangeFieldType compute_locally(const EntityType& entity)
  {
    const auto local_difference = difference_->local_function(entity);
    local_operator_.apply(*local_difference,
                          *local_difference,
                          tmp_local_matrices_.matrices()[0][0],
                          tmp_local_matrices_.matrices()[1]);
    assert(tmp_local_matrices_.matrices()[0][0].rows() >= 1);
    assert(tmp_local_matrices_.matrices()[0][0].cols() >= 1);
    return tmp_local_matrices_.matrices()[0][0][0][0];
  } // ... compute_locally(...)

  virtual void apply_local(const EntityType &entity)
  {
    result_ += compute_locally(entity);
  }

private:
  const DiscType& disc_;
  const VectorType& vector_;
  const ConstDiscreteFunctionType discrete_solution_;
  DiscreteFunctionType oswald_interpolation_;
  std::unique_ptr< const DifferenceType > difference_;
  const LocalOperatorType local_operator_;
  TmpStorageProviderType tmp_local_matrices_;
public:
  RangeFieldType result_;
}; // class LocalNonconformityESV2007< ..., ALUGrid< 2, 2, simplex, conforming >, 1 >

//#endif // HAVE_ALUGRID


class LocalResidualESV2007Base
{
public:
  static std::string id() { return "eta_R_ESV2007"; }
};


template< class DiscType, class GridType, int dimRange >
class LocalResidualESV2007
  : public LocalResidualESV2007Base
{
public:
  static const bool available = false;
};

//#if HAVE_ALUGRID

/**
 *  \brief computes the local residual estimator as defined in ESV2007
 */
template< class DiscType >
class LocalResidualESV2007< DiscType, ALUGrid< 2, 2, simplex, conforming >, 1 >
  : public LocalResidualESV2007Base
  , public GDT::Functor::Codim0< typename DiscType::GridViewType >
{
  typedef LocalResidualESV2007< DiscType, ALUGrid< 2, 2, simplex, conforming >, 1 > ThisType;
  typedef GDT::Functor::Codim0< typename DiscType::GridViewType > FunctorBaseType;
public:
  static const bool available = true;

  typedef typename FunctorBaseType::GridViewType GridViewType;
  typedef typename FunctorBaseType::EntityType   EntityType;

  typedef typename DiscType::RangeFieldType RangeFieldType;
  typedef typename DiscType::VectorType     VectorType;

private:
  typedef GDT::Spaces::FiniteVolume::Default< GridViewType, RangeFieldType, 1, 1 > P0SpaceType;
  typedef GDT::DiscreteFunction< P0SpaceType, VectorType > DiscreteFunctionType;
  typedef typename DiscreteFunctionType::DifferenceType DifferenceType;

  typedef typename DiscType::ProblemType ProblemType;
  typedef typename ProblemType::DiffusionFactorType::NonparametricType DiffusionFactorType;
  typedef typename ProblemType::DiffusionTensorType::NonparametricType DiffusionTensorType;

  typedef typename Stuff::Functions::ESV2007::Cutoff< DiffusionFactorType, DiffusionTensorType > CutoffFunctionType;
  typedef GDT::LocalOperator::Codim0Integral< GDT::LocalEvaluation::Product< CutoffFunctionType > > LocalOperatorType;
  typedef GDT::TmpStorageProvider::Matrices< RangeFieldType > TmpStorageProviderType;

  static const DiscType& assert_disc(const DiscType& disc)
  {
    if (disc.parametric())
      DUNE_THROW(NotImplemented, "Not implemented yet for parametric discretizations!");
    assert(disc.problem().diffusion_factor()->has_affine_part());
    assert(disc.problem().diffusion_tensor()->has_affine_part());
    assert(disc.problem().force()->has_affine_part());
    return disc;
  } // ... assert_disc(...)

public:
  static RangeFieldType estimate(const DiscType& disc, const VectorType& /*vector*/)
  {
    ThisType estimator(disc);
    GDT::GridWalker< GridViewType > grid_walker(*disc.grid_view());
    grid_walker.add(estimator);
    grid_walker.walk();
    return std::sqrt(estimator.result_);
  } // ... estimate(...)

  LocalResidualESV2007(const DiscType& disc)
    : disc_(assert_disc(disc))
    , p0_space_(disc_.grid_view())
    , p0_force_(p0_space_)
    , difference_(Stuff::Common::make_unique< DifferenceType >(*disc_.problem().force()->affine_part() - p0_force_))
    , cutoff_function_(*disc.problem().diffusion_factor()->affine_part(),
                       *disc.problem().diffusion_tensor()->affine_part())
    , local_operator_(Estimators::over_integrate, cutoff_function_)
    , tmp_local_matrices_({1, local_operator_.numTmpObjectsRequired()}, 1, 1)
    , result_(0.0)
  {}

  virtual void prepare()
  {
    const GDT::Operators::Projection< GridViewType > projection_operator(*disc_.grid_view());
    projection_operator.apply(*disc_.problem().force()->affine_part(), p0_force_);
    result_ = 0.0;
  } // ... prepare(...)

  RangeFieldType compute_locally(const EntityType& entity)
  {
    const auto local_difference = difference_->local_function(entity);
    local_operator_.apply(*local_difference,
                          *local_difference,
                          tmp_local_matrices_.matrices()[0][0],
                          tmp_local_matrices_.matrices()[1]);
    assert(tmp_local_matrices_.matrices()[0][0].rows() >= 1);
    assert(tmp_local_matrices_.matrices()[0][0].cols() >= 1);
    return tmp_local_matrices_.matrices()[0][0][0][0];
  } // ... compute_locally(...)

  virtual void apply_local(const EntityType &entity)
  {
    result_ += compute_locally(entity);
  }

private:
  const DiscType& disc_;
  const P0SpaceType p0_space_;
  DiscreteFunctionType p0_force_;
  std::unique_ptr< const DifferenceType > difference_;
  const CutoffFunctionType cutoff_function_;
  const LocalOperatorType local_operator_;
  TmpStorageProviderType tmp_local_matrices_;
public:
  RangeFieldType result_;
}; // class LocalResidualESV2007< ..., ALUGrid< 2, 2, simplex, conforming >, 1 >

//#endif // HAVE_ALUGRID


class LocalDiffusiveFluxESV2007Base
{
public:
  static std::string id() { return "eta_DF_ESV2007"; }
};


template< class DiscType, class GridType, int dimRange >
class LocalDiffusiveFluxESV2007
  : public LocalDiffusiveFluxESV2007Base
{
public:
  static const bool available = false;
};

//#if HAVE_ALUGRID

/**
 *  \brief computes the local diffusive flux estimator as defined in ESV2007
 */
template< class DiscType >
class LocalDiffusiveFluxESV2007< DiscType, ALUGrid< 2, 2, simplex, conforming >, 1 >
  : public LocalDiffusiveFluxESV2007Base
  , public GDT::Functor::Codim0< typename DiscType::GridViewType >
{
  typedef LocalDiffusiveFluxESV2007< DiscType, ALUGrid< 2, 2, simplex, conforming >, 1 > ThisType;
  typedef GDT::Functor::Codim0< typename DiscType::GridViewType > FunctorBaseType;
public:
  static const bool available = true;

  typedef typename FunctorBaseType::GridViewType GridViewType;
  typedef typename FunctorBaseType::EntityType   EntityType;

  typedef typename DiscType::RangeFieldType RangeFieldType;
  typedef typename DiscType::VectorType     VectorType;

  static const unsigned int dimDomain = DiscType::dimDomain;

private:
  typedef typename DiscType::AnsatzSpaceType AnsatzSpaceType;
  typedef GDT::ConstDiscreteFunction< AnsatzSpaceType, VectorType > ConstDiscreteFunctionType;
  typedef GDT::Spaces::RaviartThomas::PdelabBased< GridViewType, 0, RangeFieldType, dimDomain > RTN0SpaceType;
  typedef GDT::DiscreteFunction< RTN0SpaceType, VectorType > RTN0DiscreteFunctionType;

  typedef typename DiscType::ProblemType ProblemType;
  typedef typename ProblemType::DiffusionFactorType::NonparametricType DiffusionFactorType;
  typedef typename ProblemType::DiffusionTensorType::NonparametricType DiffusionTensorType;

  typedef GDT::LocalOperator::Codim0Integral<
      GDT::LocalEvaluation::ESV2007::DiffusiveFluxEstimate< DiffusionFactorType,
                                                            RTN0DiscreteFunctionType,
                                                            DiffusionTensorType > > LocalOperatorType;
  typedef GDT::TmpStorageProvider::Matrices< RangeFieldType > TmpStorageProviderType;

  static const DiscType& assert_disc(const DiscType& disc)
  {
    if (disc.parametric())
      DUNE_THROW(NotImplemented, "Not implemented yet for parametric discretizations!");
    assert(disc.problem().diffusion_factor()->has_affine_part());
    assert(disc.problem().diffusion_tensor()->has_affine_part());
    return disc;
  } // ... assert_disc(...)

public:
  static RangeFieldType estimate(const DiscType& disc, const VectorType& vector)
  {
    ThisType estimator(disc, vector);
    GDT::GridWalker< GridViewType > grid_walker(*disc.grid_view());
    grid_walker.add(estimator);
    grid_walker.walk();
    return std::sqrt(estimator.result_);
  } // ... estimate(...)

  LocalDiffusiveFluxESV2007(const DiscType& disc, const VectorType& vector)
    : disc_(assert_disc(disc))
    , vector_(vector)
    , discrete_solution_(*disc_.ansatz_space(), vector_)
    , rtn0_space_(disc.grid_view())
    , diffusive_flux_(rtn0_space_)
    , local_operator_(Estimators::over_integrate,
                      *disc.problem().diffusion_factor()->affine_part(),
                      *disc.problem().diffusion_tensor()->affine_part(),
                      diffusive_flux_)
    , tmp_local_matrices_({1, local_operator_.numTmpObjectsRequired()}, 1, 1)
    , result_(0.0)
  {}

  virtual void prepare()
  {
    const GDT::Operators::DiffusiveFluxReconstruction< GridViewType, DiffusionFactorType, DiffusionTensorType >
      diffusive_flux_reconstruction(*disc_.grid_view(),
                                    *disc_.problem().diffusion_factor()->affine_part(),
                                    *disc_.problem().diffusion_tensor()->affine_part());
    diffusive_flux_reconstruction.apply(discrete_solution_, diffusive_flux_);
    result_ = 0.0;
  } // ... prepare(...)

  RangeFieldType compute_locally(const EntityType& entity)
  {
    const auto local_discrete_solution = discrete_solution_.local_function(entity);
    local_operator_.apply(*local_discrete_solution,
                          *local_discrete_solution,
                          tmp_local_matrices_.matrices()[0][0],
                          tmp_local_matrices_.matrices()[1]);
    assert(tmp_local_matrices_.matrices()[0][0].rows() >= 1);
    assert(tmp_local_matrices_.matrices()[0][0].cols() >= 1);
    return tmp_local_matrices_.matrices()[0][0][0][0];
  } // ... compute_locally(...)

  virtual void apply_local(const EntityType &entity)
  {
    result_ += compute_locally(entity);
  }

private:
  const DiscType& disc_;
  const VectorType& vector_;
  const ConstDiscreteFunctionType discrete_solution_;
  const RTN0SpaceType rtn0_space_;
  RTN0DiscreteFunctionType diffusive_flux_;
  const LocalOperatorType local_operator_;
  TmpStorageProviderType tmp_local_matrices_;
public:
  RangeFieldType result_;
}; // class LocalDiffusiveFluxESV2007< ..., ALUGrid< 2, 2, simplex, conforming >, 1 >

//#endif // HAVE_ALUGRID


class ESV2007Base
{
public:
  static std::string id()
  {
    return "eta_ESV2007";
  }
};


template< class DiscType, class GridType, int dimRange >
class ESV2007
  : public ESV2007Base
{
public:
  static const bool available = false;
};


//#if HAVE_ALUGRID

template< class DiscType >
class ESV2007< DiscType, ALUGrid< 2, 2, simplex, conforming >, 1 >
  : public ESV2007Base
{
public:
  static const bool available = true;

  typedef typename DiscType::RangeFieldType RangeFieldType;
  typedef typename DiscType::VectorType     VectorType;

  static RangeFieldType estimate(const DiscType& disc, const VectorType& vector)
  {
    LocalNonconformityESV2007< DiscType, ALUGrid< 2, 2, simplex, conforming >, 1 > eta_nc(disc, vector);
    LocalResidualESV2007< DiscType, ALUGrid< 2, 2, simplex, conforming >, 1 >      eta_r(disc);
    LocalDiffusiveFluxESV2007< DiscType, ALUGrid< 2, 2, simplex, conforming >, 1 > eta_df(disc, vector);
    eta_nc.prepare();
    eta_r.prepare();
    eta_df.prepare();

    RangeFieldType eta_squared(0.0);

    const auto grid_view = *disc.grid_view();
    const auto entity_it_end = grid_view.template end< 0 >();
    for (auto entity_it = grid_view.template begin< 0 >(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      eta_squared += eta_nc.compute_locally(entity)
                   + std::pow(std::sqrt(eta_r.compute_locally(entity)) + std::sqrt(eta_df.compute_locally(entity)), 2);
    }
    return std::sqrt(eta_squared);
  } // ... estimate(...)
}; // class ESV2007< ..., ALUGrid< 2, 2, simplex, conforming >, 1 >

//#endif // HAVE_ALUGRID


class ESV2007AlternativeSummationBase
{
public:
  static std::string id()
  {
    return "eta_ESV2007_alt";
  }
};


template< class DiscType, class GridType, int dimRange >
class ESV2007AlternativeSummation
  : public ESV2007AlternativeSummationBase
{
public:
  static const bool available = false;
};


//#if HAVE_ALUGRID

template< class DiscType >
class ESV2007AlternativeSummation< DiscType, ALUGrid< 2, 2, simplex, conforming >, 1 >
  : public ESV2007AlternativeSummationBase
{
public:
  static const bool available = true;

  typedef typename DiscType::RangeFieldType RangeFieldType;
  typedef typename DiscType::VectorType     VectorType;

  static RangeFieldType estimate(const DiscType& disc, const VectorType& vector)
  {
    LocalNonconformityESV2007< DiscType, ALUGrid< 2, 2, simplex, conforming >, 1 > eta_nc(disc, vector);
    LocalResidualESV2007< DiscType, ALUGrid< 2, 2, simplex, conforming >, 1 >      eta_r(disc);
    LocalDiffusiveFluxESV2007< DiscType, ALUGrid< 2, 2, simplex, conforming >, 1 > eta_df(disc, vector);
    eta_nc.prepare();
    eta_r.prepare();
    eta_df.prepare();

    RangeFieldType eta_nc_squared(0.0);
    RangeFieldType eta_r_squared(0.0);
    RangeFieldType eta_df_squared(0.0);

    const auto grid_view = *disc.grid_view();
    const auto entity_it_end = grid_view.template end< 0 >();
    for (auto entity_it = grid_view.template begin< 0 >(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      eta_nc_squared += eta_nc.compute_locally(entity);
      eta_r_squared += eta_r.compute_locally(entity);
      eta_df_squared += eta_df.compute_locally(entity);
    }
    return std::sqrt(eta_nc_squared) + std::sqrt(eta_r_squared) + std::sqrt(eta_df_squared);
  } // ... estimate(...)
}; // class ESV2007AlternativeSummation< ..., ALUGrid< 2, 2, simplex, conforming >, 1 >

//#endif // HAVE_ALUGRID


} // namespace Estimators


template< class D, class G, int r >
class SWIPDGEstimator
{
public:
  typedef typename D::RangeFieldType RangeFieldType;
  typedef typename D::VectorType     VectorType;

private:
  template< class IndividualEstimator, bool available = false >
  class Caller
  {
  public:
    static std::vector< std::string > append(std::vector< std::string > in)
    {
      return in;
    }

    static bool equals(const std::string& /*type*/)
    {
      return false;
    }

    static RangeFieldType estimate(const D& /*disc*/, const VectorType& /*vector*/)
    {
      DUNE_THROW(Stuff::Exceptions::internal_error, "This should not happen!");
      return RangeFieldType(0);
    }
  }; // class Caller

  template< class IndividualEstimator >
  class Caller< IndividualEstimator, true >
  {
  public:
    static std::vector< std::string > append(std::vector< std::string > in)
    {
      in.push_back(IndividualEstimator::id());
      return in;
    }

    static bool equals(const std::string& type)
    {
      return IndividualEstimator::id() == type;
    }

    static RangeFieldType estimate(const D& disc, const VectorType& vector)
    {
      return IndividualEstimator::estimate(disc, vector);
    }
  }; // class Caller< ..., true >

  template< class IndividualEstimator >
  static std::vector< std::string > append(std::vector< std::string > in)
  {
    return Caller< IndividualEstimator, IndividualEstimator::available >::append(in);
  }

  template< class IndividualEstimator >
  static bool equals(const std::string& type)
  {
    return Caller< IndividualEstimator, IndividualEstimator::available >::equals(type);
  }

  template< class IndividualEstimator >
  static RangeFieldType compute_estimator(const D& disc, const VectorType& vector)
  {
    return Caller< IndividualEstimator, IndividualEstimator::available >::estimate(disc, vector);
  }

  typedef Estimators::LocalNonconformityESV2007< D, G, r >   LocalNonconformityESV2007Type;
  typedef Estimators::LocalResidualESV2007< D, G, r >        LocalResidualESV2007Type;
  typedef Estimators::LocalDiffusiveFluxESV2007< D, G, r >   LocalDiffusiveFluxESV2007Type;
  typedef Estimators::ESV2007< D, G, r >                     ESV2007Type;
  typedef Estimators::ESV2007AlternativeSummation< D, G, r > ESV2007AlternativeSummationType;

public:
  static std::vector< std::string > available()
  {
    std::vector< std::string > tmp;
    tmp = append< LocalNonconformityESV2007Type >(tmp);
    tmp = append< LocalResidualESV2007Type >(tmp);
    tmp = append< LocalDiffusiveFluxESV2007Type >(tmp);
    tmp = append< ESV2007Type >(tmp);
    tmp = append< ESV2007AlternativeSummationType >(tmp);
    return tmp;
  } // ... available(...)

  static RangeFieldType estimate(const D& disc, const VectorType& vector, const std::string type)
  {
    if (equals< LocalNonconformityESV2007Type >(type))
      return compute_estimator< LocalNonconformityESV2007Type >(disc, vector);
    else if (equals< LocalResidualESV2007Type >(type))
      return compute_estimator< LocalResidualESV2007Type >(disc, vector);
    else if (equals< LocalDiffusiveFluxESV2007Type >(type))
      return compute_estimator< LocalDiffusiveFluxESV2007Type >(disc, vector);
    else if (equals< ESV2007Type >(type))
      return compute_estimator< ESV2007Type >(disc, vector);
    else if (equals< ESV2007AlternativeSummationType >(type))
      return compute_estimator< ESV2007AlternativeSummationType >(disc, vector);
    else
      DUNE_THROW(Stuff::Exceptions::you_are_using_this_wrong,
                 "Requested type '" << type << "' is not one of available()!");
  } // ... estimate(...)
}; // class SWIPDGEstimator


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

  typedef Stuff::Grid::ConstProviderInterface< GridType > GridProviderType;
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
  static std::string static_id()
  {
    return DiscretizationInterface< Traits >::static_id() + ".swipdg";
  }

  static std::vector< std::string > available_estimators()
  {
    return internal::SWIPDGEstimator< ThisType, GridType, dimRange >::available();
  }

  SWIPDG(const GridProviderType& grid_provider,
         const Stuff::Common::Configuration& bound_inf_cfg,
         const ProblemType& prob,
         const int level_or_subdomain = 0)
    : BaseType(std::make_shared< TestSpaceType >(SpaceProvider::create(grid_provider, level_or_subdomain)),
               std::make_shared< AnsatzSpaceType >(SpaceProvider::create(grid_provider, level_or_subdomain)),
               bound_inf_cfg,
               prob)
    , beta_(GDT::LocalEvaluation::SWIPDG::internal::default_beta(dimDomain))
    , pattern_(EllipticOperatorType::pattern(*(BaseType::test_space()), *(BaseType::ansatz_space())))
  {
    // in case of parametric diffusion tensor this discretization is not affinely decomposable any more
    if (this->problem_.diffusion_tensor()->parametric())
      DUNE_THROW(NotImplemented, "The diffusion tensor must not be parametric!");
    if (!this->problem_.diffusion_tensor()->has_affine_part())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given, "The diffusion tensor must not be empty!");
  } // SWIPDG(...)

#if HAVE_DUNE_GRID_MULTISCALE
  SWIPDG(const MsGridProviderType& grid_provider,
         const Stuff::Common::Configuration& bound_inf_cfg,
         const ProblemType& prob,
         const int level_or_subdomain = 0)
    : BaseType(std::make_shared< TestSpaceType >(SpaceProvider::create(grid_provider, level_or_subdomain)),
               std::make_shared< AnsatzSpaceType >(SpaceProvider::create(grid_provider, level_or_subdomain)),
               bound_inf_cfg,
               prob)
    , beta_(GDT::LocalEvaluation::SWIPDG::internal::default_beta(dimDomain))
    , pattern_(EllipticOperatorType::pattern(*(BaseType::test_space()), *(BaseType::ansatz_space())))
  {
    // in case of parametric diffusion tensor this discretization is not affinely decomposable any more
    if (this->problem_.diffusion_tensor()->parametric())
      DUNE_THROW(NotImplemented, "The diffusion tensor must not be parametric!");
    if (!this->problem_.diffusion_tensor()->has_affine_part())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given, "The diffusion tensor must not be empty!");
  } // SWIPDG(...)
#endif // HAVE_DUNE_GRID_MULTISCALE

  const PatternType& pattern() const
  {
    return pattern_;
  }

  void init(std::ostream& out = Stuff::Common::Logger().devnull(), const std::string prefix = "")
  {
    if (!this->container_based_initialized_) {
      using namespace GDT;

      auto& matrix = *(this->matrix_);
      auto& rhs = *(this->rhs_);
      const auto& space = *(this->test_space_);
      const auto& boundary_info = *(this->boundary_info_);

      out << prefix << "assembling... " << std::flush;
      Dune::Timer timer;
      SystemAssembler< TestSpaceType > system_assembler(space);
      Functor::DirichletDetector< GridViewType > dirichlet_detector(this->boundary_info());
      system_assembler.add(dirichlet_detector);

      // lhs operator
      const auto& diffusion_factor = *(this->problem_.diffusion_factor());
      const auto& diffusion_tensor = *(this->problem_.diffusion_tensor());
      assert(!diffusion_tensor.parametric());
      assert(diffusion_tensor.has_affine_part());
      std::vector< std::unique_ptr< EllipticOperatorType > > elliptic_operators;
      for (size_t qq = 0; qq < diffusion_factor.num_components(); ++qq) {
        const size_t id = matrix.register_component(diffusion_factor.coefficient(qq),
                                                    space.mapper().size(), space.mapper().size(), pattern_);
        elliptic_operators.emplace_back(new EllipticOperatorType(
            *(diffusion_factor.component(qq)),
            *(diffusion_tensor.affine_part()),
            boundary_info,
            *(matrix.component(id)),
            space));
      }
      if (diffusion_factor.has_affine_part()) {
        if (!matrix.has_affine_part())
          matrix.register_affine_part(space.mapper().size(), space.mapper().size(), pattern_);
        elliptic_operators.emplace_back(new EllipticOperatorType(
            *(diffusion_factor.affine_part()),
            *(diffusion_tensor.affine_part()),
            boundary_info,
            *(matrix.affine_part()),
            space));
      }
      for (auto& elliptic_operator : elliptic_operators)
        system_assembler.add(*elliptic_operator);

      // rhs functional
      // * volume
      typedef typename ProblemType::FunctionType::NonparametricType FunctionType;
      const auto& force = *(this->problem_.force());
      typedef Functionals::L2Volume< FunctionType, VectorType, TestSpaceType > L2VolumeFunctionalType;
      std::vector< std::unique_ptr< L2VolumeFunctionalType > > force_functionals;
      for (size_t qq = 0; qq < force.num_components(); ++qq) {
        const size_t id = rhs.register_component(force.coefficient(qq), space.mapper().size());
        force_functionals.emplace_back(new L2VolumeFunctionalType(*(force.component(qq)),
                                                                  *(rhs.component(id)),
                                                                  space));
      }
      if (force.has_affine_part()) {
        if (!rhs.has_affine_part())
          rhs.register_affine_part(space.mapper().size());
        force_functionals.emplace_back(new L2VolumeFunctionalType(*(force.affine_part()),
                                                                  *(rhs.affine_part()),
                                                                  space));
      }
      for (auto& force_functional : force_functionals)
        system_assembler.add(*force_functional);
      // * dirichlet boundary
      const auto& dirichlet = *(this->problem_.dirichlet());
      typedef Functionals::DirichletBoundarySWIPDG< DiffusionFactorType, FunctionType, VectorType, TestSpaceType,
                                                    GridViewType, DiffusionTensorType > DirichletBoundaryFunctionalType;
      std::vector< std::unique_ptr< DirichletBoundaryFunctionalType > > dirichlet_boundary_functionals;
      if (diffusion_factor.has_affine_part() && dirichlet.has_affine_part()) {
        if (!rhs.has_affine_part())
          rhs.register_affine_part(space.mapper().size());
        dirichlet_boundary_functionals.emplace_back(new DirichletBoundaryFunctionalType(
            *(diffusion_factor.affine_part()),
            *(diffusion_tensor.affine_part()),
            *(dirichlet.affine_part()),
            boundary_info,
            *(rhs.affine_part()),
            space));
      }
      if (diffusion_factor.has_affine_part()) {
        for (size_t qq = 0; qq < dirichlet.num_components(); ++qq) {
          const size_t id = rhs.register_component(dirichlet.coefficient(qq), space.mapper().size());
          dirichlet_boundary_functionals.emplace_back(new DirichletBoundaryFunctionalType(
              *(diffusion_factor.affine_part()),
              *(diffusion_tensor.affine_part()),
              *(dirichlet.component(qq)),
              boundary_info,
              *(rhs.component(id)),
              space));
        }
      }
      if (dirichlet.has_affine_part()) {
        for (size_t qq = 0; qq < diffusion_factor.num_components(); ++qq) {
          const size_t id = rhs.register_component(diffusion_factor.coefficient(qq), space.mapper().size());
          dirichlet_boundary_functionals.emplace_back(new DirichletBoundaryFunctionalType(
              *(diffusion_factor.component(qq)),
              *(diffusion_tensor.affine_part()),
              *(dirichlet.affine_part()),
              boundary_info,
              *(rhs.component(id)),
              space));
        }
      }
      Pymor::ParameterType param;
      for (const auto& key : diffusion_factor.parameter_type().keys())
        param.set(key, diffusion_factor.parameter_type().get(key));
      for (const auto& key : dirichlet.parameter_type().keys())
        param.set(key, dirichlet.parameter_type().get(key));
      for (size_t pp = 0; pp < diffusion_factor.num_components(); ++ pp) {
        for (size_t qq = 0; qq < dirichlet.num_components(); ++qq) {
          const std::string expression = "(" + diffusion_factor.coefficient(pp)->expression()
                                             + ")*(" + dirichlet.coefficient(qq)->expression() + ")";
          const size_t id = rhs.register_component(param, expression, space.mapper().size());
          dirichlet_boundary_functionals.emplace_back(new DirichletBoundaryFunctionalType(
              *(diffusion_factor.component(pp)),
              *(diffusion_tensor.affine_part()),
              *(dirichlet.component(qq)),
              boundary_info,
              *(rhs.component(id)),
              space));
        }
      }
      for (auto& dirichlet_boundary_functional : dirichlet_boundary_functionals)
        system_assembler.add(*dirichlet_boundary_functional);

      // * neumann boundary
      const auto& neumann = *(this->problem_.neumann());
      typedef Functionals::L2Face< FunctionType, VectorType, TestSpaceType > L2FaceFunctionalType;
      std::vector< std::unique_ptr< L2FaceFunctionalType > > neumann_boundary_functionals;
      for (size_t qq = 0; qq < neumann.num_components(); ++qq) {
        const size_t id = rhs.register_component(neumann.coefficient(qq), space.mapper().size());
        neumann_boundary_functionals.emplace_back(new L2FaceFunctionalType(
            *(neumann.component(qq)),
            *(rhs.component(id)),
            space,
            new ApplyOn::NeumannIntersections< GridViewType >(boundary_info)));
      }
      if (neumann.has_affine_part()) {
        if (!rhs.has_affine_part())
          rhs.register_affine_part(space.mapper().size());
        neumann_boundary_functionals.emplace_back(new L2FaceFunctionalType(
            *(neumann.affine_part()),
            *(rhs.affine_part()),
            space,
            new ApplyOn::NeumannIntersections< GridViewType >(boundary_info)));
      }
      for (auto& neumann_boundary_functional : neumann_boundary_functionals)
        system_assembler.add(*neumann_boundary_functional);

      // products
      // * L2
      typedef Products::L2Assemblable< MatrixType, TestSpaceType > L2ProductType;
      auto l2_product_matrix = std::make_shared< AffinelyDecomposedMatrixType >();
      l2_product_matrix->register_affine_part(new MatrixType(space.mapper().size(),
                                                             space.mapper().size(),
                                                             L2ProductType::pattern(space)));
      L2ProductType l2_product(*(l2_product_matrix->affine_part()), space);
      system_assembler.add(l2_product);
      // * H1 semi
      typedef Products::H1SemiAssemblable< MatrixType, TestSpaceType > H1ProductType;
      auto h1_product_matrix = std::make_shared< AffinelyDecomposedMatrixType >();
      h1_product_matrix->register_affine_part(new MatrixType(space.mapper().size(),
                                                             space.mapper().size(),
                                                             H1ProductType::pattern(space)));
      H1ProductType h1_product(*(h1_product_matrix->affine_part()), space);
      system_assembler.add(h1_product);

      // do the actual assembling
      system_assembler.walk();
      out << "done (took " << timer.elapsed() << "s)" << std::endl;

      if (!dirichlet_detector.found())
        this->purely_neumann_ = true;

      // build parameter type
      this->inherit_parameter_type(matrix.parameter_type(), "lhs");
      this->inherit_parameter_type(rhs.parameter_type(), "rhs");

      // finalize
      this->products_.insert(std::make_pair("l2", l2_product_matrix));
      this->products_.insert(std::make_pair("h1_semi", h1_product_matrix));

      this->container_based_initialized_ = true;
    } // if (!this->container_based_initialized_)
  } // ... init(...)

  RangeFieldType estimate(const VectorType& vector, const std::string type) const
  {
    return internal::SWIPDGEstimator< ThisType, GridType, dimRange >::estimate(*this, vector, type);
  }

  friend class BlockSWIPDG< GridImp, RangeFieldImp, rangeDim, polynomialOrder, la_backend >;

  const RangeFieldType beta_;
  PatternType pattern_;
}; // class SWIPDG


} // namespace Discretizations
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_SWIPDG_HH
