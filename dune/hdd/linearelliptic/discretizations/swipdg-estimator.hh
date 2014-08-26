// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_SWIPDG_ESTIMATOR_HH
#define DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_SWIPDG_ESTIMATOR_HH

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
namespace internal {
namespace SWIPDGEstimators {


static const size_t over_integrate = 2;


class LocalNonconformityESV2007Base
{
public:
  static std::string id() { return "eta_NC_ESV2007"; }
};


template< class SpaceType, class VectorType, class ProblemType, class GridType/*, class GridViewType*/ >
class LocalNonconformityESV2007
  : public LocalNonconformityESV2007Base
{
public:
  static const bool available = false;
};

#if HAVE_ALUGRID

/**
 *  \brief computes the local nonconformity estimator as defined in ESV2007
 */
template< class SpaceType, class VectorType, class ProblemType/*, class GridViewType*/ >
class LocalNonconformityESV2007< SpaceType, VectorType, ProblemType, ALUGrid< 2, 2, simplex, conforming >/*, GridViewType*/ >
  : public LocalNonconformityESV2007Base
  , public GDT::Functor::Codim0< typename SpaceType::GridViewType >
{
  typedef LocalNonconformityESV2007< SpaceType, VectorType, ProblemType, ALUGrid< 2, 2, simplex, conforming >/*, GridViewType*/ >
    ThisType;
  typedef GDT::Functor::Codim0< typename SpaceType::GridViewType > FunctorBaseType;
public:
  static const bool available = true;

  typedef std::map< std::string, Pymor::Parameter > ParametersMapType;

  typedef typename FunctorBaseType::GridViewType GridViewType;
  typedef typename FunctorBaseType::EntityType   EntityType;

  typedef typename ProblemType::RangeFieldType   RangeFieldType;

private:
  typedef GDT::ConstDiscreteFunction< SpaceType, VectorType > ConstDiscreteFunctionType;
  typedef GDT::DiscreteFunction< SpaceType, VectorType > DiscreteFunctionType;
  typedef typename ConstDiscreteFunctionType::DifferenceType DifferenceType;

  typedef typename ProblemType::DiffusionFactorType::NonparametricType DiffusionFactorType;
  typedef typename ProblemType::DiffusionTensorType::NonparametricType DiffusionTensorType;

  typedef GDT::LocalOperator::Codim0Integral< GDT::LocalEvaluation::Elliptic< DiffusionFactorType,
                                                                              DiffusionTensorType > > LocalOperatorType;
  typedef GDT::TmpStorageProvider::Matrices< RangeFieldType > TmpStorageProviderType;

  static const ProblemType& assert_problem(const ProblemType& problem, const Pymor::Parameter& mu_bar)
  {
    if (mu_bar.type() != problem.diffusion_factor()->parameter_type())
      DUNE_THROW(Pymor::Exceptions::wrong_parameter_type,
                 "Given mu_bar is " << mu_bar.type() << " and should be " << problem.parameter_type() << "!");
    assert(!problem.diffusion_tensor()->parametric());
    return problem;
  } // ... assert_problem(...)

public:
  static RangeFieldType estimate(const SpaceType& space,
                                 const VectorType& vector,
                                 const ProblemType& problem,
                                 const ParametersMapType parameters = ParametersMapType())
  {
    const Pymor::Parameter mu_bar = problem.parametric() ? parameters.at("mu_bar") : Pymor::Parameter();
    ThisType estimator(space, vector, problem, mu_bar);
    GDT::GridWalker< GridViewType > grid_walker(*space.grid_view());
    grid_walker.add(estimator);
    grid_walker.walk();
    return std::sqrt(estimator.result_);
  } // ... estimate(...)

  LocalNonconformityESV2007(const SpaceType& space,
                            const VectorType& vector,
                            const ProblemType& problem,
                            const Pymor::Parameter mu_bar = Pymor::Parameter())
    : space_(space)
    , vector_(vector)
    , problem_(assert_problem(problem, mu_bar))
    , diffusion_factor_mu_bar_(problem_.diffusion_factor()->with_mu(mu_bar))
    , discrete_solution_(space_, vector_)
    , oswald_interpolation_(space_)
    , difference_(Stuff::Common::make_unique< DifferenceType >(discrete_solution_ - oswald_interpolation_))
    , local_operator_(over_integrate,
                      *diffusion_factor_mu_bar_,
                      *problem_.diffusion_tensor()->affine_part())
    , tmp_local_matrices_({1, local_operator_.numTmpObjectsRequired()}, 1, 1)
    , result_(0.0)
  {}

  virtual void prepare()
  {
    const GDT::Operators::OswaldInterpolation< GridViewType > oswald_interpolation_operator(*space_.grid_view());
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
  const SpaceType& space_;
  const VectorType& vector_;
  const ProblemType& problem_;
  const std::shared_ptr< const DiffusionFactorType > diffusion_factor_mu_bar_;
  const ConstDiscreteFunctionType discrete_solution_;
  DiscreteFunctionType oswald_interpolation_;
  std::unique_ptr< const DifferenceType > difference_;
  const LocalOperatorType local_operator_;
  TmpStorageProviderType tmp_local_matrices_;
public:
  RangeFieldType result_;
}; // class LocalNonconformityESV2007< ..., ALUGrid< 2, 2, simplex, conforming >, ... >

#endif // HAVE_ALUGRID


class LocalResidualESV2007Base
{
public:
  static std::string id() { return "eta_R_ESV2007"; }
};


template< class SpaceType, class VectorType, class ProblemType, class GridType/*, class GridViewType*/ >
class LocalResidualESV2007
  : public LocalResidualESV2007Base
{
public:
  static const bool available = false;
};

#if HAVE_ALUGRID

/**
 *  \brief computes the local residual estimator as defined in ESV2007
 */
template< class SpaceType, class VectorType, class ProblemType/*, class GridViewType*/ >
class LocalResidualESV2007< SpaceType, VectorType, ProblemType, ALUGrid< 2, 2, simplex, conforming >/*, GridViewType*/ >
  : public LocalResidualESV2007Base
  , public GDT::Functor::Codim0< typename SpaceType::GridViewType >
{
  typedef LocalResidualESV2007< SpaceType, VectorType, ProblemType, ALUGrid< 2, 2, simplex, conforming >/*, GridViewType*/ >
      ThisType;
  typedef GDT::Functor::Codim0< typename SpaceType::GridViewType > FunctorBaseType;
public:
  static const bool available = true;

  typedef typename FunctorBaseType::GridViewType GridViewType;
  typedef typename FunctorBaseType::EntityType   EntityType;

  typedef typename ProblemType::RangeFieldType RangeFieldType;

private:
  typedef GDT::Spaces::FiniteVolume::Default< GridViewType, RangeFieldType, 1, 1 > P0SpaceType;
  typedef GDT::DiscreteFunction< P0SpaceType, VectorType > DiscreteFunctionType;
  typedef typename DiscreteFunctionType::DifferenceType DifferenceType;

  typedef typename ProblemType::DiffusionFactorType::NonparametricType DiffusionFactorType;
  typedef typename ProblemType::DiffusionTensorType::NonparametricType DiffusionTensorType;

  typedef typename Stuff::Functions::ESV2007::Cutoff< DiffusionFactorType, DiffusionTensorType > CutoffFunctionType;
  typedef GDT::LocalOperator::Codim0Integral< GDT::LocalEvaluation::Product< CutoffFunctionType > > LocalOperatorType;
  typedef GDT::TmpStorageProvider::Matrices< RangeFieldType > TmpStorageProviderType;

  static const ProblemType& assert_problem(const ProblemType& problem)
  {
    if (problem.parametric())
      DUNE_THROW(NotImplemented, "Not implemented yet for parametric problems!");
    assert(problem.diffusion_factor()->has_affine_part());
    assert(problem.diffusion_tensor()->has_affine_part());
    assert(problem.force()->has_affine_part());
    return problem;
  } // ... assert_problem(...)

public:
  static RangeFieldType estimate(const SpaceType& space, const VectorType& /*vector*/, const ProblemType& problem)
  {
    ThisType estimator(space, problem);
    GDT::GridWalker< GridViewType > grid_walker(*space.grid_view());
    grid_walker.add(estimator);
    grid_walker.walk();
    return std::sqrt(estimator.result_);
  } // ... estimate(...)

  LocalResidualESV2007(const SpaceType& space, const ProblemType& problem)
    : space_(space)
    , problem_(assert_problem(problem))
    , p0_space_(space_.grid_view())
    , p0_force_(p0_space_)
    , difference_(Stuff::Common::make_unique< DifferenceType >(*problem_.force()->affine_part() - p0_force_))
    , cutoff_function_(*problem_.diffusion_factor()->affine_part(),
                       *problem_.diffusion_tensor()->affine_part())
    , local_operator_(over_integrate, cutoff_function_)
    , tmp_local_matrices_({1, local_operator_.numTmpObjectsRequired()}, 1, 1)
    , result_(0.0)
  {}

  virtual void prepare()
  {
    const GDT::Operators::Projection< GridViewType > projection_operator(*space_.grid_view());
    projection_operator.apply(*problem_.force()->affine_part(), p0_force_);
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
  const SpaceType& space_;
  const ProblemType& problem_;
  const P0SpaceType p0_space_;
  DiscreteFunctionType p0_force_;
  std::unique_ptr< const DifferenceType > difference_;
  const CutoffFunctionType cutoff_function_;
  const LocalOperatorType local_operator_;
  TmpStorageProviderType tmp_local_matrices_;
public:
  RangeFieldType result_;
}; // class LocalResidualESV2007< ..., ALUGrid< 2, 2, simplex, conforming >, ... >

#endif // HAVE_ALUGRID


class LocalDiffusiveFluxESV2007Base
{
public:
  static std::string id() { return "eta_DF_ESV2007"; }
};


template< class SpaceType, class VectorType, class ProblemType, class GridType/*, class GridViewType*/ >
class LocalDiffusiveFluxESV2007
  : public LocalDiffusiveFluxESV2007Base
{
public:
  static const bool available = false;
};

#if HAVE_ALUGRID

/**
 *  \brief computes the local diffusive flux estimator as defined in ESV2007
 */
template< class SpaceType, class VectorType, class ProblemType/*, class GridViewType*/ >
class LocalDiffusiveFluxESV2007< SpaceType, VectorType, ProblemType, ALUGrid< 2, 2, simplex, conforming >/*, GridViewType*/ >
  : public LocalDiffusiveFluxESV2007Base
  , public GDT::Functor::Codim0< typename SpaceType::GridViewType >
{
  typedef LocalDiffusiveFluxESV2007< SpaceType, VectorType, ProblemType, ALUGrid< 2, 2, simplex, conforming >/*, GridViewType*/ >
      ThisType;
  typedef GDT::Functor::Codim0< typename SpaceType::GridViewType > FunctorBaseType;
public:
  static const bool available = true;

  typedef typename FunctorBaseType::GridViewType GridViewType;
  typedef typename FunctorBaseType::EntityType   EntityType;

  typedef typename ProblemType::RangeFieldType RangeFieldType;

  static const unsigned int dimDomain = SpaceType::dimDomain;

private:
  typedef GDT::ConstDiscreteFunction< SpaceType, VectorType > ConstDiscreteFunctionType;
  typedef GDT::Spaces::RaviartThomas::PdelabBased< GridViewType, 0, RangeFieldType, dimDomain > RTN0SpaceType;
  typedef GDT::DiscreteFunction< RTN0SpaceType, VectorType > RTN0DiscreteFunctionType;

  typedef typename ProblemType::DiffusionFactorType::NonparametricType DiffusionFactorType;
  typedef typename ProblemType::DiffusionTensorType::NonparametricType DiffusionTensorType;

  typedef GDT::LocalOperator::Codim0Integral<
      GDT::LocalEvaluation::ESV2007::DiffusiveFluxEstimate< DiffusionFactorType,
                                                            RTN0DiscreteFunctionType,
                                                            DiffusionTensorType > > LocalOperatorType;
  typedef GDT::TmpStorageProvider::Matrices< RangeFieldType > TmpStorageProviderType;

  static const ProblemType& assert_problem(const ProblemType& problem)
  {
    if (problem.parametric())
      DUNE_THROW(NotImplemented, "Not implemented yet for parametric problems!");
    assert(problem.diffusion_factor()->has_affine_part());
    assert(problem.diffusion_tensor()->has_affine_part());
    return problem;
  } // ... assert_problem(...)

public:
  static RangeFieldType estimate(const SpaceType& space, const VectorType& vector, const ProblemType& problem)
  {
    ThisType estimator(space, vector, problem);
    GDT::GridWalker< GridViewType > grid_walker(*space.grid_view());
    grid_walker.add(estimator);
    grid_walker.walk();
    return std::sqrt(estimator.result_);
  } // ... estimate(...)

  LocalDiffusiveFluxESV2007(const SpaceType& space, const VectorType& vector, const ProblemType& problem)
    : space_(space)
    , vector_(vector)
    , problem_(assert_problem(problem))
    , discrete_solution_(space_, vector_)
    , rtn0_space_(space.grid_view())
    , diffusive_flux_(rtn0_space_)
    , local_operator_(over_integrate,
                      *problem_.diffusion_factor()->affine_part(),
                      *problem_.diffusion_tensor()->affine_part(),
                      diffusive_flux_)
    , tmp_local_matrices_({1, local_operator_.numTmpObjectsRequired()}, 1, 1)
    , result_(0.0)
  {}

  virtual void prepare()
  {
    const GDT::Operators::DiffusiveFluxReconstruction< GridViewType, DiffusionFactorType, DiffusionTensorType >
      diffusive_flux_reconstruction(*space_.grid_view(),
                                    *problem_.diffusion_factor()->affine_part(),
                                    *problem_.diffusion_tensor()->affine_part());
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
  const SpaceType& space_;
  const VectorType& vector_;
  const ProblemType& problem_;
  const ConstDiscreteFunctionType discrete_solution_;
  const RTN0SpaceType rtn0_space_;
  RTN0DiscreteFunctionType diffusive_flux_;
  const LocalOperatorType local_operator_;
  TmpStorageProviderType tmp_local_matrices_;
public:
  RangeFieldType result_;
}; // class LocalDiffusiveFluxESV2007< ..., ALUGrid< 2, 2, simplex, conforming >, ... >

#endif // HAVE_ALUGRID


class ESV2007Base
{
public:
  static std::string id()
  {
    return "eta_ESV2007";
  }
};


template< class SpaceType, class VectorType, class ProblemType, class GridType/*, class GridViewType*/ >
class ESV2007
  : public ESV2007Base
{
public:
  static const bool available = false;
};


#if HAVE_ALUGRID

template< class SpaceType, class VectorType, class ProblemType/*, class GridViewType*/ >
class ESV2007< SpaceType, VectorType, ProblemType, ALUGrid< 2, 2, simplex, conforming >/*, GridViewType*/ >
  : public ESV2007Base
{
  typedef ALUGrid< 2, 2, simplex, conforming > GridType;
public:
  static const bool available = true;

  typedef typename ProblemType::RangeFieldType RangeFieldType;

  static RangeFieldType estimate(const SpaceType& space, const VectorType& vector, const ProblemType& problem)
  {
    LocalNonconformityESV2007< SpaceType, VectorType, ProblemType, GridType > eta_nc(space, vector, problem);
    LocalResidualESV2007< SpaceType, VectorType, ProblemType, GridType >      eta_r(space, problem);
    LocalDiffusiveFluxESV2007< SpaceType, VectorType, ProblemType, GridType > eta_df(space, vector, problem);
    eta_nc.prepare();
    eta_r.prepare();
    eta_df.prepare();

    RangeFieldType eta_squared(0.0);

    const auto grid_view = space.grid_view();
    const auto entity_it_end = grid_view->template end< 0 >();
    for (auto entity_it = grid_view->template begin< 0 >(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      eta_squared += eta_nc.compute_locally(entity)
                   + std::pow(std::sqrt(eta_r.compute_locally(entity)) + std::sqrt(eta_df.compute_locally(entity)), 2);
    }
    return std::sqrt(eta_squared);
  } // ... estimate(...)
}; // class ESV2007< ..., ALUGrid< 2, 2, simplex, conforming >, ... >

#endif // HAVE_ALUGRID


class ESV2007AlternativeSummationBase
{
public:
  static std::string id()
  {
    return "eta_ESV2007_alt";
  }
};


template< class SpaceType, class VectorType, class ProblemType, class GridType/*, class GridViewType*/ >
class ESV2007AlternativeSummation
  : public ESV2007AlternativeSummationBase
{
public:
  static const bool available = false;
};


#if HAVE_ALUGRID

template< class SpaceType, class VectorType, class ProblemType/*, class GridViewType*/ >
class ESV2007AlternativeSummation< SpaceType, VectorType, ProblemType, ALUGrid< 2, 2, simplex, conforming >/*, GridViewType*/ >
  : public ESV2007AlternativeSummationBase
{
  typedef ALUGrid< 2, 2, simplex, conforming > GridType;
public:
  static const bool available = true;

  typedef typename ProblemType::RangeFieldType RangeFieldType;

  static RangeFieldType estimate(const SpaceType& space, const VectorType& vector, const ProblemType& problem)
  {
    LocalNonconformityESV2007< SpaceType, VectorType, ProblemType, GridType > eta_nc(space, vector, problem);
    LocalResidualESV2007< SpaceType, VectorType, ProblemType, GridType >      eta_r(space, problem);
    LocalDiffusiveFluxESV2007< SpaceType, VectorType, ProblemType, GridType > eta_df(space, vector, problem);
    eta_nc.prepare();
    eta_r.prepare();
    eta_df.prepare();

    RangeFieldType eta_nc_squared(0.0);
    RangeFieldType eta_r_squared(0.0);
    RangeFieldType eta_df_squared(0.0);

    const auto grid_view = space.grid_view();
    const auto entity_it_end = grid_view->template end< 0 >();
    for (auto entity_it = grid_view->template begin< 0 >(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      eta_nc_squared += eta_nc.compute_locally(entity);
      eta_r_squared += eta_r.compute_locally(entity);
      eta_df_squared += eta_df.compute_locally(entity);
    }
    return std::sqrt(eta_nc_squared) + std::sqrt(eta_r_squared) + std::sqrt(eta_df_squared);
  } // ... estimate(...)
}; // class ESV2007AlternativeSummation< ..., ALUGrid< 2, 2, simplex, conforming >, ... >

#endif // HAVE_ALUGRID


} // namespace SWIPDGEstimators
} // namespace internal


template< class SpaceType, class VectorType, class ProblemType, class GridType >
class SWIPDGEstimator
{
public:
  typedef typename ProblemType::RangeFieldType RangeFieldType;

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

    static RangeFieldType estimate(const SpaceType& /*space*/,
                                   const VectorType& /*vector*/,
                                   const ProblemType& /*problem*/)
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

    static RangeFieldType estimate(const SpaceType& space, const VectorType& vector, const ProblemType& problem)
    {
      return IndividualEstimator::estimate(space, vector, problem);
    }
  }; // class Caller< ..., true >

  template< class IndividualEstimator >
  static std::vector< std::string > call_append(std::vector< std::string > in)
  {
    return Caller< IndividualEstimator, IndividualEstimator::available >::append(in);
  }

  template< class IndividualEstimator >
  static bool call_equals(const std::string& type)
  {
    return Caller< IndividualEstimator, IndividualEstimator::available >::equals(type);
  }

  template< class IndividualEstimator >
  static RangeFieldType call_estimate(const SpaceType& space, const VectorType& vector, const ProblemType& problem)
  {
    return Caller< IndividualEstimator, IndividualEstimator::available >::estimate(space, vector, problem);
  }

  typedef internal::SWIPDGEstimators::LocalNonconformityESV2007
      < SpaceType, VectorType, ProblemType, GridType >              LocalNonconformityESV2007Type;
  typedef internal::SWIPDGEstimators::LocalResidualESV2007
      < SpaceType, VectorType, ProblemType, GridType >              LocalResidualESV2007Type;
  typedef internal::SWIPDGEstimators::LocalDiffusiveFluxESV2007
      < SpaceType, VectorType, ProblemType, GridType >              LocalDiffusiveFluxESV2007Type;
  typedef internal::SWIPDGEstimators::ESV2007
      < SpaceType, VectorType, ProblemType, GridType >              ESV2007Type;
  typedef internal::SWIPDGEstimators::ESV2007AlternativeSummation
      < SpaceType, VectorType, ProblemType, GridType >              ESV2007AlternativeSummationType;

public:
  static std::vector< std::string > available()
  {
    std::vector< std::string > tmp;
    tmp = call_append< LocalNonconformityESV2007Type >(tmp);
    tmp = call_append< LocalResidualESV2007Type >(tmp);
    tmp = call_append< LocalDiffusiveFluxESV2007Type >(tmp);
    tmp = call_append< ESV2007Type >(tmp);
    tmp = call_append< ESV2007AlternativeSummationType >(tmp);
    return tmp;
  } // ... available(...)

  static RangeFieldType estimate(const SpaceType& space,
                                 const VectorType& vector,
                                 const ProblemType& problem,
                                 const std::string type)
  {
    if (call_equals< LocalNonconformityESV2007Type >(type))
      return call_estimate< LocalNonconformityESV2007Type >(space, vector, problem);
    else if (call_equals< LocalResidualESV2007Type >(type))
      return call_estimate< LocalResidualESV2007Type >(space, vector, problem);
    else if (call_equals< LocalDiffusiveFluxESV2007Type >(type))
      return call_estimate< LocalDiffusiveFluxESV2007Type >(space, vector, problem);
    else if (call_equals< ESV2007Type >(type))
      return call_estimate< ESV2007Type >(space, vector, problem);
    else if (call_equals< ESV2007AlternativeSummationType >(type))
      return call_estimate< ESV2007AlternativeSummationType >(space, vector, problem);
    else
      DUNE_THROW(Stuff::Exceptions::you_are_using_this_wrong,
                 "Requested type '" << type << "' is not one of available()!");
  } // ... estimate(...)
}; // class SWIPDGEstimator


#if HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_EIGEN


extern template class SWIPDGEstimator< GDT::Spaces::DiscontinuousLagrange::FemBased< Fem::LeafGridPart< ALUGrid< 2, 2, simplex, conforming > >, 1, double, 1, 1 >,
                                       Dune::Stuff::LA::EigenDenseVector< double > ,
                                       ProblemInterface< typename ALUGrid< 2, 2, simplex, conforming >::template Codim< 0 >::Entity, double, 2, double, 1 >,
                                       ALUGrid< 2, 2, simplex, conforming > >;


#endif // HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_EIGEN

} // namespace Discretizations
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_SWIPDG_ESTIMATOR_HH
