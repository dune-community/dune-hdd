// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_BLOCK_SWIPDG_ESTIMATOR_HH
#define DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_BLOCK_SWIPDG_ESTIMATOR_HH

#include <memory>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <limits>
#include <type_traits>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/stuff/common/disable_warnings.hh>
# if HAVE_EIGEN
#   include <Eigen/Eigenvalues>
# endif

# include <dune/common/timer.hh>
# include <dune/common/dynmatrix.hh>

# if HAVE_ALUGRID
#   include <dune/grid/alugrid.hh>
# endif

# include <dune/geometry/quadraturerules.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/grid/multiscale/provider.hh>

#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/float_cmp.hh>
#include <dune/stuff/common/fixed_map.hh>
#include <dune/stuff/grid/layers.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/la/container.hh>
#include <dune/stuff/la/solver.hh>

#include <dune/pymor/common/exceptions.hh>

#include <dune/gdt/spaces/discontinuouslagrange.hh>
#include <dune/gdt/playground/spaces/block.hh>
#include <dune/gdt/playground/localevaluation/swipdg.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/oswaldinterpolation.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/playground/spaces/finitevolume/default.hh>
#include <dune/gdt/playground/spaces/raviartthomas/pdelab.hh>
#include <dune/gdt/playground/operators/fluxreconstruction.hh>
#include <dune/gdt/products/l2.hh>
#include <dune/gdt/products/h1.hh>
#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/assembler/gridwalker.hh>
#include <dune/gdt/playground/localevaluation/OS2014.hh>

#include <dune/hdd/linearelliptic/problems/default.hh>
#include <dune/hdd/playground/linearelliptic/problems/zero-boundary.hh>

#include "../../../linearelliptic/discretizations/base.hh"
#include "../../../linearelliptic/discretizations/swipdg.hh"
#include "../../../linearelliptic/discretizations/swipdg-estimator.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Discretizations {
namespace internal {


template< class FunctionType, class EntityType, int dimRange, int dimRangeCols >
class Minimum
{
  static_assert(AlwaysFalse< FunctionType >::value, "Not implemented for these dimensions!");
};

template< class FunctionType, class EntityType >
class Minimum< FunctionType, EntityType, 1, 1 >
{
  typedef typename FunctionType::RangeFieldType RangeFieldType;

public:
  /**
   * We try to find the minimum of a polynomial of given order by evaluating it at the points of a quadrature that
   * would integrate this polynomial exactly.
   * \todo These are just some heuristics and should be replaced by something proper.
   */
  static RangeFieldType compute(const FunctionType& function, const EntityType& entity)
  {
    typename FunctionType::RangeType tmp_value(0);
    RangeFieldType ret = std::numeric_limits< RangeFieldType >::max();
    const auto local_function = function.local_function(entity);
    const size_t ord = local_function->order();
    assert(ord < std::numeric_limits< int >::max());
    const auto& quadrature = QuadratureRules< typename FunctionType::DomainFieldType
                                            , FunctionType::dimDomain >::rule(entity.type(), int(ord));
    const auto quad_point_it_end = quadrature.end();
    for (auto quad_point_it = quadrature.begin(); quad_point_it != quad_point_it_end; ++quad_point_it) {
      local_function->evaluate(quad_point_it->position(), tmp_value);
      ret = std::min(ret, tmp_value[0]);
    }
    return ret;
  } // ... compute(...)
}; // class Minimum< ..., 1, 1 >

template< class FunctionType, class EntityType, int dimDomain >
class Minimum< FunctionType, EntityType, dimDomain, dimDomain >
{
  typedef typename FunctionType::RangeFieldType RangeFieldType;

public:
  static RangeFieldType compute(const FunctionType& function, const EntityType& entity)
  {
#if HAVE_EIGEN
    const auto local_function = function.local_function(entity);
    assert(local_function->order() == 0);
    const auto& reference_element = ReferenceElements< typename FunctionType::DomainFieldType
                                                     , FunctionType::dimDomain >::general(entity.type());
    const Stuff::LA::EigenDenseMatrix< RangeFieldType >
        tensor = local_function->evaluate(reference_element.position(0, 0));
    ::Eigen::EigenSolver< typename Stuff::LA::EigenDenseMatrix< RangeFieldType >::BackendType >
        eigen_solver(tensor.backend());
    assert(eigen_solver.info() == ::Eigen::Success);
    const auto eigenvalues = eigen_solver.eigenvalues(); // <- this should be an Eigen vector of std::complex
    RangeFieldType min_ev = std::numeric_limits< RangeFieldType >::max();
    for (size_t ii = 0; ii < boost::numeric_cast< size_t >(eigenvalues.size()); ++ii) {
      // assert this is real
      assert(std::abs(eigenvalues[ii].imag()) < 1e-15);
      // assert that this eigenvalue is positive
      const RangeFieldType eigenvalue = eigenvalues[ii].real();
      assert(eigenvalue > 1e-15);
      min_ev = std::min(min_ev, eigenvalue);
    }
    return min_ev;
#else
    static_assert(AlwaysFalse< FunctionType >::value, "You are missing eigen!");
#endif
  } // ... compute(...)
}; // class Minimum< ..., dimDomain, dimDomain >


template< class FunctionType, class EntityType >
static typename FunctionType::RangeFieldType compute_minimum(const FunctionType& function, const EntityType& entity)
{
  static const int dimRange = FunctionType::dimRange;
  static const int dimRangeCols = FunctionType::dimRangeCols;
  return Minimum< FunctionType, EntityType, dimRange, dimRangeCols >::compute(function, entity);
} // ... compute_minimum(...)


namespace BlockSWIPDGEstimators {


template< class BlockSpaceType, class VectorType, class ProblemType, class GridType >
class LocalNonconformityOS2014
  : public SWIPDGEstimators::LocalNonconformityESV2007< BlockSpaceType, VectorType, ProblemType, GridType >
{
public:
  static std::string id() { return "eta_NC_OS2014"; }
};


class LocalResidualOS2014Base
{
public:
  static std::string id() { return "eta_R_OS2014"; }
};


template< class BlockSpaceType, class VectorType, class ProblemType, class GridType >
class LocalResidualOS2014
  : public LocalResidualOS2014Base
{
public:
  static const bool available = false;
};

#if HAVE_ALUGRID

template< class BlockSpaceType, class VectorType, class ProblemType >
class LocalResidualOS2014< BlockSpaceType, VectorType, ProblemType, ALUGrid< 2, 2, simplex, conforming > >
  : public LocalResidualOS2014Base
  , public GDT::Functor::Codim0< typename BlockSpaceType::LocalSpaceType::GridViewType >
{
  typedef LocalResidualOS2014< BlockSpaceType, VectorType, ProblemType, ALUGrid< 2, 2, simplex, conforming > >
                                                                                        ThisType;
  typedef GDT::Functor::Codim0< typename BlockSpaceType::LocalSpaceType::GridViewType > FunctorBaseType;
public:
  static const bool available = true;

  typedef typename FunctorBaseType::GridViewType GridViewType;
  typedef typename FunctorBaseType::EntityType   EntityType;

  typedef typename ProblemType::RangeFieldType RangeFieldType;

  typedef std::map< std::string, Pymor::Parameter > ParametersMapType;

private:
  static const unsigned int dimDomain = GridViewType::dimension;

  typedef typename BlockSpaceType::LocalSpaceType LocalSpaceType;

  typedef GDT::Spaces::FiniteVolume::Default< GridViewType, RangeFieldType, 1, 1 > P0SpaceType;
  typedef GDT::DiscreteFunction< P0SpaceType, VectorType > DiscreteFunctionType;
  typedef typename DiscreteFunctionType::DifferenceType DifferenceType;

  typedef typename ProblemType::DiffusionFactorType::NonparametricType DiffusionFactorType;
  typedef typename ProblemType::DiffusionTensorType::NonparametricType DiffusionTensorType;

  typedef typename ProblemType::DomainFieldType DomainFieldType;
  typedef Stuff::Functions::Constant< EntityType, DomainFieldType, dimDomain, RangeFieldType, 1 >
      ConstantFunctionType;
  typedef typename ConstantFunctionType::DomainType DomainType;
  typedef GDT::LocalOperator::Codim0Integral< GDT::LocalEvaluation::Product< ConstantFunctionType > > LocalOperatorType;
  typedef GDT::TmpStorageProvider::Matrices< RangeFieldType > TmpStorageProviderType;

  static const ProblemType& assert_problem(const ProblemType& problem, const Pymor::Parameter& mu_minimizing)
  {
    if (mu_minimizing.type() != problem.parameter_type())
      DUNE_THROW(Pymor::Exceptions::wrong_parameter_type,
                 "Given mu_minimizing is of type " << mu_minimizing.type() << " and should be of type "
                 << problem.parameter_type() << "!");
    if (problem.diffusion_tensor()->parametric())
      DUNE_THROW(NotImplemented, "Not implemented for parametric diffusion_tensor!");
    if (problem.force()->parametric())
      DUNE_THROW(NotImplemented, "Not implemented for parametric force!");
    return problem;
  } // ... assert_problem(...)

public:
  static RangeFieldType estimate(const BlockSpaceType& space,
                                 const VectorType& /*vector*/,
                                 const ProblemType& problem,
                                 const ParametersMapType parameters = ParametersMapType())
  {
    if (problem.diffusion_factor()->parametric() && parameters.find("mu_minimizing") == parameters.end())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Given parameters are missing 'mu_minimizing'!");
    const Pymor::Parameter mu_minimizing = problem.parametric() ? parameters.at("mu_minimizing") : Pymor::Parameter();
    // walk the subdomains
    double eta_r_squared = 0.0;
    for (size_t subdomain = 0; subdomain < space.ms_grid()->size(); ++subdomain) {
      const auto local_space = space.local_spaces()[subdomain];
      ThisType eta_r_T(*local_space, problem, mu_minimizing);
      GDT::GridWalker< GridViewType > grid_walker(*local_space->grid_view());
      grid_walker.add(eta_r_T);
      grid_walker.walk();
      eta_r_squared += eta_r_T.result_;
    } // walk the subdomains
    return std::sqrt(eta_r_squared);
  } // ... estimate(...)

  LocalResidualOS2014(const LocalSpaceType& local_space,
                      const ProblemType& problem,
                      const Pymor::Parameter mu_minimizing = Pymor::Parameter())
    : local_space_(local_space)
    , problem_(assert_problem(problem, mu_minimizing))
    , problem_mu_minimizing_(problem_.with_mu(mu_minimizing))
    , p0_space_(local_space_.grid_view())
    , p0_force_(p0_space_)
    , difference_(Stuff::Common::make_unique< DifferenceType >(*problem_.force()->affine_part() - p0_force_))
    , constant_one_(1)
    , local_operator_(SWIPDGEstimators::over_integrate, constant_one_)
    , tmp_local_matrices_({1, local_operator_.numTmpObjectsRequired()}, 1, 1)
    , vertices_()
    , min_diffusion_value_(std::numeric_limits< RangeFieldType >::max())
    , result_(0.0)
  {}

  virtual void prepare()
  {
    const GDT::Operators::Projection< GridViewType > projection_operator(*local_space_.grid_view());
    projection_operator.apply(*problem_.force()->affine_part(), p0_force_);
    result_ = 0.0;
  } // ... prepare(...)

  virtual void apply_local(const EntityType &entity)
  {
    // collect all vertices
    for (int cc = 0; cc < entity.template count< dimDomain >(); ++cc)
      vertices_.push_back(entity.template subEntity< dimDomain >(cc)->geometry().center());
    // compute minimum diffusion
    const RangeFieldType min_diffusion_value_entity
        =   internal::compute_minimum(*problem_mu_minimizing_->diffusion_factor()->affine_part(), entity)
          * internal::compute_minimum(*problem_.diffusion_tensor()->affine_part(), entity);
    min_diffusion_value_ = std::min(min_diffusion_value_, min_diffusion_value_entity);
    // compute the local product
    const auto local_difference = difference_->local_function(entity);
    local_operator_.apply(*local_difference,
                          *local_difference,
                          tmp_local_matrices_.matrices()[0][0],
                          tmp_local_matrices_.matrices()[1]);
    assert(tmp_local_matrices_.matrices()[0][0].rows() >= 1);
    assert(tmp_local_matrices_.matrices()[0][0].cols() >= 1);
    result_ += tmp_local_matrices_.matrices()[0][0][0][0];
  } // ... apply_local(...)

  virtual void finalize()
  {
    // compute diameter
    DomainFieldType diameter(0);
    DomainType tmp_diff;
    for (size_t cc = 0; cc < vertices_.size(); ++cc) {
      const auto& vertex = vertices_[cc];
      for (size_t dd = cc + 1; dd < vertices_.size(); ++dd) {
        const auto& other_vertex = vertices_[dd];
        tmp_diff = vertex - other_vertex;
        diameter = std::max(diameter, tmp_diff.two_norm());
      }
    }
    // compute local estimator
    const RangeFieldType poincare_constant = 1.0 / (M_PIl * M_PIl);
    assert(min_diffusion_value_ > 0.0);
    result_ *= ((poincare_constant * diameter * diameter) / min_diffusion_value_);
  } // ... finalize(...)

private:
  const LocalSpaceType& local_space_;
  const ProblemType& problem_;
  const std::shared_ptr< const typename ProblemType::NonparametricType > problem_mu_minimizing_;
  const P0SpaceType p0_space_;
  DiscreteFunctionType p0_force_;
  std::unique_ptr< const DifferenceType > difference_;
  const ConstantFunctionType constant_one_;
  const LocalOperatorType local_operator_;
  TmpStorageProviderType tmp_local_matrices_;
  std::vector< DomainType > vertices_;
  RangeFieldType min_diffusion_value_;
public:
  RangeFieldType result_;
}; // class LocalResidualOS2014< ..., ALUGrid< 2, 2, simplex, conforming >, ... >

#endif // HAVE_ALUGRID


template< class BlockSpaceType, class VectorType, class ProblemType, class GridType >
class LocalDiffusiveFluxOS2014
  : public SWIPDGEstimators::LocalDiffusiveFluxESV2007< BlockSpaceType, VectorType, ProblemType, GridType >
{
public:
  static std::string id() { return "eta_DF_OS2014"; }
};


class LocalDiffusiveFluxOS2014StarBase
{
public:
  static std::string id() { return "eta_DF_OS2014_*"; }
};


template< class BlockSpaceType, class VectorType, class ProblemType, class GridType >
class LocalDiffusiveFluxOS2014Star
  : public LocalDiffusiveFluxOS2014StarBase
{
public:
  static const bool available = false;
};

#if HAVE_ALUGRID

template< class BlockSpaceType, class VectorType, class ProblemType >
class LocalDiffusiveFluxOS2014Star< BlockSpaceType, VectorType, ProblemType, ALUGrid< 2, 2, simplex, conforming > >
  : public LocalDiffusiveFluxOS2014StarBase
  , public GDT::Functor::Codim0< typename BlockSpaceType::GridViewType >
{
  typedef LocalDiffusiveFluxOS2014Star
      < BlockSpaceType, VectorType, ProblemType, ALUGrid< 2, 2, simplex, conforming > > ThisType;
  typedef GDT::Functor::Codim0< typename BlockSpaceType::GridViewType > FunctorBaseType;
public:
  static const bool available = true;

  typedef std::map< std::string, Pymor::Parameter > ParametersMapType;

  typedef typename FunctorBaseType::GridViewType GridViewType;
  typedef typename FunctorBaseType::EntityType   EntityType;

  typedef typename ProblemType::RangeFieldType RangeFieldType;

  static const unsigned int dimDomain = BlockSpaceType::dimDomain;

private:
  typedef GDT::ConstDiscreteFunction< BlockSpaceType, VectorType > ConstDiscreteFunctionType;
  typedef GDT::Spaces::RaviartThomas::PdelabBased< GridViewType, 0, RangeFieldType, dimDomain > RTN0SpaceType;
  typedef GDT::DiscreteFunction< RTN0SpaceType, VectorType > RTN0DiscreteFunctionType;

  typedef typename ProblemType::DiffusionFactorType::NonparametricType DiffusionFactorType;
  typedef typename ProblemType::DiffusionTensorType::NonparametricType DiffusionTensorType;

  typedef GDT::LocalOperator::Codim0Integral<
      GDT::LocalEvaluation::OS2014::DiffusiveFluxEstimateStar< DiffusionFactorType,
                                                               DiffusionFactorType,
                                                               DiffusionTensorType,
                                                               RTN0DiscreteFunctionType > >
                                                              LocalOperatorType;
  typedef GDT::TmpStorageProvider::Matrices< RangeFieldType > TmpStorageProviderType;

  static const ProblemType& assert_problem(const ProblemType& problem,
                                           const Pymor::Parameter& mu,
                                           const Pymor::Parameter& mu_hat)
  {
    if (mu.type() != problem.parameter_type())
      DUNE_THROW(Pymor::Exceptions::wrong_parameter_type,
                 "Given mu is of type " << mu.type() << " and should be of type " << problem.parameter_type()
                 << "!");
    if (mu_hat.type() != problem.parameter_type())
      DUNE_THROW(Pymor::Exceptions::wrong_parameter_type,
                 "Given mu_hat is of type " << mu_hat.type() << " and should be of type " << problem.parameter_type()
                 << "!");
    if (problem.diffusion_tensor()->parametric())
      DUNE_THROW(NotImplemented, "Not implemented for parametric diffusion_tensor!");
    return problem;
  } // ... assert_problem(...)

public:
  static RangeFieldType estimate(const BlockSpaceType& space,
                                 const VectorType& vector,
                                 const ProblemType& problem,
                                 const ParametersMapType parameters = ParametersMapType())
  {
    if (problem.diffusion_factor()->parametric() && parameters.find("mu") == parameters.end())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Given parameters are missing 'mu'!");
    if (problem.diffusion_factor()->parametric() && parameters.find("mu_hat") == parameters.end())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Given parameters are missing 'mu_hat'!");
    const Pymor::Parameter mu =     problem.parametric() ? parameters.at("mu")     : Pymor::Parameter();
    const Pymor::Parameter mu_hat = problem.parametric() ? parameters.at("mu_hat") : Pymor::Parameter();
    ThisType estimator(space, vector, problem, mu, mu_hat);
    GDT::GridWalker< GridViewType > grid_walker(*space.grid_view());
    grid_walker.add(estimator);
    grid_walker.walk();
    return std::sqrt(estimator.result_);
  } // ... estimate(...)

  LocalDiffusiveFluxOS2014Star(const BlockSpaceType& space,
                               const VectorType& vector,
                               const ProblemType& problem,
                               const Pymor::Parameter mu = Pymor::Parameter(),
                               const Pymor::Parameter mu_hat = Pymor::Parameter())
    : space_(space)
    , vector_(vector)
    , problem_(assert_problem(problem, mu, mu_hat))
    , problem_mu_(problem_.with_mu(mu))
    , problem_mu_hat_(problem.with_mu(mu_hat))
    , discrete_solution_(space_, vector_)
    , rtn0_space_(space.grid_view())
    , diffusive_flux_(rtn0_space_)
    , local_operator_(SWIPDGEstimators::over_integrate,
                      *problem_mu_->diffusion_factor()->affine_part(),
                      *problem_mu_hat_->diffusion_factor()->affine_part(),
                      *problem_.diffusion_tensor()->affine_part(),
                      diffusive_flux_)
    , tmp_local_matrices_({1, local_operator_.numTmpObjectsRequired()}, 1, 1)
    , result_(0.0)
  {}

  virtual void prepare()
  {
    const GDT::Operators::DiffusiveFluxReconstruction< GridViewType, DiffusionFactorType, DiffusionTensorType >
      diffusive_flux_reconstruction(*space_.grid_view(),
                                    *problem_mu_->diffusion_factor()->affine_part(),
                                    *problem_.diffusion_tensor()->affine_part());
    diffusive_flux_reconstruction.apply(discrete_solution_, diffusive_flux_);
    result_ = 0.0;
  } // ... prepare(...)

  virtual void apply_local(const EntityType &entity)
  {
    const auto local_discrete_solution = discrete_solution_.local_function(entity);
    local_operator_.apply(*local_discrete_solution,
                          *local_discrete_solution,
                          tmp_local_matrices_.matrices()[0][0],
                          tmp_local_matrices_.matrices()[1]);
    assert(tmp_local_matrices_.matrices()[0][0].rows() >= 1);
    assert(tmp_local_matrices_.matrices()[0][0].cols() >= 1);
    result_ += tmp_local_matrices_.matrices()[0][0][0][0];
  }

private:
  const BlockSpaceType& space_;
  const VectorType& vector_;
  const ProblemType& problem_;
  const std::shared_ptr< const typename ProblemType::NonparametricType > problem_mu_;
  const std::shared_ptr< const typename ProblemType::NonparametricType > problem_mu_hat_;
  const ConstDiscreteFunctionType discrete_solution_;
  const RTN0SpaceType rtn0_space_;
  RTN0DiscreteFunctionType diffusive_flux_;
  const LocalOperatorType local_operator_;
  TmpStorageProviderType tmp_local_matrices_;
public:
  RangeFieldType result_;
}; // class LocalDiffusiveFluxOS2014Star< ..., ALUGrid< 2, 2, simplex, conforming >, ... >

#endif // HAVE_ALUGRID


class OS2014Base
{
public:
  static std::string id() { return "eta_OS2014"; }
};


template< class BlockSpaceType, class VectorType, class ProblemType, class GridType >
class OS2014
  : public OS2014Base
{
public:
  static const bool available = false;
};

#if HAVE_ALUGRID

template< class BlockSpaceType, class VectorType, class ProblemType >
class OS2014< BlockSpaceType, VectorType, ProblemType, ALUGrid< 2, 2, simplex, conforming > >
  : public OS2014Base
{
  typedef ALUGrid< 2, 2, simplex, conforming > GridType;
public:
  static const bool available = true;

  typedef std::map< std::string, Pymor::Parameter > ParametersMapType;

  typedef typename ProblemType::RangeFieldType RangeFieldType;

  static RangeFieldType estimate(const BlockSpaceType& space,
                                 const VectorType& vector,
                                 const ProblemType& problem,
                                 const ParametersMapType parameters = ParametersMapType())
  {
    // check parameters
    if (problem.diffusion_factor()->parametric() && parameters.find("mu") == parameters.end())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Given parameters are missing 'mu'!");
    if (problem.diffusion_factor()->parametric() && parameters.find("mu_hat") == parameters.end())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Given parameters are missing 'mu_hat'!");
    if (problem.diffusion_factor()->parametric() && parameters.find("mu_bar") == parameters.end())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Given parameters are missing 'mu_bar'!");
    if (problem.diffusion_tensor()->parametric())
      DUNE_THROW(Stuff::Exceptions::requirements_not_met, "Not implemented for parametric diffusion_tensor!");
    if (problem.force()->parametric())
      DUNE_THROW(Stuff::Exceptions::requirements_not_met, "Not implemented for parametric force!");
    if (problem.dirichlet()->parametric())
      DUNE_THROW(Stuff::Exceptions::requirements_not_met, "Not implemented for parametric dirichlet!");
    if (problem.neumann()->parametric())
      DUNE_THROW(Stuff::Exceptions::requirements_not_met, "Not implemented for parametric neumann!");
    const Pymor::Parameter mu =     problem.parametric() ? parameters.at("mu")     : Pymor::Parameter();
    const Pymor::Parameter mu_hat = problem.parametric() ? parameters.at("mu_hat") : Pymor::Parameter();
    const Pymor::Parameter mu_bar = problem.parametric() ? parameters.at("mu_bar") : Pymor::Parameter();
    // compute parameter factors
    const double alpha_mu_mu_bar = problem.diffusion_factor()->alpha(mu, mu_bar);
    const double alpha_mu_mu_hat = problem.diffusion_factor()->alpha(mu, mu_hat);
    const double gamma_mu_mu_bar = problem.diffusion_factor()->gamma(mu, mu_bar);
    const double gamma_mu_mu_hat = problem.diffusion_factor()->gamma(mu, mu_hat);
    assert(alpha_mu_mu_bar > 0.0);
    assert(alpha_mu_mu_hat > 0.0);
    assert(gamma_mu_mu_bar > 0.0);
    assert(gamma_mu_mu_hat > 0.0);
    const double sqrt_gamma_tile = std::max(std::sqrt(gamma_mu_mu_hat), 1.0/std::sqrt(alpha_mu_mu_hat));
    // compute estimator
    typedef LocalNonconformityOS2014< BlockSpaceType, VectorType, ProblemType, GridType > LocalNonconformityOS2014Type;
    typedef LocalResidualOS2014< BlockSpaceType, VectorType, ProblemType, GridType >      LocalResidualOS2014Type;
    typedef LocalDiffusiveFluxOS2014< BlockSpaceType, VectorType, ProblemType, GridType > LocalDiffusiveFluxOS2014Type;
    return
        (1.0/std::sqrt(alpha_mu_mu_bar)) * (
            std::sqrt(gamma_mu_mu_bar) * LocalNonconformityOS2014Type::estimate(space, vector, problem, parameters)
          +                              LocalResidualOS2014Type::estimate(space, vector, problem, parameters)
          + sqrt_gamma_tile            * LocalDiffusiveFluxOS2014Type::estimate(space, vector, problem, parameters)
        );
  } // ... estimate(...)
}; // class OS2014

#endif // HAVE_ALUGRID

} // namespace BlockSWIPDGEstimators
} // namespace internal


template< class BlockSpaceType, class VectorType, class ProblemType, class GridType >
class BlockSWIPDGEstimator
{
public:
  typedef typename ProblemType::RangeFieldType RangeFieldType;
  typedef std::map< std::string, Pymor::Parameter > ParametersMapType;

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

    static RangeFieldType estimate(const BlockSpaceType& /*space*/,
                                   const VectorType& /*vector*/,
                                   const ProblemType& /*problem*/,
                                   const ParametersMapType& /*parameters*/)
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

    static RangeFieldType estimate(const BlockSpaceType& space,
                                   const VectorType& vector,
                                   const ProblemType& problem,
                                   const ParametersMapType& parameters)
    {
      return IndividualEstimator::estimate(space, vector, problem, parameters);
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
  static RangeFieldType call_estimate(const BlockSpaceType& space,
                                      const VectorType& vector,
                                      const ProblemType& problem,
                                      const ParametersMapType& parameters = ParametersMapType())
  {
    return Caller< IndividualEstimator, IndividualEstimator::available >::estimate(space, vector, problem, parameters);
  }

  typedef internal::BlockSWIPDGEstimators::LocalNonconformityOS2014
      < BlockSpaceType, VectorType, ProblemType, GridType >         LocalNonconformityOS2014Type;
  typedef internal::BlockSWIPDGEstimators::LocalResidualOS2014
      < BlockSpaceType, VectorType, ProblemType, GridType >         LocalResidualOS2014Type;
  typedef internal::BlockSWIPDGEstimators::LocalDiffusiveFluxOS2014
      < BlockSpaceType, VectorType, ProblemType, GridType >         LocalDiffusiveFluxOS2014Type;
  typedef internal::BlockSWIPDGEstimators::LocalDiffusiveFluxOS2014Star
      < BlockSpaceType, VectorType, ProblemType, GridType >         LocalDiffusiveFluxOS2014StarType;
  typedef internal::BlockSWIPDGEstimators::OS2014
      < BlockSpaceType, VectorType, ProblemType, GridType >         OS2014Type;

public:
  static std::vector< std::string > available()
  {
    std::vector< std::string > tmp;
    tmp = call_append< LocalNonconformityOS2014Type >(tmp);
    tmp = call_append< LocalResidualOS2014Type >(tmp);
    tmp = call_append< LocalDiffusiveFluxOS2014Type >(tmp);
    tmp = call_append< LocalDiffusiveFluxOS2014StarType >(tmp);
    tmp = call_append< OS2014Type >(tmp);
    return tmp;
  } // ... available(...)

  static RangeFieldType estimate(const BlockSpaceType& space,
                                 const VectorType& vector,
                                 const ProblemType& problem,
                                 const std::string type,
                                 const ParametersMapType parameters = ParametersMapType())
  {
    if (call_equals< LocalNonconformityOS2014Type >(type))
      return call_estimate< LocalNonconformityOS2014Type >(space, vector, problem, parameters);
    else if (call_equals< LocalResidualOS2014Type >(type))
      return call_estimate< LocalResidualOS2014Type >(space, vector, problem, parameters);
    else if (call_equals< LocalDiffusiveFluxOS2014Type >(type))
      return call_estimate< LocalDiffusiveFluxOS2014Type >(space, vector, problem, parameters);
    else if (call_equals< LocalDiffusiveFluxOS2014StarType >(type))
      return call_estimate< LocalDiffusiveFluxOS2014StarType >(space, vector, problem, parameters);
    else if (call_equals< OS2014Type >(type))
      return call_estimate< OS2014Type >(space, vector, problem, parameters);
    else
      DUNE_THROW(Stuff::Exceptions::you_are_using_this_wrong,
                 "Requested type '" << type << "' is not one of available()!");
  } // ... estimate(...)
}; // class BlockSWIPDGEstimator


//#if HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_EIGEN && HAVE_DUNE_GRID_MULTISCALE


//extern template class BlockSWIPDGEstimator<
//    GDT::Spaces::Block< GDT::Spaces::DiscontinuousLagrange::FemBased<
//        grid::Part::Local::IndexBased::Const< grid::Part::Leaf::Const< ALUGrid< 2, 2, simplex, conforming > > >,
//        1, double, 1, 1 > >,
//    Dune::Stuff::LA::EigenDenseVector< double >,
//    ProblemInterface< typename ALUGrid< 2, 2, simplex, conforming >::template Codim< 0 >::Entity, double, 2, double, 1 >,
//    ALUGrid< 2, 2, simplex, conforming > >;


//#endif // HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_EIGEN && HAVE_DUNE_GRID_MULTISCALE


} // namespace Discretizations
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_BLOCK_SWIPDG_ESTIMATOR_HH
