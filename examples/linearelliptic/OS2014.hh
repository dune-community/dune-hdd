// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_EXAMPLES_LINEARELLIPTIC_OS2014_HH
#define DUNE_HDD_EXAMPLES_LINEARELLIPTIC_OS2014_HH

#define DUNE_HDD_LINEARELLIPTIC_TESTCASES_BASE_DISABLE_WARNING

#include <boost/numeric/conversion/cast.hpp>

#include <dune/fem/misc/mpimanager.hh>

#include <dune/stuff/common/logging.hh>
#include <dune/stuff/functions/expression.hh>

#include <dune/pymor/parameters/base.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/operators/prolongations.hh>
#include <dune/gdt/playground/products/elliptic-swipdg.hh>

#include <dune/hdd/linearelliptic/testcases/spe10.hh>
#include <dune/hdd/linearelliptic/discretizations/block-swipdg.hh>
#include <dune/hdd/linearelliptic/estimators/block-swipdg.hh>

namespace internal {


class Initializer
{
public:
  Initializer(const ssize_t info_log_levels,
              const ssize_t debug_log_levels,
              const bool enable_warnings,
              const bool enable_colors,
              const std::string info_color,
              const std::string debug_color,
              const std::string warn_color)
  {
    try {
      int argc = 0;
      char** argv = new char* [0];
      Dune::Fem::MPIManager::initialize(argc, argv);
    } catch (...) {}
    DSC::TimedLogger().create(info_log_levels,
                              debug_log_levels,
                              enable_warnings,
                              enable_colors,
                              info_color,
                              debug_color,
                              warn_color);
    DSC::TimedLogger().get("OS2014.spe10model1example").info() << "creating grid and problem... " << std::endl;
  }
}; // class Initializer


} // namespace internal


template< class GridType >
class OS2014Spe10Model1Example
  : internal::Initializer
{
  static_assert(GridType::dimension == 2, "Only available in 2d!");
public:
  typedef Dune::HDD::LinearElliptic::TestCases::Spe10::ParametricBlockModel1< GridType > TestCaseType;
  typedef double RangeFieldType;
  typedef Dune::HDD::LinearElliptic::Discretizations::BlockSWIPDG< GridType, RangeFieldType, 1 > DiscretizationType;
  typedef typename DiscretizationType::VectorType VectorType;

private:
  typedef Dune::HDD::LinearElliptic::Estimators::BlockSWIPDG< typename DiscretizationType::AnsatzSpaceType,
                                                              VectorType,
                                                              typename DiscretizationType::ProblemType,
                                                              GridType > Estimator;

public:
  OS2014Spe10Model1Example(const std::string partitioning = "[1 1 1]",
                           const DUNE_STUFF_SSIZE_T num_refinements = 0,
                           const std::vector< std::string > products = {},
                           const ssize_t info_log_levels  = 0,
                           const ssize_t debug_log_levels = -1,
                           const bool enable_warnings = true,
                           const bool enable_colors   = true,
                           const std::string info_color  = DSC::TimedLogging::default_info_color(),
                           const std::string debug_color = DSC::TimedLogging::default_debug_color(),
                           const std::string warn_color  = DSC::TimedLogging::default_warning_color())
    : internal::Initializer(info_log_levels,
                            debug_log_levels,
                            enable_warnings,
                            enable_colors,
                            info_color,
                            debug_color,
                            warn_color)
    , test_case_({{"mu", Dune::Pymor::Parameter("mu", 1)},     // <- it does not matter which parameters we give to the
                  {"mu_hat", Dune::Pymor::Parameter("mu", 1)}, //    test case here, since we use test_case_.problem()
                  {"mu_bar", Dune::Pymor::Parameter("mu", 1)}, //    (which is the parametric problem) anyway
                  {"mu_minimizing", Dune::Pymor::Parameter("mu", 1)}},
                 partitioning,
                 boost::numeric_cast< size_t >(num_refinements))
    , reference_test_case_({{"mu", Dune::Pymor::Parameter("mu", 1)},
                            {"mu_hat", Dune::Pymor::Parameter("mu", 1)},
                            {"mu_bar", Dune::Pymor::Parameter("mu", 1)},
                            {"mu_minimizing", Dune::Pymor::Parameter("mu", 1)}},
                           partitioning,
                           boost::numeric_cast< size_t >(num_refinements + 1))
    , discretization_(*test_case_./*reference*/level_provider(0),
                      test_case_.boundary_info(),
                      test_case_.problem(),
                      products)
    , reference_discretization_(*reference_test_case_.reference_provider(),
                                reference_test_case_.boundary_info(),
                                reference_test_case_.problem(),
                                products)
  {
    auto logger = DSC::TimedLogger().get("OS2014.spe10model1example");
    logger.info() << "initializing discretization... " << std::flush;
    discretization_.init();
    reference_discretization_.init();
    logger.info() << "done (grid has " << discretization_.grid_view().indexSet().size(0)
                  << " elements, discretization has " << discretization_.ansatz_space()->mapper().size() << " DoFs)"
                  << std::endl;
  } // ... OS2014Spe10Model1Example(...)

  const TestCaseType& test_case() const
  {
    return test_case_;
  }

  DiscretizationType& discretization()
  {
    return discretization_;
  }

  DiscretizationType* discretization_and_return_ptr() const
  {
    return new DiscretizationType(discretization_);
  }

  VectorType project(const std::string expression) const
  {
    using namespace Dune;
    typedef typename DiscretizationType::AnsatzSpaceType AnsatzSpaceType;
    typedef GDT::DiscreteFunction< AnsatzSpaceType, VectorType > DiscreteFunctionType;
    DiscreteFunctionType discrete_function(*discretization_.ansatz_space());

    typedef Stuff::Functions::Expression< typename DiscreteFunctionType::EntityType,
                                          typename DiscreteFunctionType::DomainFieldType,
                                          DiscreteFunctionType::dimDomain,
                                          typename DiscreteFunctionType::RangeFieldType,
                                          DiscreteFunctionType::dimRange > ExpressionFunctionType;
    ExpressionFunctionType func("x", expression);

    GDT::Operators::Projection< typename DiscretizationType::GridViewType > projection(discretization_.grid_view());
    projection.apply(func, discrete_function);

    return discrete_function.vector();
  } // ... project(...)

  RangeFieldType compute_error(const VectorType& solution,
                               const std::string product_type = "",
                               const Dune::Pymor::Parameter mu = Dune::Pymor::Parameter(),
                               const Dune::Pymor::Parameter mu_product = Dune::Pymor::Parameter())
  {
    const std::string type = (product_type.empty() && reference_discretization_.available_products().size() == 1)
                             ? reference_discretization_.available_products()[0]
                             : product_type;
    // wrap solution into a discrete function
    using namespace Dune;
    typedef typename DiscretizationType::AnsatzSpaceType AnsatzSpaceType;
    typedef GDT::ConstDiscreteFunction< AnsatzSpaceType, VectorType > ConstDiscreteFunctionType;
    ConstDiscreteFunctionType coarse_solution(*discretization_.ansatz_space(), solution);
    // prolong to reference grid view
    typedef GDT::DiscreteFunction< AnsatzSpaceType, VectorType > DiscreteFunctionType;
    DiscreteFunctionType fine_solution(*reference_discretization_.ansatz_space());
    GDT::Operators::Prolongation< typename DiscretizationType::GridViewType >
        prolongation_operator(reference_discretization_.grid_view());
    prolongation_operator.apply(coarse_solution, fine_solution);
    // compute reference solution
    DiscreteFunctionType reference_solution(*reference_discretization_.ansatz_space());
    reference_discretization_.solve(reference_solution.vector(), mu);
    // compute error
    const auto product = reference_discretization_.get_product(type).freeze_parameter(mu_product);
    const auto difference = reference_solution.vector() - fine_solution.vector();
    return std::sqrt(product.apply2(difference, difference));
  } // ... compute_error(...)

  RangeFieldType compute_jump_norm(const VectorType& solution_vector,
                                   const Dune::Pymor::Parameter mu_product = Dune::Pymor::Parameter())
  {
    using namespace Dune;
    typedef typename DiscretizationType::AnsatzSpaceType AnsatzSpaceType;
    typedef GDT::ConstDiscreteFunction< AnsatzSpaceType, VectorType > ConstDiscreteFunctionType;
    ConstDiscreteFunctionType solution(*discretization_.ansatz_space(), solution_vector);
    const auto diffusion_factor = discretization_.problem().diffusion_factor()->with_mu(mu_product);
    const auto diffusion_tensor = discretization_.problem().diffusion_tensor()->with_mu(mu_product);
    typedef GDT::Products::EllipticSWIPDGPenaltyLocalizable
        < typename DiscretizationType::GridViewType,
          typename DiscretizationType::ProblemType::DiffusionFactorType::NonparametricType,
          ConstDiscreteFunctionType,
          ConstDiscreteFunctionType,
          RangeFieldType,
          typename DiscretizationType::ProblemType::DiffusionTensorType::NonparametricType > ProductType;
    ProductType penalty_product(discretization_.grid_view(),
                                solution,
                                solution,
                                *diffusion_factor,
                                *diffusion_tensor,
                                2);
    return std::sqrt(penalty_product.apply2());
  } // ... compute_jump_norm(...)

  std::vector< std::string > available_estimators() const
  {
    return Estimator::available();
  }

  RangeFieldType estimate(const VectorType& vector,
                          const std::string type,
                          const Dune::Pymor::Parameter mu_hat = Dune::Pymor::Parameter(),
                          const Dune::Pymor::Parameter mu_bar = Dune::Pymor::Parameter(),
                          const Dune::Pymor::Parameter mu     = Dune::Pymor::Parameter())
  {
    return Estimator::estimate(*discretization_.ansatz_space(),
                               vector,
                               discretization_.problem(),
                               type,
                               {{"mu_hat",        mu_hat},
                                {"mu_bar",        mu_bar},
                                {"mu",            mu},
                                {"mu_minimizing", Dune::Pymor::Parameter("mu", 0.1)}});
  } // ... estimate(...)

public:
  TestCaseType test_case_;
  TestCaseType reference_test_case_;
  DiscretizationType discretization_;
  DiscretizationType reference_discretization_;
}; // class OS2014Spe10Model1Example


#endif // DUNE_HDD_EXAMPLES_LINEARELLIPTIC_OS2014_HH
