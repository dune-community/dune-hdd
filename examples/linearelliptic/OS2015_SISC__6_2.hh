// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_EXAMPLES_LINEARELLIPTIC_OS2014_HH
#define DUNE_HDD_EXAMPLES_LINEARELLIPTIC_OS2014_HH

#define DUNE_HDD_LINEARELLIPTIC_TESTCASES_BASE_DISABLE_WARNING

#include <boost/numeric/conversion/cast.hpp>

#include <dune/fem/misc/mpimanager.hh>

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/functions/expression.hh>

#include <dune/pymor/parameters/base.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/operators/prolongations.hh>
#include <dune/gdt/playground/products/elliptic-swipdg.hh>
#include <dune/gdt/playground/spaces/finitevolume/default.hh>

#include <dune/hdd/linearelliptic/testcases/OS2014.hh>
#include <dune/hdd/linearelliptic/testcases/spe10.hh>
#include <dune/hdd/linearelliptic/discretizations/block-swipdg.hh>
#include <dune/hdd/linearelliptic/estimators/block-swipdg.hh>

namespace internal {


class Initializer
{
public:
  Initializer(const DUNE_STUFF_SSIZE_T info_log_levels,
              const DUNE_STUFF_SSIZE_T debug_log_levels,
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
    DSC::TimedLogger().get("OS2014.initializer").info() << "creating grid and problem... " << std::endl;
  }
}; // class Initializer


template< class TestCaseType >
class Example
  : Initializer
{
public:
  typedef typename TestCaseType::GridType GridType;
  typedef double RangeFieldType;
  typedef Dune::HDD::LinearElliptic::Discretizations::BlockSWIPDG< GridType, RangeFieldType, 1 > DiscretizationType;
  typedef typename DiscretizationType::VectorType VectorType;
  typedef typename TestCaseType::ParametersMapType ParametersMapType;

private:
  typedef Dune::HDD::LinearElliptic::Estimators::BlockSWIPDG< typename DiscretizationType::AnsatzSpaceType,
                                                              VectorType,
                                                              typename DiscretizationType::ProblemType,
                                                              GridType > Estimator;

public:
  Example(const ParametersMapType& parameter_range,
          const std::string partitioning = "[1 1 1]",
          const DUNE_STUFF_SSIZE_T num_refinements = 0,
          const DUNE_STUFF_SSIZE_T oversampling_layers = 0,
          const std::vector< std::string > products = {},
          const bool with_reference = true,
          const DUNE_STUFF_SSIZE_T info_log_levels  = 0,
          const DUNE_STUFF_SSIZE_T debug_log_levels = -1,
          const bool enable_warnings = true,
          const bool enable_colors   = true,
          const std::string info_color  = DSC::TimedLogging::default_info_color(),
          const std::string debug_color = DSC::TimedLogging::default_debug_color(),
          const std::string warn_color  = DSC::TimedLogging::default_warning_color())
    : Initializer(info_log_levels,
                  debug_log_levels,
                  enable_warnings,
                  enable_colors,
                  info_color,
                  debug_color,
                  warn_color)
    , with_reference_(with_reference),
      parameter_range_(parameter_range)
    , test_case_(merge_parameters({{"mu",     Dune::Pymor::Parameter("mu", 1)},  // <- it does not matter which parameters we give to the
                                   {"mu_hat", Dune::Pymor::Parameter("mu", 1)},  //    test case here, since we use test_case_.problem()
                                   {"mu_bar", Dune::Pymor::Parameter("mu", 1)}}, //    (which is the parametric problem) anyway
                                  parameter_range_),
                 partitioning,
                 boost::numeric_cast< size_t >(num_refinements),
                 boost::numeric_cast< size_t >(oversampling_layers))
    , reference_test_case_(with_reference_
                           ? new TestCaseType(merge_parameters({{"mu",     Dune::Pymor::Parameter("mu", 1)},
                                                                {"mu_hat", Dune::Pymor::Parameter("mu", 1)},
                                                                {"mu_bar", Dune::Pymor::Parameter("mu", 1)}},
                                                               parameter_range_),
                                              partitioning,
                                              boost::numeric_cast< size_t >(num_refinements + 1))
                           : nullptr)
    , discretization_(*test_case_.reference_provider(),
                      test_case_.boundary_info(),
                      test_case_.problem(),
                      products)
    , reference_discretization_(with_reference_
                                ? new DiscretizationType(*reference_test_case_->reference_provider(),
                                                         reference_test_case_->boundary_info(),
                                                         reference_test_case_->problem(),
                                                         products)
                                : nullptr)
  {
    auto logger = DSC::TimedLogger().get("OS2014.example");
    logger.info() << "initializing discretization... " << std::flush;
    discretization_.init();
    if (with_reference_)
      reference_discretization_->init();
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

  void visualize(const std::string& filename_prefix) const
  {
    test_case_.reference_provider()->visualize(filename_prefix + ".grid", /*coupling=*/ false);
    test_case_.problem().visualize(test_case_.reference_provider()->leaf_view(),
                                   filename_prefix + ".problem",
                                   /*subsampling=*/ false);
  } // ... visualize(...)

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
    if (!with_reference_)
      DUNE_THROW(Dune::Stuff::Exceptions::you_are_using_this_wrong,
                 "Do not call compute error() if with_reference is false!");
    const std::string type = (product_type.empty() && reference_discretization_->available_products().size() == 1)
                             ? reference_discretization_->available_products()[0]
                             : product_type;
    // wrap solution into a discrete function
    using namespace Dune;
    typedef typename DiscretizationType::AnsatzSpaceType AnsatzSpaceType;
    typedef GDT::ConstDiscreteFunction< AnsatzSpaceType, VectorType > ConstDiscreteFunctionType;
    ConstDiscreteFunctionType coarse_solution(*discretization_.ansatz_space(), solution);
    // prolong to reference grid view
    typedef GDT::DiscreteFunction< AnsatzSpaceType, VectorType > DiscreteFunctionType;
    DiscreteFunctionType fine_solution(*reference_discretization_->ansatz_space());
    GDT::Operators::Prolongation< typename DiscretizationType::GridViewType >
        prolongation_operator(reference_discretization_->grid_view());
    prolongation_operator.apply(coarse_solution, fine_solution);
    // compute reference solution
    DiscreteFunctionType reference_solution(*reference_discretization_->ansatz_space());
    reference_discretization_->solve(reference_solution.vector(), mu);
    // compute error
    const auto difference = reference_solution.vector() - fine_solution.vector();
    const auto product = reference_discretization_->get_product(type);
    return std::sqrt(product.apply2(difference, difference, mu_product));
  } // ... compute_error(...)

  VectorType* pb_project_global_to_oversampled(const VectorType& global_vector, const DUNE_STUFF_SSIZE_T subdomain) const
  {
    using namespace Dune;
    size_t ss = std::numeric_limits< size_t >::max();
    try {
      ss = boost::numeric_cast< size_t >(subdomain);
    } catch (boost::bad_numeric_cast& ee) {
      DUNE_THROW(Stuff::Exceptions::wrong_input_given,
                 "There was an error in boost converting " << subdomain << " to "
                 << Stuff::Common::Typename< size_t >::value() << ": \n\n" << ee.what());
    }
    const GDT::ConstDiscreteFunction< typename DiscretizationType::AnsatzSpaceType, VectorType >
        global_function(*discretization_.ansatz_space(), global_vector);
    const auto oversampled_discretization = discretization_.get_oversampled_discretization(subdomain, "dirichlet");
    GDT::DiscreteFunction< typename DiscretizationType::OversampledDiscretizationType::AnsatzSpaceType, VectorType >
        oversampled_function(*oversampled_discretization.ansatz_space());
    const GDT::Operators::Projection< typename DiscretizationType::OversampledDiscretizationType::GridViewType >
        projection_operator(oversampled_discretization.grid_view());
    projection_operator.apply(global_function, oversampled_function);
    return new VectorType(oversampled_function.vector());
  } // ... pb_project_global_to_oversampled(...)

  VectorType* pb_project_global_to_local(const VectorType& global_vector, const DUNE_STUFF_SSIZE_T subdomain) const
  {
    using namespace Dune;
    size_t ss = std::numeric_limits< size_t >::max();
    try {
      ss = boost::numeric_cast< size_t >(subdomain);
    } catch (boost::bad_numeric_cast& ee) {
      DUNE_THROW(Stuff::Exceptions::wrong_input_given,
                 "There was an error in boost converting " << subdomain << " to "
                 << Stuff::Common::Typename< size_t >::value() << ": \n\n" << ee.what());
    }
    const GDT::ConstDiscreteFunction< typename DiscretizationType::AnsatzSpaceType, VectorType >
        global_function(*discretization_.ansatz_space(), global_vector);
    const auto local_discretization = discretization_.get_local_discretization(subdomain);
    GDT::DiscreteFunction< typename DiscretizationType::LocalDiscretizationType::AnsatzSpaceType, VectorType >
        local_function(*local_discretization.ansatz_space());
    const GDT::Operators::Projection< typename DiscretizationType::LocalDiscretizationType::GridViewType >
        projection_operator(local_discretization.grid_view());
    projection_operator.apply(global_function, local_function);
    return new VectorType(local_function.vector());
  } // ... pb_project_global_to_local(...)

  VectorType* pb_project_oversampled_to_local(const VectorType& oversampled_vector, const DUNE_STUFF_SSIZE_T subdomain) const
  {
    using namespace Dune;
    size_t ss = std::numeric_limits< size_t >::max();
    try {
      ss = boost::numeric_cast< size_t >(subdomain);
    } catch (boost::bad_numeric_cast& ee) {
      DUNE_THROW(Stuff::Exceptions::wrong_input_given,
                 "There was an error in boost converting " << subdomain << " to "
                 << Stuff::Common::Typename< size_t >::value() << ": \n\n" << ee.what());
    }
    const auto oversampled_discretization = discretization_.get_oversampled_discretization(subdomain, "dirichlet");
    const GDT::ConstDiscreteFunction
        < typename DiscretizationType::OversampledDiscretizationType::AnsatzSpaceType, VectorType >
        oversampled_function(*oversampled_discretization.ansatz_space(), oversampled_vector);
    const auto local_discretization = discretization_.get_local_discretization(subdomain);
    const auto& local_space = *local_discretization.ansatz_space();
    GDT::DiscreteFunction< typename DiscretizationType::LocalDiscretizationType::AnsatzSpaceType, VectorType >
        local_function(local_space);
    const GDT::Operators::Projection< typename DiscretizationType::LocalDiscretizationType::GridViewType >
        projection_operator(local_space.grid_view());
    projection_operator.apply(oversampled_function, local_function);
    return new VectorType(local_function.vector());
  } // ... pb_project_oversampled_to_local(...)

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
                               merge_parameters({{"mu_hat", mu_hat},
                                                 {"mu_bar", mu_bar},
                                                 {"mu",     mu}},
                                                parameter_range_));
  } // ... estimate(...)

  std::vector< std::string > available_local_estimators() const
  {
    return Estimator::available_local();
  }

  std::vector< RangeFieldType > estimate_local(const VectorType& vector,
                                               const std::string type,
                                               const Dune::Pymor::Parameter mu_hat = Dune::Pymor::Parameter(),
                                               const Dune::Pymor::Parameter mu_bar = Dune::Pymor::Parameter(),
                                               const Dune::Pymor::Parameter mu     = Dune::Pymor::Parameter())
  {
    return Estimator::estimate_local(*discretization_.ansatz_space(),
                                     vector,
                                     discretization_.problem(),
                                     type,
                                     merge_parameters({{"mu_hat", mu_hat},
                                                       {"mu_bar", mu_bar},
                                                       {"mu",     mu}},
                                                      parameter_range_));
  } // ... estimate_local(...)

  VectorType solve_for_local_correction(const std::vector< VectorType >& local_vectors,
                                        const DUNE_STUFF_SSIZE_T subdomain,
                                        const Dune::Pymor::Parameter mu = Dune::Pymor::Parameter()) const
  {
    const size_t ss = boost::numeric_cast< size_t >(subdomain);
    return discretization_.solve_for_local_correction(local_vectors, ss, mu);
  }

  VectorType solve_oversampled(const DUNE_STUFF_SSIZE_T subdomain,
                               const std::string& boundary_value_type,
                               const VectorType& boundary_values,
                               const Dune::Pymor::Parameter mu = Dune::Pymor::Parameter()) const
  {
    using namespace Dune;
    DSC::Configuration boundary_cfg;
    if (boundary_value_type == "dirichlet")
      boundary_cfg = Stuff::Grid::BoundaryInfoConfigs::AllDirichlet::default_config();
    else if (boundary_value_type == "neumann")
      boundary_cfg = Stuff::Grid::BoundaryInfoConfigs::AllNeumann::default_config();
    else
      DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Unknown boundary_value_type given: " << boundary_value_type);
    // we need this one temporarily for the space
    const auto tmp_oversampled_discretization = discretization_.get_oversampled_discretization(boost::numeric_cast< size_t >(subdomain),
                                                                                               boundary_value_type);
    const auto oversampled_space = tmp_oversampled_discretization.ansatz_space();
    typedef typename DiscretizationType::OversampledDiscretizationType OversampledDiscretizationType;
    typedef typename HDD::LinearElliptic::Problems::ConvertToDefault< typename TestCaseType::ProblemType >::Type
                                                                       OversampledProblemType;
    typedef typename OversampledDiscretizationType::AnsatzSpaceType    OversampledAnsatzSpaceType;
    auto bv_function = std::make_shared< GDT::ConstDiscreteFunction< OversampledAnsatzSpaceType, VectorType > >(
          *oversampled_space, boundary_values, "boundary_values");
    auto nonparametric_problem = test_case_.problem().with_mu(mu);
    const OversampledProblemType oversampled_problem(nonparametric_problem->diffusion_factor()->affine_part(),
                                                     nonparametric_problem->diffusion_tensor()->affine_part(),
                                                     nonparametric_problem->force()->affine_part(),
                                                     bv_function,
                                                     bv_function);
    OversampledDiscretizationType discretization(*test_case_.reference_provider(),
                                                 boundary_cfg,
                                                 oversampled_problem,
                                                 boost::numeric_cast< int >(subdomain));
    discretization.init();
    VectorType solution = discretization.create_vector();
    discretization.solve(solution);
    return solution;
  } // ... solve_oversampled(...)

  void visualize_on_coarse_grid(const std::vector< double >& vector,
                                const std::string& filename,
                                const std::string& name)
  {
    using namespace Dune;
    const auto ms_grid = discretization_.ansatz_space()->ms_grid();
    const auto& grid_view = discretization_.grid_view();
    if (vector.size() != ms_grid->size())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given,
                 "Given vector has wrong lenght!\n"
                 << "  expected: " << ms_grid->size() << "\n"
                 << "  actual:   " << vector.size() << "\n");
    if (filename.empty())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Given filename must not be empty!");
    if (name.empty())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Given name must not be empty!");
    GDT::Spaces::FiniteVolume::Default< typename std::remove_reference< decltype(grid_view) >::type, RangeFieldType, 1 >
        fv_space(grid_view);
    auto visualization = GDT::make_discrete_function< VectorType >(fv_space, name);
    for (const auto& entity : Stuff::Common::entityRange(grid_view))
      visualization.vector().set_entry(grid_view.indexSet().index(entity), vector[ms_grid->subdomainOf(entity)]);
    visualization.visualize(filename);
  } // ... visualize_on_coarse_grid(...)

  RangeFieldType alpha(const Dune::Pymor::Parameter& mu_1, const Dune::Pymor::Parameter& mu_2)
  {
    return test_case_.problem().diffusion_factor()->alpha(mu_1, mu_2);
  }

  RangeFieldType gamma(const Dune::Pymor::Parameter& mu_1, const Dune::Pymor::Parameter& mu_2)
  {
    return test_case_.problem().diffusion_factor()->gamma(mu_1, mu_2);
  }

private:
  static ParametersMapType merge_parameters(const ParametersMapType& first,
                                            const ParametersMapType& second)
  {
    ParametersMapType ret = first;
    for (const auto& element : second)
      ret.insert(element);
    return ret;
  }

  const bool with_reference_;
  const ParametersMapType parameter_range_;
  TestCaseType test_case_;
  std::unique_ptr< TestCaseType > reference_test_case_;
  DiscretizationType discretization_;
  std::unique_ptr< DiscretizationType > reference_discretization_;
}; // class Example


} // namespace internal


template< class GridImp >
class Spe10Model1Example
  : public internal::Example< Dune::HDD::LinearElliptic::TestCases::Spe10::ParametricBlockModel1< GridImp > >
{
  static_assert(GridImp::dimension == 2, "Only available in 2d!");
  typedef internal::Example< Dune::HDD::LinearElliptic::TestCases::Spe10::ParametricBlockModel1< GridImp > > BaseType;

public:
  template< class... Args >
  Spe10Model1Example(Args&& ...args)
    : BaseType({{"parameter_range_min", Dune::Pymor::Parameter("mu", 0.1)},
                {"parameter_range_max", Dune::Pymor::Parameter("mu", 1.0)}},
               std::forward< Args >(args)...)
  {}
}; // class Spe10Model1Example


template< class GridImp >
class OS2014Example
  : public internal::Example< Dune::HDD::LinearElliptic::TestCases::OS2014::ParametricBlockConvergence< GridImp > >
{
  static_assert(GridImp::dimension == 2, "Only available in 2d!");
  typedef internal::Example
      < Dune::HDD::LinearElliptic::TestCases::OS2014::ParametricBlockConvergence< GridImp > > BaseType;

public:
  template< class... Args >
  OS2014Example(Args&& ...args)
    : BaseType({{"parameter_range_min", Dune::Pymor::Parameter("mu", 0.1)},
                {"parameter_range_max", Dune::Pymor::Parameter("mu", 1.0)}},
               std::forward< Args >(args)...)
  {}
}; // class OS2014Example


template< class GridImp >
class FiveSpotExample
  : public internal::Example< Dune::HDD::LinearElliptic::TestCases::OS2014::FiveSpotBlock< GridImp > >
{
  static_assert(GridImp::dimension == 2, "Only available in 2d!");
  typedef internal::Example
      < Dune::HDD::LinearElliptic::TestCases::OS2014::FiveSpotBlock< GridImp > > BaseType;

public:
  template< class... Args >
  FiveSpotExample(Args&& ...args)
    : BaseType({{"parameter_range_min", Dune::Pymor::Parameter("mu", 0.1)},
                {"parameter_range_max", Dune::Pymor::Parameter("mu", 0.9)}},
               std::forward< Args >(args)...)
  {}
}; // class FiveSpotExample


template< class GridImp >
class LocalThermalblockExample
  : public internal::Example< Dune::HDD::LinearElliptic::TestCases::OS2014::LocalThermalblockBlock< GridImp > >
{
  static_assert(GridImp::dimension == 2, "Only available in 2d!");
  typedef internal::Example
      < Dune::HDD::LinearElliptic::TestCases::OS2014::LocalThermalblockBlock< GridImp > > BaseType;

public:
  template< class... Args >
  LocalThermalblockExample(Args&& ...args)
    : BaseType({{"parameter_range_min", Dune::Pymor::Parameter("mu", 0.1)},
                {"parameter_range_max", Dune::Pymor::Parameter("mu", 1.0)}},
               std::forward< Args >(args)...)
  {}
}; // class LocalThermalblockExample


#endif // DUNE_HDD_EXAMPLES_LINEARELLIPTIC_OS2014_HH
