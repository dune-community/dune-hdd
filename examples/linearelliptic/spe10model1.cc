// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include <dune/stuff/common/disable_warnings.hh>
# if HAVE_ALUGRID
#   include <dune/grid/alugrid.hh>
# else
#   error This example requires alugrid!
# endif

# include <dune/common/timer.hh>
# include <dune/common/fvector.hh>

# if HAVE_DUNE_FEM
#   include <dune/fem/misc/mpimanager.hh>
# else
#   include <dune/common/parallel/mpihelper.hh>
# endif

# include <dune/grid/io/file/dgfparser.hh>
# include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/grid/multiscale/provider/cube.hh>

#include <dune/stuff/functions/spe10.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/common/float_cmp.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/common/string.hh>

#include <dune/gdt/operators/prolongations.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/playground/products/elliptic.hh>
#include <dune/gdt/operators/fluxreconstruction.hh>

#include <dune/hdd/linearelliptic/problems/default.hh>
#include <dune/hdd/playground/linearelliptic/discretizations/block-swipdg.hh>

using namespace Dune;
using namespace HDD;
using namespace GDT;


template< class E, class D, int d, class R >
class IndicatorFunction
  : public Stuff::LocalizableFunctionInterface< E, D, d, R, 1 >
{
  typedef Stuff::LocalizableFunctionInterface< E, D, d, R, 1 > BaseType;
  typedef IndicatorFunction< E, D, d, R > ThisType;

  class Localfunction
    : public Stuff::LocalfunctionInterface< E, D, d, R, 1 >
  {
    typedef Stuff::LocalfunctionInterface< E, D, d, R, 1 > InterfaceType;
  public:
    using typename InterfaceType::EntityType;
    using typename InterfaceType::DomainType;
    using typename InterfaceType::RangeType;
    using typename InterfaceType::JacobianRangeType;

    Localfunction(const EntityType& entity, const RangeType& value)
      : InterfaceType(entity)
      , value_(value)
    {}

    virtual size_t order() const DS_OVERRIDE DS_FINAL
    {
      return 0;
    }

    virtual void evaluate(const DomainType& xx, RangeType& ret) const DS_OVERRIDE DS_FINAL
    {
      assert(this->is_a_valid_point(xx));
      ret = value_;
    }

    virtual void jacobian(const DomainType& xx, JacobianRangeType& ret) const DS_OVERRIDE DS_FINAL
    {
      assert(this->is_a_valid_point(xx));
      ret *= 0.0;
    }

  private:
    const RangeType value_;
  }; // class Localfunction

public:
  using typename BaseType::EntityType;
  typedef Stuff::Common::FieldVector< D, d > DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::LocalfunctionType;

  IndicatorFunction(std::vector< std::pair< std::pair< DomainType, DomainType >, R > > values,
                    const std::string name = "indicator")
    : values_(values)
    , name_(name)
  {}

  virtual ~IndicatorFunction() {}

  virtual ThisType* copy() const DS_OVERRIDE DS_FINAL
  {
    DUNE_THROW(NotImplemented, "");
  }

  virtual std::string name() const DS_OVERRIDE DS_FINAL
  {
    return name_;
  }

  virtual std::unique_ptr< LocalfunctionType > local_function(const EntityType& entity) const DS_OVERRIDE DS_FINAL
  {
    const auto center = entity.geometry().center();
    for (const auto& element : values_) {
      const auto& domain = element.first;
      const auto& lower_left = domain.first;
      const auto& upper_right = domain.second;
      if (Stuff::Common::FloatCmp::lt(lower_left, center) && Stuff::Common::FloatCmp::lt(center, upper_right)) {
        const auto& value = element.second;
        return Stuff::Common::make_unique< Localfunction >(entity, value);
      }
    }
    return Stuff::Common::make_unique< Localfunction >(entity, 0.0);
  } // ... local_function(...)

private:
  const std::vector< std::pair< std::pair< DomainType, DomainType >, R > > values_;
  const std::string name_;
}; // class IndicatorFunction


int main(int argc, char** argv)
{
  try {
#if HAVE_DUNE_FEM
    Fem::MPIManager::initialize(argc, argv);
#else
    MPIHelper::instance(argc, argv);
#endif

    typedef ALUGrid< 2, 2, simplex, conforming > GridType;
    typedef GridType::Codim< 0 >::Entity         EntityType;
    typedef GridType::ctype   DomainFieldType;
    static const unsigned int dimDomain = GridType::dimension;
    typedef Dune::FieldVector< DomainFieldType, dimDomain > DomainType;
    typedef double            RangeFieldType;
    static const unsigned int dimRange = 1;

    typedef Stuff::Functions::Spe10Model1
        < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimDomain, dimDomain > Spe10Model1FunctionType;
    Stuff::Common::Configuration config = Spe10Model1FunctionType::default_config();
    const std::string x_elements = "5";
    const std::string y_elements = "1";
    config["upper_right"] = "[5 1]";
    config["num_elements"] = "[100 20]";
    config["num_partitions"] = "[" + x_elements + " " + y_elements + "]";
    std::shared_ptr< Spe10Model1FunctionType > spe10_model1_function(Spe10Model1FunctionType::create(config));
    typedef Stuff::Functions::Constant< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >
        ConstantFunctionType;
    auto one = std::make_shared< ConstantFunctionType >(1.0, "diffusion_factor");
    auto dirichlet = std::make_shared< ConstantFunctionType >(0.0, "dirichlet");
    auto neumann = std::make_shared< ConstantFunctionType >(0.0, "neumann");
    typedef IndicatorFunction< EntityType, DomainFieldType, dimDomain, RangeFieldType > IndicatorFunctionType;
    const RangeFieldType scale = 1.0;
    auto force = std::shared_ptr< IndicatorFunctionType >(
                   new IndicatorFunctionType({{{{0.55, 0.70}, {0.70, 0.85}},  1.0 * scale},
                                              {{{0.45, 0.15}, {0.60, 0.30}},  1.0 * scale},
                                              {{{3.00, 0.75}, {3.15, 0.90}}, -1.0 * scale},
                                              {{{4.30, 0.50}, {4.45, 0.65}}, -1.0 * scale}}, "force"));
    typedef LinearElliptic::Problems::Default< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >
        ProblemType;
    const ProblemType problem(one, spe10_model1_function, force, dirichlet, neumann);

    typedef Stuff::Grid::Providers::Cube< GridType > GridProvider;
    auto level_grid_provider = GridProvider::create(config);
    level_grid_provider->grid().globalRefine(1);

    typedef grid::Multiscale::Providers::Cube< GridType > MsGridProviderType;
    MsGridProviderType level_provider(level_grid_provider->grid_ptr(),
                                      config.get< DomainType >("lower_left"),
                                      config.get< DomainType >("upper_right"),
                                      config.get< std::vector< size_t >  >("num_partitions", dimDomain));

    typedef Stuff::Grid::Providers::Cube< GridType > GridProvider;
    auto reference_grid_provider = GridProvider::create(config);
    reference_grid_provider->grid().globalRefine(1);
    reference_grid_provider->grid().globalRefine(2 * DGFGridInfo< GridType >::refineStepsForHalf());
    MsGridProviderType reference_level_provider(reference_grid_provider->grid_ptr(),
                                                config.get< DomainType >("lower_left"),
                                                config.get< DomainType >("upper_right"),
                                                config.get< std::vector< size_t >  >("num_partitions", dimDomain));

    problem.visualize(*level_grid_provider->leaf_view(), "problem");

    std::cout << "solving... " << std::flush;
    Dune::Timer timer;
    typedef LinearElliptic::Discretizations::BlockSWIPDG < GridType, RangeFieldType, dimRange, 1 > DiscretizationType;
    DiscretizationType discretization(level_provider,
                                      Stuff::Grid::BoundaryInfoConfigs::AllDirichlet::default_config(),
                                      problem);
    discretization.init();
    auto solution = discretization.create_vector();
    discretization.solve(solution);
    std::cout << "done (took " << timer.elapsed() << "s)" << std::endl;

    std::cout << "computing reference solution... " << std::flush;
    timer.reset();
    DiscretizationType reference_discretization(reference_level_provider,
                                                Stuff::Grid::BoundaryInfoConfigs::AllDirichlet::default_config(),
                                                problem);
    reference_discretization.init();
    auto reference_solution = reference_discretization.create_vector();
    reference_discretization.solve(reference_solution);
    std::cout << "done (took " << timer.elapsed() << "s)" << std::endl;

    std::cout << "prolonging solution... " << std::flush;
    timer.reset();
    typedef typename DiscretizationType::AnsatzSpaceType SpaceType;
    typedef typename SpaceType::GridViewType GridViewType;
    typedef typename DiscretizationType::VectorType VectorType;
    ConstDiscreteFunction< SpaceType, VectorType > discrete_solution_on_level(*discretization.ansatz_space(),
                                                                              solution);
    discrete_solution_on_level.visualize("solution_grid_"
                                         + DSC::toString(discrete_solution_on_level.space().grid_view()->indexSet().size(0))
                                         + "_msgrid_" + x_elements + "x" + y_elements);

    typedef Spaces::RaviartThomas::PdelabBased< typename SpaceType::GridViewType, 0, RangeFieldType, dimDomain >
        RTN0SpaceType;
    const RTN0SpaceType rtn0_space(discretization.grid_view());
    VectorType diffusive_flux_vector(rtn0_space.mapper().size());
    DiscreteFunction< RTN0SpaceType, VectorType > diffusive_flux(rtn0_space, diffusive_flux_vector);
    const Operators::DiffusiveFluxReconstruction< typename SpaceType::GridViewType, ConstantFunctionType, Spe10Model1FunctionType >
      diffusive_flux_reconstruction(*discretization.grid_view(),
                                    *one,
                                    *spe10_model1_function);
    diffusive_flux_reconstruction.apply(discrete_solution_on_level, diffusive_flux);
    diffusive_flux.visualize("diffusive_flux_grid_"
                             + DSC::toString(discrete_solution_on_level.space().grid_view()->indexSet().size(0))
                             + "_msgrid_" + x_elements + "x" + y_elements);

    auto solution_on_reference_level = reference_discretization.create_vector();
    const auto& reference_grid_view = *reference_discretization.ansatz_space()->grid_view();
    DiscreteFunction< SpaceType, VectorType > discrete_solution_on_reference_level(*reference_discretization.ansatz_space(),
                                                                                   solution_on_reference_level);
    Operators::Prolongation< GridViewType >(reference_grid_view).apply(discrete_solution_on_level,
                                                                       discrete_solution_on_reference_level);
    std::cout << "done (took " << timer.elapsed() << "s)" << std::endl;

    ConstDiscreteFunction< SpaceType, VectorType > discrete_reference_solution(*reference_discretization.ansatz_space(),
                                                                               reference_solution);
    auto difference_vector = discrete_solution_on_reference_level.vector() - discrete_reference_solution.vector();
    for (auto& element : difference_vector)
      element = std::abs(element);
    ConstDiscreteFunction< SpaceType, VectorType > difference(*reference_discretization.ansatz_space(),
                                                              difference_vector);
    difference.visualize("difference_grid_"
                         + DSC::toString(discrete_reference_solution.space().grid_view()->indexSet().size(0))
                         + "_msgrid_" + x_elements + "x" + y_elements);

    std::cout << "computing energy error... " << std::flush;
    timer.reset();
    Products::Elliptic< ConstantFunctionType, GridViewType, RangeFieldType, Spe10Model1FunctionType >
        elliptic_product(*one, *spe10_model1_function, reference_grid_view);
    const RangeFieldType energy_error = std::sqrt(elliptic_product.apply2(difference, difference));
    std::cout << "done (is " << energy_error << ", took " << timer.elapsed() << "s)" << std::endl;

    std::cout << "estimating error... " << std::flush;
    timer.reset();
    const RangeFieldType estimate = discretization.estimate(solution);
    std::cout << "done (is " << estimate << ", took " << timer.elapsed() << "s)" << std::endl;

    std::cout << "computing local error/estimator contributions... " << std::flush;
    timer.reset();
    // walk the multiscale grid
    std::vector< RangeFieldType > local_error_contributions(reference_level_provider.ms_grid()->size());
    std::vector< RangeFieldType > local_estimator_contributions(reference_level_provider.ms_grid()->size());
    typedef typename DiscretizationType::LocalDiscretizationType LocalDiscretizationType;
    const auto& local_reference_discretizations = reference_discretization.local_discretizations();
    assert(local_reference_discretizations.size() == reference_level_provider.ms_grid()->size());
    for (size_t ss = 0; ss < reference_level_provider.ms_grid()->size(); ++ss) {
      const auto& local_reference_discretization = *local_reference_discretizations[ss];
      typedef typename LocalDiscretizationType::AnsatzSpaceType LocalAnsatzSpaceType;
      typedef typename LocalAnsatzSpaceType::GridViewType LocalGridViewType;
      const auto local_grid_view = local_reference_discretization.ansatz_space()->grid_view();
      // project level and reference solutions to local reference grid parts
      auto local_level_solution_vector = local_reference_discretization.create_vector();
      DiscreteFunction< LocalAnsatzSpaceType, VectorType >
          local_level_solution(*local_reference_discretization.ansatz_space(), local_level_solution_vector);
      Operators::Prolongation< LocalGridViewType >(*local_grid_view).apply(discrete_solution_on_level,
                                                                           local_level_solution);
      auto local_reference_solution_vector = local_reference_discretization.create_vector();
      DiscreteFunction< LocalAnsatzSpaceType, VectorType >
          local_reference_solution(*local_reference_discretization.ansatz_space(), local_reference_solution_vector);
      Operators::Prolongation< LocalGridViewType >(*local_grid_view).apply(discrete_reference_solution,
                                                                           local_reference_solution);
      // compute energy error
      Products::Elliptic< ConstantFunctionType, LocalGridViewType, RangeFieldType, Spe10Model1FunctionType >
          local_elliptic_product(*one, *spe10_model1_function, *local_grid_view);
      const auto local_difference = local_level_solution - local_reference_solution;
      local_error_contributions[ss] = local_elliptic_product.apply2(local_difference, local_difference);

      // compute estimator
      local_estimator_contributions[ss] = discretization.estimate_local(solution, ss);
    } // walk the multiscale grid
    std::cout << "done (took " << timer.elapsed() << "s):" << std::endl;
    RangeFieldType sum_local_errors = 0.0;
    RangeFieldType sum_local_estimators = 0.0;
    for (size_t ss = 0; ss < reference_level_provider.ms_grid()->size(); ++ss) {
      sum_local_errors += local_error_contributions[ss];
      sum_local_estimators += local_estimator_contributions[ss];
    }
    std::cout << "  energy error: " << std::sqrt(sum_local_errors) << std::endl;
    std::cout << "  estimator: " << std::sqrt(3.0) * std::sqrt(sum_local_estimators) << std::endl;

    std::map< RangeFieldType, size_t > subdomains_by_error;
    std::map< RangeFieldType, size_t > subdomains_by_estimator;
    for (size_t ss = 0; ss < reference_level_provider.ms_grid()->size(); ++ss) {
      subdomains_by_error[local_error_contributions[ss] / sum_local_errors] = ss;
      subdomains_by_estimator[local_estimator_contributions[ss] / sum_local_estimators] = ss;
    }

    if (reference_level_provider.ms_grid()->size() <= 20) {
      std::cout << "local error contributions:" << std::endl;
      for (auto element : subdomains_by_error)
        std::cout << "  " << element.second << ": " << std::floor(100.0 * element.first) << " %" << std::endl;
      std::cout << "local estimator contributions:" << std::endl;
      for (auto element : subdomains_by_estimator)
        std::cout << "  " << element.second << ": " << std::floor(100.0 * element.first) << " %" << std::endl;
    }

    std::cout << "visualizing local errors/estimators... " << std::flush;
    timer.reset();
    const auto& global_level_grid_view = level_provider.ms_grid()->globalGridPart()->gridView();
    std::vector< double > local_error_visualization(global_level_grid_view.indexSet().size(0));
    std::vector< double > local_estimator_visualization(global_level_grid_view.indexSet().size(0));
    const auto entity_it_end = global_level_grid_view.template end< 0 >();
    for (auto entity_it = global_level_grid_view.template begin< 0 >(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      const size_t index = global_level_grid_view.indexSet().index(entity);
      const size_t subdomain = level_provider.ms_grid()->subdomainOf(entity);
      assert(subdomain < level_provider.ms_grid()->size());
      local_error_visualization[index] = local_error_contributions[subdomain] / sum_local_errors;
      local_estimator_visualization[index] = local_estimator_contributions[subdomain] / sum_local_estimators;
    }
    VTKWriter< typename MsGridProviderType::MsGridType::GlobalGridPartType::GridViewType >
        vtk_writer(global_level_grid_view);
    vtk_writer.addCellData(local_error_visualization, "local error contribution");
    vtk_writer.addCellData(local_estimator_visualization, "local estimator contribution");
    vtk_writer.write("local_contributions_msgrid_" + x_elements + "x" + y_elements, VTK::appendedraw);
    std::cout << "done (took " << timer.elapsed() << "s):" << std::endl;

    // if we came that far we can as well be happy about it
    return 0;
  } catch (Dune::Exception& e) {
    std::cerr << "\ndune reported error: " << e.what() << std::endl;
    std::abort();
  } catch (std::exception& e) {
    std::cerr << "\n" << e.what() << std::endl;
    std::abort();
  } catch (...) {
    std::cerr << "\nUnknown exception thrown!" << std::endl;
    std::abort();
  } // try
} // ... main(...)
