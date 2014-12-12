// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

/**
  * This file is intended to serve as a starting point for quick tests.
  */

#include "config.h"

#ifndef DUNE_GDT_LOCALEVALUATION_SWIPDG_DISABLE_WARNINGS
# define DUNE_GDT_LOCALEVALUATION_SWIPDG_DISABLE_WARNINGS
#endif
#ifndef DUNE_HDD_LINEARELLIPTIC_TESTCASES_BASE_DISABLE_WARNING
# define DUNE_HDD_LINEARELLIPTIC_TESTCASES_BASE_DISABLE_WARNING
#endif

#if !HAVE_ALUGRID
# error Alugrid is required!
#endif

#if HAVE_DUNE_FEM
# include <dune/fem/misc/mpimanager.hh>
#endif

#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/common/timedlogging.hh>
#include <dune/stuff/grid/search.hh>

#include <dune/gdt/operators/prolongations.hh>
#include <dune/gdt/playground/operators/fluxreconstruction.hh>
#include <dune/gdt/playground/spaces/finitevolume/default.hh>
#include <dune/gdt/playground/spaces/raviartthomas/pdelab.hh>
#include <dune/gdt/products/elliptic.hh>

#include <dune/hdd/linearelliptic/testcases/spe10.hh>
#include <dune/hdd/linearelliptic/discretizations/block-swipdg.hh>
#include <dune/hdd/linearelliptic/estimators/block-swipdg.hh>

using namespace Dune;
using namespace Dune::Stuff;
using namespace Dune::Pymor;
using namespace Dune::GDT;
using namespace Dune::HDD;


int main(int argc, char** argv)
{
  try {
#if HAVE_DUNE_FEM
    Dune::Fem::MPIManager::initialize(argc, argv);
#endif

    DSC::TimedLogger().create(0, 0);
    auto logger = DSC::TimedLogger().get("main");
    logger.info() << "creating grid and test case ..." << std::endl;

    typedef ALUGrid< 2, 2, simplex, conforming > GridType;
    typedef LinearElliptic::TestCases::Spe10::ParametricBlockModel1< GridType > TestCaseType;
    const TestCaseType test_case({{"mu",     Parameter("mu", 1)},
                                  {"mu_bar", Parameter("mu", 0.1)},
                                  {"mu_hat", Parameter("mu", 0.1)},
                                  {"parameter_range_min", Parameter("mu", 0.1)},
                                  {"parameter_range_max", Parameter("mu", 1.0)}},
                                 "[25 5 1 ]");
    const auto& problem = test_case.problem();
    auto problem_mu_min = problem.with_mu(Parameter("mu", 0.1));
    auto problem_mu_max = problem.with_mu(Parameter("mu", 1.0));
    const size_t over_integrate = 2;

    logger.info() << "discretizing ..." << std::endl;
    auto disc = LinearElliptic::Discretizations::make_block_swipdg< Stuff::LA::ChooseBackend::eigen_sparse >(
                  *test_case.level_provider(0),
                  test_case.boundary_info(),
                  problem);
    disc.init();

    Parameter mu_min("mu", 0.1);
    Parameter mu_max("mu", 1.0);
    logger.info() << "solving for " << mu_min << " ..." << std::endl;
    typedef typename decltype(disc)::VectorType VectorType;
    auto pressure_mu_min = make_discrete_function< VectorType >(*disc.ansatz_space(), "pressure");
    disc.solve(pressure_mu_min.vector(), mu_min);
    logger.info() << "solving for " << mu_max << " ..." << std::endl;
    auto pressure_mu_max = make_discrete_function< VectorType >(*disc.ansatz_space(), "pressure");
    disc.solve(pressure_mu_max.vector(), mu_max);

    logger.info() << "visualizing solutions ..." << std::endl;
    pressure_mu_min.visualize("pressure_mu_min");
    pressure_mu_max.visualize("pressure_mu_max");

    logger.info() << "reconstructing diffusive flux ..." << std::endl;
    Spaces::RaviartThomas::PdelabBased< typename std::remove_reference< decltype(disc) >::type::GridViewType,
                                        0, double, 2 > rt_space(disc.grid_view());
    auto reconstruction_mu_min = GDT::Operators::make_diffusive_flux_reconstruction(disc.grid_view(),
                                                                                    *problem_mu_min->diffusion_factor()->affine_part(),
                                                                                    *problem_mu_min->diffusion_tensor()->affine_part(),
                                                                                    over_integrate);
    auto velocity_mu_min = make_discrete_function< VectorType >(rt_space, "velocity");
    reconstruction_mu_min.apply(pressure_mu_min, velocity_mu_min);

    auto reconstruction_mu_max = GDT::Operators::make_diffusive_flux_reconstruction(disc.grid_view(),
                                                                                    *problem_mu_max->diffusion_factor()->affine_part(),
                                                                                    *problem_mu_max->diffusion_tensor()->affine_part(),
                                                                                    over_integrate);
    auto velocity_mu_max = make_discrete_function< VectorType >(rt_space, "velocity");
    reconstruction_mu_max.apply(pressure_mu_max, velocity_mu_max);

    logger.info() << "visualizing diffusive flux ..." << std::endl;
    velocity_mu_min.visualize("velocity_mu_min");
    velocity_mu_max.visualize("velocity_mu_max");

    logger.info() << "computing local estimator contributions ..." << std::endl;
    typedef LinearElliptic::Estimators::BlockSWIPDG< typename decltype(disc)::AnsatzSpaceType,
                                                     VectorType,
                                                     std::remove_reference< decltype(problem) >::type,
                                                     GridType > Estimator;
    auto eta_Ts = Estimator::estimate_local(*disc.ansatz_space(),
                                            pressure_mu_max.vector(),
                                            problem,
                                            "eta_OS2014_*",
                                            {{"mu",     mu_max},
                                             {"mu_bar", mu_min},
                                             {"mu_hat", mu_min},
                                             {"parameter_range_min", mu_min},
                                             {"parameter_range_max", mu_max}});
    double eta = 0.0;
    for (const auto& eta_T : eta_Ts)
      eta += std::pow(eta_T, 2);
    eta = std::sqrt(eta);
    for (auto& eta_T : eta_Ts)
      eta_T /= eta;

    logger.info() << "computing reference solution ..." << std::endl;
    auto reference_disc = LinearElliptic::Discretizations::make_block_swipdg< Stuff::LA::ChooseBackend::eigen_sparse >(
                  *test_case.reference_provider(),
                  test_case.boundary_info(),
                  problem);
    reference_disc.init();
    auto reference_solution = make_discrete_function< VectorType >(*reference_disc.ansatz_space());
    reference_disc.solve(reference_solution.vector(), Parameter("mu", 1.0));

    logger.info() << "prolonging solution to reference grid ..." << std::endl;
    auto coarse_solution = make_discrete_function< VectorType >(*reference_disc.ansatz_space());
    GDT::Operators::prolong(pressure_mu_max, coarse_solution);

    logger.info() << "computing local error contributions ..." << std::endl;
    auto difference = reference_solution - coarse_solution;
    auto error_product = Products::make_elliptic_localizable(reference_disc.grid_view(),
                                                             difference,
                                                             *problem_mu_min->diffusion_factor()->affine_part(),
                                                             *problem.diffusion_tensor()->affine_part(),
                                                             2);
    error_product.prepare();
    auto error_distribution = disc.create_vector();
    auto entity_search = Stuff::Grid::make_entity_in_level_search(disc.grid_view());
    for (const auto& entity : DSC::entityRange(reference_disc.grid_view())) {
      // search for the father entity in the (coarser) grid view
      const auto points = {entity.geometry().center()};
      const auto father_entity_ptr_ptrs = entity_search(points);
      if (father_entity_ptr_ptrs.size() != 1)
        DUNE_THROW(Stuff::Exceptions::internal_error, father_entity_ptr_ptrs.size());
      const auto& father_entity_ptr_ptr = father_entity_ptr_ptrs[0];
      assert(father_entity_ptr_ptr);
      const auto father_entity_ptr = *father_entity_ptr_ptr;
      const auto& father_entity = *father_entity_ptr;
      assert(disc.grid_view().contains(father_entity));
      // compute local error
      error_distribution.add_to_entry(disc.grid_view().indexSet().index(father_entity),
                                      error_product.compute_locally(entity));
    }
    // sum local errors on the subdomain level
    const auto ms_grid = disc.ansatz_space()->ms_grid();
    auto error_distribution_subdomains = eta_Ts.copy();
    error_distribution_subdomains *= 0.0;
    for (size_t ss = 0; ss < ms_grid->size(); ++ss)
      for (const auto& entity : DSC::entityRange(ms_grid->localGridPart(ss)))
        error_distribution_subdomains[ss] += error_distribution[disc.grid_view().indexSet().index(entity)];
    double error = 0.0;
    for (const auto& element : error_distribution_subdomains)
      error += element;
    error = std::sqrt(error);
    for (auto& element : error_distribution_subdomains)
      element = std::sqrt(element) / error;

    logger.info() << "visualizing local errors and indicators ..." << std::endl;
    Spaces::FiniteVolume::Default< typename std::remove_reference< decltype(disc) >::type::GridViewType,
                                   double, 1 > fv_space(disc.grid_view());
    auto indicator_visualization = make_discrete_function< VectorType >(fv_space, "distribution");
    auto error_visualization = make_discrete_function< VectorType >(fv_space, "distribution");
    for (size_t ss = 0; ss < ms_grid->size(); ++ss)
      for (const auto& entity : DSC::entityRange(ms_grid->localGridPart(ss))) {
        const auto index = disc.grid_view().indexSet().index(entity);
        indicator_visualization.vector().set_entry(index, eta_Ts[ss]);
        error_visualization.vector().set_entry(disc.grid_view().indexSet().index(entity),
                                               error_distribution_subdomains[ss]);
      }
    indicator_visualization.visualize("indicator_distribution");
    error_visualization.visualize("error_distribution");

    logger.info() << "finished" << std::endl;

  } catch (Dune::Exception& ee) {
    std::cout << "Dune reported error: " << ee.what() << std::endl;
  } catch (std::exception& ee) {
    std::cout << "stl reported error:\n" << ee.what() << std::endl;
  }

  return 0;
} // ... main(...)
