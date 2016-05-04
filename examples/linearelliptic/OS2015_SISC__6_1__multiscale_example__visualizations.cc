// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#if HAVE_DUNE_FEM && HAVE_ALUGRID

# ifndef DUNE_GDT_LOCALEVALUATION_SWIPDG_DISABLE_WARNINGS
#   define DUNE_GDT_LOCALEVALUATION_SWIPDG_DISABLE_WARNINGS
# endif
# ifndef DUNE_HDD_LINEARELLIPTIC_TESTCASES_BASE_DISABLE_WARNING
#   define DUNE_HDD_LINEARELLIPTIC_TESTCASES_BASE_DISABLE_WARNING
# endif

# include <dune/fem/misc/mpimanager.hh>

# include <dune/stuff/common/ranges.hh>
# include <dune/stuff/common/timedlogging.hh>
# include <dune/stuff/grid/search.hh>

# include <dune/gdt/operators/prolongations.hh>
# include <dune/gdt/playground/operators/fluxreconstruction.hh>
# include <dune/gdt/playground/spaces/finitevolume/default.hh>
# include <dune/gdt/playground/spaces/raviartthomas/pdelab.hh>
# include <dune/gdt/products/elliptic.hh>

# include <dune/hdd/linearelliptic/testcases/OS2015.hh>
# include <dune/hdd/linearelliptic/discretizations/block-swipdg.hh>
# include <dune/hdd/linearelliptic/estimators/block-swipdg.hh>

using namespace Dune;
using namespace Dune::Stuff;
using namespace Dune::Pymor;
using namespace Dune::GDT;
using namespace Dune::HDD;


static const Parameter mu_min("mu", 0.1);
static const Parameter mu_max("mu", 1.0);
static const Parameter mu_bar("mu", 0.1);
static const Parameter mu_hat("mu", 0.1);

static const size_t over_integrate = 2;


template< class SpaceType, class VectorType, class ProblemType >
Stuff::LA::CommonDenseVector< double > compute_indicator_distribution(const SpaceType& space,
                                                                      const VectorType& vector,
                                                                      const ProblemType& problem,
                                                                      const Parameter& mu)
{
  typedef LinearElliptic::Estimators::BlockSWIPDG< SpaceType,
                                                   VectorType,
                                                   ProblemType,
                                                   typename SpaceType::GridViewType::Grid > Estimator;
  auto eta_Ts = Estimator::estimate_local(space,
                                          vector,
                                          problem,
                                          "eta_OS2014_*",
                                          {{"mu",     mu},
                                           {"mu_bar", mu_bar},
                                           {"mu_hat", mu_hat},
                                           {"parameter_range_min", mu_min},
                                           {"parameter_range_max", mu_max}});
  double eta = 0.0;
  for (const auto& eta_T : eta_Ts)
    eta += std::pow(eta_T, 2);
  eta = std::sqrt(eta);
  for (auto& eta_T : eta_Ts)
    eta_T /= eta;
  return eta_Ts;
} // ... compute_indicator_distribution(...)


template< class DiscretizationType, class DifferenceType, class DiffusionFactorType, class DiffusionTensorType >
Stuff::LA::CommonDenseVector< double > compute_error_distribution(const DiscretizationType& discretization,
                                                                  const DiscretizationType& reference_discretization,
                                                                  const DifferenceType& difference,
                                                                  const DiffusionFactorType& diffusion_factor,
                                                                  const DiffusionTensorType& diffusion_tensor)
{
  auto error_product = Products::make_elliptic_localizable(reference_discretization.grid_view(),
                                                           difference,
                                                           diffusion_factor,
                                                           diffusion_tensor,
                                                           over_integrate);
  error_product.prepare();
  auto error_distribution = discretization.create_vector();
  auto entity_search = Stuff::Grid::make_entity_in_level_search(discretization.grid_view());
  for (const auto& entity : DSC::entityRange(reference_discretization.grid_view())) {
    // search for the father entity in the (coarser) grid view
    const auto points = {entity.geometry().center()};
    const auto father_entity_ptr_ptrs = entity_search(points);
    if (father_entity_ptr_ptrs.size() != 1)
      DUNE_THROW(Stuff::Exceptions::internal_error, father_entity_ptr_ptrs.size());
    const auto& father_entity_ptr_ptr = father_entity_ptr_ptrs[0];
    assert(father_entity_ptr_ptr);
    const auto father_entity_ptr = *father_entity_ptr_ptr;
    const auto& father_entity = *father_entity_ptr;
    assert(discretization.grid_view().contains(father_entity));
    // compute local error
    error_distribution.add_to_entry(discretization.grid_view().indexSet().index(father_entity),
                                    error_product.compute_locally(entity));
  }
  // sum local errors on the subdomain level
  const auto ms_grid = discretization.ansatz_space()->ms_grid();
  auto error_distribution_on_subdomains = Stuff::LA::CommonDenseVector< double >(ms_grid->size(), 0.0);
  error_distribution_on_subdomains *= 0.0;
  for (size_t ss = 0; ss < ms_grid->size(); ++ss)
    for (const auto& entity : DSC::entityRange(ms_grid->localGridPart(ss)))
      error_distribution_on_subdomains[ss]
          += error_distribution[discretization.grid_view().indexSet().index(entity)];
  double error = 0.0;
  for (const auto& element : error_distribution_on_subdomains)
    error += element;
  error = std::sqrt(error);
  for (auto& element : error_distribution_on_subdomains)
    element = std::sqrt(element) / error;
  return error_distribution_on_subdomains;
} // ... compute_error_distribution(...)


int main(int argc, char** argv)
{
  try {
    Dune::Fem::MPIManager::initialize(argc, argv);

    DSC::TimedLogger().create(0, 0);
    auto logger = DSC::TimedLogger().get("main");
    logger.info() << "creating grid and test case ..." << std::endl;

    typedef ALUGrid< 2, 2, simplex, conforming > GridType;
    typedef LinearElliptic::TestCases::OS2015::Multiscale< GridType > TestCaseType;
    const TestCaseType test_case({{"mu",     mu_max},
                                  {"mu_bar", mu_bar},
                                  {"mu_hat", mu_hat},
                                  {"parameter_range_min", mu_min},
                                  {"parameter_range_max", mu_max}},
                                 "[25 5 1 ]");
    const auto& problem = test_case.problem();
    auto problem_mu_min = problem.with_mu(mu_min);
    auto problem_mu_max = problem.with_mu(mu_max);
    auto problem_mu_bar = problem.with_mu(mu_bar);

    logger.info() << "discretizing ..." << std::endl;
    auto disc = LinearElliptic::Discretizations::make_block_swipdg< Stuff::LA::ChooseBackend::eigen_sparse >(
                  *test_case.level_provider(0),
                  test_case.boundary_info(),
                  problem);
    disc.init();

    logger.info() << "visualizing problems ..." << std::endl;
    problem.visualize(disc.grid_view(), "parametric_problem");
    problem_mu_min->visualize(disc.grid_view(), "problem_mu_min");
    problem_mu_max->visualize(disc.grid_view(), "problem_mu_max");

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

    logger.info() << "reconstructing diffusive fluxes ..." << std::endl;
    Spaces::RaviartThomas::PdelabBased< typename std::remove_reference< decltype(disc) >::type::GridViewType,
                                        0, double, 2 > rt_space(disc.grid_view());
    auto reconstruction_mu_min = GDT::Operators::make_diffusive_flux_reconstruction(disc.grid_view(),
                                                                                    *problem_mu_min->diffusion_factor()->affine_part(),
                                                                                    *problem_mu_min->diffusion_tensor()->affine_part(),
                                                                                    over_integrate);
    auto reconstruction_mu_max = GDT::Operators::make_diffusive_flux_reconstruction(disc.grid_view(),
                                                                                    *problem_mu_max->diffusion_factor()->affine_part(),
                                                                                    *problem_mu_max->diffusion_tensor()->affine_part(),
                                                                                    over_integrate);
    auto velocity_mu_min = make_discrete_function< VectorType >(rt_space, "velocity");
    auto velocity_mu_max = make_discrete_function< VectorType >(rt_space, "velocity");
    reconstruction_mu_min.apply(pressure_mu_min, velocity_mu_min);
    reconstruction_mu_max.apply(pressure_mu_max, velocity_mu_max);

    logger.info() << "visualizing diffusive fluxes ..." << std::endl;
    velocity_mu_min.visualize("velocity_mu_min");
    velocity_mu_max.visualize("velocity_mu_max");

    logger.info() << "computing local estimator contributions ..." << std::endl;
    const auto eta_Ts_mu_min = compute_indicator_distribution(*disc.ansatz_space(),
                                                                    pressure_mu_min.vector(),
                                                                    problem,
                                                                    mu_min);
    const auto eta_Ts_mu_max = compute_indicator_distribution(*disc.ansatz_space(),
                                                                    pressure_mu_max.vector(),
                                                                    problem,
                                                                    mu_max);

    logger.info() << "computing reference solutions ..." << std::endl;
    auto reference_disc = LinearElliptic::Discretizations::make_block_swipdg< Stuff::LA::ChooseBackend::eigen_sparse >(
                  *test_case.reference_provider(),
                  test_case.boundary_info(),
                  problem);
    reference_disc.init();
    auto reference_pressure_mu_min = make_discrete_function< VectorType >(*reference_disc.ansatz_space());
    auto reference_pressure_mu_max = make_discrete_function< VectorType >(*reference_disc.ansatz_space());
    reference_disc.solve(reference_pressure_mu_min.vector(), mu_min);
    reference_disc.solve(reference_pressure_mu_max.vector(), mu_max);

    logger.info() << "prolonging solutions to reference grid ..." << std::endl;
    auto pressure_mu_min_on_reference = make_discrete_function< VectorType >(*reference_disc.ansatz_space());
    auto pressure_mu_max_on_reference = make_discrete_function< VectorType >(*reference_disc.ansatz_space());
    GDT::Operators::prolong(pressure_mu_min, pressure_mu_min_on_reference);
    GDT::Operators::prolong(pressure_mu_max, pressure_mu_max_on_reference);

    logger.info() << "computing local error contributions ..." << std::endl;
    auto error_distribution_mu_min_on_subdomains
        = compute_error_distribution(disc,
                                     reference_disc,
                                     reference_pressure_mu_min - pressure_mu_min_on_reference,
                                     *problem_mu_bar->diffusion_factor()->affine_part(),
                                     *problem_mu_bar->diffusion_tensor()->affine_part());
    auto error_distribution_mu_max_on_subdomains
        = compute_error_distribution(disc,
                                     reference_disc,
                                     reference_pressure_mu_max - pressure_mu_max_on_reference,
                                     *problem_mu_bar->diffusion_factor()->affine_part(),
                                     *problem_mu_bar->diffusion_tensor()->affine_part());

    logger.info() << "visualizing local errors and indicators ..." << std::endl;
    Spaces::FiniteVolume::Default< typename std::remove_reference< decltype(disc) >::type::GridViewType,
                                   double, 1 > fv_space(disc.grid_view());
    auto indicator_visualization_mu_min = make_discrete_function< VectorType >(fv_space, "distribution");
    auto indicator_visualization_mu_max = make_discrete_function< VectorType >(fv_space, "distribution");
    auto error_visualization_mu_min = make_discrete_function< VectorType >(fv_space, "distribution");
    auto error_visualization_mu_max = make_discrete_function< VectorType >(fv_space, "distribution");
    const auto ms_grid = disc.ansatz_space()->ms_grid();
    for (size_t ss = 0; ss < ms_grid->size(); ++ss)
      for (const auto& entity : DSC::entityRange(ms_grid->localGridPart(ss))) {
        const auto index = disc.grid_view().indexSet().index(entity);
        indicator_visualization_mu_min.vector().set_entry(index, eta_Ts_mu_min[ss]);
        indicator_visualization_mu_max.vector().set_entry(index, eta_Ts_mu_max[ss]);
        error_visualization_mu_min.vector().set_entry(disc.grid_view().indexSet().index(entity),
                                                      error_distribution_mu_min_on_subdomains[ss]);
        error_visualization_mu_max.vector().set_entry(disc.grid_view().indexSet().index(entity),
                                                      error_distribution_mu_max_on_subdomains[ss]);
      }
    indicator_visualization_mu_min.visualize("indicator_distribution_mu_min");
    indicator_visualization_mu_max.visualize("indicator_distribution_mu_max");
    error_visualization_mu_min.visualize("error_distribution_mu_min");
    error_visualization_mu_max.visualize("error_distribution_mu_max");

    logger.info() << "finished" << std::endl;

  } catch (Dune::Exception& ee) {
    std::cout << "Dune reported error: " << ee.what() << std::endl;
  } catch (std::exception& ee) {
    std::cout << "stl reported error:\n" << ee.what() << std::endl;
  }

  return 0;
} // ... main(...)


#elif HAVE_DUNE_FEM // HAVE_DUNE_FEM && HAVE_ALUGRID


int main(int /*argc*/, char** /*argv*/)
{
  std::cerr << "Error: you are missing ALUGrid!" << std::endl;
}


#elif HAVE_ALUGRID // HAVE_DUNE_FEM && HAVE_ALUGRID


int main(int /*argc*/, char** /*argv*/)
{
  std::cerr << "Error: you are missing dune-fem!" << std::endl;
}


#else // HAVE_DUNE_FEM && HAVE_ALUGRID


int main(int /*argc*/, char** /*argv*/)
{
  std::cerr << "Error: you are missing dune-fem and ALUGrid!" << std::endl;
}


#endif // HAVE_DUNE_FEM && HAVE_ALUGRID
