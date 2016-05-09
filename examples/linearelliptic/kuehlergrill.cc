// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include <dune/grid/alugrid.hh>
#include <dune/grid/yaspgrid.hh>

#if HAVE_DUNE_FEM
# include <dune/fem/misc/mpimanager.hh>
#endif

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/timedlogging.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/functions.hh>
#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/grid/provider/dgf.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/spaces/fv/default.hh>
#include <dune/gdt/operators/projections.hh>

#include <dune/hdd/linearelliptic/problems.hh>
#include <dune/hdd/linearelliptic/discretizations/cg.hh>

using namespace Dune;


int main(int argc, char** argv)
{
  try {
    Stuff::Common::TimedLogger().create(0);
    auto logger = Stuff::Common::TimedLogger().get("main");
    logger.info() << "creating grid... " << std::flush;

#if HAVE_DUNE_FEM
    Fem::MPIManager::initialize(argc, argv);
#endif
    typedef ALUGrid< 3, 3, simplex, conforming > GridType;
    typedef Stuff::Grid::Providers::DGF< GridType > GridProviderType;
    GridProviderType grid_provider("meshExport10.dgf");
    logger.info() << "done (has " << grid_provider.leaf_view().size(0) << " elements)" << std::endl;

    logger.info() << "visualizing grid... " << std::flush;
    typedef HDD::LinearElliptic::ProblemsProvider< typename GridType::Codim< 0 >::Entity,
                                                   typename GridType::ctype, 3, double, 1 > ProblemProvider;
    Stuff::Common::Configuration problem_cfg;
    problem_cfg["type"] = "hdd.linearelliptic.problem.default";
    problem_cfg["diffusion_factor.type"] = "stuff.function.constant";
    problem_cfg["diffusion_factor.name"] = "diffusion_factor";
    problem_cfg["diffusion_factor.value"] = "1";
    problem_cfg["diffusion_tensor.type"] = "stuff.function.constant";
    problem_cfg["diffusion_tensor.name"] = "diffusion_tensor";
    problem_cfg["diffusion_tensor.value"] = "[1 0 0; 0 1 0; 0 0 1]";
    problem_cfg["force.type"] = "stuff.function.constant";
    problem_cfg["force.name"] = "force";
    problem_cfg["force.value"] = "0";
    problem_cfg["dirichlet.type"] = "stuff.function.constant";
    problem_cfg["dirichlet.name"] = "dirichlet";
    problem_cfg["dirichlet.value"] = "1";
    problem_cfg["neumann.type"] = "stuff.function.indicator";
    problem_cfg["neumann.0.domain"] = "[-1.05 1.05; -2.5 2.5; -15.25 -0.75]";
    problem_cfg["neumann.0.value"] = "-1";
    auto problem = ProblemProvider::create(problem_cfg);

    Stuff::Common::Configuration boundary_info_cfg;
    boundary_info_cfg["type"] = Stuff::Grid::BoundaryInfoConfigs::FunctionBased::static_id();
    boundary_info_cfg["default"] = "neumann";
    boundary_info_cfg["dirichlet.value_range"] = "[0.9 1.1]";
    boundary_info_cfg["dirichlet.function.type"]     = "stuff.function.indicator";
    boundary_info_cfg["dirichlet.function.0.value"]  = "1";
    boundary_info_cfg["dirichlet.function.0.domain"] = "[-1.05 1.05; -4 -3; -16 0]";
    boundary_info_cfg["dirichlet.function.1.value"]  = "1";
    boundary_info_cfg["dirichlet.function.1.domain"] = "[-1.05 1.05; 3 4; -16 0]";

    grid_provider.visualize("grid", boundary_info_cfg);
    problem->visualize(grid_provider.leaf_view(), "problem");
    logger.info() << "done" << std::endl;

    logger.info() << "discretizting... " << std::flush;
    typedef HDD::LinearElliptic::Discretizations::CG< GridType, Stuff::Grid::ChooseLayer::leaf, double, 1, 1 >
        DiscretizationType;
    DiscretizationType discretization(grid_provider, boundary_info_cfg, *problem);
    discretization.init();
    logger.info() << "done" << std::endl;

    logger.info() << "solving... " << std::flush;
    auto solution = discretization.create_vector();
    discretization.solve(solution);
    logger.info() << "done" << std::endl;

    logger.info() << "projecting onto FV space... " << std::flush; // otherwise paraview will crash
    typedef GDT::Spaces::FV::Default< typename DiscretizationType::GridViewType, double, 1 > FvSpaceType;
    FvSpaceType fv_space(discretization.grid_view());
    const auto cg_solution = GDT::make_const_discrete_function(discretization.ansatz_space(), solution, "cg");
    auto fv_solution = GDT::make_discrete_function< typename DiscretizationType::VectorType >(fv_space, "temperature");
    GDT::project(cg_solution, fv_solution);
    auto neumann_values = GDT::make_discrete_function< typename DiscretizationType::VectorType >(fv_space,
                                                                                                 "neumann");
    GDT::project(*problem->neumann()->affine_part(), neumann_values);
    logger.info() << "done" << std::endl;

    logger.info() << "visualizing solution... " << std::flush;
    fv_solution.visualize("solution");
    neumann_values.visualize("neumann");
    logger.info() << "done" << std::endl;

    logger.info() << std::endl;
    return EXIT_SUCCESS;
  } catch (Dune::Exception& e) {
    std::cerr << "\n" << e.what() << std::endl;
    std::abort();
  } catch (std::exception& e) {
    std::cerr << "\n" << e.what() << std::endl;
    std::abort();
  } catch (...) {
    std::cerr << "\nUnknown exception thrown!" << std::endl;
    std::abort();
  } // try
} // ... main(...)
