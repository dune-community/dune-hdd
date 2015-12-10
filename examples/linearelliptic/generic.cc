// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifdef ENABLE_MPI
# undef ENABLE_MPI
#endif
#define ENABLE_MPI 0

#ifdef ENABLE_PARMETIS
# undef ENABLE_PARMETIS
#endif
#define ENABLE_PARMETIS 0

#include "config.h"

#if !(HAVE_ALUGRID || HAVE_DUNE_ALUGRID)
# error Missing ALUGrid!
#endif

#include <dune/grid/alugrid.hh>

#include "generic.hh"

using namespace Dune;


int main(int /*argc*/, char** /*argv*/)
{
  try {
    typedef GenericLinearellipticExample< ALUGrid< 3, 3, simplex, nonconforming >,
                                          GDT::ChooseSpaceBackend::fem,
                                          Stuff::LA::ChooseBackend::istl_sparse > ExampleType;

    auto logger_cfg = ExampleType::logger_options();
    logger_cfg["info"] = "99";
    logger_cfg["info_color"] = "blue";

    Stuff::Common::Configuration battery_geometry;
    battery_geometry["lower_left"]   = "[0      0     0]";
    battery_geometry["upper_right"]  = "[0.0184 0.008 0.008]";
    battery_geometry["num_elements"] = "[46     20    20]";
    battery_geometry["separator"]    = "[0.0084 0.01; 0 0.008; 0 0.008]";
    battery_geometry["filename"]     = "geometry__46x20x20_h4e-6m";

    auto grid_cfg = ExampleType::grid_options("stuff.grid.provider.cube");
    grid_cfg.add(battery_geometry, "", true);

    Stuff::Common::Configuration boundary_cfg;
    boundary_cfg["type"]        = "stuff.grid.boundaryinfo.normalbased";
    boundary_cfg["default"]     = "neumann";
    boundary_cfg["dirichlet.0"] = "[-1 0 0]";
    boundary_cfg["dirichlet.1"] = "[1 0 0]";

    auto problem_cfg = ExampleType::problem_options("hdd.linearelliptic.problem.battery");
    problem_cfg.add(battery_geometry, "diffusion_factor", true);
    problem_cfg["dirichlet.expression"] = "0";
    problem_cfg["neumann.expression"] = "0";
    problem_cfg["force.value"] = "1000";

    auto solver_options = ExampleType::solver_options("bicgstab.amg.ilu0");

    ExampleType example(logger_cfg, grid_cfg, boundary_cfg, problem_cfg);
    auto& disc = example.discretization();
    auto U = disc.create_vector();
    disc.solve(solver_options, U, Pymor::Parameter("ELECTROLYTE", 0.6));
    disc.visualize(U, "solution", "solution");

    DSC::TimedLogger().get("main").info() << "finished!" << std::endl;

    // if we came that far we can as well be happy about it
    return EXIT_SUCCESS;
  } catch (Dune::Exception& e) {
    std::cerr << "\n" << e.what() << std::endl;
  } catch (std::exception& e) {
    std::cerr << "\n" << e.what() << std::endl;
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  } // try
  return EXIT_FAILURE;
} // ... main(...)
