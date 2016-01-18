// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include <dune/grid/yaspgrid.hh>

#include "generic.hh"

using namespace Dune;
using namespace Dune::Stuff;
using namespace Dune::GDT;


int main(int argc, char** argv)
{
  try {
    typedef GenericLinearellipticExample<YaspGrid<3>, ChooseSpaceBackend::pdelab, LA::ChooseBackend::istl_sparse>
        ExampleType;

    auto logger_cfg = ExampleType::logger_options();
    logger_cfg["info_color"] = "blue";
    logger_cfg["info"] = "99";

    auto grid_cfg = ExampleType::grid_options("stuff.grid.provider.cube");
    grid_cfg["type"] = "stuff.grid.provider.cube";
    grid_cfg["num_elements"] = "[60 220 85]";
    grid_cfg["upper_right"] = "365.76 670.56 51.816";

    auto boundary_cfg = ExampleType::boundary_options("stuff.grid.boundaryinfo.alldirichlet");

    Common::Configuration problem_cfg;
    problem_cfg["type"] = "hdd.linearelliptic.problem.default";
    problem_cfg["diffusion_factor.type"] = "stuff.function.constant";
    problem_cfg["diffusion_factor.name"] = "diffusion_factor";
    problem_cfg["diffusion_factor.value"] = "1";
    problem_cfg["diffusion_tensor.type"] = "stuff.function.spe10.model2";
    problem_cfg["force.type"] = "stuff.function.constant";
    problem_cfg["force.name"] = "force";
    problem_cfg["force.value"] = "1";
    problem_cfg["dirichlet.type"] = "stuff.function.constant";
    problem_cfg["dirichlet.name"] = "dirichlet";
    problem_cfg["dirichlet.value"] = "0";
    problem_cfg["neumann.type"] = "stuff.function.constant";
    problem_cfg["neumann.name"] = "neumann";
    problem_cfg["neumann.value"] = "0";

    ExampleType example(logger_cfg, grid_cfg, boundary_cfg, problem_cfg);


    // if we came that far we can as well be happy about it
    return EXIT_SUCCESS;
  } catch (Dune::Exception& e) {
    std::cerr << "\nDune reported error: " << e.what() << std::endl;
  } catch (std::exception& e) {
    std::cerr << "\n" << e.what() << std::endl;
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  } // try
  return EXIT_FAILURE;
} // ... main(...)
