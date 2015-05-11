// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include <string>
#include <vector>

#include "mpi_cg_example.hh"

#include <boost/filesystem.hpp>

#include <dune/stuff/test/common.hh>
#include <dune/stuff/common/convergence-study.hh>
#include <dune/hdd/linearelliptic/testcases/ESV2007.hh>
#include <dune/hdd/linearelliptic/testcases/spe10model2.hh>
#include <dune/hdd/linearelliptic/testcases/OS2014.hh>
#include <dune/hdd/test/linearelliptic_cg.hh>

void run_example(int argc, char** argv)
{
  // create empty example
  typedef MpiCGExample ExampleType;


  // read or write config file
  const std::string config_file_name = ExampleType::static_id() + ".cfg";
  if (!boost::filesystem::exists(config_file_name)) {
    std::cout << "Writing default configuration to '" << config_file_name << "'... " << std::flush;
    ExampleType::write_config_file(config_file_name);
    std::cout << "done.\n"
              << "Please review the configuration and start me again!" << std::endl;
  } else {

    // init discrete problem and discretization
    ExampleType example;
    example.initialize(std::vector< std::string >(argv, argv + argc));
    const auto& discreteProblem = example.discrete_problem();
    auto& info = DSC_LOG_INFO;
    Dune::Timer timer;

    const auto& discretization = example.discretization();
    auto solution = discretization.create_vector();

    // solve
    if (discretization.parametric()) {
      info << "discretization is parametric with parameter_type: " << discretization.parameter_type() << std::endl;
      const auto& config = discreteProblem.config();
      if (config.has_sub("parameter")) {
        const auto parameters = config.sub("parameter");
        size_t pp = 0;
        while (parameters.has_sub(Dune::Stuff::Common::toString(pp))) {
          const auto parameter = parameters.sub(Dune::Stuff::Common::toString(pp));
          Dune::Pymor::Parameter mu;
          for (std::string key : parameter.getValueKeys())
            mu.set(key, parameter.get< std::vector< double > >(key));
          info << "solving for mu = " << mu << "... " << std::flush;
          timer.reset();
          discretization.solve(solution, mu);
          info << " done (took " << timer.elapsed() << "s)" << std::endl;
          discretization.visualize(solution,
                                   example.static_id() + ".solution_to_parameter_" + Dune::Stuff::Common::toString(pp),
                                   "solution to parameter " + Dune::Stuff::Common::toString(pp));
          ++pp;
        }
      } else
        info << "doing nothing, since there is no 'parameter' specified in the config!" << std::endl;
    } else {
      info << "discretization is not parametric, solving... " << std::flush;
      timer.reset();
      discretization.solve(solution, Dune::Pymor::Parameter());
      info << " done (took " << timer.elapsed() << "s)" << std::endl;
      discretization.visualize(solution, example.static_id() + ".solution", "solution");
    }

  } // read or write config file
}

void run_eoc_study(DSC::Configuration& config)
{
  using namespace Dune;
  using namespace Dune::HDD;
//  typedef Dune::SPGrid< double, 3 > SPG3;
  typedef Dune::SPGrid< double, 2 > SPG2;
  typedef LinearElliptic::TestCases::ESV2007< SPG2 > TestCase;

//  typedef LinearElliptic::TestCases::Spe10::Model2< SPG3 > TestCase;
  TestCase test_case(0, config.get<size_t>("grids.refinements", 4u));
  test_case.print_header(DSC_LOG_INFO_0);
  DSC_LOG_INFO << std::endl;
  LinearElliptic::Tests::CGStudy< TestCase, 1, GDT::ChooseSpaceBackend::pdelab, Stuff::LA::ChooseBackend::istl_sparse >
      ::DiscretizationType  disc(test_case,
            test_case.boundary_info(),
            test_case.problem(),
            test_case.level_of(0));
//  test_case.problem().visualize(test_case.grid().leafGridView(), "problem");
  disc.init();
  auto options = disc.solver_options();
  config.add(options, "solver");
  if(!config.has_key("solver.verbose")) config["verbose"] = "6";
  if(!config.has_key("solver.precision")) config["precision"] = "1e-1";
  if(!config.has_key("solver.type")) config["type"] = "bicgstab.amg.ssor";
  if(!config.has_key("solver.max_iter")) config["max_iter"] = "700";
  auto solution = disc.create_vector();
//  test_case.problem().neumann()->affine_part()->visualize(disc.grid_view(), "neumann");
//  test_case.visualize(test_case.boundary_info());
  try {
    disc.solve(config.sub("solver"), solution);
    double min = std::numeric_limits< double >::max();
    for (const auto& element : solution)
      min = std::min(min, element);

    solution += std::abs(min) + 0.1;
    DSC_LOG_INFO_0 << "Min " << min << " \n";
  }
   catch (Dune::Stuff::Exceptions::linear_solver_failed) {

  }
  disc.visualize(solution, "spe10_solution", "solution");
}

int main(int argc, char** argv)
{
  auto& helper = Dune::MPIHelper::instance(argc, argv);
  try {
    DSC::Configuration config(false, true, true);
    config.read_command_line(argc, argv);
    config.set_record_defaults(true);
    DSC::Logger().create(DSC::LOG_CONSOLE | DSC::LOG_INFO | DSC::LOG_DEBUG | DSC::LOG_ERROR , "", "", "");
    DSC::testCreateDirectory(DSC_CONFIG_GET("global.datadir", "data/"));

    DSC::TimedLogger().create(-1, -1);
    DS::threadManager().set_max_threads(1u);
    run_eoc_study(config);
    run_example(argc, argv);
    std::unique_ptr<boost::filesystem::ofstream> of(DSC::make_ofstream(
                           std::string(DSC_CONFIG_GET("global.datadir", "data/")) + std::string("/env.txt")));
    config.report(*of);
    return 0;
  } catch (Dune::Exception& e) {
    std::cerr << "\nDune reported error: " << e.what() << std::endl;
  } catch (std::exception& e) {
    std::cerr << "\n" << e.what() << std::endl;
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
  return 1;
} // ... main(...)
