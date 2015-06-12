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
#include <dune/hdd/linearelliptic/testcases/thermalblock.hh>
#include <dune/hdd/linearelliptic/testcases/random_block_testcase.hh>
#include <dune/hdd/test/linearelliptic_cg.hh>
#include <dune/grid/alugrid.hh>
#include <dune/grid/uggrid.hh>

//void run_example(int argc, char** argv)
//{
//  // create empty example
//  typedef MpiCGExample ExampleType;


//  // read or write config file
//  const std::string config_file_name = ExampleType::static_id() + ".cfg";
//  if (!boost::filesystem::exists(config_file_name)) {
//    std::cout << "Writing default configuration to '" << config_file_name << "'... " << std::flush;
//    ExampleType::write_config_file(config_file_name);
//    std::cout << "done.\n"
//              << "Please review the configuration and start me again!" << std::endl;
//  } else {

//    // init discrete problem and discretization
//    ExampleType example;
//    example.initialize(std::vector< std::string >(argv, argv + argc));
//    const auto& discreteProblem = example.discrete_problem();
//    auto& info = DSC_LOG_INFO;
//    Dune::Timer timer;

//    const auto& discretization = example.discretization();
//    auto solution = discretization.create_vector();

//    // solve
//    if (discretization.parametric()) {
//      info << "discretization is parametric with parameter_type: " << discretization.parameter_type() << std::endl;
//      const auto& config = discreteProblem.config();
//      if (config.has_sub("parameter")) {
//        const auto parameters = config.sub("parameter");
//        size_t pp = 0;
//        while (parameters.has_sub(Dune::Stuff::Common::toString(pp))) {
//          const auto parameter = parameters.sub(Dune::Stuff::Common::toString(pp));
//          Dune::Pymor::Parameter mu;
//          for (std::string key : parameter.getValueKeys())
//            mu.set(key, parameter.get< std::vector< double > >(key));
//          info << "solving for mu = " << mu << "... " << std::flush;
//          timer.reset();
//          discretization.solve(solution, mu);
//          info << " done (took " << timer.elapsed() << "s)" << std::endl;
//          discretization.visualize(solution,
//                                   example.static_id() + ".solution_to_parameter_" + Dune::Stuff::Common::toString(pp),
//                                   "solution to parameter " + Dune::Stuff::Common::toString(pp));
//          ++pp;
//        }
//      } else
//        info << "doing nothing, since there is no 'parameter' specified in the config!" << std::endl;
//    } else {
//      info << "discretization is not parametric, solving... " << std::flush;
//      timer.reset();
//      discretization.solve(solution, Dune::Pymor::Parameter());
//      info << " done (took " << timer.elapsed() << "s)" << std::endl;
//      discretization.visualize(solution, example.static_id() + ".solution", "solution");
//    }

//  } // read or write config file
//}

template <class F>
void log_shift(F& f)
{
  //    double min = std::numeric_limits< double >::max();
  //    for (const auto& element : solution)
  //      min = std::min(min, element);

  //    const double all_min = test_case.grid().comm().min(min);
  //    constexpr double target = 1.0;
  //    const double shift = target - all_min;
  //    solution += std::abs(all_min) + 0.1;
  //    solution += shift;
  //    DSC_LOG_INFO_0 << "Min " << all_min << "Shift " << shift << " \n";
}

void run_eoc_study(DSC::Configuration& config)
{
  using namespace Dune;
  using namespace Dune::HDD;
  constexpr size_t dim = 3;
  typedef Dune::SPGrid< double, dim > SPG2;
//  typedef Dune::YaspGrid< 2 > SPG2;
//  typedef Dune::UGGrid< 2 > SPG2;
//  typedef Dune::ALUGrid< 2, 2, simplex, conforming, MPI_Comm > SPG2;
//  typedef LinearElliptic::TestCases::ESV2007< SPG2 > TestCase;
  //  typedef LinearElliptic::TestCases::Spe10::Model1< SPG2 > TestCase;
    typedef LinearElliptic::TestCases::RandomBlockTestcase< SPG2 > TestCase;
//  typedef LinearElliptic::TestCases::Spe10::Model2< SPG3 > TestCase;

  const DSC::ValueInitFieldVector<size_t, dim, 2u> blocks;
  const unsigned int overlap_size = config.get<size_t>("grids.overlap", 4u);
  TestCase test_case(config.get<size_t>("grids.refinements", 4u), blocks, overlap_size);
  test_case.print_header(DSC_LOG_INFO_0);
  DSC_LOG_INFO << std::endl;
  LinearElliptic::Tests::CGStudy< TestCase, 1, GDT::ChooseSpaceBackend::pdelab, Stuff::LA::ChooseBackend::istl_sparse >
      ::DiscretizationType  disc(test_case,
            test_case.boundary_info(),
            test_case.problem());
//  test_case.problem().visualize(test_case.grid().leafGridView(), "problem");
  disc.init();
  auto solution = disc.create_vector();
  test_case.problem().visualize(disc.grid_view(), "foo");

//  test_case.visualize(test_case.boundary_info());
  try {
//    {0.1,1,1,1}
    const auto mu = Dune::Pymor::Parameter("diffusion", {0.15227525});//, 0.87955853, 0.24041678, 0.24039507, 1, 1, 1, 1  });
//    const auto mu = Dune::Pymor::Parameter("diffusion", {1, 1, 1, 1 , 1, 1, 1, 1  });
    const auto sub = config.sub("solver");
//    disc.solve(sub, solution, mu);
    disc.solve(sub, solution, mu);

    DSC_LOG_DEBUG << test_case.grid().comm().rank() << " | "<< test_case.grid().comm().size() << std::endl;
   }
   catch (Dune::Stuff::Exceptions::linear_solver_failed& e) {
    DSC_LOG_ERROR <<  e.what();
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
//    run_example(argc, argv);
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
