// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include <string>
#include <vector>

#include "cg_mpi.hh"

#include <boost/filesystem.hpp>

#include <dune/stuff/test/common.hh>
#include <dune/stuff/common/convergence-study.hh>
#include <dune/hdd/linearelliptic/testcases/ESV2007.hh>
#include <dune/hdd/linearelliptic/testcases/spe10.hh>
#include <dune/hdd/linearelliptic/testcases/OS2014.hh>
#include <dune/hdd/test/linearelliptic_cg.hh>

void run_example(int argc, char** argv)
{
  // create empty example
  typedef LinearellipticExampleCG ExampleType;


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

void run_eoc_study(int argc, char** argv)
{
  using namespace Dune;
  using namespace Dune::HDD;
  typedef Dune::SPGrid< double, 3 > SPG3;
  typedef Dune::SPGrid< double, 2 > SPG2;
//  typedef LinearElliptic::TestCases::ESV2007< SPG2 > TT; TT test_case;

  typedef LinearElliptic::TestCases::Spe10::Model2< SPG3 > TT; TT test_case;
  test_case.print_header(DSC_LOG_INFO_0);
  DSC_LOG_INFO << std::endl;
  LinearElliptic::Tests::CGStudy< TT, 1, GDT::ChooseSpaceBackend::pdelab, Stuff::LA::ChooseBackend::istl_sparse >
      eoc_study(test_case, {},  {}, "cg_study_");
  const auto fr = eoc_study.run_eoc(DSC_LOG_INFO_0);
  Stuff::Test::print_collected_eoc_study_results(fr, std::cout);
}

int main(int argc, char** argv)
{
  auto& helper = Dune::MPIHelper::instance(argc, argv);
  try {
    DSC::Logger().create(DSC::LOG_CONSOLE | DSC::LOG_INFO | DSC::LOG_DEBUG | DSC::LOG_ERROR , "", "", "");

    DSC::TimedLogger().create(-1, -1);
    DS::threadManager().set_max_threads(1u);
    run_eoc_study(argc, argv);
//    run_example(argc, argv);

    // if we came that far we can as well be happy about it
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
