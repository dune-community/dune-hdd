// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#ifdef HAVE_DUNE_GRID_MULTISCALE
#undef HAVE_DUNE_GRID_MULTISCALE
#endif

#include <string>
#include <vector>

#if HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
# define ENABLE_ALUGRID 1
# include <dune/grid/alugrid.hh>
#else
# error This example requires alugrid!
#endif

#include "swipdg.hh"

#include <boost/filesystem.hpp>

int main(int argc, char** argv)
{
  try {
    // create empty example
    typedef LinearellipticExampleSWIPDG< Dune::ALUConformGrid< 2, 2 > > ExampleType;

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
        info << "discretization is parametric with parameter_type:" << std::endl;
        info << "  " << discretization.parameter_type() << std::endl;
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

    // if we came that far we can as well be happy about it
    return 0;
  } catch (Dune::Exception& e) {
    std::cerr << "\ndune reported error: " << e.what() << std::endl;
    std::abort();
  } catch (std::exception& e) {
    std::cerr << "\n" << e.what() << std::endl;
    std::abort();
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
    std::abort();
  } // try
} // ... main(...)
