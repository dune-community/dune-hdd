// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include "config.h"
#endif

#include "cg_gdt.hh"

#include <boost/filesystem.hpp>

int main(int argc, char** argv)
{
  try {
    const std::string settingsFilename = LinearellipticExampleCG::static_id() + ".settings";
    if (!boost::filesystem::exists(settingsFilename))
      LinearellipticExampleCG::writeSettingsFile(settingsFilename);

    LinearellipticExampleCG example(std::vector< std::string >(argv, argv + argc));

    if (example.parametric()) {

    } else {
      const Dune::Pymor::LA::EigenDenseVector* solution = example.solve();
    }

    // if we came that far we can as well be happy about it
    return 0;
  } catch (Dune::PymorException &e) {
    std::cerr << Dune::Stuff::Common::colorStringRed("dune-pymor reported: ") << e << std::endl;
  } catch (Dune::Exception &e) {
    std::cerr << Dune::Stuff::Common::colorStringRed("Dune reported error: ") << e << std::endl;
  } catch (...) {
    std::cerr << Dune::Stuff::Common::colorStringRed("Unknown exception thrown!") << std::endl;
  }
}
