// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include "config.h"
#endif

#include <string>
#include <vector>

#include <dune/pymor/common/exceptions.hh>

#include "swipdg-with-dune-gdt.hh"

#include <boost/filesystem.hpp>

int main(int argc, char** argv)
{
  try {
    // create empty example
    typedef LinearellipticExampleSWIPDG< Dune::SGrid< 2, 2 > > LinearellipticExample;
    LinearellipticExample example;

    // write settings file
    const std::string settingsFilename = example.static_id() + ".settings";
    if (!boost::filesystem::exists(settingsFilename))
      example.writeSettingsFile(settingsFilename);

    // init problem
    example.init_discrete_probem(std::vector< std::string >(argv, argv + argc));

    // if we came that far we can as well be happy about it
    return 0;
  } catch (Dune::PymorException& e) {
    std::cerr << e << std::endl;
  } catch (Dune::Exception& e) {
    std::cerr << Dune::Stuff::Common::colorStringRed("Dune reported error: ") << e << std::endl;
  } catch (...) {
    std::cerr << Dune::Stuff::Common::colorStringRed("Unknown exception thrown!") << std::endl;
  }
}
