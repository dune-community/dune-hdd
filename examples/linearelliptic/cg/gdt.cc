// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include <config.h>
#endif

#include "gdt.hh"

#include <boost/filesystem.hpp>

int main(int argc, char** argv)
{
  const std::string settingsFilename = id() + ".settings";
  if (!boost::filesystem::exists(settingsFilename))
    ProblemType::writeSettingsFile(settingsFilename);

  return run(argc, argv);
} // main
