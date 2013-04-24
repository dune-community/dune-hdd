#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include <config.h>
#endif

#include "detailed_discretizations.hh"

#include <boost/filesystem.hpp>

int main(int argc, char** argv)
{
  const std::string descriptionFilename = id() + ".description";
  if (!boost::filesystem::exists(descriptionFilename))
    ProblemType::writeDescriptionFile(descriptionFilename);

  return run(argc, argv);
} // main
