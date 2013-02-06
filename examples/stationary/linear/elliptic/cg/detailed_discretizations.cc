#include "detailed_discretizations.hh"

#include <boost/filesystem.hpp>

int main(int argc, char** argv)
{
  const std::string descriptionFilename = id() + ".description";
  if (!boost::filesystem::exists(descriptionFilename))
    writeDescriptionFile(descriptionFilename);

  return run(argc, argv);
} // main
