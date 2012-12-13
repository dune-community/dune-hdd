#include "detailed_discretizations.hh"

#include <boost/filesystem.hpp>

int main(int argc, char** argv)
{
  const std::string paramFilename = id + ".param";
  if (!boost::filesystem::exists(paramFilename))
    writeParamFile(paramFilename);

  run(argc, argv);

  // if we came that far we can as well be happy about it
  return 0;
} // main
