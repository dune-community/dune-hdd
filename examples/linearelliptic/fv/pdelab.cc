#include "pdelab.hh"

#include <boost/filesystem.hpp>

int main(int argc, char** argv)
{
  const std::string paramFilename = id + ".param";
  if (!boost::filesystem::exists(paramFilename))
    writeParamFile(paramFilename);

  return run(argc, argv);
} // main
