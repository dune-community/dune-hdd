// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

/**
  * This file is intended to serve as a starting point for quick tests.
  */

#include "config.h"

#include <exception>

#include <boost/exception/diagnostic_information.hpp>
#include <boost/exception/exception.hpp>

#if HAVE_DUNE_FEM
# include <dune/fem/misc/mpimanager.hh>
#endif

#include <dune/stuff/common/exceptions.hh>

using namespace Dune;


int main(int argc, char** argv)
{
  try {
#if HAVE_DUNE_FEM
    Fem::MPIManager::initialize(argc, argv);
#endif

  } catch (Dune::Exception& ee) {
    std::cerr << ee.what() << std::endl;
    return EXIT_FAILURE;
  } catch (boost::exception& ee) {
    std::cerr << boost::diagnostic_information(ee) << std::endl;
    return EXIT_FAILURE;
  } catch (std::exception& ee) {
    std::cerr << ee.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
} // ... main(...)
