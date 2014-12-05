// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

/**
  * This file is intended to serve as a starting point for quick tests.
  */

#include "config.h"

#if HAVE_DUNE_FEM
# include <dune/fem/misc/mpimanager.hh>
#endif

using namespace Dune;


int main(int argc, char** argv)
{
#if HAVE_DUNE_FEM
  Fem::MPIManager::initialize(argc, argv);
#endif

  return 0;
} // ... main(...)
