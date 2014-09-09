// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include "swipdg.hh"


template class LinearellipticExampleSWIPDG< Dune::SGrid< 1, 1 > >;

#if HAVE_ALUGRID

template class LinearellipticExampleSWIPDG< Dune::ALUConformGrid< 2, 2 > >;

#endif // HAVE_ALUGRID
