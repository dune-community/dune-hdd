// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include "config.h"
#endif

#include <dune/stuff/common/disable_warnings.hh>
  #include <dune/grid/sgrid.hh>
#include <dune/stuff/common/reenable_warnings.hh>
#include <dune/stuff/common/string.hh>


#include "swipdg-with-dune-gdt.hh"

template< class G >
std::string LinearellipticExampleSWIPDG< G >::static_id()
{
  return "linearelliptic.swipdg-with-dune-gdt";
}

template< class G >
void LinearellipticExampleSWIPDG< G >::writeSettingsFile(const std::string filename)
{
  DiscreteProblemType::writeSettingsFile(filename, static_id());
}

template< class G >
LinearellipticExampleSWIPDG< G >::LinearellipticExampleSWIPDG()
{}

template< class G >
void LinearellipticExampleSWIPDG< G >::initialize(const std::vector< std::string >& arguments)
{
  discreteProblem_ = std::make_shared< const DiscreteProblemType >(static_id(), arguments);
}

//template< class G >
//void LinearellipticExampleSWIPDG< G >::init_discretization()
//{

//}

template class LinearellipticExampleSWIPDG< Dune::SGrid< 2, 2 > >;
