// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_EXAMPLES_LINEARELLIPTIC_SWIPDG_WITH_DUNE_GDT_HH
#define DUNE_HDD_EXAMPLES_LINEARELLIPTIC_SWIPDG_WITH_DUNE_GDT_HH

#include <memory>

#include <dune/pymor/common/exceptions.hh>

#include "discreteproblem.hh"


template< class GridImp >
class LinearellipticExampleSWIPDG
{
  typedef DiscreteProblem< GridImp >                    DiscreteProblemType;
//  typedef typename DiscreteProblemType::SettingsType    SettingsType;
//  typedef typename DiscreteProblemType::GridPartType    GridPartType;
//  typedef typename DiscreteProblemType::RangeFieldType  RangeFieldType;
//  static const int DUNE_UNUSED(                         dimDomain) = DiscreteProblemType::dimDomain;
//  static const int DUNE_UNUSED(                         dimRange) = DiscreteProblemType::dimRange;

public:
  static std::string static_id();

  static void writeSettingsFile(const std::string filename);

  LinearellipticExampleSWIPDG();

  void initialize(const std::vector< std::string >& arguments);

  bool initialized();

private:
  std::shared_ptr< const DiscreteProblemType > discreteProblem_;
//  std::shared_ptr< const DiscreteProblem
}; // class LinearellipticExampleCG


#endif // DUNE_HDD_EXAMPLES_LINEARELLIPTIC_SWIPDG_WITH_DUNE_GDT_HH
