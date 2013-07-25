// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_EXAMPLES_LINEARELLIPTIC_CG_WITH_DUNE_GDT_HH
#define DUNE_HDD_EXAMPLES_LINEARELLIPTIC_CG_WITH_DUNE_GDT_HH

#include <memory>

#include <dune/pymor/common/exceptions.hh>

#include <dune/hdd/linearelliptic/discretizations/cg-with-dune-gdt.hh>

#include "discreteproblem.hh"

namespace Elliptic = Dune::HDD::LinearElliptic;

template< class GridImp, int polynomialOrder = 1 >
class LinearellipticExampleCG
{
public:
  typedef DiscreteProblem< GridImp >                    DiscreteProblemType;
private:
  typedef typename DiscreteProblemType::GridPartType    GridPartType;
  typedef typename DiscreteProblemType::RangeFieldType  RangeFieldType;
  static const unsigned int DUNE_UNUSED(                dimRange) = DiscreteProblemType::dimRange;
  static const unsigned int DUNE_UNUSED(                polOrder) = polynomialOrder;
public:
  typedef Elliptic::Discretization::ContinuousGalerkinWithDuneGDT<  GridPartType, RangeFieldType, dimRange,
                                                                    polOrder > DiscretizationType;

  static std::string static_id();

  static void write_settings_file(const std::string filename);

  LinearellipticExampleCG();

  void initialize(const std::vector< std::string >& arguments);

  bool initialized() const;

  DiscreteProblemType discrete_problem() const;

  DiscretizationType discretization() const;

private:
  bool initialized_;
  std::shared_ptr< const DiscreteProblemType > discreteProblem_;
  std::shared_ptr< DiscretizationType > discretization_;
}; // class LinearellipticExampleCG


#endif // DUNE_HDD_EXAMPLES_LINEARELLIPTIC_CG_WITH_DUNE_GDT_HH
