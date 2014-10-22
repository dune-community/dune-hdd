// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_PROBLEMS_ZERO_BOUNDARY_HH
#define DUNE_HDD_LINEARELLIPTIC_PROBLEMS_ZERO_BOUNDARY_HH

#include <dune/stuff/functions/constant.hh>

#include <dune/pymor/functions/default.hh>

#include <dune/hdd/linearelliptic/problems/default.hh>

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Problems {


template< class ProblemImp >
class ZeroBoundary
  : public Problems::Default< typename ProblemImp::EntityType
                            , typename ProblemImp::DomainFieldType, ProblemImp::dimDomain
                            , typename ProblemImp::RangeFieldType, ProblemImp::dimRange >
{
  typedef Problems::Default< typename ProblemImp::EntityType
                           , typename ProblemImp::DomainFieldType, ProblemImp::dimDomain
                           , typename ProblemImp::RangeFieldType, ProblemImp::dimRange > BaseType;
public:
  using typename BaseType::EntityType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::RangeFieldType;
  static const unsigned int dimDomain = BaseType::dimDomain;
  static const unsigned int dimRange = BaseType::dimRange;

private:
  typedef Stuff::Functions::Constant< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >
    ConstantFunctionType;
  typedef Pymor::Functions::NonparametricDefault< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >
    NonparametricFunctionType;

public:
  static std::string static_id()
  {
    return BaseType::static_id() + ".zero-boundary";
  }

  virtual std::string type() const
  {
    return BaseType::static_id() + ".zero-boundary";
  }

  ZeroBoundary(const ProblemImp& problem)
    : BaseType(problem.diffusion_factor(),
               problem.diffusion_tensor(),
               problem.force(),
               std::make_shared< NonparametricFunctionType >(new ConstantFunctionType(0)),
               std::make_shared< NonparametricFunctionType >(new ConstantFunctionType(0)))
  {}
}; // class ZeroBoundary


} // namespace Problems
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_PROBLEMS_ZERO_BOUNDARY_HH
