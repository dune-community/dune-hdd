// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_PROBLEMS_ZERO_BOUNDARY_HH
#define DUNE_HDD_LINEARELLIPTIC_PROBLEMS_ZERO_BOUNDARY_HH

#include <dune/stuff/functions/constant.hh>

#include <dune/pymor/functions/default.hh>

#include <dune/hdd/linearelliptic/problems/interfaces.hh>

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Problems {


template< class ProblemImp >
class ZeroBoundary
  : public ProblemInterface< typename ProblemImp::EntityType
                           , typename ProblemImp::DomainFieldType, ProblemImp::dimDomain
                           , typename ProblemImp::RangeFieldType, ProblemImp::dimRange >
{
  typedef ProblemInterface< typename ProblemImp::EntityType
                          , typename ProblemImp::DomainFieldType, ProblemImp::dimDomain
                          , typename ProblemImp::RangeFieldType, ProblemImp::dimRange > BaseType;
public:
  using typename BaseType::EntityType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::RangeFieldType;
  static const unsigned int dimDomain = BaseType::dimDomain;
  static const unsigned int dimRange = BaseType::dimRange;

  using typename BaseType::DiffusionFactorType;
  using typename BaseType::DiffusionTensorType;
  using typename BaseType::FunctionType;
private:
  typedef Stuff::Functions::Constant< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >
    ConstantFunctionType;
  typedef Pymor::Function::NonparametricDefault< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >
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
    : problem_(problem)
    , zero_(new ConstantFunctionType(0))
  {}

  virtual const DiffusionFactorType& diffusion_factor() const DS_OVERRIDE DS_FINAL
  {
    return problem_.diffusion_factor();
  }

  virtual const DiffusionTensorType& diffusion_tensor() const DS_OVERRIDE DS_FINAL
  {
    return problem_.diffusion_tensor();
  }

  virtual const FunctionType& force() const DS_OVERRIDE DS_FINAL
  {
    return problem_.force();
  }

  virtual const FunctionType& dirichlet() const DS_OVERRIDE DS_FINAL
  {
    return zero_;
  }

  virtual const FunctionType& neumann() const DS_OVERRIDE DS_FINAL
  {
    return zero_;
  }

private:
  const ProblemImp& problem_;
  const NonparametricFunctionType zero_;
}; // class ZeroBoundary


} // namespace Problems
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_PROBLEMS_ZERO_BOUNDARY_HH
