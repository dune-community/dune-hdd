// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_INTERFACES_HH
#define DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_INTERFACES_HH

#include <memory>

#include <dune/common/exceptions.hh>

#include <dune/grid/part/interface.hh>

#include <dune/stuff/grid/boundaryinfo.hh>
//#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/parameter/tree.hh>

#include <dune/pymor/common/crtp.hh>
#include <dune/pymor/parameters/base.hh>
#include <dune/pymor/discretizations/interfaces.hh>

#include "../problems/interfaces.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {


template< class Traits >
class DiscretizationInterface
  : public Pymor::CRTPInterface< DiscretizationInterface< Traits >, Traits >
  , public Pymor::StationaryDiscretizationInterface< Traits >
{
  typedef Pymor::CRTPInterface< DiscretizationInterface< Traits >, Traits > CRTP;
  typedef Pymor::StationaryDiscretizationInterface< Traits >                BaseType;
public:
  typedef typename Traits::derived_type derived_type;

  typedef Dune::grid::Part::Interface< typename Traits::GridPartTraits > GridPartType;

  static const int polOrder = Traits::polOrder;

  typedef typename GridPartType::ctype    DomainFieldType;
  static const unsigned int               dimDomain = GridPartType::dimension;
  typedef typename Traits::RangeFieldType RangeFieldType;
  static const unsigned int               dimRange = Traits::dimRange;

  typedef typename Traits::BoundaryInfoType BoundaryInfoType;
  typedef typename Traits::ProblemType      ProblemType;

  typedef Dune::Stuff::Common::ExtendedParameterTree SettingsType;

  DiscretizationInterface(const Pymor::ParameterType tt = Pymor::ParameterType())
    : BaseType(tt)
  {}

  DiscretizationInterface(const Pymor::Parametric& other)
    : BaseType(other)
  {}

  static const std::string static_id()
  {
    return "dune.hdd.linearelliptic.discretization";
  }

  std::shared_ptr< const GridPartType > gridPart() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(CRTP::as_imp(*this).gridPart());
    return CRTP::as_imp(*this).gridPart();
  }

  std::shared_ptr< const BoundaryInfoType > boundaryInfo() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(CRTP::as_imp(*this).boundaryInfo());
    return CRTP::as_imp(*this).boundaryInfo();
  }

  const ProblemType& problem() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(CRTP::as_imp(*this).model());
    return CRTP::as_imp(*this).model();
  }

  void init(std::ostream& out = Dune::Stuff::Common::Logger().devnull(),
            const std::string prefix = "")
  {
    CHECK_INTERFACE_IMPLEMENTATION(CRTP::as_imp(*this).init(out, prefix));
    CRTP::as_imp(*this).init(out, prefix);
  }

  bool initialized() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(CRTP::as_imp(*this).initialized());
    return CRTP::as_imp(*this).initialized();
  }
}; // class DiscretizationInterface


//template< class Traits >
//class SolverParametricInterface
//{
//public:
//  typedef SolverParametricInterface< Traits > ThisType;
//  typedef typename Traits::derived_type       derived_type;

//  typedef Dune::grid::Part::Interface< typename Traits::GridPartTraits > GridPartType;

//  static const int polOrder = Traits::polOrder;

//  typedef typename GridPartType::ctype    DomainFieldType;
//  static const int                        dimDomain = GridPartType::dimension;
//  typedef typename Traits::RangeFieldType RangeFieldType;
//  static const int                        dimRange = Traits::dimRange;

//  typedef Dune::Stuff::GridboundaryInterface< typename GridPartType::GridViewType > BoundaryInfoType;
//  typedef ModelInterface< DomainFieldType, dimDomain, RangeFieldType, dimRange >            ModelType;

//  typedef typename ModelType::ParamType ParamType;

//  typedef typename Traits::VectorType VectorType;

//  static const std::string id()
//  {
//    return "solver.linearelliptic";
//  }

//  std::shared_ptr< const GridPartType > gridPart() const
//  {
//    CHECK_INTERFACE_IMPLEMENTATION(CRTP::as_imp(*this).gridPart());
//    return CRTP::as_imp(*this).gridPart();
//  }

//  std::shared_ptr< const BoundaryInfoType > boundaryInfo() const
//  {
//    CHECK_INTERFACE_IMPLEMENTATION(CRTP::as_imp(*this).boundaryInfo());
//    return CRTP::as_imp(*this).boundaryInfo();
//  }

//  std::shared_ptr< const ModelType > model() const
//  {
//    CHECK_INTERFACE_IMPLEMENTATION(CRTP::as_imp(*this).model());
//    return CRTP::as_imp(*this).model();
//  }

//  std::shared_ptr< VectorType > createVector() const
//  {
//    CHECK_INTERFACE_IMPLEMENTATION(CRTP::as_imp(*this).createVector());
//    return CRTP::as_imp(*this).createVector();
//  }

//  void visualize(const std::shared_ptr< const VectorType > vector,
//                 const std::string filename,
//                 const std::string name,
//                 std::ostream& out = Dune::Stuff::Common::Logger().devnull(),
//                 const std::string prefix = "") const
//  {
//    CHECK_INTERFACE_IMPLEMENTATION(CRTP::as_imp(*this).visualize(vector, filename, name, out, prefix));
//    CRTP::as_imp(*this).visualize(vector, filename, name, out, prefix);
//  }

//  void init(std::ostream& out = Dune::Stuff::Common::Logger().devnull(),
//            const std::string prefix = "")
//  {
//    CHECK_INTERFACE_IMPLEMENTATION(CRTP::as_imp(*this).init(out, prefix));
//    CRTP::as_imp(*this).init(out, prefix);
//  }

//  bool initialized() const
//  {
//    CHECK_INTERFACE_IMPLEMENTATION(CRTP::as_imp(*this).initialized());
//    return CRTP::as_imp(*this).initialized();
//  }

//  void solve(std::shared_ptr< VectorType > solution,
//             const ParamType& mu,
//             const std::string linearSolverType,
//             const double linearSolverPrecision,
//             const size_t linearSolverMaxIterations,
//             std::ostream& out = Dune::Stuff::Common::Logger().devnull(),
//             const std::string prefix = "") const
//  {
//    if (!initialized())
//      DUNE_THROW(Dune::InvalidStateException,
//                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
//                 << " call init() before calling solve()!");
//    CHECK_INTERFACE_IMPLEMENTATION(CRTP::as_imp(*this).solve(solution,
//                                                 mu,
//                                                 linearSolverType,
//                                                 linearSolverPrecision,
//                                                 linearSolverMaxIterations,
//                                                 out,
//                                                 prefix));
//    CRTP::as_imp(*this).solve(solution, mu, linearSolverType, linearSolverPrecision, linearSolverMaxIterations, out, prefix);
//  }

//  derived_type& CRTP::as_imp(*this)
//  {
//    return static_cast< derived_type& >(*this);
//  }

//  const derived_type& CRTP::as_imp(*this) const
//  {
//    return static_cast< const derived_type& >(*this);
//  }
//}; // class SolverParametricInterface< DomainFieldType, dimDomain, RangeFieldType, 1 >


} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_INTERFACES_HH
