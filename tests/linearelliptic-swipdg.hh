// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_TEST_LINEARELLIPTIC_SWIPDG_HH
#define DUNE_HDD_TEST_LINEARELLIPTIC_SWIPDG_HH

#include <dune/stuff/common/exceptions.hh>

#include <dune/gdt/products/l2.hh>
#include <dune/gdt/products/h1.hh>

#include <dune/hdd/linearelliptic/discretizations/swipdg.hh>
#include <dune/hdd/playground/linearelliptic/testcases/ESV2007.hh>

#include "linearelliptic.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Tests {
namespace internal {


template< class TestCaseType, int polOrder, GDT::ChooseSpaceBackend space_backend, Stuff::LA::ChooseBackend la_backend >
class Discretization
{
  typedef typename TestCaseType::GridType GridType;
  typedef typename TestCaseType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange = TestCaseType::dimRange;
public:
  typedef Discretizations::SWIPDG< GridType, Stuff::Grid::ChooseLayer::level, RangeFieldType, dimRange, polOrder > Type;
}; // class Discretization


} // namespace internal


template< class TestCaseType, int polOrder, GDT::ChooseSpaceBackend space_backend, Stuff::LA::ChooseBackend la_backend >
class EocStudySWIPDG
  : public EocStudyBase< TestCaseType, typename internal::Discretization< TestCaseType, polOrder, space_backend, la_backend >::Type >
{
  typedef EocStudyBase
      < TestCaseType, typename internal::Discretization< TestCaseType, polOrder, space_backend, la_backend >::Type >
    BaseType;

  typedef typename BaseType::DiscretizationType DiscretizationType;
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::FunctionType FunctionType;
  typedef typename BaseType::VectorType   VectorType;

public:
  EocStudySWIPDG(const TestCaseType& test_case)
    : BaseType(test_case)
  {}

  virtual ~EocStudySWIPDG() {}

  virtual std::string identifier() const DS_OVERRIDE DS_FINAL
  {
    return DiscretizationType::static_id() + " (polorder " + Stuff::Common::toString(polOrder) + ")";
  }

  virtual std::vector< std::string > available_norms() const DS_OVERRIDE DS_FINAL
  {
    return {"L2", "H1_semi"};
  }

  virtual std::vector< std::string > available_estimators() const DS_OVERRIDE DS_FINAL
  {
    return {};
  }

  virtual size_t expected_rate(const std::string type) const DS_OVERRIDE DS_FINAL
  {
    if (type.compare("L2") == 0)
      return polOrder + 1;
    else if (type.compare("H1_semi") == 0)
      return polOrder;
    else
      DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Wrong type '" << type << "' requested!");
  } // ... expected_rate(...)

  std::vector< double > expected_results(const std::string type) const
  {
    if (std::is_same< TestCaseType, TestCases::ESV2007< ALUConformGrid< 2, 2 > > >::value
        || std::is_same< TestCaseType, TestCases::ESV2007< ALUGrid< 2, 2, simplex, conforming > > >::value) {
      if (polOrder == 1) {
        if (type.compare("L2") == 0)
          return {1.15e-01, 3.04e-02, 7.51e-03, 1.86e-03};
        else if (type.compare("H1_semi") == 0)
          return {3.79e-01, 1.90e-01, 9.38e-02, 4.67e-02};
        else
          DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Wrong type '" << type << "' requested!");
      } else if (polOrder == 2) {
        if (type.compare("L2") == 0)
          return {1.25e-02, 1.42e-03, 1.69e-04, 2.08e-05};
        else if (type.compare("H1_semi") == 0)
          return {7.84e-02, 2.01e-02, 5.02e-03, 1.26e-03};
        else
          DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Wrong type '" << type << "' requested!");
      } else
        DUNE_THROW(NotImplemented, "Please record the expected results for this polOrder!");
//    } else if (std::is_same< TestCaseType, EllipticTestCaseType::LocalThermalBlock< Dune::ALUConformGrid< 2, 2 > > >::value) {
//      if (polOrder == 1) {
//        if (type.compare("L2") == 0)
//          return {5.85e-02, 2.00e-02, 5.55e-03, 1.30e-03};
//        else if (type.compare("H1_semi") == 0)
//          return {4.33e-01, 2.94e-01, 1.51e-01, 6.55e-02};
//        else
//          DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Wrong type '" << type << "' requested!");
//      } else if (polOrder == 2) {
//        if (type.compare("L2") == 0)
//          return {1.19e-02, 2.12e-03, 3.90e-04, 7.77e-05};
//        else if (type.compare("H1_semi") == 0)
//          return {1.70e-01, 5.97e-02, 1.95e-02, 6.05e-03};
//        else
//          DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Wrong type '" << type << "' requested!");
//      } else
//        DUNE_THROW(NotImplemented, "Please record the expected results for this polOrder!");
//    } else if (std::is_same< TestCaseType, EllipticTestCaseType::ER07< Dune::ALUConformGrid< 2, 2 > > >::value) {
//      if (polOrder == 1) {
//        if (type.compare("L2") == 0)
//          return {6.10e-02, 1.66e-02, 4.23e-03};
//        else if (type.compare("H1_semi") == 0)
//          return {2.99e-01, 1.47e-01, 7.26e-02};
//        else
//          DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Wrong type '" << type << "' requested!");
//      } else if (polOrder == 2) {
//        if (type.compare("L2") == 0)
//          return {6.43e-03, 8.24e-04, 1.05e-04};
//        else if (type.compare("H1_semi") == 0)
//          return {5.41e-02, 1.42e-02, 3.56e-03};
//        else
//          DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Wrong type '" << type << "' requested!");
//      } else
//        DUNE_THROW(NotImplemented, "Please record the expected results for this polOrder!");
//    } else if (std::is_same< TestCaseType, EllipticTestCaseType::MixedBoundaryTypes< Dune::ALUConformGrid< 2, 2 > > >::value) {
//      if (polOrder == 1) {
//        if (type.compare("L2") == 0)
//          return {4.03e-02, 1.13e-02, 2.84e-03, 6.34e-04};
//        else if (type.compare("H1_semi") == 0)
//          return {2.70e-01, 1.40e-01, 6.88e-02, 3.09e-02};
//        else
//          DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Wrong type '" << type << "' requested!");
//      } else if (polOrder == 2) {
//        if (type.compare("L2") == 0)
//          return {3.59e-03, 6.26e-04, 1.22e-04, 2.69e-05};
//        else if (type.compare("H1_semi") == 0)
//          return {4.82e-02, 1.80e-02, 7.20e-03, 2.86e-03};
//        else
//          DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Wrong type '" << type << "' requested!");
//      } else
//        DUNE_THROW(NotImplemented, "Please record the expected results for this polOrder!");
//    } else if (std::is_same< TestCaseType, EllipticTestCaseType::Spe10Model1< Dune::ALUConformGrid< 2, 2 > > >::value) {
//      if (polOrder == 1) {
//        if (type.compare("L2") == 0)
//          return {7.23e-02, 2.60e-02};
//        else if (type.compare("H1_semi") == 0)
//          return {5.29e-01, 3.49e-01};
//        else
//          DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Wrong type '" << type << "' requested!");
//      } else if (polOrder == 2) {
//        if (type.compare("L2") == 0)
//          return {2.09e-02, 3.75e-03};
//        else if (type.compare("H1_semi") == 0)
//          return {2.57e-01, 8.48e-02};
//        else
//          DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Wrong type '" << type << "' requested!");
//      } else
//        DUNE_THROW(NotImplemented, "Please record the expected results for this polOrder!");
    } else
      DUNE_THROW(NotImplemented, "Please record the expected results for this TestCaseType/GridType combination!");
  } // ... expected_results(...)

private:
  virtual double compute_norm_(const GridViewType& grid_view,
                               const FunctionType& function,
                               const std::string type) const
  {
    using namespace GDT;
    if (type == "L2") {
      return Products::L2< GridViewType >(grid_view).induced_norm(function);
    } else if (type == "H1_semi") {
      return Products::H1SemiGeneric< GridViewType >(grid_view).induced_norm(function);
    } else
      DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Wrong type '" << type << "' requested!");
  } // ... compute_norm_(...)

  virtual double estimate_(const VectorType& /*vector*/, const std::string /*type*/) const
  {
    DUNE_THROW(NotImplemented, "");
  }
}; // class EocStudySWIPDG


} // namespace Tests
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_TEST_LINEARELLIPTIC_SWIPDG_HH
