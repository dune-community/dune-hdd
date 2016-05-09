// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_TEST_LINEARELLIPTIC_SWIPDG_EXPECTATIONS_HH
#define DUNE_HDD_TEST_LINEARELLIPTIC_SWIPDG_EXPECTATIONS_HH

#include <dune/stuff/common/disable_warnings.hh>
# if HAVE_ALUGRID
#   include <dune/grid/alugrid.hh>
# endif
# include <dune/grid/sgrid.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/stuff/test/gtest/gtest.h>
#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/type_utils.hh>

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Tests {
namespace internal {


template< class TestCaseType, int polOrder >
class SWIPDGStudyExpectationsBase
{
public:
  static size_t rate(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "L2")
      return polOrder + 1;
    else if (type == "H1_semi")
      return polOrder;
    else if (type == "energy")
      return polOrder;
    else if (type == "eta_NC_ESV2007")
      return polOrder;
    else if (type.substr(0, 14) == "eta_R_ESV2007")
      return polOrder + 1;
    else if (type == "eta_DF_ESV2007")
      return polOrder;
    else if (type == "eta_ESV2007")
      return polOrder;
    else if (type == "eff_ESV2007")
      return 0;
    else if (type == "eta_ESV2007_alt")
      return polOrder;
    else if (type == "eff_ESV2007_alt")
      return 0;
    else
      EXPECT_TRUE(false) << "expected rate missing for type: " << type;
    return 0;
  } // ... rate(...)
}; // class SWIPDGStudyExpectationsBase


} // namespace internal


template< class TestCaseType, int polOrder, bool anything = true >
class SWIPDGStudyExpectations
  : public internal::SWIPDGStudyExpectationsBase< TestCaseType, polOrder >
{
public:
  static std::vector< double > results(const TestCaseType& /*test_case*/, const std::string type)
  {
    EXPECT_TRUE(false) << "Please record the expected results for\n"
                       << "TestCaseType: " << Stuff::Common::Typename< TestCaseType >::value() << "\n"
                       << "polOrder: " << polOrder << "\n"
                       << "type: " << type << "\n"
                       << "Please put an appropriate specialiaztion of SWIPDGStudyExpectations for this TestCaseType "
                       << "in a separate object file (see examples below) or add\n"
                       << "  'template class SWIPDGStudyExpectations< TestCaseType, " << polOrder << " >'\n"
                       << "for this polOrder in the appropriate object file!\n\n"
                       << "Oh: and do not forget to add\n"
                       << "  'extern template class SWIPDGStudyExpectations< ... >'\n"
                       << "to each test source using these results!";
    return {};
  } // ... results(...)
}; // SWIPDGStudyExpectations


} // namespace Tests
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_TEST_LINEARELLIPTIC_SWIPDG_EXPECTATIONS_HH
