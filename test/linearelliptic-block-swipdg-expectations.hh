// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_TEST_LINEARELLIPTIC_BLOCK_SWIPDG_EXPECTATIONS_HH
#define DUNE_HDD_TEST_LINEARELLIPTIC_BLOCK_SWIPDG_EXPECTATIONS_HH

#include <dune/stuff/common/disable_warnings.hh>
# if HAVE_ALUGRID
#   include <dune/grid/alugrid.hh>
# endif
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/stuff/test/gtest/gtest.h>
#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/type_utils.hh>

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace TestCases {


// forwards
template< class GridType >
class ESV2007;


namespace OS2014 {


template< class GridType >
class ParametricConvergence;


} // namespace OS2014


#if HAVE_DUNE_GRID_MULTISCALE


template< class GridType >
class ESV2007Multiscale;


namespace OS2014 {


template< class GridType >
class ParametricBlockConvergence;


} // namespace OS2014


#endif // HAVE_DUNE_GRID_MULTISCALE


} // namespace TestCases
namespace Tests {
namespace internal {


template< class TestCaseType, int polOrder >
class BlockSWIPDGStudyExpectationsBase
{
public:
  static size_t rate(const TestCaseType& test_case, const std::string type)
  {
    const auto partitioning = test_case.partitioning();
    if (type == "L2")
      return polOrder + 1;
    else if (type == "H1_semi" || type.substr(0, 6) == "energy")
      return polOrder;
    else if (type == "eta_NC_OS2014")
      return polOrder;
    else if (type.substr(0, 12) == "eta_R_OS2014") {
      if (partitioning.size() >= 8 && partitioning.substr(partitioning.size() - 8) == "H_with_h")
        return polOrder + 1;
      else
        return polOrder;
    } else if (type.substr(0, 13) == "eta_DF_OS2014")
      return polOrder;
    else if (type.substr(0, 10) == "eta_OS2014")
      return polOrder;
    else if (type.substr(0, 10) == "eff_OS2014")
      return 0;
    else
      EXPECT_TRUE(false) << "expected rate missing for type: " << type;
    return 0;
  } // ... rate(...)
}; // class BlockSWIPDGStudyExpectationsBase


} // namespace internal


template< class TestCaseType, int polOrder = 1, bool anything = true >
class BlockSWIPDGStudyExpectations
  : public internal::BlockSWIPDGStudyExpectationsBase< TestCaseType, polOrder >
{
public:
  static std::vector< double > results(const TestCaseType& /*test_case*/, const std::string type)
  {
    EXPECT_TRUE(false) << "Please record the expected results for\n"
                       << "TestCaseType: " << Stuff::Common::Typename< TestCaseType >::value() << "\n"
                       << "polOrder: " << polOrder << "\n"
                       << "type: " << type << "\n"
                       << "Please put an appropriate specialiaztion of BlockSWIPDGStudyExpectations for this TestCaseType "
                       << "in a separate object file (see examples below) or add\n"
                       << "  'template class BlockSWIPDGStudyExpectations< TestCaseType, " << polOrder << " >'\n"
                       << "for this polOrder in the appropriate object file!\n\n"
                       << "Oh: and do not forget to add\n"
                       << "  'extern template class BlockSWIPDGStudyExpectations< ... >'\n"
                       << "to each test source using these results!";
    return {};
  } // ... results(...)
}; // BlockSWIPDGStudyExpectations


} // namespace Tests
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_TEST_LINEARELLIPTIC_BLOCK_SWIPDG_EXPECTATIONS_HH
