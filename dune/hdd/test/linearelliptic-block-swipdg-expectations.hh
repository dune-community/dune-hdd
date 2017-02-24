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
class OS2014;


#if HAVE_DUNE_GRID_MULTISCALE


template< class GridType >
class OS2014Multiscale;


namespace OS2015 {


template< class GridType >
class Academic;


} // namespace OS2015


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

using SPGRIDGridType = SPGrid<double, 2, SPIsotropicRefinement, MPI_Comm>;
typedef TestCases::OS2014Multiscale< SPGRIDGridType > SPGRIDGridTypeTestCaseType;
template< bool anything >
class BlockSWIPDGStudyExpectations< SPGRIDGridTypeTestCaseType, 1, anything >
    : public internal::BlockSWIPDGStudyExpectationsBase< SPGRIDGridTypeTestCaseType, 1 >
{
public:
  static std::vector< double > results(const SPGRIDGridTypeTestCaseType& test_case, const std::string type)
  {
    if (test_case.partitioning() == "[1 1 1]") {
      if (type == "L2")
        return {1.83e-02, 4.53e-03, 1.12e-03, 2.78e-04};
      else if (type == "H1_semi")
        return {2.77e-01, 1.39e-01, 6.98e-02, 3.50e-02};
      else if (type == "energy")
        return {2.77e-01, 1.39e-01, 6.98e-02, 3.50e-02};
      else if (type == "eta_NC_OS2014")
        return {1.66e-01, 7.89e-02, 3.91e-02, 1.95e-02};
      else if (type.substr(0, 12) == "eta_R_OS2014")
        return {5.79e-01, 2.90e-01, 1.45e-01, 7.27e-02};
      else if (type == "eta_DF_OS2014_*")
        return {3.55e-01, 1.76e-01, 8.73e-02, 4.35e-02};
      else if (type == "eta_DF_OS2014")
        return {3.55e-01, 1.76e-01, 8.73e-02, 4.35e-02};
      else if (type == "eta_OS2014")
        return {1.10e+00, 5.45e-01, 2.72e-01, 1.36e-01};
      else if (type == "eta_OS2014_*")
        return {1.10e+00, 5.45e-01, 2.72e-01, 1.36e-01};
      else if (type == "eff_OS2014")
        return {3.35e+00, 3.37e+00, 3.38e+00, 3.39e+00};
      else if (type == "eff_OS2014_*")
        return {3.35e+00, 3.37e+00, 3.38e+00, 3.39e+00};
      else
        EXPECT_TRUE(false) << "test results missing for type: " << type;
    } else if (test_case.partitioning() == "[2 2 1]") {
      if (type == "L2")
        return {1.83e-02, 4.53e-03, 1.12e-03, 2.78e-04};
      else if (type == "H1_semi")
        return {3.28e-01, 1.62e-01, 8.04e-02, 4.01e-02};
      else if (type == "energy")
        return {3.28e-01, 1.62e-01, 8.04e-02, 4.01e-02};
      else if (type == "eta_NC_OS2014")
        return {1.66e-01, 7.89e-02, 3.91e-02, 1.95e-02};
      else if (type.substr(0, 12) == "eta_R_OS2014")
        return {2.89e-01, 1.45e-01, 7.27e-02, 3.63e-02};
      else if (type == "eta_DF_OS2014")
        return {3.55e-01, 1.76e-01, 8.73e-02, 4.35e-02};
      else if (type == "eta_DF_OS2014_*")
        return {3.55e-01, 1.76e-01, 8.73e-02, 4.35e-02};
      else if (type == "eta_OS2014")
        return {8.10e-01, 4.00e-01, 1.99e-01, 9.94e-02};
      else if (type == "eta_OS2014_*")
        return {8.10e-01, 4.00e-01, 1.99e-01, 9.94e-02};
      else if (type == "eff_OS2014")
        return {2.47, 2.47, 2.48, 2.48};
      else if (type == "eff_OS2014_*")
        return {2.47, 2.47, 2.48, 2.48};
      else
        EXPECT_TRUE(false) << "test results missing for type: " << type;
    } else if (test_case.partitioning() == "[4 4 1]") {
      if (type == "L2")
        return {1.83e-02, 4.53e-03, 1.12e-03, 2.78e-04};
      else if (type == "H1_semi")
        return {3.28e-01, 1.62e-01, 8.04e-02, 4.01e-02};
      else if (type == "energy")
        return {3.28e-01, 1.62e-01, 8.04e-02, 4.01e-02};
      else if (type == "eta_NC_OS2014")
        return {1.66e-01, 7.89e-02, 3.91e-02, 1.95e-02};
      else if (type.substr(0, 12) == "eta_R_OS2014")
        return {1.45e-01, 7.26e-02, 3.63e-02, 1.82e-02};
      else if (type == "eta_DF_OS2014")
        return {3.55e-01, 1.76e-01, 8.73e-02, 4.35e-02};
      else if (type == "eta_DF_OS2014_*")
        return {3.55e-01, 1.76e-01, 8.73e-02, 4.35e-02};
      else if (type == "eta_OS2014")
        return {6.65e-01, 3.27e-01, 1.63e-01, 8.12e-02};
      else if (type == "eta_OS2014_*")
        return {6.65e-01, 3.27e-01, 1.63e-01, 8.12e-02};
      else if (type == "eff_OS2014")
        return {2.03, 2.02, 2.02, 2.03};
      else if (type == "eff_OS2014_*")
        return {2.03, 2.02, 2.02, 2.03};
      else
        EXPECT_TRUE(false) << "test results missing for type: " << type;
    } else if (test_case.partitioning() == "[8 8 1]") {
      if (type == "L2")
        return {1.83e-02, 4.53e-03, 1.12e-03, 2.78e-04};
      else if (type == "H1_semi")
        return {3.28e-01, 1.62e-01, 8.04e-02, 4.01e-02};
      else if (type == "energy")
        return {3.28e-01, 1.62e-01, 8.04e-02, 4.01e-02};
      else if (type == "eta_NC_OS2014")
        return {1.66e-01, 7.89e-02, 3.91e-02, 1.95e-02};
      else if (type.substr(0, 12) == "eta_R_OS2014")
        return {7.23e-02, 3.63e-02, 1.82e-02, 9.09e-03};
      else if (type == "eta_DF_OS2014")
        return {3.55e-01, 1.76e-01, 8.73e-02, 4.35e-02};
      else if (type == "eta_DF_OS2014_*")
        return {3.55e-01, 1.76e-01, 8.73e-02, 4.35e-02};
      else if (type == "eta_OS2014")
        return {5.93e-01, 2.91e-01, 1.45e-01, 7.21e-02};
      else if (type == "eta_OS2014_*")
        return {5.93e-01, 2.91e-01, 1.45e-01, 7.21e-02};
      else if (type == "eff_OS2014")
        return {1.81, 1.80, 1.80, 1.80};
      else if (type == "eff_OS2014_*")
        return {1.81, 1.80, 1.80, 1.80};
      else
        EXPECT_TRUE(false) << "test results missing for type: " << type;
    } else
      EXPECT_TRUE(false) << "test results missing for partitioning: " << test_case.partitioning();
    return {};
  } // ... results(...)
}; // BlockSWIPDGStudyExpectations
template class BlockSWIPDGStudyExpectations< SPGRIDGridTypeTestCaseType, 1 >;
} // namespace Tests
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_TEST_LINEARELLIPTIC_BLOCK_SWIPDG_EXPECTATIONS_HH
