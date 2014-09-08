// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include <dune/stuff/test/gtest/gtest.h>

#if HAVE_ALUGRID

#include <dune/hdd/playground/linearelliptic/testcases/ESV2007.hh>

#include "linearelliptic-block-swipdg-expectations.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Tests {


template< bool anything >
class BlockSWIPDGStudyExpectations< TestCases::ESV2007Multiscale< ALUGrid< 2, 2, simplex, conforming > >, 1, anything >
  : public internal::BlockSWIPDGStudyExpectationsBase
        < TestCases::ESV2007Multiscale< ALUGrid< 2, 2, simplex, conforming > >, 1 >
{
  typedef TestCases::ESV2007Multiscale< ALUGrid< 2, 2, simplex, conforming > > TestCaseType;

public:
  static std::vector< double > results(const TestCaseType& test_case, const std::string type)
  {
    if (test_case.partitioning() == "[1 1 1]") {
      if (type == "energy")
        return {3.28e-01, 1.62e-01, 8.04e-02, 4.01e-02};
      else if (type == "eta_NC_OS2014")
        return {1.66e-01, 7.89e-02, 3.91e-02, 1.95e-02};
      else if (type == "eta_R_OS2014")
        return {5.79e-01, 2.90e-01, 1.45e-01, 7.27e-02};
      else if (type == "eta_DF_OS2014")
        return {3.55e-01, 1.76e-01, 8.73e-02, 4.35e-02};
      else if (type == "eta_OS2014")
        return {1.10e+00, 5.45e-01, 2.72e-01, 1.36e-01};
      else if (type == "eff_OS2014")
        return {3.35e+00, 3.37e+00, 3.38e+00, 3.39e+00};
      else
        EXPECT_TRUE(false) << "test results missing for type: " << type;
    } else if (test_case.partitioning() == "[2 2 1]") {
      if (type == "energy")
        return {3.28e-01, 1.62e-01, 8.04e-02, 4.01e-02};
      else if (type == "eta_NC_OS2014")
        return {1.66e-01, 7.89e-02, 3.91e-02, 1.95e-02};
      else if (type == "eta_R_OS2014")
        return {2.89e-01, 1.45e-01, 7.27e-02, 3.63e-02};
      else if (type == "eta_DF_OS2014")
        return {3.55e-01, 1.76e-01, 8.73e-02, 4.35e-02};
      else if (type == "eta_OS2014")
        return {8.10e-01, 4.00e-01, 1.99e-01, 9.94e-02};
      else if (type == "eff_OS2014")
        return {2.47, 2.47, 2.48, 2.48};
      else
        EXPECT_TRUE(false) << "test results missing for type: " << type;
    } else if (test_case.partitioning() == "[4 4 1]") {
      if (type == "energy")
        return {3.28e-01, 1.62e-01, 8.04e-02, 4.01e-02};
      else if (type == "eta_NC_OS2014")
        return {1.66e-01, 7.89e-02, 3.91e-02, 1.95e-02};
      else if (type == "eta_R_OS2014")
        return {1.45e-01, 7.26e-02, 3.63e-02, 1.82e-02};
      else if (type == "eta_DF_OS2014")
        return {3.55e-01, 1.76e-01, 8.73e-02, 4.35e-02};
      else if (type == "eta_OS2014")
        return {6.65e-01, 3.27e-01, 1.63e-01, 8.12e-02};
      else if (type == "eff_OS2014")
        return {2.03, 2.02, 2.02, 2.03};
      else
        EXPECT_TRUE(false) << "test results missing for type: " << type;
    } else if (test_case.partitioning() == "[8 8 1]") {
      if (type == "energy")
        return {3.28e-01, 1.62e-01, 8.04e-02, 4.01e-02};
      else if (type == "eta_NC_OS2014")
        return {1.66e-01, 7.89e-02, 3.91e-02, 1.95e-02};
      else if (type == "eta_R_OS2014")
        return {7.23e-02, 3.63e-02, 1.82e-02, 9.09e-03};
      else if (type == "eta_DF_OS2014")
        return {3.55e-01, 1.76e-01, 8.73e-02, 4.35e-02};
      else if (type == "eta_OS2014")
        return {5.93e-01, 2.91e-01, 1.45e-01, 7.21e-02};
      else if (type == "eff_OS2014")
        return {1.81, 1.80, 1.80, 1.80};
      else
        EXPECT_TRUE(false) << "test results missing for type: " << type;
    } else
      EXPECT_TRUE(false) << "test results missing for partitioning: " << test_case.partitioning();
    return {};
  } // ... results(...)
}; // BlockSWIPDGStudyExpectations


template class BlockSWIPDGStudyExpectations< TestCases::ESV2007Multiscale< ALUGrid< 2, 2, simplex, conforming > >, 1 >;


} // namespace Tests
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // HAVE_ALUGRID
