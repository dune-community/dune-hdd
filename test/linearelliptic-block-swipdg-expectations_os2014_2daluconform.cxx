// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#if HAVE_ALUGRID
# include <dune/grid/alugrid.hh>

# include <dune/hdd/linearelliptic/testcases/OS2014.hh>

# include "linearelliptic-block-swipdg-expectations.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {


typedef TestCases::OS2014::ParametricBlockConvergence< ALUGrid< 2, 2, simplex, conforming > > TestCaseType;

namespace Tests {


template< bool anything >
class BlockSWIPDGStudyExpectations< TestCaseType, 1, anything >
  : public internal::BlockSWIPDGStudyExpectationsBase< TestCaseType, 1 >
{
public:
  static std::vector< double > results(const TestCaseType& test_case, const std::string type)
  {
    if (test_case.num_refinements() != 3) {
      EXPECT_TRUE(false) << "test results missing for num_refinements: " << test_case.num_refinements();
      return {};
    }
    const auto mu = test_case.parameters().at("mu");
    const auto mu_bar = test_case.parameters().at("mu_bar");
    const auto mu_hat = test_case.parameters().at("mu_hat");
    if (test_case.partitioning() == "[1 1 1]") {
      if (mu == 0.1 && mu_bar == 0.1 && mu_hat == 1) {
        if (type == "eta_DF_OS2014")
          return {};
        else if (type == "eta_DF_OS2014_*")
          return {};
        else if (type == "eta_OS2014")
          return {};
        else if (type == "eta_OS2014_*")
          return {};
        else if (type == "eff_OS2014_*_mu")
          return {};
        else if (type == "eff_OS2014_mu")
          return {};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else if (mu == 1 && mu_bar == 1 && mu_hat == 1) {
        if (type == "energy_mu")
          return {3.28e-01, 1.60e-01, 7.78e-02, 3.47e-02};
        else if (type == "eta_NC_OS2014")
          return {1.66e-01, 7.89e-02, 3.91e-02, 1.95e-02};
        else if (type == "eta_R_OS2014_*")
          return {5.79e-01, 2.90e-01, 1.45e-01, 7.27e-02};
        else if (type == "eta_DF_OS2014_*")
          return {3.55e-01, 1.76e-01, 8.73e-02, 4.35e-02};
        else if (type == "eta_OS2014_*")
          return {1.10e+00, 5.45e-01, 2.72e-01, 1.36e-01};
        else if (type == "eff_OS2014_*_mu")
          return {3.36e+00, 3.40e+00, 3.49e+00, 3.91e+00};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else
        EXPECT_TRUE(false) << "test results missing for parameters: mu     = " << mu << "\n"
                           << "                                     mu_bar = " << mu_bar << "\n"
                           << "                                     mu_hat = " << mu_hat;
    } else if (test_case.partitioning() == "[2 2 1]") {
      if (mu == 1 && mu_bar == 1 && mu_hat == 1) {
        if (type == "eta_R_OS2014_*")
          return {};
        else if (type == "eta_OS2014_*")
          return {};
        else if (type == "eff_OS2014_*_mu")
          return {};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else
        EXPECT_TRUE(false) << "test results missing for parameters: mu     = " << mu << "\n"
                           << "                                     mu_bar = " << mu_bar << "\n"
                           << "                                     mu_hat = " << mu_hat;
    } else if (test_case.partitioning() == "[2 2 1]_H_with_h") {
      if (mu == 1 && mu_bar == 0.1 && mu_hat == 0.1) {
        if (type == "energy_mu_bar")
          return {3.81e-01, 1.87e-01, 9.08e-02, 4.05e-02};
        else if (type == "eta_NC_OS2014")
          return {1.82e-01, 8.57e-02, 4.22e-02, 2.11e-02};
        else if (type == "eta_R_OS2014_*")
          return {2.89e-01, 7.26e-02, 1.82e-02, 4.54e-03};
        else if (type == "eta_DF_OS2014_*")
          return {3.16e-01, 1.56e-01, 7.74e-02, 3.85e-02};
        else if (type == "eta_OS2014_*")
          return {1.18e+00, 5.00e-01, 2.29e-01, 1.10e-01};
        else if (type == "eff_OS2014_*_mu_bar")
          return {3.10e+00, 2.67e+00, 2.52e+00, 2.71e+00};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else if (mu == 1 && mu_bar == 1 && mu_hat == 0.1) {
        if (type == "eta_DF_OS2014_*")
          return {3.16e-01, 1.56e-01, 7.74e-02, 3.85e-02};
        else if (type == "eta_OS2014_*")
          return {7.71e-01, 3.08e-01, 1.35e-01, 6.26e-02};
        else if (type == "eff_OS2014_*_mu")
          return {2.35e+00, 1.92e+00, 1.73e+00, 1.80e+00};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else if (mu == 1 && mu_bar == 1 && mu_hat == 1) {
        if (type == "eta_R_OS2014_*")
          return {2.89e-01, 7.26e-02, 1.82e-02, 4.54e-03};
        else if (type == "eta_OS2014_*")
          return {8.10e-01, 3.27e-01, 1.45e-01, 6.76e-02};
        else if (type == "eff_OS2014_*_mu")
          return {2.47e+00, 2.04e+00, 1.86e+00, 1.95e+00};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else
        EXPECT_TRUE(false) << "test results missing for parameters: mu     = " << mu << "\n"
                           << "                                     mu_bar = " << mu_bar << "\n"
                           << "                                     mu_hat = " << mu_hat;
    } else if (test_case.partitioning() == "[4 4 1]") {
      if (mu == 0.1 && mu_bar == 0.1 && mu_hat == 0.1) {
        if (type == "eff_OS2014_*_mu")
          return {};
        else if (type == "eff_OS2014_mu")
          return {};
        else if (type == "eta_DF_OS2014")
          return {};
        else if (type == "eta_DF_OS2014_*")
          return {};
        else if (type == "eta_OS2014")
          return {};
        else if (type == "eta_OS2014_*")
          return {};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else if (mu == 1 && mu_bar == 1 && mu_hat == 0.1) {
        if (type == "eff_OS2014_*_mu")
          return {};
        else if (type == "eff_OS2014_mu")
          return {};
        else if (type == "eta_DF_OS2014")
          return {};
        else if (type == "eta_DF_OS2014_*")
          return {};
        else if (type == "eta_OS2014")
          return {};
        else if (type == "eta_OS2014_*")
          return {};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else if (mu == 0.1 && mu_bar == 0.1 && mu_hat == 1) {
        if (type == "eff_OS2014_*_mu")
          return {};
        else if (type == "eff_OS2014_mu")
          return {};
        else if (type == "eta_DF_OS2014")
          return {};
        else if (type == "eta_DF_OS2014_*")
          return {};
        else if (type == "eta_OS2014")
          return {};
        else if (type == "eta_OS2014_*")
          return {};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else if (mu == 1 && mu_bar == 1 && mu_hat == 1) {
        if (type == "eta_R_OS2014_*")
          return {};
        else if (type == "eta_OS2014_*")
          return {};
        else if (type == "eff_OS2014_*_mu")
          return {};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else
        EXPECT_TRUE(false) << "test results missing for parameters: mu     = " << mu << "\n"
                           << "                                     mu_bar = " << mu_bar << "\n"
                           << "                                     mu_hat = " << mu_hat;
    } else if (test_case.partitioning() == "[8 8 1]") {
      if (mu == 1 && mu_bar == 1 && mu_hat == 1) {
        if (type == "eta_R_OS2014_*")
          return {};
        else if (type == "eta_OS2014_*")
          return {};
        else if (type == "eff_OS2014_*_mu")
          return {};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else
        EXPECT_TRUE(false) << "test results missing for parameters: mu     = " << mu << "\n"
                           << "                                     mu_bar = " << mu_bar << "\n"
                           << "                                     mu_hat = " << mu_hat;
    } else
      EXPECT_TRUE(false) << "test results missing for partitioning: " << test_case.partitioning();
    return {};
  } // ... results(...)
}; // BlockSWIPDGStudyExpectations


template class BlockSWIPDGStudyExpectations< TestCaseType, 1 >;


} // namespace Tests
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // HAVE_ALUGRID
