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
          return {};
        else if (type == "eta_NC_OS2014")
          return {};
        else if (type == "eta_R_OS2014_*")
          return {};
        else if (type == "eta_DF_OS2014_*")
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
      if (mu == 1 && mu_bar == 1 && mu_hat == 0.1) {
        if (type == "eta_OS2014_*_mu")
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
