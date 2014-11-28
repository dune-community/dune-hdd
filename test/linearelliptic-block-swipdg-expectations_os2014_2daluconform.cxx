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
      if (   mu            == 0.1
          && mu_bar        == 0.1
          && mu_hat        == 0.1) {
        if (type == "energy_mu")
          return {};
        else if (type == "eta_DF_OS2014")
          return {};
        else if (type == "eta_DF_OS2014_*")
          return {};
        else if (type == "eta_OS2014")
          return {};
        else if (type == "eta_OS2014_*")
          return {};
        else if (type == "eff_OS2014_mu")
          return {};
        else if (type == "eff_OS2014_*_mu")
          return {};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else if (   mu            == 0.3
                 && mu_bar        == 0.3
                 && mu_hat        == 0.1) {
        if (type == "energy_mu")
          return {};
        else if (type == "eta_DF_OS2014")
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
      } else if (   mu            == 0.5
                 && mu_bar        == 0.5
                 && mu_hat        == 0.1) {
        if (type == "energy_mu")
          return {};
        else if (type == "eta_DF_OS2014")
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
      } else if (   mu            == 0.1
                 && mu_bar        == 0.1
                 && mu_hat        == 1) {
        if (type == "eta_DF_OS2014")
          return {1.01e+00, 1.21e+00, 1.35e+00, 1.41e+00};
        else if (type == "eta_DF_OS2014_*")
          return {1.16e+00, 6.90e-01, 3.34e-01, 1.62e-01};
        else if (type == "eta_OS2014")
          return {4.67e+00, 4.65e+00, 4.67e+00, 4.64e+00};
        else if (type == "eta_OS2014_*")
          return {5.15e+00, 3.01e+00, 1.45e+00, 6.97e-01};
        else if (type == "eff_OS2014_*_mu")
          return {5.86e+00, 5.65e+00, 5.77e+00, 6.41e+00};
        else if (type == "eff_OS2014_mu")
          return {5.31e+00, 8.74e+00, 1.86e+01, 4.27e+01};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
//      } else if (   mu            == Parameter("mu", )
//                 && mu_bar        == Parameter("mu", )
//                 && mu_hat        == Parameter("mu", )) {
//        if (type == "energy_mu")
//          return {};
//        else if (type == "eta_DF_OS2014")
//          return {};
//        else if (type == "eta_DF_OS2014_*")
//          return {};
//        else if (type == "eta_OS2014")
//          return {};
//        else if (type == "eta_OS2014_*")
//          return {};
//        else if (type == "eff_OS2014_*_mu")
//          return {};
//        else if (type == "eff_OS2014_mu")
//          return {};
//        else
//          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else if (   mu            == 1
                 && mu_bar        == 1
                 && mu_hat        == 1) {
        if (type == "energy_mu")
          return {};
        else if (type == "eta_NC_OS2014")
          return {};
        else if (type == "eta_R_OS2014")
          return {};
        else if (type == "eta_DF_OS2014")
          return {};
        else if (type == "eta_OS2014")
          return {};
        else if (type == "eff_OS2014_mu")
          return {};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else
        EXPECT_TRUE(false) << "test results missing for parameters: mu = " << mu << "\n"
                           << "                                     mu_bar = " << mu_bar << "\n"
                           << "                                     mu_hat = " << mu_hat;
    } else if (test_case.partitioning() == "[4 4 1]") {
      if (mu == 0.1 && mu_bar == 0.1 && mu_hat == 0.1) {
        if (type == "eff_OS2014_*_mu")
          return {2.24e+00, 2.22e+00, 2.27e+00, 2.49e+00};
        else if (type == "eff_OS2014_mu")
          return {2.24e+00, 2.22e+00, 2.27e+00, 2.49e+00};
        else if (type == "eta_DF_OS2014")
          return {1.25e+00, 7.37e-01, 3.69e-01, 1.83e-01};
        else if (type == "eta_DF_OS2014_*")
          return {1.25e+00, 7.37e-01, 3.69e-01, 1.83e-01};
        else if (type == "eta_OS2014")
          return {1.97e+00, 1.18e+00, 5.71e-01, 2.71e-01};
        else if (type == "eta_OS2014_*")
          return {1.97e+00, 1.18e+00, 5.71e-01, 2.71e-01};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else if (mu == 1 && mu_bar == 1 && mu_hat == 0.1) {
        if (type == "eff_OS2014_*_mu")
          return {1.68e+00, 1.69e+00, 1.73e+00, 1.94e+00};
        else if (type == "eff_OS2014_mu")
          return {1.44e+01, 2.75e+01, 5.52e+01, 1.22e+02};
        else if (type == "eta_DF_OS2014")
          return {1.36e+00, 1.33e+00, 1.33e+00, 1.32e+00};
        else if (type == "eta_DF_OS2014_*")
          return {4.13e-01, 2.05e-01, 1.02e-01, 5.06e-02};
        else if (type == "eta_OS2014")
          return {4.71e+00, 4.42e+00, 4.30e+00, 4.24e+00};
        else if (type == "eta_OS2014_*")
          return {5.50e-01, 2.71e-01, 1.35e-01, 6.74e-02};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else if (mu == 0.1 && mu_bar == 0.1 && mu_hat == 1) {
        if (type == "eff_OS2014_*_mu")
          return {4.99e+00, 4.94e+00, 5.01e+00, 5.53e+00};
        else if (type == "eff_OS2014_mu")
          return {4.44e+00, 8.02e+00, 1.78e+01, 4.18e+01};
        else if (type == "eta_DF_OS2014")
          return {1.01e+00, 1.21e+00, 1.35e+00, 1.41e+00};
        else if (type == "eta_DF_OS2014_*")
          return {1.16e+00, 6.90e-01, 3.34e-01, 1.62e-01};
        else if (type == "eta_OS2014")
          return {3.91e+00, 4.27e+00, 4.48e+00, 4.55e+00};
        else if (type == "eta_OS2014_*")
          return {4.39e+00, 2.63e+00, 1.26e+00, 6.01e-01};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else if (mu == 1 && mu_bar == 1 && mu_hat == 1) {
        if (type == "eff_OS2014_*_mu")
          return {2.36e+00, 2.38e+00, 2.44e+00, 2.73e+00};
        else if (type == "eff_OS2014_mu")
          return {2.36e+00, 2.38e+00, 2.44e+00, 2.73e+00};
        else if (type == "eta_DF_OS2014")
          return {3.55e-01, 1.76e-01, 8.73e-02, 4.35e-02};
        else if (type == "eta_DF_OS2014_*")
          return {3.55e-01, 1.76e-01, 8.73e-02, 4.35e-02};
        else if (type == "eta_OS2014")
          return {7.74e-01, 3.82e-01, 1.90e-01, 9.49e-02};
        else if (type == "eta_OS2014_*")
          return {7.74e-01, 3.82e-01, 1.90e-01, 9.49e-02};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else
        EXPECT_TRUE(false) << "test results missing for parameters: mu = " << mu << "\n"
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
