// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

//#if HAVE_ALUGRID

#include <dune/hdd/playground/linearelliptic/testcases/OS2014.hh>

#include "linearelliptic-block-swipdg-expectations.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Tests {


template< bool anything >
class BlockSWIPDGStudyExpectations< TestCases::OS2014Multiscale< ALUGrid< 2, 2, simplex, conforming > >, 1, anything >
  : public internal::BlockSWIPDGStudyExpectationsBase
        < TestCases::OS2014Multiscale< ALUGrid< 2, 2, simplex, conforming > >, 1 >
{
  typedef TestCases::OS2014Multiscale< ALUGrid< 2, 2, simplex, conforming > > TestCaseType;

public:
  static std::vector< double > results(const TestCaseType& test_case, const std::string type)
  {
    using Parameter = Pymor::Parameter;
    const auto mu = test_case.parameters().at("mu");
    const auto mu_bar = test_case.parameters().at("mu_bar");
    const auto mu_hat = test_case.parameters().at("mu_hat");
    const auto mu_minimizing = test_case.parameters().at("mu_minimizing");
    if (test_case.partitioning() == "[1 1 1]") {
      if (   mu            == Parameter("mu", 0.1)
          && mu_bar        == Parameter("mu", 0.1)
          && mu_hat        == Parameter("mu", 0.1)
          && mu_minimizing == Parameter("mu", 0.1)) {
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
      } else if (   mu            == Parameter("mu", 0.3)
                 && mu_bar        == Parameter("mu", 0.3)
                 && mu_hat        == Parameter("mu", 0.1)
                 && mu_minimizing == Parameter("mu", 0.1)) {
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
      } else if (   mu            == Parameter("mu", 0.5)
                 && mu_bar        == Parameter("mu", 0.5)
                 && mu_hat        == Parameter("mu", 0.1)
                 && mu_minimizing == Parameter("mu", 0.1)) {
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
      } else if (   mu            == Parameter("mu", 0.1)
                 && mu_bar        == Parameter("mu", 0.1)
                 && mu_hat        == Parameter("mu", 1.0)
                 && mu_minimizing == Parameter("mu", 0.1)) {
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
//                 && mu_hat        == Parameter("mu", )
//                 && mu_minimizing == Parameter("mu", )) {
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
//          DUNE_THROW(Stuff::Exceptions::test_results_missing, type);
      } else if (   mu            == Parameter("mu", 1)
                 && mu_bar        == Parameter("mu", 1)
                 && mu_hat        == Parameter("mu", 1)
                 && mu_minimizing == Parameter("mu", 1)) {
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
                           << "                                     mu_hat = " << mu_hat << "\n"
                           << "                                     mu_minimizing = " << mu_minimizing;
    } else
      EXPECT_TRUE(false) << "test results missing for partitioning: " << test_case.partitioning();
  } // ... results(...)
}; // BlockSWIPDGStudyExpectations


template class BlockSWIPDGStudyExpectations< TestCases::OS2014Multiscale< ALUGrid< 2, 2, simplex, conforming > >, 1 >;


} // namespace Tests
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

//#endif // HAVE_ALUGRID
