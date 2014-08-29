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


template< bool implemented >
class EocStudyBlockSWIPDGExpectations< TestCases::OS2014Multiscale< ALUGrid< 2, 2, simplex, conforming > >,
                                       1, implemented >
  : public internal::EocStudyBlockSWIPDGExpectationsBase
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
          return {8.80e-01, 5.32e-01, 2.51e-01, 1.09e-01};
        else if (type == "eta_DF_OS2014")
          return {1.25e+00, 7.37e-01, 3.69e-01, 1.83e-01};
        else if (type == "eta_DF_OS2014_*")
          return {1.25e+00, 7.37e-01, 3.69e-01, 1.83e-01};
        else if (type == "eta_OS2014")
          return {2.73e+00, 1.56e+00, 7.62e-01, 3.67e-01};
        else if (type == "eta_OS2014_*")
          return {2.73e+00, 1.56e+00, 7.62e-01, 3.67e-01};
        else if (type == "eff_OS2014_mu")
          return {3.11e+00, 2.94e+00, 3.03e+00, 3.37e+00};
        else if (type == "eff_OS2014_*_mu")
          return {3.11e+00, 2.94e+00, 3.03e+00, 3.37e+00};
        else
          DUNE_THROW(Stuff::Exceptions::test_results_missing, type);
      } else if (   mu            == Parameter("mu", 0.3)
                 && mu_bar        == Parameter("mu", 0.3)
                 && mu_hat        == Parameter("mu", 0.1)
                 && mu_minimizing == Parameter("mu", 0.1)) {
        if (type == "energy_mu")
          return {6.28e-01, 3.64e-01, 1.74e-01, 7.64e-02};
        else if (type == "eta_DF_OS2014")
          return {1.19e+00, 7.72e-01, 5.24e-01, 4.31e-01};
        else if (type == "eta_DF_OS2014_*")
          return {9.41e-01, 5.40e-01, 2.72e-01, 1.36e-01};
        else if (type == "eta_OS2014")
          return {3.39e+00, 2.05e+00, 1.26e+00, 9.15e-01};
        else if (type == "eta_OS2014_*")
          return {1.88e+00, 1.03e+00, 5.07e-01, 2.47e-01};
        else if (type == "eff_OS2014_*_mu")
          return {2.99e+00, 2.82e+00, 2.91e+00, 3.23e+00};
        else if (type == "eff_OS2014_mu")
          return {5.40e+00, 5.64e+00, 7.21e+00, 1.20e+01};
        else
          DUNE_THROW(Stuff::Exceptions::test_results_missing, type);
      } else if (   mu            == Parameter("mu", 0.5)
                 && mu_bar        == Parameter("mu", 0.5)
                 && mu_hat        == Parameter("mu", 0.1)
                 && mu_minimizing == Parameter("mu", 0.1)) {
        if (type == "energy_mu")
          return {4.74e-01, 2.62e-01, 1.27e-01, 5.61e-02};
        else if (type == "eta_DF_OS2014")
          return {1.20e+00, 9.06e-01, 7.61e-01, 7.13e-01};
        else if (type == "eta_DF_OS2014_*")
          return {7.14e-01, 3.95e-01, 1.99e-01, 9.91e-02};
        else if (type == "eta_OS2014")
          return {3.94e+00, 2.68e+00, 2.02e+00, 1.75e+00};
        else if (type == "eta_OS2014_*")
          return {1.57e+00, 8.28e-01, 4.12e-01, 2.03e-01};
        else if (type == "eff_OS2014_*_mu")
          return {3.32e+00, 3.17e+00, 3.25e+00, 3.62e+00};
        else if (type == "eff_OS2014_mu")
          return {8.32e+00, 1.02e+01, 1.60e+01, 3.13e+01};
        else
          DUNE_THROW(Stuff::Exceptions::test_results_missing, type);
      } else if (   mu            == Parameter("mu", 0.1)
                 && mu_bar        == Parameter("mu", 0.1)
                 && mu_hat        == Parameter("mu", 1.0)
                 && mu_minimizing == Parameter("mu", 0.1)) {
        if (type == "eta_DF_OS2014")
          return {1.36e+00, 1.33e+00, 1.33e+00, 1.32e+00};
        else if (type == "eta_DF_OS2014_*")
          return {4.13e-01, 2.05e-01, 1.02e-01, 5.06e-02};
        else if (type == "eta_OS2014")
          return {5.47e+00, 4.80e+00, 4.49e+00, 4.34e+00};
        else if (type == "eta_OS2014_*")
          return {1.31e+00, 6.53e-01, 3.26e-01, 1.63e-01};
        else if (type == "eff_OS2014_*_mu")
          return {4.00e+00, 4.07e+00, 4.19e+00, 4.69e+00};
        else if (type == "eff_OS2014_mu")
          return {1.67e+01, 2.99e+01, 5.77e+01, 1.25e+02};
        else
          DUNE_THROW(Stuff::Exceptions::test_results_missing, type);
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
          return {3.29e-01, 1.61e-01, 7.79e-02, 3.48e-02};
        else if (type == "eta_NC_OS2014")
          return {1.67e-01, 7.90e-02, 3.92e-02, 1.96e-02};
        else if (type == "eta_R_OS2014")
          return {5.80e-01, 2.91e-01, 1.46e-01, 7.28e-02};
        else if (type == "eta_DF_OS2014")
          return {3.56e-01, 1.77e-01, 8.74e-02, 4.36e-02};
        else if (type == "eta_OS2014")
          return {1.11e+00, 5.46e-01, 2.73e-01, 1.37e-01};
        else if (type == "eff_OS2014_mu")
          return {3.37e+00, 3.41e+00, 3.50e+00, 3.92e+00};
        else
          DUNE_THROW(Stuff::Exceptions::test_results_missing, type);
      } else
          DUNE_THROW(Stuff::Exceptions::test_results_missing,
                     "mu = " << mu << "\nmu_bar = " << mu_bar << "\nmu_hat = " << mu_hat
                     << "\nmu_minimizing = " << mu_minimizing);
    } else
      DUNE_THROW(Stuff::Exceptions::test_results_missing, test_case.partitioning());
  } // ... results(...)
}; // EocStudyBlockSWIPDGExpectations


template class EocStudyBlockSWIPDGExpectations< TestCases::OS2014Multiscale< ALUGrid< 2, 2, simplex, conforming > >, 1, true >;


} // namespace Tests
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

//#endif // HAVE_ALUGRID
