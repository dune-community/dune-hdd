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
class EocStudyBlockSWIPDGExpectations< TestCases::OS2014Multiscale< ALUGrid< 2, 2, simplex, conforming > >, 1, implemented >
{
  typedef TestCases::OS2014Multiscale< ALUGrid< 2, 2, simplex, conforming > > TestCaseType;

public:
  static std::vector< double > results(const TestCaseType& test_case, const std::string type)
  {
    const auto mu = test_case.parameters().at("mu");
    const auto mu_bar = test_case.parameters().at("mu_bar");
    const auto mu_hat = test_case.parameters().at("mu_hat");
    const auto mu_minimizing = test_case.parameters().at("mu_minimizing");
    if (test_case.partitioning() == "[1 1 1]") {
      if (   mu ==            Pymor::Parameter("mu", 1)
          && mu_bar ==        Pymor::Parameter("mu", 1)
          && mu_hat ==        Pymor::Parameter("mu", 1)
          && mu_minimizing == Pymor::Parameter("mu", 1)) {
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
