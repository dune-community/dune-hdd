// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#if HAVE_ALUGRID

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
      if (type == "energy_mu") {
        if (   mu ==            Pymor::Parameter("mu", 1)
            && mu_bar ==        Pymor::Parameter("mu", 1)
            && mu_hat ==        Pymor::Parameter("mu", 1)
            && mu_minimizing == Pymor::Parameter("mu", 1))
          return {3.29e-01, 1.61e-01, 7.79e-02, 3.48e-02};
        else
          DUNE_THROW(Stuff::Exceptions::test_results_missing,
                     "mu = " << mu << "\nmu_bar = " << mu_bar << "\nmu_hat = " << mu_hat
                     << "\nmu_minimizing = " << mu_minimizing);
      } else
        DUNE_THROW(Stuff::Exceptions::test_results_missing, type);
    } else
      DUNE_THROW(Stuff::Exceptions::test_results_missing, test_case.partitioning());
  } // ... results(...)
}; // EocStudyBlockSWIPDGExpectations


template class EocStudyBlockSWIPDGExpectations< TestCases::OS2014Multiscale< ALUGrid< 2, 2, simplex, conforming > >, 1, true >;


} // namespace Tests
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // HAVE_ALUGRID
