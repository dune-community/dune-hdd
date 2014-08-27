// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#if HAVE_ALUGRID

#include <dune/hdd/playground/linearelliptic/testcases/ESV2007.hh>

#include "linearelliptic-block-swipdg-expectations.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Tests {


template< bool implemented >
class EocStudyBlockSWIPDGExpectations< TestCases::ESV2007Multiscale< ALUGrid< 2, 2, simplex, conforming > >, 1, implemented >
{
  typedef TestCases::ESV2007Multiscale< ALUGrid< 2, 2, simplex, conforming > > TestCaseType;

public:
  static std::vector< double > results(const TestCaseType& test_case, const std::string type)
  {
    if (test_case.partitioning() == "[1 1 1]") {
      if (type == "energy")
        return {3.29e-01, 1.63e-01, 8.05e-02, 4.02e-02};
      else if (type == "eta_NC_OS2014")
        return {1.67e-01, 7.90e-02, 3.92e-02, 1.96e-02};
      else if (type == "eta_R_OS2014")
        return {5.80e-01, 2.91e-01, 1.46e-01, 7.28e-02};
      else if (type == "eta_DF_OS2014")
        return {3.56e-01, 1.77e-01, 8.74e-02, 4.36e-02};
      else if (type == "eta_OS2014")
        return {1.11e+00, 5.46e-01, 2.73e-01, 1.37e-01};
      else if (type == "eff_OS2014")
        return {3.36, 3.38, 3.39, 3.40};
      else
        DUNE_THROW(Stuff::Exceptions::test_results_missing, type);
    } else if (test_case.partitioning() == "[2 2 1]") {
      if (type == "energy")
        return {3.29e-01, 1.63e-01, 8.05e-02, 4.02e-02};
      else if (type == "eta_NC_OS2014")
        return {1.67e-01, 7.90e-02, 3.92e-02, 1.96e-02};
      else if (type == "eta_R_OS2014")
        return {2.90e-01, 1.46e-01, 7.28e-02, 3.64e-02};
      else if (type == "eta_DF_OS2014")
        return {3.56e-01, 1.77e-01, 8.74e-02, 4.36e-02};
      else if (type == "eta_OS2014")
        return {1.11e+00, 5.46e-01, 2.73e-01, 1.37e-01};
      else if (type == "eff_OS2014")
        return {2.48, 2.48, 2.49, 2.49};
      else
        DUNE_THROW(Stuff::Exceptions::test_results_missing, type);
    } else if (test_case.partitioning() == "[4 4 1]") {
      if (type == "energy")
        return {3.29e-01, 1.63e-01, 8.05e-02, 4.02e-02};
      else if (type == "eta_NC_OS2014")
        return {1.67e-01, 7.90e-02, 3.92e-02, 1.96e-02};
      else if (type == "eta_R_OS2014")
        return {1.46e-01, 7.27e-02, 3.64e-02, 1.82e-02};
      else if (type == "eta_DF_OS2014")
        return {3.56e-01, 1.77e-01, 8.74e-02, 4.36e-02};
      else if (type == "eta_OS2014")
        return {1.11e+00, 5.46e-01, 2.73e-01, 1.37e-01};
      else if (type == "eff_OS2014")
        return {2.04, 2.03, 2.03, 2.04};
      else
        DUNE_THROW(Stuff::Exceptions::test_results_missing, type);
    } else if (test_case.partitioning() == "[8 8 1]") {
      if (type == "energy")
        return {3.29e-01, 1.63e-01, 8.05e-02, 4.02e-02};
      else if (type == "eta_NC_OS2014")
        return {1.67e-01, 7.90e-02, 3.92e-02, 1.96e-02};
      else if (type == "eta_R_OS2014")
        return {7.24e-02, 3.64e-02, 1.83e-02, 9.10e-03};
      else if (type == "eta_DF_OS2014")
        return {3.56e-01, 1.77e-01, 8.74e-02, 4.36e-02};
      else if (type == "eta_OS2014")
        return {1.11e+00, 5.46e-01, 2.73e-01, 1.37e-01};
      else if (type == "eff_OS2014")
        return {1.82, 1.81, 1.81, 1.81};
      else
        DUNE_THROW(Stuff::Exceptions::test_results_missing, type);
    } else
      DUNE_THROW(Stuff::Exceptions::test_results_missing, test_case.partitioning());
  } // ... results(...)
}; // EocStudyBlockSWIPDGExpectations


template class EocStudyBlockSWIPDGExpectations< TestCases::ESV2007Multiscale< ALUGrid< 2, 2, simplex, conforming > >, 1, true >;


} // namespace Tests
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // HAVE_ALUGRID
