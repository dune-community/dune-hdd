// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#if HAVE_ALUGRID
# include <dune/grid/alugrid.hh>

# include <dune/hdd/linearelliptic/testcases/spe10.hh>

# include "linearelliptic-block-swipdg-expectations.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {


typedef TestCases::OS2015::Multiscale< ALUGrid< 2, 2, simplex, conforming > > TestCaseType;

namespace Tests {


template< bool anything >
class BlockSWIPDGStudyExpectations< TestCaseType, 1, anything >
  : public internal::BlockSWIPDGStudyExpectationsBase< TestCaseType, 1 >
{
public:
  static std::vector< double > results(const TestCaseType& test_case, const std::string type)
  {
    const auto mu     = test_case.parameters().at("mu");
    const auto mu_bar = test_case.parameters().at("mu_bar");
    const auto mu_hat = test_case.parameters().at("mu_hat");
    if (test_case.num_refinements() == 3 && test_case.partitioning() == "[25 5 1]_H_with_h") {
      // fine grid sizes: 16000, 64000, 256000, 1024000
      // fine grid width: 6.25e-03, 3.12e-03, 1.56e-03, 7.81e-04
      if (mu == 1 && mu_bar == 1 && mu_hat == 1) {
        if (type == "energy_mu")
          return {7.49e-01, 4.52e-01, 2.58e-01, 1.26e-01}; // average order 0.86
        else if (type == "eta_NC_OS2014")
          return {2.13e+00, 1.46e+00, 1.02e+00, 7.20e-01}; // average order 0.52
        else if (type == "eta_R_OS2014_*")
          return {1.88e-09, 7.05e-10, 7.44e-11, 2.00e-10}; // average order 1.07 [1.41, 3.24, -1.43]
        else if (type == "eta_DF_OS2014_*")
          return {9.66e-01, 6.05e-01, 3.85e-01, 2.49e-01}; // average order 0.65
        else if (type == "eta_OS2014_*")
          return {3.10e+00, 2.07e+00, 1.40e+00,9.69e-01 }; // average order 0.56
        else if (type == "eff_OS2014_*_mu")
          return {4.14e+00, 4.58e+00, 5.44e+00, 7.70e+00};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else if (mu == 1 && mu_bar == 1 && mu_hat == 0.1) {
        if (type == "eta_DF_OS2014_*")
          return {1.50e+00, 9.28e-01, 5.80e-01, 3.68e-01}; // average order 0.68
        else if (type == "eta_OS2014_*")
          return {3.63e+00, 2.39e+00, 1.60e+00, 1.09e+00}; // average order 0.58
        else if (type == "eff_OS2014_*_mu")
          return {4.85e+00, 5.29e+00, 6.20e+00, 8.64e+00};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else if (mu == 1 && mu_bar == 0.1 && mu_hat == 0.1) {
        if (type == "energy_mu_bar")
          return {7.43e-01, 4.48e-01, 2.56e-01, 1.25e-01}; // average order 0.86
        else if (type == "eta_NC_OS2014")
          return {1.51e+00, 9.64e-01, 6.29e-01, 4.17e-01}; // average order 0.62
        else if (type == "eta_DF_OS2014_*")
          return {1.50e+00, 9.28e-01, 5.80e-01, 3.68e-01}; // average order 0.68
        else if (type == "eta_OS2014_*")
          return {6.27e+00, 3.97e+00, 2.57e+00, 1.69e+00}; // average order 0.63
        else if (type == "eff_OS2014_*_mu_bar")
          return {8.44e+00, 8.88e+00, 1.00e+01, 1.35e+01};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else
        EXPECT_TRUE(false) << "test results missing for parameters: mu = " << mu << "\n"
                           << "                                     mu_bar = " << mu_bar << "\n"
                           << "                                     mu_hat = " << mu_hat;
    } else
      EXPECT_TRUE(false) << "test results missing for partitioning: " << test_case.partitioning() << "\n"
                         << "                  and num_refinements: " << test_case.num_refinements();
    return {};
  } // ... results(...)
}; // BlockSWIPDGStudyExpectations


template class BlockSWIPDGStudyExpectations< TestCaseType, 1 >;


} // namespace Tests
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // HAVE_ALUGRID
