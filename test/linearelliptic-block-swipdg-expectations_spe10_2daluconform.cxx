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


typedef TestCases::Spe10::ParametricBlockModel1< ALUGrid< 2, 2, simplex, conforming > > TestCaseType;

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
    if (test_case.num_refinements() == 2 && test_case.partitioning() == "[10 2 1]") {
      // fine grid sizes: 8000, 32000, 128000
      // fine grid width: 8.84e-03, 4.42e-03, 2.21e-03
      if (mu == 1 && mu_bar == 1 && mu_hat == 1) {
        if (type == "energy_mu")
          return {9.24e-01, 5.21e-01, 2.51e-01};
        else if (type == "eta_NC_OS2014")
          return {2.74e+00, 1.84e+00, 1.27e+00};
        else if (type == "eta_R_OS2014_*")
          return {1.51e-09, 1.98e-08, 4.15e-08};
        else if (type == "eta_DF_OS2014_*")
          return {1.22e+00, 7.62e-01, 4.81e-01};
        else if (type == "eta_OS2014_*")
          return {3.97e+00, 2.60e+00, 1.75e+00};
        else if (type == "eff_OS2014_*_mu")
          return {4.29e+00, 4.99e+00, 6.98e+00};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else
        EXPECT_TRUE(false) << "test results missing for parameters: mu = " << mu << "\n"
                           << "                                     mu_bar = " << mu_bar << "\n"
                           << "                                     mu_hat = " << mu_hat;
    } else if (test_case.num_refinements() == 2 && test_case.partitioning() == "[25 5 1]_H_with_h") {
      // fine grid sizes: 8000, 32000, 128000
      // fine grid width: 8.84e-03, 4.42e-03, 2.21e-03
      if (mu == 1 && mu_bar == 1 && mu_hat == 1) {
        if (type == "energy_mu")
          return {9.24e-01, 5.21e-01, 2.51e-01}; // average order 0.94
        else if (type == "eta_NC_OS2014")
          return {2.74e+00, 1.84e+00, 1.27e+00}; // average order 0.56
        else if (type == "eta_R_OS2014_*")
          return {3.61e-09, 4.65e-10, 1.03e-10}; // average order 2.57
        else if (type == "eta_DF_OS2014_*")
          return {1.22e+00, 7.62e-01, 4.81e-01}; // average order 0.67
        else if (type == "eta_OS2014_*")
          return {3.97e+00, 2.60e+00, 1.75e+00}; // average order 0.59
        else if (type == "eff_OS2014_*_mu")
          return {4.29e+00, 4.99e+00, 6.98e+00};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else if (mu == 1 && mu_bar == 1 && mu_hat == 0.1) {
        if (type == "eta_DF_OS2014_*")
          return {1.89e+00, 1.16e+00, 7.26e-01}; // average order 0.69
        else if (type == "eta_OS2014_*")
          return {4.63e+00, 3.00e+00, 2.00e+00}; // average order 0.61
        else if (type == "eff_OS2014_*_mu")
          return {5.01e+00, 5.76e+00, 7.95e+00};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else if (mu == 1 && mu_bar == 0.1 && mu_hat == 0.1) {
        if (type == "energy_mu_bar")
          return {9.16e-01, 5.16e-01, 2.49e-01}; // average order 0.94
        else if (type == "eta_NC_OS2014")
          return {2.01e+00, 1.25e+00, 8.06e-01}; // average order 0.66
        else if (type == "eta_DF_OS2014_*")
          return {1.89e+00, 1.16e+00, 7.26e-01}; // average order 0.69
        else if (type == "eta_OS2014_*")
          return {8.24e+00, 5.12e+00, 3.28e+00}; // average order 0.67
        else if (type == "eff_OS2014_*_mu_bar")
          return {9.00e+00, 9.92e+00, 1.32e+01};
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
