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
    const auto mu            = test_case.parameters().at("mu");
    const auto mu_bar        = test_case.parameters().at("mu_bar");
    const auto mu_hat        = test_case.parameters().at("mu_hat");
    const auto mu_minimizing = test_case.parameters().at("mu_minimizing");
    if (test_case.num_refinements() == 1) {
      if (test_case.partitioning() == "[20 4 1]") {
        if (mu == 0.1 && mu_bar == 0.1 && mu_hat == 0.1 && mu_minimizing == 0.1) {
          if (type == "energy_mu")
            return {9.13e-01, 4.40e-01};
          else if (type == "eta_OS2014")
            return {3.67e+00, 2.29e+00};
          else if (type == "eta_OS2014_*")
            return {3.67e+00, 2.29e+00};
          else
            EXPECT_TRUE(false) << "test results missing for type: " << type;
        } else if (mu == 1 && mu_bar == 1 && mu_hat == 0.1 && mu_minimizing == 0.1) {
          if (type == "energy_mu")
            return {8.38e-01, 4.02e-01};
          else if (type == "eta_OS2014")
            return {3.54e+01, 3.43e+01};
          else if (type == "eta_OS2014_*")
            return {3.34e+00, 2.21e+00};
          else
            EXPECT_TRUE(false) << "test results missing for type: " << type;
        } else if (mu == 0.1 && mu_bar == 0.1 && mu_hat == 1 && mu_minimizing == 0.1) {
          if (type == "energy_mu")
            return {9.13e-01, 4.40e-01};
          else if (type == "eta_OS2014")
            return {2.55e+01, 2.43e+01};
          else if (type == "eta_OS2014_*")
            return {6.52e+00, 4.06e+00};
          else
            EXPECT_TRUE(false) << "test results missing for type: " << type;
        } else if (mu == 1 && mu_bar == 1 && mu_hat == 1 && mu_minimizing == 0.1) {
          if (type == "energy_mu")
            return {8.38e-01, 4.02e-01};
          else if (type == "eta_OS2014")
            return {3.97e+00, 2.60e+00};
          else if (type == "eta_OS2014_*")
            return {3.97e+00, 2.60e+00};
          else
            EXPECT_TRUE(false) << "test results missing for type: " << type;
        } else
          EXPECT_TRUE(false) << "test results missing for parameters: mu = " << mu << "\n"
                             << "                                     mu_bar = " << mu_bar << "\n"
                             << "                                     mu_hat = " << mu_hat << "\n"
                             << "                                     mu_minimizing = " << mu_minimizing;
      } else
        EXPECT_TRUE(false) << "test results missing for partitioning: " << test_case.partitioning();
    } else
      EXPECT_TRUE(false) << "test results missing for num_refinements: " << test_case.num_refinements();
    return {};
  } // ... results(...)
}; // BlockSWIPDGStudyExpectations


template class BlockSWIPDGStudyExpectations< TestCaseType, 1 >;


} // namespace Tests
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // HAVE_ALUGRID
