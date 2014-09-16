// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#if HAVE_ALUGRID
# include <dune/grid/alugrid.hh>

# include <dune/stuff/test/gtest/gtest.h>

# include <dune/hdd/linearelliptic/testcases/spe10.hh>

# include "linearelliptic-swipdg-expectations.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Tests {


template< bool anything >
class SWIPDGStudyExpectations< TestCases::Spe10::Model1< ALUGrid< 2, 2, simplex, conforming > >, 1, anything >
  : public internal::SWIPDGStudyExpectationsBase< TestCases::Spe10::Model1< ALUGrid< 2, 2, simplex, conforming > >, 1 >
{
  typedef TestCases::Spe10::Model1< ALUGrid< 2, 2, simplex, conforming > > TestCaseType;

public:
  static std::vector< double > results(const TestCaseType& test_case, const std::string type)
  {
    if (test_case.num_refinements() == 1) {
      if (type == "L2")
        return {3.81e-03, 1.06e-03};
      else if (type == "H1_semi")
        return {3.38e-01, 1.67e-01};
      else if (type == "energy")
        return {8.38e-01, 4.02e-01};
      else if (type == "eta_NC_ESV2007")
        return {2.74e+00, 1.84e+00};
      else if (type == "eta_R_ESV2007")
        return {3.43e-17, 1.64e-17};
      else if (type == "eta_DF_ESV2007")
        return {1.22e+00, 7.62e-01};
      else if (type == "eta_ESV2007")
        return {3.00e+00, 1.99e+00};
      else if (type == "eta_ESV2007_alt")
        return {3.97e+00, 2.60e+00};
      else if (type == "eff_ESV2007")
        return {3.59e+00, 4.95e+00};
      else if (type == "eff_ESV2007_alt")
        return {4.74e+00, 6.46e+00};
      else
        EXPECT_TRUE(false) << "test results missing for type: " << type;
    } else if (test_case.num_refinements() == 2) {
      if (type == "L2")
        return {4.21e-03, 1.49e-03, 4.58e-04};
      else if (type == "H1_semi")
        return {3.76e-01, 2.18e-01, 1.07e-01};
      else if (type == "energy")
        return {9.24e-01, 5.21e-01, 2.51e-01};
      else
        EXPECT_TRUE(false) << "test results missing for type: " << type;
    } else
      EXPECT_TRUE(false) << "test results missing for num_refinements: " << test_case.num_refinements();
    return {};
  } // ... results(...)
}; // SWIPDGStudyExpectations


template< bool anything >
class SWIPDGStudyExpectations< TestCases::Spe10::ParametricModel1< ALUGrid< 2, 2, simplex, conforming > >, 1, anything >
  : public internal::SWIPDGStudyExpectationsBase
        < TestCases::Spe10::ParametricModel1< ALUGrid< 2, 2, simplex, conforming > >, 1 >
{
  typedef TestCases::Spe10::ParametricModel1< ALUGrid< 2, 2, simplex, conforming > > TestCaseType;

public:
  static std::vector< double > results(const TestCaseType& test_case, const std::string type)
  {
    const auto mu = test_case.parameters().at("mu");
    const auto mu_bar = test_case.parameters().at("mu_bar");
    const auto mu_hat = test_case.parameters().at("mu_hat");
    const auto mu_minimizing = test_case.parameters().at("mu_minimizing");
    if (test_case.num_refinements() == 2) {
      if (mu == 0.1) {
        if (type == "L2")
          return {5.33e-03, 1.85e-03, 5.72e-04};
        else if (type == "H1_semi")
          return {4.62e-01, 2.67e-01, 1.30e-01};
        else if (type == "energy")
          return {9.50e-01, 5.33e-01, 2.55e-01};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else if (mu == 0.2) {
        if (type == "L2")
          return {4.63e-03, 1.62e-03, 4.99e-04};
        else if (type == "H1_semi")
          return {4.30e-01, 2.49e-01, 1.21e-01};
        else if (type == "energy")
          return {9.58e-01, 5.39e-01, 2.59e-01};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else if (mu == 0.3) {
        if (type == "L2")
          return {4.40e-03, 1.54e-03, 4.78e-04};
        else if (type == "H1_semi")
          return {4.13e-01, 2.39e-01, 1.17e-01};
        else if (type == "energy")
          return {9.44e-01, 5.31e-01, 2.55e-01};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else if (mu == 0.4) {
        if (type == "L2")
          return {4.30e-03, 1.51e-03, 4.69e-04};
        else if (type == "H1_semi")
          return {4.02e-01, 2.33e-01, 1.14e-01};
        else if (type == "energy")
          return {9.35e-01, 5.26e-01, 2.53e-01};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else if (mu == 0.5) {
        if (type == "L2")
          return {4.25e-03, 1.50e-03, 4.64e-04};
        else if (type == "H1_semi")
          return {3.94e-01, 2.28e-01, 1.12e-01};
        else if (type == "energy")
          return {9.28e-01, 5.23e-01, 2.52e-01};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else if (mu == 0.6) {
        if (type == "L2")
          return {4.22e-03, 1.49e-03, 4.61e-04};
        else if (type == "H1_semi")
          return {3.89e-01, 2.25e-01, 1.10e-01};
        else if (type == "energy")
          return {9.24e-01, 5.20e-01, 2.50e-01};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else if (mu == 0.7) {
        if (type == "L2")
          return {4.21e-03, 1.49e-03, 4.60e-04};
        else if (type == "H1_semi")
          return {3.84e-01, 2.23e-01, 1.09e-01};
        else if (type == "energy")
          return {9.21e-01, 5.19e-01, 2.50e-01};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else if (mu == 0.8) {
        if (type == "L2")
          return {4.21e-03, 1.49e-03, 4.59e-04};
        else if (type == "H1_semi")
          return {3.81e-01, 2.21e-01, 1.08e-01};
        else if (type == "energy")
          return {9.18e-01, 5.17e-01, 2.49e-01};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else if (mu == 0.9) {
        if (type == "L2")
          return {4.21e-03, 1.49e-03, 4.58e-04};
        else if (type == "H1_semi")
          return {3.78e-01, 2.19e-01, 1.08e-01};
        else if (type == "energy")
          return {9.16e-01, 5.16e-01, 2.49e-01};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else if (mu == 1.0) {
        if (type == "L2")
          return {4.21e-03, 1.49e-03, 4.58e-04};
        else if (type == "H1_semi")
          return {3.76e-01, 2.18e-01, 1.07e-01};
        else if (type == "energy")
          return {9.15e-01, 5.16e-01, 2.48e-01};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else
        EXPECT_TRUE(false) << "test results missing for parameters: mu = " << mu << "\n"
                           << "                                     mu_bar = " << mu_bar << "\n"
                           << "                                     mu_hat = " << mu_hat << "\n"
                           << "                                     mu_minimizing = " << mu_minimizing;
    } else
      EXPECT_TRUE(false) << "test results missing for num_refinements: " << test_case.num_refinements();
    return {};
  } // ... results(...)
}; // SWIPDGStudyExpectations


template class SWIPDGStudyExpectations< TestCases::Spe10::Model1< ALUGrid< 2, 2, simplex, conforming > >, 1 >;

template class SWIPDGStudyExpectations< TestCases::Spe10::ParametricModel1< ALUGrid< 2, 2, simplex, conforming > >, 1 >;


} // namespace Tests
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // HAVE_ALUGRID
