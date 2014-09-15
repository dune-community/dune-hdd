// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include <dune/grid/sgrid.hh>

#include <dune/stuff/test/gtest/gtest.h>

#include <dune/hdd/playground/linearelliptic/testcases/spe10.hh>

#include "linearelliptic-swipdg-expectations.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Tests {


template< bool anything >
class SWIPDGStudyExpectations< TestCases::Spe10::Model1< SGrid< 2, 2 > >, 1, anything >
  : public internal::SWIPDGStudyExpectationsBase< TestCases::Spe10::Model1< SGrid< 2, 2 > >, 1 >
{
  typedef TestCases::Spe10::Model1< SGrid< 2, 2 > > TestCaseType;

public:
  static std::vector< double > results(const TestCaseType& test_case, const std::string type)
  {
    if (test_case.num_refinements() == 1) {
      if (type == "L2")
        return {1.10e-02, 9.46e-03};
      else if (type == "H1_semi")
        return {7.91e-01, 1.46e+00};
      else if (type == "energy")
        return {7.39e+00, 1.82e+01};
      else
        EXPECT_TRUE(false) << "test results missing for type: " << type;
    } else if (test_case.num_refinements() == 2) {
      if (type == "L2")
        return {1.10e-02, 9.16e-03, 3.86e-03};
      else if (type == "H1_semi")
        return {8.10e-01, 1.63e+00, 1.19e+00};
      else if (type == "energy")
        return {7.58e+00, 2.04e+01, 1.48e+01};
      else
        EXPECT_TRUE(false) << "test results missing for type: " << type;
    } else
      EXPECT_TRUE(false) << "test results missing for num_refinements: " << test_case.num_refinements();
    return {};
  } // ... results(...)
}; // SWIPDGStudyExpectations


template< bool anything >
class SWIPDGStudyExpectations< TestCases::Spe10::ParametricModel1< SGrid< 2, 2 > >, 1, anything >
  : public internal::SWIPDGStudyExpectationsBase< TestCases::Spe10::ParametricModel1< SGrid< 2, 2 > >, 1 >
{
  typedef TestCases::Spe10::ParametricModel1< SGrid< 2, 2 > > TestCaseType;

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
          return {1.45e-02, 9.39e-03, 3.78e-03};
        else if (type == "H1_semi")
          return {1.05e+00, 1.60e+00, 1.06e+00};
        else if (type == "energy")
          return {7.04e+00, 1.88e+01, 1.17e+01};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else if (mu == 0.2) {
        if (type == "L2")
          return {1.21e-02, 1.02e-02, 3.48e-03};
        else if (type == "H1_semi")
          return {8.53e-01, 1.81e+00, 1.06e+00};
        else if (type == "energy")
          return {6.84e+00, 2.11e+01, 1.27e+01};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else if (mu == 0.3) {
        if (type == "L2")
          return {1.16e-02, 8.79e-03, 3.43e-03};
        else if (type == "H1_semi")
          return {8.26e-01, 1.54e+00, 1.03e+00};
        else if (type == "energy")
          return {6.87e+00, 1.96e+01, 1.22e+01};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else if (mu == 0.4) {
        if (type == "L2")
          return {1.16e-02, 8.33e-03, 4.05e-03};
        else if (type == "H1_semi")
          return {8.52e-01, 1.46e+00, 1.27e+00};
        else if (type == "energy")
          return {7.50e+00, 1.80e+01, 1.51e+01};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else if (mu == 0.5) {
        if (type == "L2")
          return {1.15e-02, 8.26e-03, 4.01e-03};
        else if (type == "H1_semi")
          return {8.53e-01, 1.45e+00, 1.26e+00};
        else if (type == "energy")
          return {7.69e+00, 1.78e+01, 1.66e+01};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else if (mu == 0.6) {
        if (type == "L2")
          return {1.10e-02, 7.85e-03, 5.81e-03};
        else if (type == "H1_semi")
          return {7.96e-01, 1.37e+00, 1.90e+00};
        else if (type == "energy")
          return {6.89e+00, 1.63e+01, 1.89e+01};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else if (mu == 0.7) {
        if (type == "L2")
          return {1.09e-02, 4.08e-02, 3.39e-03};
        else if (type == "H1_semi")
          return {7.95e-01, 7.68e+00, 1.02e+00};
        else if (type == "energy")
          return {7.04e+00, 7.27e+01, 1.25e+01};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else if (mu == 0.8) {
        if (type == "L2")
          return {1.29e-02, 1.00e-02, 9.75e-03};
        else if (type == "H1_semi")
          return {1.04e+00, 1.80e+00, 3.27e+00};
        else if (type == "energy")
          return {9.84e+00, 2.20e+01, 3.30e+01};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else if (mu == 0.9) {
        if (type == "L2")
          return {1.11e-02, 9.47e-03, 4.44e-03};
        else if (type == "H1_semi")
          return {8.16e-01, 1.68e+00, 1.37e+00};
        else if (type == "energy")
          return {7.53e+00, 2.08e+01, 1.63e+01};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else if (mu == 1.0) {
        if (type == "L2")
          return {1.10e-02, 9.16e-03, 3.86e-03};
        else if (type == "H1_semi")
          return {8.10e-01, 1.63e+00, 1.19e+00};
        else if (type == "energy")
          return {7.48e+00, 2.02e+01, 2.02e+01};
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


template class SWIPDGStudyExpectations< TestCases::Spe10::Model1< SGrid< 2, 2 > >, 1 >;

template class SWIPDGStudyExpectations< TestCases::Spe10::ParametricModel1< SGrid< 2, 2 > >, 1 >;


} // namespace Tests
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune
