// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_TEST_MAIN_CATCH_EXCEPTIONS
# define DUNE_STUFF_TEST_MAIN_CATCH_EXCEPTIONS
#endif
#ifndef DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING
# define DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING
#endif

#define DUNE_HDD_LINEARELLIPTIC_TESTCASES_BASE_DISABLE_WARNING

// This one has to come first (includes the config.h)!
#include <dune/stuff/test/main.hxx>

#include <vector>
#include <string>
#include <map>

#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>

#include <dune/pymor/parameters/base.hh>

#include <dune/hdd/playground/linearelliptic/testcases/OS2014.hh>

#include "linearelliptic-block-swipdg.hh"
#include "linearelliptic-block-swipdg-expectations.hh"

using namespace Dune;
using namespace HDD;


typedef ALUGrid< 2, 2, simplex, conforming > GridType;

typedef LinearElliptic::TestCases::OS2014::ParametricBlockConvergence< GridType > SmoothTestCaseType;
typedef LinearElliptic::Tests::BlockSWIPDGStudy< SmoothTestCaseType >             SmoothStudyType;


template< class TestCaseType >
void print_parameter_information(const TestCaseType& test_case)
{
  const auto& parameters = test_case.parameters();
  const auto& problem = test_case.problem();
  for (auto parameter : parameters)
    EXPECT_EQ(parameter.second.type(), problem.parameter_type())
        << "          id: " << parameter.first << ", parameter: " << parameter.second;
  const auto& diffusion_factor = *problem.diffusion_factor();
  const auto& diffusion_tensor = *problem.diffusion_tensor();
  const auto& force = *problem.force();
  const auto& dirichlet = *problem.dirichlet();
  const auto& neumann = *problem.neumann();
  EXPECT_TRUE(diffusion_factor.parametric());
  EXPECT_FALSE(diffusion_tensor.parametric());
  EXPECT_FALSE(force.parametric());
  EXPECT_FALSE(dirichlet.parametric());
  EXPECT_FALSE(neumann.parametric());
  DSC_LOG_INFO << "| mu            = " << parameters.at("mu") << "\n"
               << "| mu_bar        = " << parameters.at("mu_bar") << "\n"
               << "| mu_hat        = " << parameters.at("mu_hat") << "\n"
               << "| mu_minimizing = " << parameters.at("mu_minimizing") << "\n";
  const double alpha = diffusion_factor.alpha(parameters.at("mu"), parameters.at("mu_hat"));
  const double gamma = diffusion_factor.gamma(parameters.at("mu"), parameters.at("mu_hat"));
  DSC_LOG_INFO << "| alpha(mu, mu_hat)^-1/2    = " << std::setprecision(2) << std::scientific
                                                   << 1.0/std::sqrt(alpha) << "\n"
               << "| gamma_tilde(mu, mu_hat)^2 = " << std::setprecision(2) << std::scientific
                                                   << std::max(std::sqrt(gamma), 1.0/std::sqrt(alpha)) << "\n"
               << "+==================================================================+\n";
} // ... print_parameter_information(...)


template< class TestCaseType, class StudyType >
void run_eoc_study(const std::string partitioning,
                   const std::vector< std::string >& only_these_norms,
                   const std::map< std::string, Pymor::Parameter >& parameters,
                   const bool print_header,
                   const std::string visualization)
{
  const TestCaseType test_case(parameters, partitioning);
  if (print_header)
    test_case.print_header(DSC_LOG_INFO);
  print_parameter_information(test_case);
  StudyType study(test_case, only_these_norms, {}, visualization);
  Stuff::Test::check_eoc_study_for_success(study, study.run_eoc(DSC_LOG_INFO));
} // ... parametric_block_convergence_study(...)


TEST(OS2014_parametric_convergence_study, eta_DF_comparison)
{
  const std::string partitioning = "[4 4 1]";
  const std::vector< std::string > only_these_norms = {"eta_DF_OS2014", "eta_DF_OS2014_*", "eta_OS2014",
                                                       "eta_OS2014_*", "eff_OS2014_mu", "eff_OS2014_*_mu"};
  const std::string visualization_prefix = "parametric_block_convergence_study_eta_DF_comparison";
  bool print_header = true;
  for (auto mu_hat_value : {0.1, 0.5, 1.0}) {
    const auto mu_hat = Pymor::Parameter("mu", mu_hat_value);
    for (auto mu_value : {0.1, 0.3, 0.5, 0.75, 1.0}) {
      const auto mu = Pymor::Parameter("mu", mu_value);
      const auto mu_bar = mu;
      run_eoc_study< SmoothTestCaseType, SmoothStudyType >(partitioning,
                                                           only_these_norms,
                                                           {{"mu_hat",        mu_hat},
                                                            {"mu_bar",        mu_bar},
                                                            {"mu",            mu},
                                                            {"mu_minimizing", Pymor::Parameter("mu", 0.1)}},
                                                           print_header,
                                                           visualization_prefix);
      if (print_header)
        print_header = false;
    }
  }
} // TEST(OS2014_parametric_convergence_study, eta_DF_comparison)


extern template class Dune::HDD::LinearElliptic::Tests::BlockSWIPDGStudyExpectations< SmoothTestCaseType >;


#else // HAVE_ALUGRID


TEST(DISABLED_OS2014_parametric_convergence_study, eta_DF_comparison) {}


#endif // HAVE_ALUGRID
