// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include <dune/stuff/test/gtest/gtest.h>
#include <dune/stuff/test/common.hh>

#include "OS2014.hh"

#if HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_GRID_MULTISCALE


void print_parameter_information(const ParametricBlockTestCaseType& parametric_test_case)
{
  const auto& parameters = parametric_test_case.parameters();
  const auto& parametric_problem = parametric_test_case.problem();
  for (auto parameter : parameters)
    EXPECT_EQ(parameter.second.type(), parametric_problem.parameter_type())
        << "          id: " << parameter.first << ", parameter: " << parameter.second;
  const auto& diffusion_factor = *parametric_problem.diffusion_factor();
  const auto& diffusion_tensor = *parametric_problem.diffusion_tensor();
  const auto& force = *parametric_problem.force();
  const auto& dirichlet = *parametric_problem.dirichlet();
  const auto& neumann = *parametric_problem.neumann();
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


void OS2014_nonparametric_convergence_study__SWIPDG_fine_triangulation()
{
  const NonparametricTestCaseType test_case;
  test_case.print_header(DSC_LOG_INFO);
  DSC_LOG_INFO << std::endl;
  NonparametricEocStudyType eoc_study(test_case,
                                      {"energy", "eta_NC_ESV2007", "eta_R_ESV2007", "eta_DF_ESV2007", "eta_ESV2007",
                                       "eff_ESV2007"});
  Stuff::Test::check_for_success(eoc_study, eoc_study.run(false, DSC_LOG_INFO));
} // ... OS2014_nonparametric_convergence_study__SWIPDG_fine_triangulation(...)


void OS2014_nonparametric_convergence_study__SWIPDG_fine_triangulation_alternative_summation()
{
  const NonparametricTestCaseType test_case;
  NonparametricEocStudyType eoc_study(test_case,
                                      {"energy", "eta_ESV2007", "eff_ESV2007", "eta_ESV2007_alt", "eff_ESV2007_alt"});
  Stuff::Test::check_for_success(eoc_study, eoc_study.run(false, DSC_LOG_INFO));
} // ... OS2014_nonparametric_convergence_study__SWIPDG_fine_triangulation_alternative_summation(...)


void nonparametric_block_convergence_study(const std::string& partitioning)
{
  const NonparametricBlockTestCaseType test_case(partitioning);
  NonparametricBlockEocStudyType eoc_study(test_case,
                                           {"energy", "eta_NC_OS2014", "eta_R_OS2014", "eta_DF_OS2014", "eta_OS2014",
                                            "eff_OS2014"});
  Dune::Stuff::Test::check_for_success(eoc_study, eoc_study.run(false, DSC_LOG_INFO));
} // ... nonparametric_block_convergence_study(...)


{


#endif // HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_GRID_MULTISCALE
