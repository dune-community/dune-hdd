// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#define DUNE_STUFF_FUNCTIONS_DISABLE_CHECKS 1
#define DUNE_HDD_LINEARELLIPTIC_TESTCASES_BASE_DISABLE_WARNING 1

#ifndef DUNE_STUFF_TEST_MAIN_CATCH_EXCEPTIONS
# define DUNE_STUFF_TEST_MAIN_CATCH_EXCEPTIONS 1
#endif
#ifndef DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING
# define DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING 1
#endif

// This one has to come first (includes the config.h)!
#include <dune/stuff/test/main.hh>

#if HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_GRID_MULTISCALE

#include <dune/stuff/common/disable_warnings.hh>
# include <dune/grid/alugrid.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/hdd/playground/linearelliptic/testcases/ESV2007.hh>
#include <dune/hdd/playground/linearelliptic/testcases/OS2014.hh>

#include "linearelliptic-swipdg.hh"
#include "linearelliptic-block-swipdg.hh"

using namespace Dune;
using namespace HDD;
using Parameter = Pymor::Parameter;


static const GDT::ChooseSpaceBackend  space_backend =
#if HAVE_DUNE_FEM
                                                      GDT::ChooseSpaceBackend::fem;
#else
# error This test requires dune-fem!
#endif

static const Stuff::LA::ChooseBackend la_backend    =
#if HAVE_EIGEN
                                                      Stuff::LA::ChooseBackend::eigen_sparse;
#else
                                                      Stuff::LA::default_sparse_backend;
#endif

typedef ALUGrid< 2, 2, simplex, conforming > GridType;

typedef LinearElliptic::TestCases::ESV2007< GridType >           NonparametricTestCaseType;
typedef LinearElliptic::TestCases::ESV2007Multiscale< GridType > NonparametricBlockTestCaseType;
typedef LinearElliptic::TestCases::OS2014Multiscale< GridType >  ParametricBlockTestCaseType;

typedef LinearElliptic::Tests::EocStudySWIPDG< NonparametricTestCaseType, 1, space_backend, la_backend >
                                                                 NonparametricEocStudyType;
typedef LinearElliptic::Tests::EocStudyBlockSWIPDG< NonparametricBlockTestCaseType, 1, la_backend >
                                                                 NonparametricBlockEocStudyType;
typedef LinearElliptic::Tests::EocStudyBlockSWIPDG< ParametricBlockTestCaseType, 1, la_backend >
                                                                 ParametricBlockEocStudyType;


void nonparametric_block_convergence_study(const std::string& partitioning)
{
  const NonparametricBlockTestCaseType test_case(partitioning);
  NonparametricBlockEocStudyType eoc_study(test_case,
                                           {"energy", "eta_NC_OS2014", "eta_R_OS2014", "eta_DF_OS2014", "eta_OS2014",
                                            "eff_OS2014"});
  check_for_success(eoc_study, eoc_study.run(false, DSC_LOG_INFO));
} // ... nonparametric_block_convergence_study(...)


template< class PT >
void print_parameter_information(const PT& parametric_test_case)
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


TEST(OS2014_nonparametric_convergence_study, SWIPDG_fine_triangulation)
{
  const NonparametricTestCaseType test_case;
  test_case.print_header(DSC_LOG_INFO);
  DSC_LOG_INFO << std::endl;
  NonparametricEocStudyType eoc_study(test_case,
                                      {"energy", "eta_NC_ESV2007", "eta_R_ESV2007", "eta_DF_ESV2007", "eta_ESV2007",
                                       "eff_ESV2007"});
  check_for_success(eoc_study, eoc_study.run(false, DSC_LOG_INFO));
}
TEST(OS2014_nonparametric_convergence_study, SWIPDG_fine_triangulation_alternative_summation)
{
  const NonparametricTestCaseType test_case;
  NonparametricEocStudyType eoc_study(test_case,
                                      {"energy", "eta_ESV2007", "eff_ESV2007", "eta_ESV2007_alt", "eff_ESV2007_alt"});
  check_for_success(eoc_study, eoc_study.run(false, DSC_LOG_INFO));
}
TEST(OS2014_nonparametric_convergence_study, Block_SWIPDG_01_subdomain) {
  nonparametric_block_convergence_study("[1 1 1]");
}
TEST(OS2014_nonparametric_convergence_study, Block_SWIPDG_04_subdomain) {
  nonparametric_block_convergence_study("[2 2 1]");
}
TEST(OS2014_nonparametric_convergence_study, Block_SWIPDG_16_subdomain) {
  nonparametric_block_convergence_study("[4 4 1]");
}
TEST(OS2014_nonparametric_convergence_study, Block_SWIPDG_64_subdomain) {
  nonparametric_block_convergence_study("[8 8 1]");
}

TEST(OS2014_parametric_convergence_study, Block_SWIPDG_01_subdomain)
{
  const ParametricBlockTestCaseType test_case({{"mu_hat", Parameter("mu", 1)},
                                               {"mu_bar", Parameter("mu", 1)},
                                               {"mu", Parameter("mu", 1)},
                                               {"mu_minimizing", Parameter("mu", 1)}},
                                              "[1 1 1]");
  test_case.print_header(DSC_LOG_INFO);
  print_parameter_information(test_case);
  DSC_LOG_INFO << std::endl;
  ParametricBlockEocStudyType eoc_study(test_case,
                                        {"energy_mu", "eta_NC_OS2014", "eta_R_OS2014", "eta_DF_OS2014", "eta_OS2014",
                                         "eff_OS2014_mu"});
  check_for_success(eoc_study, eoc_study.run(false, DSC_LOG_INFO));
}


#else // HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_GRID_MULTISCALE


TEST(DISABLED_OS2014_nonparametric_convergence_study, SWIPDG_fine_triangulation) {}
TEST(DISABLED_OS2014_nonparametric_convergence_study, SWIPDG_fine_triangulation_alternative_summation) {}
TEST(DISABLED_OS2014_nonparametric_convergence_study, Block_SWIPDG_01_subdomain) {}
TEST(DISABLED_OS2014_nonparametric_convergence_study, Block_SWIPDG_04_subdomain) {}
TEST(DISABLED_OS2014_nonparametric_convergence_study, Block_SWIPDG_16_subdomain) {}
TEST(DISABLED_OS2014_nonparametric_convergence_study, Block_SWIPDG_64_subdomain) {}


#endif // HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_GRID_MULTISCALE
