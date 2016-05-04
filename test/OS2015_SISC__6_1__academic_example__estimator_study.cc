// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING
# define DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING 1
#endif
#ifndef DUNE_HDD_LINEARELLIPTIC_TESTCASES_BASE_DISABLE_WARNING
# define DUNE_HDD_LINEARELLIPTIC_TESTCASES_BASE_DISABLE_WARNING
#endif

// This one has to come first (includes the config.h)!
#include <dune/stuff/test/main.hxx>

#if HAVE_DUNE_GRID_MULTISCALE && HAVE_DUNE_FEM && HAVE_DUNE_ISTL && HAVE_ALUGRID
# include <vector>
# include <string>
# include <map>

# include <dune/grid/alugrid.hh>

# include <dune/pymor/parameters/base.hh>

# include <dune/hdd/linearelliptic/testcases/OS2015.hh>

# include "linearelliptic-block-swipdg.hh"
# include "linearelliptic-block-swipdg-expectations.hh"

using namespace Dune;
using namespace HDD;


typedef ALUGrid< 2, 2, simplex, conforming > GridType;

typedef LinearElliptic::TestCases::OS2015::Academic< GridType > TestCaseType;
typedef LinearElliptic::Tests::BlockSWIPDGStudy< TestCaseType >                   StudyType;


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
               << "| mu_hat        = " << parameters.at("mu_hat") << "\n";
  const double alpha_bar = diffusion_factor.alpha(parameters.at("mu"), parameters.at("mu_bar"));
  const double alpha_hat = diffusion_factor.alpha(parameters.at("mu"), parameters.at("mu_hat"));
  const double gamma_bar = diffusion_factor.gamma(parameters.at("mu"), parameters.at("mu_bar"));
  DSC_LOG_INFO << "| alpha(mu, mu_bar)^-1/2 = " << std::setprecision(2) << std::scientific
                                                << 1.0/std::sqrt(alpha_bar) << "\n"
               << "| alpha(mu, mu_hat)^-1/2 = " << std::setprecision(2) << std::scientific
                                                << 1.0/std::sqrt(alpha_hat) << "\n"
               << "| gamma(mu, mu_bar)^1/2  = " << std::setprecision(2) << std::scientific
                                                << std::sqrt(gamma_bar) << "\n"
               << "+-------------------------\n"
               << "| ";
} // ... print_parameter_information(...)


void run_eoc_study(const std::string partitioning,
                   const std::vector< std::string >& only_these_norms,
                   const std::map< std::string, Pymor::Parameter >& parameters,
                   const bool print_header,
                   const std::string visualization,
                   const bool H_with_h)
{
  const TestCaseType test_case(parameters,
                               partitioning,
                               3,
                               0,
                               H_with_h);
  if (print_header)
    test_case.print_header(DSC_LOG_INFO);
  print_parameter_information(test_case);
  StudyType study(test_case, only_these_norms, {}, visualization);
  Stuff::Test::check_eoc_study_for_success(study, study.run_eoc(DSC_LOG_INFO));
} // ... run_eoc_study(...)


TEST(OS2015_SISC__6_1__academic_example__estimator_study, table_1)
{
  using Pymor::Parameter;
  run_eoc_study("[1 1 1]",
                {"energy_mu", "eta_NC_OS2014", "eta_R_OS2014_*", "eta_DF_OS2014_*", "eta_OS2014_*", "eff_OS2014_*_mu"},
                {{"mu_hat", Parameter("mu", 1)},
                 {"mu_bar", Parameter("mu", 1)},
                 {"mu",     Parameter("mu", 1)}},
                true,
                "",
                false);
} // TEST(OS2015_SISC__6_1__academic_example__estimator_study, table_1)


TEST(OS2015_SISC__6_1__academic_example__estimator_study, table_2_columns_2_3_4)
{
  using Pymor::Parameter;
  run_eoc_study("[2 2 1]",
                {"eta_R_OS2014_*", "eta_OS2014_*", "eff_OS2014_*_mu"},
                {{"mu_hat", Parameter("mu", 1)},
                 {"mu_bar", Parameter("mu", 1)},
                 {"mu",     Parameter("mu", 1)}},
                true,
                "",
                true);
} // TEST(OS2015_SISC__6_1__academic_example__estimator_study, table_2_columns_2_3_4)


TEST(OS2015_SISC__6_1__academic_example__estimator_study, table_2_columns_5_6_7)
{
  using Pymor::Parameter;
  run_eoc_study("[2 2 1]",
                {"eta_DF_OS2014_*", "eta_OS2014_*", "eff_OS2014_*_mu"},
                {{"mu_hat", Parameter("mu", 0.1)},
                 {"mu_bar", Parameter("mu", 1)},
                 {"mu",     Parameter("mu", 1)}},
                true,
                "",
                true);
} // TEST(OS2015_SISC__6_1__academic_example__estimator_study, table_2_columns_5_6_7)


TEST(OS2015_SISC__6_1__academic_example__estimator_study, table_3)
{
  using Pymor::Parameter;
  run_eoc_study("[2 2 1]",
                {"energy_mu_bar", "eta_NC_OS2014", "eta_DF_OS2014_*", "eta_OS2014_*", "eff_OS2014_*_mu_bar"},
                {{"mu_hat", Parameter("mu", 0.1)},
                 {"mu_bar", Parameter("mu", 0.1)},
                 {"mu",     Parameter("mu", 1)}},
                true,
                "",
                true);
} // TEST(OS2015_SISC__6_1__academic_example__estimator_study, table_3)


extern template class Dune::HDD::LinearElliptic::Tests::BlockSWIPDGStudyExpectations< TestCaseType >;


#else // HAVE_DUNE_GRID_MULTISCALE && HAVE_DUNE_FEM && HAVE_DUNE_ISTL && HAVE_ALUGRID


TEST(DISABLED_OS2015_SISC__6_1__academic_example__estimator_study, table_1)
{
  std::cerr << "You are missing dune-fem or dune-grid-multiscale or alugrid!" << std::endl;
}
TEST(DISABLED_OS2015_SISC__6_1__academic_example__estimator_study, table_2__columns_2_3_4)
{
  std::cerr << "You are missing dune-fem or dune-grid-multiscale or alugrid!" << std::endl;
}
TEST(DISABLED_OS2015_SISC__6_1__academic_example__estimator_study, table_2__columns_5_6_7)
{
  std::cerr << "You are missing dune-fem or dune-grid-multiscale or alugrid!" << std::endl;
}
TEST(DISABLED_OS2015_SISC__6_1__academic_example__estimator_study, table_3)
{
  std::cerr << "You are missing dune-fem or dune-grid-multiscale or alugrid!" << std::endl;
}


#endif // HAVE_DUNE_GRID_MULTISCALE && HAVE_DUNE_FEM && HAVE_DUNE_ISTL && HAVE_ALUGRID
