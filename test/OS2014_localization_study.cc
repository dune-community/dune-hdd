// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

//#ifndef DUNE_STUFF_TEST_MAIN_CATCH_EXCEPTIONS
//# define DUNE_STUFF_TEST_MAIN_CATCH_EXCEPTIONS 1
//#endif
#ifndef DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING
# define DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING 1
#endif

#define DUNE_HDD_LINEARELLIPTIC_TESTCASES_BASE_DISABLE_WARNING

// This one has to come first (includes the config.h)!
#include <dune/stuff/test/main.hxx>

#include <vector>
#include <string>
#include <map>

#if HAVE_ALUGRID
# include <dune/grid/alugrid.hh>

# include <dune/pymor/parameters/base.hh>

# include <dune/hdd/linearelliptic/testcases/spe10.hh>

# include "linearelliptic-block-swipdg.hh"
# include "linearelliptic-swipdg.hh"

using namespace Dune;
using namespace HDD;


typedef ALUGrid< 2, 2, simplex, conforming > GridType;

typedef LinearElliptic::TestCases::Spe10::Model1< GridType >            NonparametricTestCaseType;
typedef LinearElliptic::Tests::SWIPDGStudy< NonparametricTestCaseType > NonparametricStudyType;

typedef LinearElliptic::TestCases::Spe10::BlockModel1< GridType >                 NonparametricBlockTestCaseType;
typedef LinearElliptic::Tests::BlockSWIPDGStudy< NonparametricBlockTestCaseType > NonparametricBlockStudyType;

typedef LinearElliptic::TestCases::Spe10::ParametricBlockModel1< GridType > ParametricTestCaseType;
typedef LinearElliptic::Tests::BlockSWIPDGStudy< ParametricTestCaseType >   ParametricStudyType;


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
               << "| mu_minimizing = " << parameters.at("mu_minimizing") << "\n"
               << "+==========================================================+\n";
} // ... print_parameter_information(...)


template< class TestCaseType, class StudyType >
void run_parametric_localization_study(const std::string partitioning,
                                       const std::vector< std::string >& only_these_indicators,
                                       const std::map< std::string, Pymor::Parameter >& parameters,
                                       const bool print_header,
                                       const std::string visualization)
{
  const TestCaseType test_case(parameters, partitioning);
  if (print_header)
    test_case.print_header(DSC_LOG_INFO);
  print_parameter_information(test_case);
  StudyType study(test_case, {}, only_these_indicators, visualization);
  study.run_localization(DSC_LOG_INFO);
} // ... run_parametric_localization_study(...)


TEST(OS2014_nonparametric_localization_study, SWIPDG_fine_triangulation) {
  const NonparametricTestCaseType test_case;
  test_case.print_header(DSC_LOG_INFO);
  DSC_LOG_INFO << std::endl;
  NonparametricStudyType study(test_case,
                               {},
                               {"eta_ESV2007", "eta_ESV2007_alt"},
                               "OS2014_nonparametric_localization_study_swipdg_fine_triangulation");
  study.run_localization(DSC_LOG_INFO);
} // TEST(OS2014_nonparametric_localization_study, SWIPDG_fine_triangulation)


TEST(OS2014_nonparametric_localization_study, Block_SWIPDG_80_subdomains) {
  const NonparametricBlockTestCaseType test_case("[20 4 1]");
  NonparametricBlockStudyType study(test_case,
                                    {},
                                    {"eta_OS2014"},
                                    "OS2014_nonparametric_localization_study_block_swipdg_80_subdomains");
  study.run_localization(DSC_LOG_INFO);
} // TEST(OS2014_nonparametric_localization_study, Block_SWIPDG_80_subdomains)


TEST(OS2014_parametric_localization_study, Block_SWIPDG_80_subdomains)
{
  const std::string partitioning = "[20 4 1]";
  const std::vector< std::string > only_these_indicators = {/*"eta_OS2014, eta_OS2014_*"*/};
  const std::string visualization_prefix = "OS2014_parametric_localization_study_block_swipdg_80_subdomain";
  bool print_header = true;
  for (auto mu_hat_value : {/*0.1, 0.5,*/ 1.0}) {
    const auto mu_hat = Pymor::Parameter("mu", mu_hat_value);
    for (auto mu_value : {0.1/*, 0.3, 0.5, 0.75, 1.0*/}) {
      const auto mu = Pymor::Parameter("mu", mu_value);
      const auto mu_bar = mu;
      run_parametric_localization_study< ParametricTestCaseType,
                                         ParametricStudyType >(partitioning,
                                                               only_these_indicators,
                                                               {{"mu_hat",        mu_hat},
                                                                {"mu_bar",        mu_bar},
                                                                {"mu",            mu},
                                                                {"mu_minimizing", Pymor::Parameter("mu", 0.1)}},
                                                               print_header,
                                                               visualization_prefix + "_" + DSC::toString(mu_value));
      if (print_header)
        print_header = false;
    }
  }
} // TEST(OS2014_parametric_localization_study, Block_SWIPDG_80_subdomains)


#else // HAVE_ALUGRID


TEST(DISABLED_OS2014_nonparametric_localization_study, SWIPDG_fine_triangulation) {}
TEST(DISABLED_OS2014_nonparametric_localization_study, Block_SWIPDG_80_subdomains) {}
TEST(DISABLED_OS2014_parametric_localization_study, Block_SWIPDG_80_subdomains) {}


#endif // HAVE_ALUGRID
