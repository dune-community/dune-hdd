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

// This one has to come first (includes the config.h)!
#include <dune/stuff/test/main.hxx>

#if HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_GRID_MULTISCALE

#include "OS2014.hh"


TEST(OS2014_nonparametric_convergence_study, SWIPDG_fine_triangulation) {
  OS2014_nonparametric_convergence_study__SWIPDG_fine_triangulation();
}
TEST(OS2014_nonparametric_convergence_study, SWIPDG_fine_triangulation_alternative_summation) {
  OS2014_nonparametric_convergence_study__SWIPDG_fine_triangulation_alternative_summation();
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

TEST(OS2014_parametric_convergence_study, eta_DF_comparison_01_subdomain)
{
  const std::string partitioning = "[1 1 1]";
  const std::vector< std::string > only_these_norms = {"eta_DF_OS2014", "eta_DF_OS2014_*", "eta_OS2014",
                                                       "eta_OS2014_*", "eff_OS2014_mu", "eff_OS2014_*_mu"};
  bool print_header = true;
  for (auto mu_hat_value : {0.1, 0.5, 1.0}) {
    const auto mu_hat = Parameter("mu", mu_hat_value);
    for (auto mu_value : {0.1, 0.3, 0.5, 0.75, 1.0}) {
      const auto mu = Parameter("mu", mu_value);
      const auto mu_bar = mu;
      parametric_convergence_study(partitioning,
                                   only_these_norms,
                                   {{"mu_hat",        mu_hat},
                                    {"mu_bar",        mu_bar},
                                    {"mu",            mu},
                                    {"mu_minimizing", Parameter("mu", 0.1)}},
                                   print_header);
      if (print_header)
        print_header = false;
    }
  }
} // OS2014_parametric_convergence_study, eta_DF_comparison_01_subdomain

TEST(OS2014_nonparametric_localization_study, SWIPDG_fine_triangulation) {
  nonparametric_localization_study__SWIPDG_fine_triangulation();
}
TEST(OS2014_nonparametric_localization_study, Block_SWIPDG_20_subdomain) {
  nonparametric_localization_study__Block_SWIPDG("[20 4 1]");
}


#else // HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_GRID_MULTISCALE


TEST(DISABLED_OS2014_nonparametric_convergence_study, SWIPDG_fine_triangulation) {}
TEST(DISABLED_OS2014_nonparametric_convergence_study, SWIPDG_fine_triangulation_alternative_summation) {}
TEST(DISABLED_OS2014_nonparametric_convergence_study, Block_SWIPDG_01_subdomain) {}
TEST(DISABLED_OS2014_nonparametric_convergence_study, Block_SWIPDG_04_subdomain) {}
TEST(DISABLED_OS2014_nonparametric_convergence_study, Block_SWIPDG_16_subdomain) {}
TEST(DISABLED_OS2014_nonparametric_convergence_study, Block_SWIPDG_64_subdomain) {}

TEST(DISABLED_OS2014_parametric_convergence_study, eta_DF_comparison_01_subdomain) {}

TEST(DISABLED_OS2014_nonparametric_localization_study, SWIPDG_fine_triangulation) {}
TEST(DISABLED_OS2014_nonparametric_localization_study, Block_SWIPDG_20_subdomain) {}


#endif // HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_GRID_MULTISCALE
