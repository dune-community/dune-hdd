// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_TEST_MAIN_CATCH_EXCEPTIONS
# define DUNE_STUFF_TEST_MAIN_CATCH_EXCEPTIONS 1
#endif
#ifndef DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING
# define DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING 1
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

TEST(OS2014_parametric_convergence_study, Block_SWIPDG_01_subdomain)
{
  parametric_convergence_study("[1 1 1]",
                               {"energy_mu", "eta_NC_OS2014", "eta_R_OS2014", "eta_DF_OS2014", "eta_OS2014",
                                "eff_OS2014_mu"},
                               {{"mu_hat", Parameter("mu", 1)},
                                {"mu_bar", Parameter("mu", 1)},
                                {"mu", Parameter("mu", 1)},
                                {"mu_minimizing", Parameter("mu", 1)}});
}


#else // HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_GRID_MULTISCALE


TEST(DISABLED_OS2014_nonparametric_convergence_study, SWIPDG_fine_triangulation) {}
TEST(DISABLED_OS2014_nonparametric_convergence_study, SWIPDG_fine_triangulation_alternative_summation) {}
TEST(DISABLED_OS2014_nonparametric_convergence_study, Block_SWIPDG_01_subdomain) {}
TEST(DISABLED_OS2014_nonparametric_convergence_study, Block_SWIPDG_04_subdomain) {}
TEST(DISABLED_OS2014_nonparametric_convergence_study, Block_SWIPDG_16_subdomain) {}
TEST(DISABLED_OS2014_nonparametric_convergence_study, Block_SWIPDG_64_subdomain) {}


#endif // HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_GRID_MULTISCALE
