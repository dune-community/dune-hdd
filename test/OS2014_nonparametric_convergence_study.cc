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

// This one has to come first (includes the config.h)!
#include <dune/stuff/test/main.hxx>

#if HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_GRID_MULTISCALE

#include "OS2014_nonparametric_convergence_study.hh"


TEST(OS2014_nonparametric_convergence_study, SWIPDG_fine_triangulation) {
  nonparametric_convergence_study(/*"nonparametric_convergence_study"*/);
}
TEST(OS2014_nonparametric_convergence_study, SWIPDG_fine_triangulation_alternative_summation) {
  nonparametric_convergence_study_alternative_summation();
}
TEST(OS2014_nonparametric_convergence_study, Block_SWIPDG_01_subdomain) {
  nonparametric_block_convergence_study("[1 1 1]");
}
TEST(OS2014_nonparametric_convergence_study, Block_SWIPDG_04_subdomains) {
  nonparametric_block_convergence_study("[2 2 1]");
}
TEST(OS2014_nonparametric_convergence_study, Block_SWIPDG_16_subdomains) {
  nonparametric_block_convergence_study("[4 4 1]");
}
TEST(OS2014_nonparametric_convergence_study, Block_SWIPDG_64_subdomains) {
  nonparametric_block_convergence_study("[8 8 1]");
}


#else // HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_GRID_MULTISCALE


TEST(DISABLED_OS2014_nonparametric_convergence_study, SWIPDG_fine_triangulation) {}
TEST(DISABLED_OS2014_nonparametric_convergence_study, SWIPDG_fine_triangulation_alternative_summation) {}
TEST(DISABLED_OS2014_nonparametric_convergence_study, Block_SWIPDG_01_subdomain) {}
TEST(DISABLED_OS2014_nonparametric_convergence_study, Block_SWIPDG_04_subdomains) {}
TEST(DISABLED_OS2014_nonparametric_convergence_study, Block_SWIPDG_16_subdomains) {}
TEST(DISABLED_OS2014_nonparametric_convergence_study, Block_SWIPDG_64_subdomains) {}


#endif // HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_GRID_MULTISCALE
