// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#define DUNE_STUFF_FUNCTIONS_DISABLE_CHECKS 1
#ifndef DUNE_STUFF_TEST_MAIN_CATCH_EXCEPTIONS
# define DUNE_STUFF_TEST_MAIN_CATCH_EXCEPTIONS 1
#endif
#ifndef DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING
# define DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING 1
#endif

// This one has to come first (includes the config.h)!
#include <dune/stuff/test/main.hxx>

#if HAVE_DUNE_GRID_MULTISCALE && HAVE_DUNE_FEM && HAVE_ALUGRID

# include <dune/stuff/common/disable_warnings.hh>
#   include <dune/grid/alugrid.hh>
# include <dune/stuff/common/reenable_warnings.hh>

# include <dune/stuff/common/exceptions.hh>

# include <dune/stuff/common/print.hh>
# include <dune/stuff/common/float_cmp.hh>

# include <dune/hdd/playground/linearelliptic/testcases/ESV2007.hh>

# include "linearelliptic-block-swipdg.hh"

using namespace Dune;
using namespace HDD;


typedef ALUGrid< 2, 2, simplex, conforming > AluConform2dGridType;

typedef testing::Types< LinearElliptic::TestCases::ESV2007Multiscale< AluConform2dGridType >
                      > AluConform2dTestCases;


template< class TestCaseType >
struct linearelliptic_SWIPDG_discretization
  : public ::testing::Test
{
  template< Stuff::LA::ChooseBackend la_backend >
  static void eoc_study(const std::string partitioning = "[1 1 1]")
  {
    const TestCaseType test_case(partitioning);
    test_case.print_header(DSC_LOG_INFO);
    DSC_LOG_INFO << std::endl;
    LinearElliptic::Tests::BlockSWIPDGStudy< TestCaseType, 1, la_backend > eoc_study(test_case);
    Stuff::Test::check_eoc_study_for_success(eoc_study, eoc_study.run_eoc(DSC_LOG_INFO));
  } // ... eoc_study()
}; // linearelliptic_SWIPDG_discretization


TYPED_TEST_CASE(linearelliptic_SWIPDG_discretization, AluConform2dTestCases);

# if HAVE_DUNE_ISTL
TYPED_TEST(linearelliptic_SWIPDG_discretization, eoc_study_using_istl_and_1_subdomain) {
  this->template eoc_study< Stuff::LA::ChooseBackend::istl_sparse >("[1 1 1]");
}
TYPED_TEST(linearelliptic_SWIPDG_discretization, eoc_study_using_istl_and_4_subdomain) {
  this->template eoc_study< Stuff::LA::ChooseBackend::istl_sparse >("[2 2 1]");
}
TYPED_TEST(linearelliptic_SWIPDG_discretization, eoc_study_using_istl_and_16_subdomain) {
  this->template eoc_study< Stuff::LA::ChooseBackend::istl_sparse >("[4 4 1]");
}
TYPED_TEST(linearelliptic_SWIPDG_discretization, eoc_study_using_istl_and_64_subdomain) {
  this->template eoc_study< Stuff::LA::ChooseBackend::istl_sparse >("[8 8 1]");
}
# elif HAVE_EIGEN // HAVE_DUNE_ISTL
TYPED_TEST(linearelliptic_SWIPDG_discretization, eoc_study_using_istl_and_1_subdomain) {
  this->template eoc_study< Stuff::LA::ChooseBackend::eigen_sparse >("[1 1 1]");
}
TYPED_TEST(linearelliptic_SWIPDG_discretization, eoc_study_using_istl_and_4_subdomain) {
  this->template eoc_study< Stuff::LA::ChooseBackend::eigen_sparse >("[2 2 1]");
}
TYPED_TEST(linearelliptic_SWIPDG_discretization, eoc_study_using_istl_and_16_subdomain) {
  this->template eoc_study< Stuff::LA::ChooseBackend::eigen_sparse >("[4 4 1]");
}
TYPED_TEST(linearelliptic_SWIPDG_discretization, eoc_study_using_istl_and_64_subdomain) {
  this->template eoc_study< Stuff::LA::ChooseBackend::eigen_sparse >("[8 8 1]");
}
#else // HAVE_EIGEN // HAVE_DUNE_ISTL
TEST(DISABLED_linearelliptic_SWIPDG_discretization, eoc_study_using_istl_and_1_subdomain) {}
TEST(DISABLED_linearelliptic_SWIPDG_discretization, eoc_study_using_istl_and_4_subdomain) {}
TEST(DISABLED_linearelliptic_SWIPDG_discretization, eoc_study_using_istl_and_16_subdomain) {}
TEST(DISABLED_linearelliptic_SWIPDG_discretization, eoc_study_using_istl_and_64_subdomain) {}
# endif // HAVE_EIGEN // HAVE_DUNE_ISTL


#else // HAVE_DUNE_GRID_MULTISCALE && HAVE_DUNE_FEM && HAVE_ALUGRID


TEST(DISABLED_linearelliptic_SWIPDG_discretization, eoc_study_using_istl_and_1_subdomain) {}
TEST(DISABLED_linearelliptic_SWIPDG_discretization, eoc_study_using_istl_and_4_subdomain) {}
TEST(DISABLED_linearelliptic_SWIPDG_discretization, eoc_study_using_istl_and_16_subdomain) {}
TEST(DISABLED_linearelliptic_SWIPDG_discretization, eoc_study_using_istl_and_64_subdomain) {}


#endif // HAVE_DUNE_GRID_MULTISCALE && HAVE_DUNE_FEM && HAVE_ALUGRID


extern template class LinearElliptic::Tests::BlockSWIPDGStudyExpectations
    < LinearElliptic::TestCases::ESV2007Multiscale< ALUGrid< 2, 2, simplex, conforming > >, 1 >;
