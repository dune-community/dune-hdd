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

#if HAVE_DUNE_GRID_MULTISCALE && HAVE_DUNE_FEM

# include <dune/stuff/common/disable_warnings.hh>
# if HAVE_ALUGRID
#   include <dune/grid/alugrid.hh>
# endif
# include <dune/grid/sgrid.hh>
# include <dune/stuff/common/reenable_warnings.hh>

# include <dune/stuff/common/exceptions.hh>

# include <dune/hdd/linearelliptic/testcases/ESV2007.hh>

# include "linearelliptic-block-swipdg.hh"

using namespace Dune;
using namespace HDD;


template< class TestCaseType >
struct linearelliptic_BlockSWIPDG_discretization
  : public ::testing::Test
{
  template< Stuff::LA::ChooseBackend la_backend >
  static void eoc_study(const std::string partitioning = "[1 1 1]")
  {
    const TestCaseType test_case(partitioning, 2);
    test_case.print_header(DSC_LOG_INFO);
    DSC_LOG_INFO << std::endl;
    LinearElliptic::Tests::BlockSWIPDGStudy< TestCaseType, 1, la_backend > eoc_study(test_case);
    Stuff::Test::check_eoc_study_for_success(eoc_study, eoc_study.run_eoc(DSC_LOG_INFO));
  } // ... eoc_study()
}; // linearelliptic_BlockSWIPDG_discretization


typedef SGrid< 2, 2 > SGridType;
# if HAVE_ALUGRID
typedef ALUGrid< 2, 2, simplex, conforming > AluConform2dGridType;
# endif // HAVE_ALUGRID

typedef testing::Types< LinearElliptic::TestCases::ESV2007Multiscale< SGridType >
# if HAVE_ALUGRID
                      , LinearElliptic::TestCases::ESV2007Multiscale< AluConform2dGridType >
# endif
                      > TestCases;
TYPED_TEST_CASE(linearelliptic_BlockSWIPDG_discretization, TestCases);


# if HAVE_EIGEN
TYPED_TEST(linearelliptic_BlockSWIPDG_discretization, eoc_study_using_eigen_on_1_subdomain) {
  this->template eoc_study< Stuff::LA::ChooseBackend::eigen_sparse >("[1 1 1]");
}
TYPED_TEST(linearelliptic_BlockSWIPDG_discretization, eoc_study_using_eigen_on_4_subdomains) {
  this->template eoc_study< Stuff::LA::ChooseBackend::eigen_sparse >("[2 2 1]");
}
TYPED_TEST(linearelliptic_BlockSWIPDG_discretization, eoc_study_using_eigen_on_16_subdomains) {
  this->template eoc_study< Stuff::LA::ChooseBackend::eigen_sparse >("[4 4 1]");
}
TYPED_TEST(linearelliptic_BlockSWIPDG_discretization, eoc_study_using_eigen_on_64_subdomains) {
  this->template eoc_study< Stuff::LA::ChooseBackend::eigen_sparse >("[8 8 1]");
}
# else // HAVE_EIGEN
TEST(DISABLED_linearelliptic_BlockSWIPDG_discretization, eoc_study_using_eigen_on_1_subdomain) {}
TEST(DISABLED_linearelliptic_BlockSWIPDG_discretization, eoc_study_using_eigen_on_4_subdomains) {}
TEST(DISABLED_linearelliptic_BlockSWIPDG_discretization, eoc_study_using_eigen_on_16_subdomains) {}
TEST(DISABLED_linearelliptic_BlockSWIPDG_discretization, eoc_study_using_eigen_on_64_subdomains) {}
# endif // HAVE_EIGEN

# if HAVE_DUNE_ISTL
TYPED_TEST(linearelliptic_BlockSWIPDG_discretization, eoc_study_using_istl_on_1_subdomain) {
  this->template eoc_study< Stuff::LA::ChooseBackend::istl_sparse >("[1 1 1]");
}
TYPED_TEST(linearelliptic_BlockSWIPDG_discretization, eoc_study_using_istl_on_4_subdomains) {
  this->template eoc_study< Stuff::LA::ChooseBackend::istl_sparse >("[2 2 1]");
}
TYPED_TEST(linearelliptic_BlockSWIPDG_discretization, eoc_study_using_istl_on_16_subdomains) {
  this->template eoc_study< Stuff::LA::ChooseBackend::istl_sparse >("[4 4 1]");
}
TYPED_TEST(linearelliptic_BlockSWIPDG_discretization, eoc_study_using_istl_on_64_subdomains) {
  this->template eoc_study< Stuff::LA::ChooseBackend::istl_sparse >("[8 8 1]");
}
# else // HAVE_DUNE_ISTL
TEST(DISABLED_linearelliptic_BlockSWIPDG_discretization, eoc_study_using_istl_on_1_subdomains) {}
TEST(DISABLED_linearelliptic_BlockSWIPDG_discretization, eoc_study_using_istl_on_4_subdomains) {}
TEST(DISABLED_linearelliptic_BlockSWIPDG_discretization, eoc_study_using_istl_on_16_subdomains) {}
TEST(DISABLED_linearelliptic_BlockSWIPDG_discretization, eoc_study_using_istl_on_64_subdomains) {}
# endif // HAVE_DUNE_ISTL


namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Tests {

//extern template class LinearElliptic::Tests::BlockSWIPDGStudyExpectations
//    < LinearElliptic::TestCases::ESV2007Multiscale< SGridType >, 1 >;

# if HAVE_ALUGRID

extern template class LinearElliptic::Tests::BlockSWIPDGStudyExpectations
    < LinearElliptic::TestCases::ESV2007Multiscale< AluConform2dGridType >, 1 >;

# endif // HAVE_ALUGRID

} // namespace Tests
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#else // HAVE_DUNE_GRID_MULTISCALE && HAVE_DUNE_FEM


TEST(DISABLED_linearelliptic_BlockSWIPDG_discretization, eoc_study_using_eigen_on_1_subdomain) {}
TEST(DISABLED_linearelliptic_BlockSWIPDG_discretization, eoc_study_using_eigen_on_4_subdomains) {}
TEST(DISABLED_linearelliptic_BlockSWIPDG_discretization, eoc_study_using_eigen_on_16_subdomains) {}
TEST(DISABLED_linearelliptic_BlockSWIPDG_discretization, eoc_study_using_eigen_on_64_subdomains) {}

TEST(DISABLED_linearelliptic_BlockSWIPDG_discretization, eoc_study_using_istl_on_1_subdomains) {}
TEST(DISABLED_linearelliptic_BlockSWIPDG_discretization, eoc_study_using_istl_on_4_subdomains) {}
TEST(DISABLED_linearelliptic_BlockSWIPDG_discretization, eoc_study_using_istl_on_16_subdomains) {}
TEST(DISABLED_linearelliptic_BlockSWIPDG_discretization, eoc_study_using_istl_on_64_subdomains) {}


#endif // HAVE_DUNE_GRID_MULTISCALE && HAVE_DUNE_FEM
