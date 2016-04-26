// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#define DUNE_STUFF_FUNCTIONS_DISABLE_CHECKS 1
#define DUNE_HDD_LINEARELLIPTIC_TESTCASES_BASE_DISABLE_WARNING 1
#ifndef DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING
# define DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING 1
#endif

// This one has to come first (includes the config.h)!
#include <dune/stuff/test/main.hxx>

#include <sstream>

#if HAVE_DUNE_GRID_MULTISCALE && HAVE_DUNE_FEM && HAVE_DUNE_ISTL && HAVE_ALUGRID

# include <dune/stuff/common/disable_warnings.hh>
#   include <dune/grid/alugrid.hh>
# include <dune/stuff/common/reenable_warnings.hh>

# include <dune/stuff/common/exceptions.hh>

# include <dune/stuff/common/print.hh>
# include <dune/stuff/common/float_cmp.hh>

# include <dune/hdd/linearelliptic/testcases/ESV2007.hh>

# include "linearelliptic-swipdg.hh"
# include "linearelliptic-block-swipdg.hh"


using namespace Dune;
using namespace HDD;


class OS2014_FVCA7_Poster_convergence_study
  : public ::testing::Test
{
  typedef ALUGrid< 2, 2, simplex, conforming > GridType;

  typedef LinearElliptic::TestCases::ESV2007< GridType >           TestCaseType;
  typedef LinearElliptic::TestCases::ESV2007Multiscale< GridType > BlockTestCaseType;

  static const GDT::ChooseSpaceBackend  space_backend = GDT::ChooseSpaceBackend::fem;
  static const Stuff::LA::ChooseBackend la_backend    = Stuff::LA::ChooseBackend::istl_sparse;

  typedef LinearElliptic::Tests::SWIPDGStudy< TestCaseType, 1, space_backend, la_backend > EocStudyType;
  typedef LinearElliptic::Tests::BlockSWIPDGStudy< BlockTestCaseType, 1, la_backend >      BlockEocStudyType;

protected:
  void ESV2007_fine_triangulation() const
  {
    const TestCaseType test_case;
    test_case.print_header(DSC_LOG_INFO);
    DSC_LOG_INFO << std::endl;
    EocStudyType eoc_study(test_case, {"energy", "eta_ESV2007"});
    Stuff::Test::check_eoc_study_for_success(eoc_study, eoc_study.run_eoc(DSC_LOG_INFO));
  } // ... ESV2007_fine_triangulation(...)

  void BlockSWIPDG_coarse_triangulation(const std::string partitioning) const
  {
    const BlockTestCaseType test_case(partitioning);
    BlockEocStudyType eoc_study(test_case, {"energy", "eta_OS2014"});
    Stuff::Test::check_eoc_study_for_success(eoc_study, eoc_study.run_eoc(DSC_LOG_INFO));
  } // ... BlockSWIPDG_coarse_triangulation(...)
}; // class OS2014_FVCA7_Poster_convergence_study


TEST_F(OS2014_FVCA7_Poster_convergence_study, ESV2007_fine_triangulation) {
  ESV2007_fine_triangulation();
}
TEST_F(OS2014_FVCA7_Poster_convergence_study, BlockSWIPDG_01_subdomain) {
  BlockSWIPDG_coarse_triangulation("[1 1 1]");
}
TEST_F(OS2014_FVCA7_Poster_convergence_study, BlockSWIPDG_04_subdomain) {
  BlockSWIPDG_coarse_triangulation("[2 2 1]");
}
TEST_F(OS2014_FVCA7_Poster_convergence_study, BlockSWIPDG_16_subdomain) {
  BlockSWIPDG_coarse_triangulation("[4 4 1]");
}
TEST_F(OS2014_FVCA7_Poster_convergence_study, BlockSWIPDG_64_subdomain) {
  BlockSWIPDG_coarse_triangulation("[8 8 1]");
}


extern template class LinearElliptic::Tests::SWIPDGStudyExpectations
    < LinearElliptic::TestCases::ESV2007< ALUGrid< 2, 2, simplex, conforming > >, 1 >;

extern template class LinearElliptic::Tests::BlockSWIPDGStudyExpectations
    < LinearElliptic::TestCases::ESV2007Multiscale< ALUGrid< 2, 2, simplex, conforming > >, 1 >;


#else // HAVE_DUNE_GRID_MULTISCALE && HAVE_DUNE_FEM && HAVE_DUNE_ISTL && HAVE_ALUGRID


TEST(DISABLED_OS2014_FVCA7_Poster_convergence_study, ESV2007_fine_triangulation) {}
TEST(DISABLED_OS2014_FVCA7_Poster_convergence_study, BlockSWIPDG_01_subdomain) {}
TEST(DISABLED_OS2014_FVCA7_Poster_convergence_study, BlockSWIPDG_04_subdomain) {}
TEST(DISABLED_OS2014_FVCA7_Poster_convergence_study, BlockSWIPDG_16_subdomain) {}
TEST(DISABLED_OS2014_FVCA7_Poster_convergence_study, BlockSWIPDG_64_subdomain) {}


#endif // HAVE_DUNE_GRID_MULTISCALE && HAVE_DUNE_FEM && HAVE_DUNE_ISTL && HAVE_ALUGRID
