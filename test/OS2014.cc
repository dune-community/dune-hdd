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

#include <sstream>

#if HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
# define ENABLE_ALUGRID 1
# include <dune/grid/alugrid.hh>
#else
# error This test requires ALUGrid!
#endif

#include <dune/stuff/common/exceptions.hh>

#include <dune/stuff/common/print.hh>
#include <dune/stuff/common/float_cmp.hh>

#include <dune/hdd/playground/linearelliptic/testcases/OS2014.hh>

#include "linearelliptic-swipdg.hh"
#include "linearelliptic-block-swipdg.hh"

using namespace Dune;
using namespace HDD;


class OS2014_nonparametric_convergence_study
  : public ::testing::Test
{
protected:
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

  typedef LinearElliptic::TestCases::ESV2007< GridType >           TestCaseType;
  typedef LinearElliptic::TestCases::ESV2007Multiscale< GridType > BlockTestCaseType;

  typedef LinearElliptic::Tests::EocStudySWIPDG< TestCaseType, 1, space_backend, la_backend > EocStudyType;
  typedef LinearElliptic::Tests::EocStudyBlockSWIPDG< BlockTestCaseType, 1, la_backend >      BlockEocStudyType;

  template< class StudyType >
  static void check_for_success(const StudyType& study,
                                const std::map< std::string, std::vector< double > >& errors_map)
  {
    for (const auto& norm : study.used_norms()) {
      const auto expected_results = study.expected_results(norm);
      const auto errors_search = errors_map.find(norm);
      EXPECT_NE(errors_search, errors_map.end())
          << "          norm = " << norm;
      const auto& errors = errors_search->second;
      EXPECT_LE(errors.size(), expected_results.size())
          << "          norm = " << norm;
      for (size_t ii = 0; ii < errors.size(); ++ii)
        EXPECT_LE(errors[ii], expected_results[ii])
            << "          norm = " << norm << ", level = " << ii;
    }
  } // ... check_for_success(...)
}; // class OS2014_nonparametric_convergence_study


TEST_F(OS2014_nonparametric_convergence_study, SWIPDG_fine_triangulation)
{
  const TestCaseType test_case;
  test_case.print_header(DSC_LOG_INFO);
  DSC_LOG_INFO << std::endl;
  EocStudyType eoc_study(test_case,
                         {"energy",
                          "eta_NC_ESV2007", "eta_R_ESV2007", "eta_DF_ESV2007", "eta_ESV2007", "eff_ESV2007"});
  check_for_success(eoc_study, eoc_study.run(false, DSC_LOG_INFO));
}
TEST_F(OS2014_nonparametric_convergence_study, SWIPDG_fine_triangulation_alternate_summation)
{
  const TestCaseType test_case;
  EocStudyType eoc_study(test_case,
                         {"energy",
                          "eta_ESV2007", "eff_ESV2007", "eta_ESV2007_alt", "eff_ESV2007_alt"});
  check_for_success(eoc_study, eoc_study.run(false, DSC_LOG_INFO));
}
TEST_F(OS2014_nonparametric_convergence_study, Block_SWIPDG_01_subdomain) {
  const BlockTestCaseType test_case("[1 1 1]");
  BlockEocStudyType eoc_study(test_case,
                              {"energy", "eta_NC_OS2014", "eta_R_OS2014", "eta_DF_OS2014", "eta_OS2014", "eff_OS2014"});
  check_for_success(eoc_study, eoc_study.run(false, DSC_LOG_INFO));
}
TEST_F(OS2014_nonparametric_convergence_study, Block_SWIPDG_04_subdomain) {
  const BlockTestCaseType test_case("[2 2 1]");
  BlockEocStudyType eoc_study(test_case,
                              {"energy", "eta_NC_OS2014", "eta_R_OS2014", "eta_DF_OS2014", "eta_OS2014", "eff_OS2014"});
  check_for_success(eoc_study, eoc_study.run(false, DSC_LOG_INFO));
}
TEST_F(OS2014_nonparametric_convergence_study, Block_SWIPDG_16_subdomain) {
  const BlockTestCaseType test_case("[4 4 1]");
  BlockEocStudyType eoc_study(test_case,
                              {"energy", "eta_NC_OS2014", "eta_R_OS2014", "eta_DF_OS2014", "eta_OS2014", "eff_OS2014"});
  check_for_success(eoc_study, eoc_study.run(false, DSC_LOG_INFO));
}
TEST_F(OS2014_nonparametric_convergence_study, Block_SWIPDG_64_subdomain) {
  const BlockTestCaseType test_case("[8 8 1]");
  BlockEocStudyType eoc_study(test_case,
                              {"energy", "eta_NC_OS2014", "eta_R_OS2014", "eta_DF_OS2014", "eta_OS2014", "eff_OS2014"});
  check_for_success(eoc_study, eoc_study.run(false, DSC_LOG_INFO));
}
