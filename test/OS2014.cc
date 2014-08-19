// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#define DUNE_STUFF_FUNCTIONS_DISABLE_CHECKS 1
#define DUNE_HDD_LINEARELLIPTIC_TESTCASES_BASE_DISABLE_WARNING 1

// This one has to come first (includes the config.h)!
#include <dune/stuff/test/test_common.hh>

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

// change this to toggle output
std::ostream& test_out = std::cout;
//std::ostream& test_out = DSC_LOG.devnull();


using namespace Dune;
using namespace HDD;


class OS2014_nonparametric_convergence_study
  : public ::testing::Test
{
  template< class StudyType >
  static void check_for_success(const StudyType& study, std::map< std::string, std::vector< double > >& errors_map)
  {
    for (const auto& norm : study.provided_norms())
    {
      const auto errors = errors_map[norm];
      const auto expected_results = study.expected_results(norm);
      assert(expected_results.size() <= errors.size());
      for (size_t ii = 0; ii < errors.size(); ++ii)
        EXPECT_LE(errors[ii], expected_results[ii]) << "          'norm' = " << norm << ", level = " << ii;
    }
  } // ... check_for_success(...)

  static const GDT::ChooseSpaceBackend  space_backend = GDT::ChooseSpaceBackend::fem;
  static const Stuff::LA::ChooseBackend la_backend    = Stuff::LA::ChooseBackend::eigen_sparse;

  typedef ALUGrid< 2, 2, simplex, conforming > GridType;

  typedef LinearElliptic::TestCases::OS2014< GridType >           TestCaseType;
//  typedef LinearElliptic::TestCases::ESV2007Multiscale< GridType > BlockTestCaseType;

  typedef LinearElliptic::Tests::EocStudySWIPDG< TestCaseType, 1, space_backend, la_backend > EocStudyType;
//  typedef LinearElliptic::Tests::EocStudyBlockSWIPDG< BlockTestCaseType, 1, la_backend >      BlockEocStudyType;

protected:
  void SWIPDG_fine_triangulation() const
  {
    const TestCaseType test_case;
    test_case.print_header(test_out);
    test_out << std::endl;
    EocStudyType eoc_study(test_case);
    auto errors = eoc_study.run(false, test_out);
    check_for_success(eoc_study, errors);
  } // ... SWIPDG_fine_triangulation(...)

//  bool BlockSWIPDG_coarse_triangulation(const std::string partitioning) const
//  {
//    const BlockTestCaseType test_case(partitioning);
//    BlockEocStudyType eoc_study(test_case);
//    auto errors = eoc_study.run(false, test_out);
//    return check_for_success(eoc_study, errors);
//  } // ... BlockSWIPDG_coarse_triangulation(...)
}; // class OS2014_nonparametric_convergence_study


TEST_F(OS2014_nonparametric_convergence_study, SWIPDG_fine_triangulation) {
  SWIPDG_fine_triangulation();
}


#include <dune/stuff/test/test_main.cxx>
