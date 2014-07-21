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
  static std::vector< double > truncate_vector_(const std::vector< double >& in, const size_t size)
  {
    assert(size <= in.size());
    if (size == in.size())
      return in;
    else {
      std::vector< double > ret(size);
      for (size_t ii = 0; ii < size; ++ii)
        ret[ii] = in[ii];
      return ret;
    }
  } // ... truncate_vector(...)

  template< class StudyType >
  static bool check_for_success_(const StudyType& study, std::map< std::string, std::vector< double > >& errors)
  {
    size_t failures = 0;
    std::stringstream ss;
    ss << "\n";
    for (const auto& norm : study.provided_norms())
      if (!Dune::Stuff::Common::FloatCmp::lt(errors[norm],
                                             truncate_vector_(study.expected_results(norm), errors[norm].size()))) {
        ++failures;
        Dune::Stuff::Common::print(errors[norm],                 "errors           (" + norm + ")", ss);
        Dune::Stuff::Common::print(study.expected_results(norm), "expected results (" + norm + ")", ss);
      }
    if (failures)
      DUNE_THROW(errors_are_not_as_expected, ss.str());
    return true;
  } // ... check_for_success_(...)

  static const GDT::ChooseSpaceBackend  space_backend = GDT::ChooseSpaceBackend::fem;
  static const Stuff::LA::ChooseBackend la_backend    = Stuff::LA::ChooseBackend::eigen_sparse;

  typedef ALUGrid< 2, 2, simplex, conforming > GridType;

  typedef LinearElliptic::TestCases::OS2014< GridType >           TestCaseType;
//  typedef LinearElliptic::TestCases::ESV2007Multiscale< GridType > BlockTestCaseType;

  typedef LinearElliptic::Tests::EocStudySWIPDG< TestCaseType, 1, space_backend, la_backend > EocStudyType;
//  typedef LinearElliptic::Tests::EocStudyBlockSWIPDG< BlockTestCaseType, 1, la_backend >      BlockEocStudyType;

protected:
  bool SWIPDG_fine_triangulation() const
  {
    const TestCaseType test_case;
    test_case.print_header(test_out);
    test_out << std::endl;
    EocStudyType eoc_study(test_case);
    auto errors = eoc_study.run(false, test_out);
    return check_for_success_(eoc_study, errors);
  } // ... SWIPDG_fine_triangulation(...)

//  bool BlockSWIPDG_coarse_triangulation(const std::string partitioning) const
//  {
//    const BlockTestCaseType test_case(partitioning);
//    BlockEocStudyType eoc_study(test_case);
//    auto errors = eoc_study.run(false, test_out);
//    return check_for_success_(eoc_study, errors);
//  } // ... BlockSWIPDG_coarse_triangulation(...)
}; // class OS2014_nonparametric_convergence_study


TEST_F(OS2014_nonparametric_convergence_study, SWIPDG_fine_triangulation) {
  EXPECT_TRUE(SWIPDG_fine_triangulation());
}


int main(int argc, char** argv)
{
  try {
    test_init(argc, argv);
    return RUN_ALL_TESTS();
  } catch (Dune::Exception& e) {
    std::cerr << "\nDune reported error: " << e.what() << std::endl;
    std::abort();
  } catch (std::exception& e) {
    std::cerr << "\n" << e.what() << std::endl;
    std::abort();
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
    std::abort();
  } // try
}
