// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#define DUNE_STUFF_FUNCTIONS_DISABLE_CHECKS

// This one has to come first (includes the config.h)!
#include <dune/stuff/test/main.hh>

#if HAVE_ALUGRID

#include <dune/grid/alugrid.hh>

#include <dune/stuff/common/exceptions.hh>

#include <dune/stuff/common/print.hh>
#include <dune/stuff/common/float_cmp.hh>

#include <dune/hdd/playground/linearelliptic/testcases/ESV2007.hh>

#include "linearelliptic-block-swipdg.hh"

using namespace Dune;
using namespace HDD;


// change this to toggle output
std::ostream& test_out = std::cout;
//std::ostream& test_out = DSC_LOG.devnull();


typedef Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming > AluConform2dGridType;

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
    test_case.print_header(test_out);
    test_out << std::endl;
    LinearElliptic::Tests::EocStudyBlockSWIPDG< TestCaseType, 1, la_backend > eoc_study(test_case);
    auto errors = eoc_study.run(false, test_out);
    for (const auto& norm : eoc_study.provided_norms())
      if (!Dune::Stuff::Common::FloatCmp::lt(errors[norm],
                                             truncate_vector(eoc_study.expected_results(norm), errors[norm].size()))) {
        std::stringstream ss;
        Dune::Stuff::Common::print(errors[norm],                        "errors           (" + norm + ")", ss);
        Dune::Stuff::Common::print(eoc_study.expected_results(norm), "   expected results (" + norm + ")", ss);
        DUNE_THROW(errors_are_not_as_expected, ss.str());
      }
  } // ... eoc_study()
}; // linearelliptic_SWIPDG_discretization


TYPED_TEST_CASE(linearelliptic_SWIPDG_discretization, AluConform2dTestCases);
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


#else // HAVE_ALUGRID


TEST(DISABLED_linearelliptic_SWIPDG_discretization, eoc_study_using_istl_and_1_subdomain);
TEST(DISABLED_linearelliptic_SWIPDG_discretization, eoc_study_using_istl_and_4_subdomain);
TEST(DISABLED_linearelliptic_SWIPDG_discretization, eoc_study_using_istl_and_16_subdomain);
TEST(DISABLED_linearelliptic_SWIPDG_discretization, eoc_study_using_istl_and_64_subdomain);


#endif // HAVE_ALUGRID
