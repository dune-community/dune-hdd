// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#define DUNE_STUFF_FUNCTIONS_DISABLE_CHECKS 1

#ifndef DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING
# define DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING 1
#endif

// This one has to come first (includes the config.h)!
#include <dune/stuff/test/main.hxx>

#if HAVE_ALUGRID
# include <dune/stuff/common/disable_warnings.hh>
#   include <dune/grid/alugrid.hh>
# include <dune/stuff/common/reenable_warnings.hh>

# include <dune/stuff/common/exceptions.hh>

# include <dune/stuff/common/print.hh>
# include <dune/stuff/common/float_cmp.hh>

# include <dune/hdd/playground/linearelliptic/testcases/ESV2007.hh>

# include "linearelliptic-swipdg.hh"

using namespace Dune;
using namespace HDD;


typedef ALUGrid< 2, 2, simplex, conforming > AluConform2dGridType;

typedef testing::Types< LinearElliptic::TestCases::ESV2007< AluConform2dGridType >
                      > AluConform2dTestCases;


template< class TestCaseType >
struct linearelliptic_SWIPDG_discretization
  : public ::testing::Test
{
  template< GDT::ChooseSpaceBackend space_backend, Stuff::LA::ChooseBackend la_backend >
  static void eoc_study()
  {
    const TestCaseType test_case;
    test_case.print_header(DSC_LOG_INFO);
    DSC_LOG_INFO << std::endl;
    LinearElliptic::Tests::SWIPDGStudy< TestCaseType, 1, space_backend, la_backend > eoc_study(test_case);
    Stuff::Test::check_eoc_study_for_success(eoc_study, eoc_study.run_eoc(DSC_LOG_INFO));
  } // ... eoc_study()

# if HAVE_DUNE_FEM
#   if HAVE_DUNE_ISTL
  static void eoc_study_using_fem_and_istl()
  {
    eoc_study< GDT::ChooseSpaceBackend::fem, Stuff::LA::ChooseBackend::istl_sparse >();
  }
#   endif // HAVE_DUNE_ISTL

#   if HAVE_EIGEN
  static void eoc_study_using_fem_and_eigen_sparse()
  {
    eoc_study< GDT::ChooseSpaceBackend::fem, Stuff::LA::ChooseBackend::eigen_sparse >();
  }
#   endif // HAVE_EIGEN
# endif // HAVE_DUNE_FEM

//  static void eoc_study_using_pdelab_and_istl()
//  {
//    eoc_study< GDT::ChooseSpaceBackend::pdelab, Stuff::LA::ChooseBackend::istl_sparse >();
//  }

//  static void eoc_study_using_pdelab_and_eigen_sparse()
//  {
//    eoc_study< GDT::ChooseSpaceBackend::pdelab, Stuff::LA::ChooseBackend::eigen_sparse >();
//  }
}; // linearelliptic_SWIPDG_discretization


# if HAVE_DUNE_FEM && (HAVE_DUNE_ISTL || HAVE_EIGEN)

TYPED_TEST_CASE(linearelliptic_SWIPDG_discretization, AluConform2dTestCases);

#   if HAVE_DUNE_ISTL
TYPED_TEST(linearelliptic_SWIPDG_discretization, eoc_study_using_fem_and_istl) {
  this->eoc_study_using_fem_and_istl();
}
#   else // HAVE_DUNE_ISTL
TYPED_TEST(DISABLED_linearelliptic_SWIPDG_discretization, eoc_study_using_fem_and_istl) {}
#   endif

#   if HAVE_EIGEN
TYPED_TEST(linearelliptic_SWIPDG_discretization, eoc_study_using_fem_and_eigen_sparse) {
  this->eoc_study_using_fem_and_eigen_sparse();
}
#   else // HAVE_EIGEN
TYPED_TEST(DISABLED_linearelliptic_SWIPDG_discretization, eoc_study_using_fem_and_eigen_sparse) {}
#   endif // HAVE_EIGEN

//TYPED_TEST(linearelliptic_SWIPDG_discretization, eoc_study_using_pdelab_and_istl) {
//  this->eoc_study_using_pdelab_and_istl();
//}

//TYPED_TEST(linearelliptic_SWIPDG_discretization, eoc_study_using_pdelab_and_eigen_sparse) {
//  this->eoc_study_using_pdelab_and_eigen_sparse();
//}

# endif // HAVE_DUNE_FEM && (HAVE_DUNE_ISTL || HAVE_EIGEN)
#else // HAVE_ALUGRID


TYPED_TEST(DISABLED_linearelliptic_SWIPDG_discretization, eoc_study_using_fem_and_istl) {}
TYPED_TEST(DISABLED_linearelliptic_SWIPDG_discretization, eoc_study_using_fem_and_eigen_sparse) {}


#endif // HAVE_ALUGRID
