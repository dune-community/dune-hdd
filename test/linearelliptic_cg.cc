// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING
# define DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING 1
#endif

// This one has to come first (includes the config.h)!
#include <dune/stuff/test/main.hxx>

#include <dune/stuff/common/disable_warnings.hh>
# if HAVE_ALUGRID
#   include <dune/grid/alugrid.hh>
# endif
# include <dune/grid/spgrid.hh>
#include <dune/stuff/common/reenable_warnings.hh>

# include <dune/stuff/common/exceptions.hh>

# include <dune/stuff/common/print.hh>
# include <dune/stuff/common/float_cmp.hh>

# include <dune/hdd/linearelliptic/testcases/ESV2007.hh>
# include <dune/hdd/linearelliptic/testcases/OS2014.hh>
# include <dune/hdd/linearelliptic/testcases/spe10.hh>

# include "linearelliptic_cg.hh"

using namespace Dune;
using namespace HDD;


template< class TestCaseType >
struct linearelliptic_CG_discretization
  : public ::testing::Test
{
  template< GDT::ChooseSpaceBackend space_backend, Stuff::LA::ChooseBackend la_backend >
  static void eoc_study()
  {
    TestCaseType test_case;
    test_case.print_header(DSC_LOG_INFO);
    DSC_LOG_INFO << std::endl;
    LinearElliptic::Tests::CGStudy< TestCaseType, 1, space_backend, la_backend > eoc_study(test_case);
    Stuff::Test::check_eoc_study_for_success(eoc_study, eoc_study.run_eoc(DSC_LOG_INFO));
  } // ... eoc_study()

#if HAVE_DUNE_FEM
# if HAVE_DUNE_ISTL
  static void eoc_study_using_fem_and_istl()
  {
    eoc_study< GDT::ChooseSpaceBackend::fem, Stuff::LA::ChooseBackend::istl_sparse >();
  }
# endif // HAVE_DUNE_ISTL

# if HAVE_EIGEN
  static void eoc_study_using_fem_and_eigen()
  {
    eoc_study< GDT::ChooseSpaceBackend::fem, Stuff::LA::ChooseBackend::eigen_sparse >();
  }
# endif // HAVE_EIGEN
#endif // HAVE_DUNE_FEM

#if HAVE_DUNE_PDELAB
# if HAVE_DUNE_ISTL
  static void eoc_study_using_pdelab_and_istl()
  {
    eoc_study< GDT::ChooseSpaceBackend::pdelab, Stuff::LA::ChooseBackend::istl_sparse >();
  }
# endif // HAVE_DUNE_ISTL

# if HAVE_EIGEN
  static void eoc_study_using_pdelab_and_eigen()
  {
    eoc_study< GDT::ChooseSpaceBackend::pdelab, Stuff::LA::ChooseBackend::eigen_sparse >();
  }
# endif // HAVE_EIGEN
#endif // HAVE_DUNE_PDELAB
}; // linearelliptic_CG_discretization

typedef SPGrid< double, 2 > SPG;
typedef testing::Types< LinearElliptic::TestCases::ESV2007< SPG >
//                      , LinearElliptic::TestCases::Spe10::Model1< SPG >
#if HAVE_ALUGRID
                      , LinearElliptic::TestCases::ESV2007< ALUGrid< 2, 2, simplex, conforming > >
//                      , LinearElliptic::TestCases::Spe10::Model1< ALUGrid< 2, 2, simplex, conforming > >
#endif
                      > TestCases;

TYPED_TEST_CASE(linearelliptic_CG_discretization, TestCases);

#if HAVE_DUNE_FEM && (HAVE_DUNE_ISTL || HAVE_EIGEN)


# if HAVE_DUNE_ISTL
TYPED_TEST(linearelliptic_CG_discretization, eoc_study_using_fem_and_istl) {
  this->eoc_study_using_fem_and_istl();
}
# else // HAVE_DUNE_ISTL
TEST(DISABLED_linearelliptic_CG_discretization, eoc_study_using_fem_and_istl) {}
# endif

# if HAVE_EIGEN
TYPED_TEST(linearelliptic_CG_discretization, eoc_study_using_fem_and_eigen) {
  this->eoc_study_using_fem_and_eigen();
}
# else // HAVE_EIGEN
TEST(DISABLED_linearelliptic_CG_discretization, eoc_study_using_fem_and_eigen) {}
# endif // HAVE_EIGEN


//TYPED_TEST(linearelliptic_CG_discretization, eoc_study_using_pdelab_and_eigen) {
//  this->eoc_study_using_pdelab_and_eigen();
//}


#else // HAVE_DUNE_FEM && (HAVE_DUNE_ISTL || HAVE_EIGEN)

TYPED_TEST(linearelliptic_CG_discretization, eoc_study_using_pdelab_and_istl) {
  this->eoc_study_using_pdelab_and_istl();
}

TEST(DISABLED_linearelliptic_CG_discretization, eoc_study_using_fem_and_istl) {}
TEST(DISABLED_linearelliptic_CG_discretization, eoc_study_using_fem_and_eigen) {}


#endif // HAVE_DUNE_FEM && (HAVE_DUNE_ISTL || HAVE_EIGEN)


namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Tests {

#if HAVE_ALUGRID


extern template class CGStudyExpectations< TestCases::ESV2007< ALUGrid< 2, 2, simplex, conforming > >, 1 >;

extern template class CGStudyExpectations< TestCases::Spe10::Model1< ALUGrid< 2, 2, simplex, conforming > >, 1 >;


#endif // HAVE_ALUGRID

} // namespace Tests
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

