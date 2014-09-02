// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_TEST_OS2014_HH
#define DUNE_HDD_TEST_OS2014_HH

#define DUNE_STUFF_FUNCTIONS_DISABLE_CHECKS 1
#define DUNE_HDD_LINEARELLIPTIC_TESTCASES_BASE_DISABLE_WARNING 1

#if HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_GRID_MULTISCALE

#include <string>
#include <vector>
#include <map>

#include <dune/stuff/common/disable_warnings.hh>
# include <dune/grid/alugrid.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/stuff/la/container.hh>

#include <dune/pymor/parameters/base.hh>

#include <dune/hdd/playground/linearelliptic/testcases/ESV2007.hh>
#include <dune/hdd/playground/linearelliptic/testcases/OS2014.hh>
#include <dune/hdd/playground/linearelliptic/testcases/spe10.hh>

#include "linearelliptic-swipdg.hh"
#include "linearelliptic-block-swipdg.hh"

using namespace Dune;
using namespace HDD;
using Parameter = Pymor::Parameter;


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


typedef LinearElliptic::TestCases::ESV2007< GridType >             NonparametricEocTestCaseType;
typedef LinearElliptic::Tests::SWIPDGStudy
    < NonparametricEocTestCaseType, 1, space_backend, la_backend > NonparametricEocStudyType;

void OS2014_nonparametric_convergence_study__SWIPDG_fine_triangulation();

void OS2014_nonparametric_convergence_study__SWIPDG_fine_triangulation_alternative_summation();


typedef LinearElliptic::TestCases::ESV2007Multiscale< GridType > NonparametricBlockEocTestCaseType;
typedef LinearElliptic::Tests::BlockSWIPDGStudy
    < NonparametricBlockEocTestCaseType, 1, la_backend >         NonparametricBlockEocStudyType;

void nonparametric_block_convergence_study(const std::string& partitioning);


typedef LinearElliptic::TestCases::OS2014Multiscale< GridType > ParametricBlockEocTestCaseType;
typedef LinearElliptic::Tests::BlockSWIPDGStudy
    < ParametricBlockEocTestCaseType, 1, la_backend >           ParametricBlockEocStudyType;

void print_parameter_information(const ParametricBlockEocTestCaseType& parametric_test_case);

void parametric_convergence_study(const std::string partitioning,
                                  const std::vector< std::string >& only_these_norms,
                                  const std::map< std::string, Pymor::Parameter >& parameters,
                                  const bool print_header = false);


typedef LinearElliptic::TestCases::Spe10Model1< GridType >                  NonparametricLocalizationTestCaseType;
typedef LinearElliptic::Tests::SWIPDGStudy
    < NonparametricLocalizationTestCaseType, 1, space_backend, la_backend > NonparametricLocalizationStudyType;

void nonparametric_localization_study__SWIPDG_fine_triangulation();


typedef LinearElliptic::TestCases::Spe10Model1Multiscale< GridType > NonparametricBlockLocalizationTestCaseType;
typedef LinearElliptic::Tests::BlockSWIPDGStudy
    < NonparametricBlockLocalizationTestCaseType, 1, la_backend >    NonparametricBlockLocalizationStudyType;

extern template class LinearElliptic::Tests::BlockSWIPDGStudyExpectations
    < LinearElliptic::TestCases::ESV2007Multiscale< ALUGrid< 2, 2, simplex, conforming > >, 1, true >;

extern template class LinearElliptic::Tests::BlockSWIPDGStudyExpectations
    < LinearElliptic::TestCases::OS2014Multiscale< ALUGrid< 2, 2, simplex, conforming > >, 1, true >;


#endif // HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_GRID_MULTISCALE

#endif // DUNE_HDD_TEST_OS2014_HH
