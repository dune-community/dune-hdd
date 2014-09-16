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

#include <dune/hdd/linearelliptic/testcases/ESV2007.hh>
#include <dune/hdd/linearelliptic/testcases/OS2014.hh>
#include <dune/hdd/linearelliptic/testcases/spe10.hh>

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

void nonparametric_convergence_study(const std::string visualization = "");

void nonparametric_convergence_study_alternative_summation();


typedef LinearElliptic::TestCases::ESV2007Multiscale< GridType > NonparametricBlockEocTestCaseType;
typedef LinearElliptic::Tests::BlockSWIPDGStudy
    < NonparametricBlockEocTestCaseType, 1, la_backend >         NonparametricBlockEocStudyType;

void nonparametric_block_convergence_study(const std::string& partitioning);


extern template class LinearElliptic::Tests::SWIPDGStudyExpectations
    < LinearElliptic::TestCases::ESV2007< ALUGrid< 2, 2, simplex, conforming > >, 1 >;

extern template class LinearElliptic::Tests::BlockSWIPDGStudyExpectations
    < LinearElliptic::TestCases::ESV2007Multiscale< ALUGrid< 2, 2, simplex, conforming > >, 1 >;

extern template class LinearElliptic::Tests::BlockSWIPDGStudyExpectations
    < LinearElliptic::TestCases::OS2014::ParametricBlockConvergence< ALUGrid< 2, 2, simplex, conforming > >, 1 >;


#endif // HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_GRID_MULTISCALE

#endif // DUNE_HDD_TEST_OS2014_HH
