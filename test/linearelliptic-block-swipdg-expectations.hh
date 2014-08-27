// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_TEST_LINEARELLIPTIC_BLOCK_SWIPDG_EXPECTATIONS_HH
#define DUNE_HDD_TEST_LINEARELLIPTIC_BLOCK_SWIPDG_EXPECTATIONS_HH

#if HAVE_ALUGRID
# include <dune/grid/alugrid.hh>
#endif

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/type_utils.hh>

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace TestCases {


// forwards
template< class GridType >
class ESV2007;


template< class GridType >
class OS2014;


#if HAVE_DUNE_GRID_MULTISCALE


template< class GridType >
class ESV2007Multiscale;


template< class GridType >
class OS2014Multiscale;


#endif // HAVE_DUNE_GRID_MULTISCALE


} // namespace TestCases
namespace Tests {


template< class TestCaseType, int polOrder, bool implemented = true >
class EocStudyBlockSWIPDGExpectations
{
public:
  static std::vector< double > results(const TestCaseType& /*test_case*/, const std::string type)
  {
    DUNE_THROW(Stuff::Exceptions::test_results_missing,
               "Please record the expected results for\n"
               << "TestCaseType: " << Stuff::Common::Typename< TestCaseType >::value() << "\n"
               << "polOrder: " << polOrder << "\n"
               << "type: " << type << "\n"
               << "Please put an appropriate specialiaztion of EocStudyBlockSWIPDGExpectations for this TestCaseType "
               << "in a separate object file (see examples below) or add\n"
               << "  'template class EocStudyBlockSWIPDGExpectations< TestCaseType, " << polOrder << ", true >'\n"
               << "for this polOrder in the appropriate object file!\n\n"
               << "Oh: and do not forget to add\n"
               << "  'template class EocStudyBlockSWIPDGExpectations...'\n"
               << "to each test source using these results!");
    return {};
  } // ... results(...)
}; // EocStudyBlockSWIPDGExpectations


} // namespace Tests
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_TEST_LINEARELLIPTIC_BLOCK_SWIPDG_EXPECTATIONS_HH
