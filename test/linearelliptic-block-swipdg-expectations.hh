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
namespace internal {


template< class TestCaseType, int polOrder >
class BlockSWIPDGStudyExpectationsBase
{
public:
  static size_t rate(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "L2")
      return polOrder + 1;
    else if (type == "H1_semi")
      return polOrder;
    else if (type.substr(0, 6) == "energy")
      return polOrder;
    else if (type == "eta_NC_OS2014")
      return polOrder;
    else if (type == "eta_R_OS2014")
      return polOrder + 1;
    else if (type.substr(0, 13) == "eta_DF_OS2014")
      return polOrder;
    else if (type == "eta_OS2014")
      return polOrder;
    else if (type == "eta_OS2014_*")
      return polOrder;
    else if (type.substr(0, 10) == "eff_OS2014")
      return 0;
    else if (type.substr(0, 12) == "eff_OS2014_*")
      return 0;
    else
      DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Wrong type '" << type << "' requested!");
  } // ... rate(...)
}; // class BlockSWIPDGStudyExpectationsBase


} // namespace internal


template< class TestCaseType, int polOrder, bool implemented = true >
class BlockSWIPDGStudyExpectations
  : public internal::BlockSWIPDGStudyExpectationsBase< TestCaseType, polOrder >
{
public:
  static std::vector< double > results(const TestCaseType& /*test_case*/, const std::string type)
  {
    DUNE_THROW(Stuff::Exceptions::test_results_missing,
               "Please record the expected results for\n"
               << "TestCaseType: " << Stuff::Common::Typename< TestCaseType >::value() << "\n"
               << "polOrder: " << polOrder << "\n"
               << "type: " << type << "\n"
               << "Please put an appropriate specialiaztion of BlockSWIPDGStudyExpectations for this TestCaseType "
               << "in a separate object file (see examples below) or add\n"
               << "  'template class BlockSWIPDGStudyExpectations< TestCaseType, " << polOrder << ", true >'\n"
               << "for this polOrder in the appropriate object file!\n\n"
               << "Oh: and do not forget to add\n"
               << "  'extern template class BlockSWIPDGStudyExpectations...'\n"
               << "to each test source using these results!");
    return {};
  } // ... results(...)
}; // BlockSWIPDGStudyExpectations


} // namespace Tests
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_TEST_LINEARELLIPTIC_BLOCK_SWIPDG_EXPECTATIONS_HH
