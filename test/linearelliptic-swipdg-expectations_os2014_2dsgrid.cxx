// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include <dune/grid/sgrid.hh>

#include <dune/stuff/test/gtest/gtest.h>

#include <dune/hdd/linearelliptic/testcases/ESV2007.hh>

#include "linearelliptic-swipdg-expectations.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Tests {


template< bool anything >
class SWIPDGStudyExpectations< TestCases::ESV2007< SGrid< 2, 2 > >, 1, anything >
  : public internal::SWIPDGStudyExpectationsBase< TestCases::ESV2007< SGrid< 2, 2 > >, 1 >
{
  typedef TestCases::ESV2007< SGrid< 2, 2 > > TestCaseType;

public:
  static std::vector< double > results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "L2")
      return {1.13e-02, 2.90e-03, 7.41e-04, 1.88e-04};
    else if (type == "H1_semi")
      return {2.77e-01, 1.39e-01, 6.98e-02, 3.50e-02};
    else if (type == "energy")
      return {2.77e-01, 1.39e-01, 6.98e-02, 3.50e-02};
    else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // SWIPDGStudyExpectations


template class SWIPDGStudyExpectations< TestCases::ESV2007< SGrid< 2, 2 > >, 1 >;


} // namespace Tests
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune
