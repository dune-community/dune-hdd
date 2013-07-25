// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_PROBLEMS_HH
#define DUNE_HDD_LINEARELLIPTIC_PROBLEMS_HH

#include <string>
#include <vector>

#include <dune/stuff/common/parameter/tree.hh>

#include <dune/pymor/common/exceptions.hh>

#include "problems/interfaces.hh"
#include "problems/default.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {


template< class D, int d, class R, int r, bool s = true >
class Problems
{
public:
  static std::vector< std::string > available()
  {
    return {
        Problem::Default< D, d, R, r, s >::static_id()
    };
  }

  static Dune::ParameterTree defaultSettings(const std::string type = available()[0],
                                             const std::string subname = "")
  {
  if (type == Problem::Default< D, d, R, r, s >::static_id())
    return Problem::Default< D, d, R, r, s >::defaultSettings(subname);
  else
    DUNE_PYMOR_THROW(Pymor::Exception::wrong_input,
                     "unknown Problem '" << type << "' requested!");
  }

  static ProblemInterface< D, d, R, r, s >* create(const std::string type = available()[0],
                                                   const Dune::ParameterTree settings = defaultSettings())
  {
    if (type == Problem::Default< D, d, R, r, s >::static_id())
      return Problem::Default< D, d, R, r, s >::create(settings);
    else
      DUNE_PYMOR_THROW(Pymor::Exception::wrong_input,
                       "unknown Problem '" << type << "' requested!");
  }
};


} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_PROBLEMS_HH
