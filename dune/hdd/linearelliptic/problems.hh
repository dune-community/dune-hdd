// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_PROBLEMS_HH
#define DUNE_HDD_LINEARELLIPTIC_PROBLEMS_HH

#include <string>
#include <vector>
#include <memory>

#include <dune/stuff/common/configtree.hh>
#include <dune/stuff/common/exceptions.hh>

#include "problems/interfaces.hh"
#include "../playground/linearelliptic/problems/default.hh"
#include "../playground/linearelliptic/problems/ESV2007.hh"
#include "../playground/linearelliptic/problems/thermalblock.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {


template< class E, class D, int d, class R, int r = 1 >
class ProblemProvider
{
public:
  typedef ProblemInterface< E, D, d, R, r > InterfaceType;

  static std::vector< std::string > available()
  {
    return {
        Problems::Default< E, D, d, R, r >::static_id()
      , Problems::ESV2007< E, D, d, R, r >::static_id()
      , Problems::Thermalblock< E, D, d, R, r >::static_id()
    };
  }

  static Stuff::Common::ConfigTree default_config(const std::string type = available()[0],
                                                  const std::string sub_name = "")
  {
  if (type == Problems::Default< E, D, d, R, r >::static_id())
    return Problems::Default< E, D, d, R, r >::default_config(sub_name);
  else if (type == Problems::ESV2007< E, D, d, R, r >::static_id())
    return Problems::ESV2007< E, D, d, R, r >::default_config(sub_name);
  else if (type == Problems::Thermalblock< E, D, d, R, r >::static_id())
    return Problems::Thermalblock< E, D, d, R, r >::default_config(sub_name);
  else
    DUNE_THROW_COLORFULLY(Stuff::Exceptions::wrong_input_given,
                          "'" << type << "' is not a valid " << InterfaceType::static_id() << "!");
  }

  static std::unique_ptr< InterfaceType > create(const std::string type = available()[0],
                                                 const Stuff::Common::ConfigTree config = default_config(available()[0]))
  {
    if (type == Problems::Default< E, D, d, R, r >::static_id())
      return Problems::Default< E, D, d, R, r >::create(config);
    else if (type == Problems::ESV2007< E, D, d, R, r >::static_id())
      return Problems::ESV2007< E, D, d, R, r >::create(config);
    else if (type == Problems::Thermalblock< E, D, d, R, r >::static_id())
      return Problems::Thermalblock< E, D, d, R, r >::create(config);
    else
      DUNE_THROW_COLORFULLY(Stuff::Exceptions::wrong_input_given,
                            "'" << type << "' is not a valid " << InterfaceType::static_id() << "!");
  }
}; // clas ProblemProvider


} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_PROBLEMS_HH
