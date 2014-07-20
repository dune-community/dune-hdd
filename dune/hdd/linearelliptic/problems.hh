// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_PROBLEMS_HH
#define DUNE_HDD_LINEARELLIPTIC_PROBLEMS_HH

#include <string>
#include <vector>
#include <memory>

#if HAVE_ALUGRID
# include <dune/grid/alugrid.hh>
#endif

#include <dune/stuff/common/configtree.hh>
#include <dune/stuff/common/exceptions.hh>

#include "problems/interfaces.hh"
#include "problems/default.hh"
#include "problems/ESV2007.hh"
#include "../playground/linearelliptic/problems/mixed-boundaries.hh"
#include "../playground/linearelliptic/problems/thermalblock.hh"
#include "../playground/linearelliptic/problems/OS2014.hh"

namespace Dune {
namespace HDD {
namespace internal {


void lib_exists();


} // namespace internal
namespace LinearElliptic {
namespace internal {


template< class E, class D, int d, class R, int r = 1 >
class ProblemProviderBase
{
public:
  typedef ProblemInterface< E, D, d, R, r > InterfaceType;

protected:
  template< class ProblemType >
  static std::unique_ptr< ProblemType > call_create(const Stuff::Common::ConfigTree& config)
  {
    if (config.empty())
      return ProblemType::create();
    else
      return ProblemType::create(config);
  } // ... call_create(...)

public:
  static std::vector< std::string > available()
  {
    return {
        Problems::Default< E, D, d, R, r >::static_id()
      , Problems::MixedBoundaries< E, D, d, R, r >::static_id()
      , Problems::Thermalblock< E, D, d, R, r >::static_id()
    };
  } // ... available(...)

  static Stuff::Common::ConfigTree default_config(const std::string type = available()[0],
                                                  const std::string sub_name = "")
  {
  if (type == Problems::Default< E, D, d, R, r >::static_id())
    return Problems::Default< E, D, d, R, r >::default_config(sub_name);
  else if (type == Problems::MixedBoundaries< E, D, d, R, r >::static_id())
    return Problems::MixedBoundaries< E, D, d, R, r >::default_config(sub_name);
  else if (type == Problems::Thermalblock< E, D, d, R, r >::static_id())
    return Problems::Thermalblock< E, D, d, R, r >::default_config(sub_name);
  else
    DUNE_THROW_COLORFULLY(Stuff::Exceptions::wrong_input_given,
                          "'" << type << "' is not a valid " << InterfaceType::static_id() << "!");
  } // ... default_config(...)

  static std::unique_ptr< InterfaceType > create(const std::string type = available()[0],
                                                 const Stuff::Common::ConfigTree config = default_config(available()[0]))
  {
    if (type == Problems::Default< E, D, d, R, r >::static_id())
      return call_create< Problems::Default< E, D, d, R, r > >(config);
    else if (type == Problems::MixedBoundaries< E, D, d, R, r >::static_id())
      return call_create< Problems::MixedBoundaries< E, D, d, R, r > >(config);
    else if (type == Problems::Thermalblock< E, D, d, R, r >::static_id())
      return call_create< Problems::Thermalblock< E, D, d, R, r > >(config);
    else
      DUNE_THROW_COLORFULLY(Stuff::Exceptions::wrong_input_given,
                            "'" << type << "' is not a valid " << InterfaceType::static_id() << "!");
  } // ... create(...)
}; // clas ProblemProviderBase


} // namespace internal


template< class E, class D, int d, class R, int r = 1 >
class ProblemProvider
  : public internal::ProblemProviderBase< E, D, d, R, r >
{};


template< class E, class D, class R >
class ProblemProvider< E, D, 2, R, 1 >
  : public internal::ProblemProviderBase< E, D, 2, R, 1 >
{
  static const unsigned int d = 2;
  static const unsigned int r = 1;
  typedef internal::ProblemProviderBase< E, D, d, R, r > BaseType;
public:
  using typename BaseType::InterfaceType;

  static std::vector< std::string > available()
  {
    auto base = BaseType::available();
    base.push_back(Problems::ESV2007< E, D, d, R, r >::static_id());
    base.push_back(Problems::OS2014< E, D, d, R, r >::static_id());
    return base;
  } // ... available(...)

  static Stuff::Common::ConfigTree default_config(const std::string type = available()[0],
                                                  const std::string sub_name = "")
  {
  if (type == Problems::ESV2007< E, D, d, R, r >::static_id())
    return Problems::ESV2007< E, D, d, R, r >::default_config(sub_name);
  else if (type == Problems::OS2014< E, D, d, R, r >::static_id())
    return Problems::OS2014< E, D, d, R, r >::default_config(sub_name);
  else
    return BaseType::default_config(type, sub_name);
  } // ... default_config(...)

  static std::unique_ptr< InterfaceType > create(const std::string type = available()[0],
                                                 const Stuff::Common::ConfigTree config = Stuff::Common::ConfigTree())
  {
    if (type == Problems::ESV2007< E, D, d, R, r >::static_id())
      return BaseType::template call_create< Problems::ESV2007< E, D, d, R, r > >(config);
    else if (type == Problems::OS2014< E, D, d, R, r >::static_id())
      return BaseType::template call_create< Problems::OS2014< E, D, d, R, r > >(config);
    else
      return BaseType::create(type, config);
  } // ... create(...)
}; // clas ProblemProvider< ..., 2, ... 1 >


#if HAVE_ALUGRID


extern template class ProblemProvider< typename ALUConformGrid< 2, 2 >::template Codim< 0 >::Entity,
                                       double, 2, double, 1 >;


#endif // HAVE_ALUGRID

} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_PROBLEMS_HH
