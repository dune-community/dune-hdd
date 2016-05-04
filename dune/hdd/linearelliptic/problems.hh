// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_PROBLEMS_HH
#define DUNE_HDD_LINEARELLIPTIC_PROBLEMS_HH

#include <string>
#include <vector>
#include <memory>

#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/exceptions.hh>

#include "problems/interfaces.hh"
#include "problems/default.hh"
#include "problems/ORS2016.hh"
#include "problems/ESV2007.hh"
#include "problems/mixed-boundaries.hh"
#include "problems/OS2014.hh"
#include "problems/spe10.hh"
#include "problems/thermalblock.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {


/**
 * \note If you want to add a new function FooBar, do the following: provide a definition that is available for all
 *       template arguments, like:
\code
template< class E, class D, int d, class R, int r = 1 >
class FooBar
  : public ProblemInterface< E, D, d, R, r >
{
  FooBar() { static_assert(AlwaysFalse< E >::value, "Not available for these dimensions!"); }
};
\endcode
 *       Every specialization that can be provided by the provider then has to define:
\code
static const bool available = true;
\endcode
 *       This is all you have to do when implementing the function. In addition you have to add the appropriate include
 *       in this file (of course) and the appropriate type below (just like the rest, should be obvious).
 */
template< class E, class D, int d, class R, int r >
class ProblemsProvider
{
public:
  typedef ProblemInterface< E, D, d, R, r > InterfaceType;

private:
  template< class P, bool available = false >
  struct Call
  {
    static std::vector< std::string > append(std::vector< std::string > in)
    {
      return in;
    }

    static bool compare(const std::string& /*type*/)
    {
      return false;
    }

    static Stuff::Common::Configuration default_config(const std::string /*sub_name*/)
    {
      DUNE_THROW(Stuff::Exceptions::internal_error, "This should not happen!");
      return Stuff::Common::Configuration(0);
    }

    static std::unique_ptr< P > create(const Stuff::Common::Configuration& /*cfg*/)
    {
      DUNE_THROW(Stuff::Exceptions::internal_error, "This should not happen!");
      return std::unique_ptr< P >(nullptr);
    }
  }; // struct Call

  template< class P >
  struct Call< P, true >
  {
    static std::vector< std::string > append(std::vector< std::string > in)
    {
      in.push_back(P::static_id());
      return in;
    }

    static bool compare(const std::string& type)
    {
      return type == P::static_id();
    }

    static Stuff::Common::Configuration default_config(const std::string sub_name)
    {
      return P::default_config(sub_name);
    }

    static std::unique_ptr< P > create(const Stuff::Common::Configuration& cfg)
    {
      if (cfg.empty())
        return P::create();
      else
        return P::create(cfg);
    }
  }; // struct Call< ..., true >

  template< class P >
  static std::vector< std::string > call_append(std::vector< std::string > in)
  {
    return Call< P, P::available >::append(in);
  }

  template< class P >
  static bool call_compare(const std::string& type)
  {
    return Call< P, P::available >::compare(type);
  }

  template< class P >
  static Stuff::Common::Configuration call_default_config(const std::string sub_name)
  {
    return Call< P, P::available >::default_config(sub_name);
  }

  template< class P >
  static std::unique_ptr< P > call_create(const Stuff::Common::Configuration& cfg)
  {
    return Call< P, P::available >::create(cfg);
  }

  static std::string available_as_str()
  {
    std::string ret = "";
    const auto vals = available();
    if (vals.size() > 0) {
      ret += vals[0];
      for (size_t ii = 1; ii < vals.size(); ++ii)
        ret += "\n   " + vals[ii];
    }
    return ret;
  } // ... available_as_str(...)

  typedef Problems::Default< E, D, d, R, r >                   DefaultType;
  typedef Problems::ESV2007< E, D, d, R, r >                   ESV2007Type;
  typedef Problems::ORS2016< E, D, d, R, r >                   ORS2016Type;
  typedef Problems::MixedBoundaries< E, D, d, R, r >           MixedBoundariesType;
  typedef Problems::OS2014::ParametricESV2007< E, D, d, R, r > OS2014ParametricESV2007Type;
  typedef Problems::Spe10::Model1< E, D, d, R, r >             Spe10Model1Type;
  typedef Problems::Thermalblock< E, D, d, R, r >              ThermalblockType;

public:
  static std::vector< std::string > available()
  {
    std::vector< std::string > ret;
    ret = call_append< DefaultType >(ret);
    ret = call_append< ESV2007Type >(ret);
    ret = call_append< ORS2016Type >(ret);
    ret = call_append< MixedBoundariesType >(ret);
    ret = call_append< OS2014ParametricESV2007Type >(ret);
    ret = call_append< Spe10Model1Type >(ret);
    ret = call_append< ThermalblockType >(ret);
    return ret;
  } // ... available(...)

  static Stuff::Common::Configuration default_config(const std::string type, const std::string sub_name = "")
  {
    if (call_compare< DefaultType >(type))
      return call_default_config< DefaultType >(sub_name);
    else if (call_compare< ESV2007Type >(type))
      return call_default_config< ESV2007Type >(sub_name);
    else if (call_compare< ORS2016Type >(type))
      return call_default_config< ORS2016Type >(sub_name);
    else if (call_compare< MixedBoundariesType >(type))
      return call_default_config< MixedBoundariesType >(sub_name);
    else if (call_compare< OS2014ParametricESV2007Type >(type))
      return call_default_config< OS2014ParametricESV2007Type >(sub_name);
    else if (call_compare< Spe10Model1Type >(type))
      return call_default_config< Spe10Model1Type >(sub_name);
    else if (call_compare< ThermalblockType >(type))
      return call_default_config< ThermalblockType >(sub_name);
    else if (available().empty())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given,
                 "There is no " << InterfaceType::static_id() << " available for dimensions " << int(d) << " -> "
                 << int(r) << "!");
    else
      DUNE_THROW(Stuff::Exceptions::wrong_input_given,
                 "Requested type '" << type << "' is not one of those avaible for dimensions " << int(d) << " -> "
                 << int(r) << ":\n" << available_as_str());  } // ... default_config(...)

  static std::unique_ptr< InterfaceType > create(const std::string type = available()[0],
                                                 const Stuff::Common::Configuration cfg = Stuff::Common::Configuration())
  {
    if (call_compare< DefaultType >(type))
      return call_create< DefaultType >(cfg);
    else if (call_compare< ESV2007Type >(type))
      return call_create< ESV2007Type >(cfg);
    else if (call_compare< ORS2016Type >(type))
      return call_create< ORS2016Type >(cfg);
    else if (call_compare< MixedBoundariesType >(type))
      return call_create< MixedBoundariesType >(cfg);
    else if (call_compare< OS2014ParametricESV2007Type >(type))
      return call_create< OS2014ParametricESV2007Type >(cfg);
    else if (call_compare< Spe10Model1Type >(type))
      return call_create< Spe10Model1Type >(cfg);
    else if (call_compare< ThermalblockType >(type))
      return call_create< ThermalblockType >(cfg);
    else if (available().empty())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given,
                 "There is no " << InterfaceType::static_id() << " available for dimensions " << int(d) << " -> "
                 << int(r) << "!");
    else
      DUNE_THROW(Stuff::Exceptions::wrong_input_given,
                 "Requested type '" << type << "' is not one of those avaible for dimensions " << int(d) << " -> "
                 << int(r) << ":\n" << available_as_str());
  } // ... create(...)
}; // clas ProblemsProvider


} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_PROBLEMS_HH
