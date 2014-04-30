// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_PROBLEMS_ESV2007_HH
#define DUNE_HDD_LINEARELLIPTIC_PROBLEMS_ESV2007_HH

#include <memory>

#include <dune/common/static_assert.hh>

#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/expression.hh>
#include <dune/stuff/common/memory.hh>

#include <dune/pymor/functions/default.hh>

#include "default.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Problems {


template< class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim = 1 >
class ESV2007
{
  static_assert(AlwaysFalse< EntityImp >::value, "Not available for dimRange > 1!");
};


template< class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp >
class ESV2007< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 >
  : public Default< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 >
{
  typedef Default< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 > BaseType;
  typedef ESV2007< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 > ThisType;

  typedef Stuff::Functions::Constant< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 > ConstantFunctionType;
  typedef Stuff::Functions::Expression< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 > ExpressionFunctionType;
  typedef Pymor::Function::NonparametricDefault< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 > FunctionType;

public:
  static std::string static_id()
  {
    return BaseType::BaseType::static_id() + ".ESV2007";
  }

  static Stuff::Common::ConfigTree default_config(const std::string sub_name = "")
  {
    Stuff::Common::ConfigTree config("integration_order", "3");
    if (sub_name.empty())
      return config;
    else {
      Stuff::Common::ConfigTree tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  static std::unique_ptr< ThisType > create(const Stuff::Common::ConfigTree config = default_config(),
                                            const std::string sub_name = static_id())
  {
    const Stuff::Common::ConfigTree cfg = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    return Stuff::Common::make_unique< ThisType >(cfg.get("integration_order",
                                                          default_config().get< size_t >("integration_order")));
  } // ... create(...)

  ESV2007(const size_t integration_order = default_config().get< size_t >("integration_order"))
    : BaseType(std::make_shared< FunctionType >(new ConstantFunctionType(1)),
               std::make_shared< FunctionType >(new ExpressionFunctionType("x", "64.0 * pi * pi * (cos(8.0 * pi * x[0]) + cos(8.0 * pi * x[1]))", integration_order)),
               std::make_shared< FunctionType >(new ExpressionFunctionType("x", "cos(8.0 * pi * x[0]) + cos(8.0 * pi * x[1])", integration_order)),
               std::make_shared< FunctionType >(new ConstantFunctionType(0)))
  {}

}; // class ESV2007< ..., 1 >


} // namespace Problems
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_PROBLEMS_ESV2007_HH
