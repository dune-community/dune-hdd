// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_PROBLEMS_THERMALBLOCK_HH
#define DUNE_HDD_LINEARELLIPTIC_PROBLEMS_THERMALBLOCK_HH

#include <memory>

#include <dune/common/static_assert.hh>

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/expression.hh>

#include <dune/pymor/functions/default.hh>
#include <dune/pymor/functions/checkerboard.hh>

#include "default.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Problems {


template< class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim = 1 >
class Thermalblock
{
  static_assert(AlwaysFalse< EntityImp >::value, "Not available for dimRange > 1!");
};


template< class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp >
class Thermalblock< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 >
  : public Default< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 >
{
  typedef Default< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 > BaseType;
  typedef Thermalblock< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 > ThisType;

public:
  typedef Pymor::Function::Checkerboard< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 >
      CheckerboardFunctionType;
  using typename BaseType::FunctionType;

  static std::string static_id()
  {
    return BaseType::BaseType::static_id() + ".thermalblock";
  }

  static Stuff::Common::ConfigTree default_config(const std::string sub_name = "")
  {
    typedef Stuff::Functions::Constant< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 > ConstantFunctionType;
    Stuff::Common::ConfigTree config;
    Stuff::Common::ConfigTree checkerboard_config = CheckerboardFunctionType::defaultSettings();
    checkerboard_config["name"] = "checkerboard_diffusion";
    checkerboard_config["parameterName"] = "diffusion_factor";
    config.add(checkerboard_config, "diffusion_factor");
    Stuff::Common::ConfigTree constant_config = ConstantFunctionType::defaultSettings();
    constant_config["type"] = ConstantFunctionType::static_id();
    constant_config["name"] = "force";
    constant_config["value"] = "1.0";
    config.add(constant_config, "force");
    constant_config["name"] = "dirichlet";
    constant_config["value"] = "0.0";
    config.add(constant_config, "dirichlet");
    constant_config["name"] = "neumann";
    config.add(constant_config, "neumann");
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
    std::shared_ptr< CheckerboardFunctionType >
        checkerboard_function(CheckerboardFunctionType::create(cfg.sub("diffusion_factor")));
    return Stuff::Common::make_unique< ThisType >(checkerboard_function,
                                                  BaseType::template create_function< 1 >("force", cfg),
                                                  BaseType::template create_function< 1 >("dirichlet", cfg),
                                                  BaseType::template create_function< 1 >("neumann", cfg));
  } // ... create(...)

  Thermalblock(const std::shared_ptr< const CheckerboardFunctionType >& checkerboard_function,
               const std::shared_ptr< const FunctionType >& force,
               const std::shared_ptr< const FunctionType >& diffusion,
               const std::shared_ptr< const FunctionType >& neumann)
    : BaseType(checkerboard_function, force, diffusion, neumann)
  {}

}; // class Thermalblock< ..., 1 >


} // namespace Problems
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_PROBLEMS_THERMALBLOCK_HH
