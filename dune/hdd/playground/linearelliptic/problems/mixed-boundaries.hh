// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_PROBLEMS_MIXED_BOUNDARIES_HH
#define DUNE_HDD_LINEARELLIPTIC_PROBLEMS_MIXED_BOUNDARIES_HH

#include <memory>

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/common/string.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/expression.hh>

#include "default.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Problems {


template< class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim = 1 >
class MixedBoundaries
  : public Default< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim >
{
  typedef Default< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim >         BaseType;
  typedef MixedBoundaries< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim > ThisType;

  typedef Stuff::Functions::Constant< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim >
      ConstantFunctionType;
  typedef Stuff::Functions::Expression< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim >
      ExpressionFunctionType;

public:
  using typename BaseType::DiffusionFactorType;
  using typename BaseType::FunctionType;

  static std::string static_id()
  {
    return BaseType::BaseType::static_id() + ".mixedboundaries";
  }

  static Stuff::Common::ConfigTree default_config(const std::string sub_name = "")
  {
    Stuff::Common::ConfigTree config;
    for (const std::string& type : {"diffusion_factor", "force", "neumann"}) {
      config.add(ConstantFunctionType::default_config(), type);
      config[type + ".type"] = ConstantFunctionType::static_id();
      config[type + ".value"] = "1";
      config[type + ".name"] = type;
    }
    config["neumann.value"] = "0.1";
    config.add(ExpressionFunctionType::default_config(), "dirichlet");
    config["dirichlet.type"] = ExpressionFunctionType::static_id();
    config["dirichlet.name"] = "dirichlet";
    config["dirichlet.expression"] = "0.25";
    for (size_t dd = 0; dd < domainDim; ++dd)
      config["dirichlet.expression"] += "*x[" + Stuff::Common::toString(dd) + "]";
    config["dirichlet.order"] = Stuff::Common::toString(domainDim);
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
    Stuff::Common::ConfigTree cfg = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Stuff::Common::ConfigTree default_cfg = default_config();
    if (!cfg.has_sub("diffusion_factor"))
      cfg.add(default_cfg.sub("diffusion_factor"), "diffusion_factor");
    if (!cfg.has_sub("force"))
      cfg.add(default_cfg.sub("force"), "force");
    if (!cfg.has_sub("dirichlet"))
      cfg.add(default_cfg.sub("dirichlet"), "dirichlet");
    if (!cfg.has_sub("neumann"))
      cfg.add(default_cfg.sub("neumann"), "neumann");
    return Stuff::Common::make_unique< ThisType >(BaseType::template create_function< 1 >("diffusion_factor", cfg),
                                                  BaseType::template create_function< 1 >("force", cfg),
                                                  BaseType::template create_function< 1 >("dirichlet", cfg),
                                                  BaseType::template create_function< 1 >("neumann", cfg));

  } // ... create(...)

  MixedBoundaries(const std::shared_ptr< const DiffusionFactorType > diffusion_factor
                      = BaseType::template create_function< 1 >("diffusion_factor", default_config()),
                  const std::shared_ptr< const FunctionType > force
                      = BaseType::template create_function< 1 >("force", default_config()),
                  const std::shared_ptr< const FunctionType > dirichlet
                      = BaseType::template create_function< 1 >("dirichlet", default_config()),
                  const std::shared_ptr< const FunctionType > neumann
                      = BaseType::template create_function< 1 >("neumann", default_config()))
    : BaseType(diffusion_factor, force, dirichlet, neumann)
  {}

  virtual std::string type() const DS_OVERRIDE
  {
    return BaseType::BaseType::static_id() + ".mixedboundaries";
  }
}; // class MixedBoundaries


} // namespace Problems
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_PROBLEMS_MIXED_BOUNDARIES_HH
