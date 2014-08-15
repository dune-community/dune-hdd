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

#include "../../../linearelliptic/problems/default.hh"

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

  typedef Stuff::Functions::Constant< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 >
      ConstantScalarFunctionType;
  typedef Stuff::Functions::Constant< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim >
      ConstantVectorFunctionType;
  typedef Stuff::Functions::Constant< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, domainDim, domainDim >
      ConstantMatrixFunctionType;
  typedef Stuff::Functions::Expression< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim >
      ExpressionFunctionType;

public:
  using typename BaseType::DiffusionFactorType;
  using typename BaseType::DiffusionTensorType;
  using typename BaseType::FunctionType;

  static std::string static_id()
  {
    return BaseType::BaseType::static_id() + ".mixedboundaries";
  }

  static Stuff::Common::Configuration default_config(const std::string sub_name = "")
  {
    Stuff::Common::Configuration config;
    config.add(ConstantScalarFunctionType::default_config(), "diffusion_factor");
    config["diffusion_factor.type"] = ConstantScalarFunctionType::static_id();
    config["diffusion_factor.name"] = "diffusion_factor";
    config.add(ConstantMatrixFunctionType::default_config(), "diffusion_tensor");
    config["diffusion_tensor.type"] = ConstantMatrixFunctionType::static_id();
    config["diffusion_tensor.name"] = "diffusion_tensor";
    for (const std::string& type : {"force", "neumann"}) {
      config.add(ConstantVectorFunctionType::default_config(), type);
      config[type + ".type"] = ConstantVectorFunctionType::static_id();
      config[type + ".name"] = type;
    }
    config["force.value"] = "1";
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
      Stuff::Common::Configuration tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  static std::unique_ptr< ThisType > create(const Stuff::Common::Configuration config = default_config(),
                                            const std::string sub_name = static_id())
  {
    Stuff::Common::Configuration cfg = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Stuff::Common::Configuration default_cfg = default_config();
    if (!cfg.has_sub("diffusion_factor"))
      cfg.add(default_cfg.sub("diffusion_factor"), "diffusion_factor");
    if (!cfg.has_sub("diffusion_tensor"))
      cfg.add(default_cfg.sub("diffusion_tensor"), "diffusion_tensor");
    if (!cfg.has_sub("force"))
      cfg.add(default_cfg.sub("force"), "force");
    if (!cfg.has_sub("dirichlet"))
      cfg.add(default_cfg.sub("dirichlet"), "dirichlet");
    if (!cfg.has_sub("neumann"))
      cfg.add(default_cfg.sub("neumann"), "neumann");
    return Stuff::Common::make_unique< ThisType >(BaseType::create_scalar_function("diffusion_factor", cfg),
                                                  BaseType::create_matrix_function("diffusion_tensor", cfg),
                                                  BaseType::create_vector_function("force", cfg),
                                                  BaseType::create_vector_function("dirichlet", cfg),
                                                  BaseType::create_vector_function("neumann", cfg));

  } // ... create(...)

  MixedBoundaries(const std::shared_ptr< const DiffusionFactorType > diffusion_factor
                      = BaseType::create_scalar_function("diffusion_factor", default_config()),
                  const std::shared_ptr< const DiffusionTensorType > diffusion_tensor
                      = BaseType::create_matrix_function("diffusion_tensor", default_config()),
                  const std::shared_ptr< const FunctionType > force
                      = BaseType::create_vector_function("force", default_config()),
                  const std::shared_ptr< const FunctionType > dirichlet
                      = BaseType::create_vector_function("dirichlet", default_config()),
                  const std::shared_ptr< const FunctionType > neumann
                      = BaseType::create_vector_function("neumann", default_config()))
    : BaseType(diffusion_factor, diffusion_tensor, force, dirichlet, neumann)
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
