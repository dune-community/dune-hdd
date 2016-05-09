// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_PROBLEMS_THERMALBLOCK_HH
#define DUNE_HDD_LINEARELLIPTIC_PROBLEMS_THERMALBLOCK_HH

#include <memory>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/static_assert.hh>
#include <dune/common/timer.hh>
#include <dune/common/typetraits.hh>

#if HAVE_DUNE_GRID_MULTISCALE
# include <dune/grid/multiscale/provider/cube.hh>
#endif

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/common/string.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/expression.hh>
#include <dune/stuff/playground/functions/indicator.hh>
#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/grid/boundaryinfo.hh>

#include <dune/pymor/functions/default.hh>
#include <dune/pymor/functions/checkerboard.hh>

#include "default.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Problems {


template< class E, class D, int d, class R, int r = 1 >
class Thermalblock
  : public ProblemInterface< E, D, d, R, r >
{
  Thermalblock() { static_assert(AlwaysFalse< E >::value, "Not available for these dimensions!"); }
};


template< class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp >
class Thermalblock< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 >
  : public Default< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 >
{
  typedef Default< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 > BaseType;
  typedef Thermalblock< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 > ThisType;

public:
  typedef Pymor::Functions::Checkerboard< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 >
      CheckerboardFunctionType;
  using typename BaseType::DiffusionTensorType;
  using typename BaseType::FunctionType;

  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::BaseType::static_id() + ".thermalblock";
  }

  static Stuff::Common::Configuration default_config(const std::string sub_name = "")
  {
    typedef Stuff::Functions::Constant< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 >
        ConstantFunctionType;
    typedef Stuff::Functions::Constant< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, domainDim, domainDim >
        ConstantMatrixFunctionType;
    Stuff::Common::Configuration config;
    config["type"] = static_id();
    Stuff::Common::Configuration checkerboard_config = CheckerboardFunctionType::default_config();
    checkerboard_config["name"] = "diffusion_factor";
    checkerboard_config["type"] = CheckerboardFunctionType::static_id();
    checkerboard_config["num_elements"] = "[4 4 4]";
    checkerboard_config["parameter_name"] = "diffusion";
    config.add(checkerboard_config, "diffusion_factor");
    Stuff::Common::Configuration diffusion_tensor_config = ConstantMatrixFunctionType::default_config();
    diffusion_tensor_config["name"] = "diffusion_tensor";
    diffusion_tensor_config["type"] = ConstantMatrixFunctionType::static_id();
    config.add(diffusion_tensor_config, "diffusion_tensor");
    Stuff::Common::Configuration constant_config = ConstantFunctionType::default_config();
    constant_config["type"] = ConstantFunctionType::static_id();
    constant_config["name"] = "force";
    constant_config["value"] = "1";
    config.add(constant_config, "force");
    constant_config["name"] = "dirichlet";
    constant_config["value"] = "0";
    config.add(constant_config, "dirichlet");
    constant_config["name"] = "neumann";
    config.add(constant_config, "neumann");
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
    const Stuff::Common::Configuration cfg = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    std::shared_ptr< CheckerboardFunctionType >
        checkerboard_function(CheckerboardFunctionType::create(cfg.sub("diffusion_factor")));
    return Stuff::Common::make_unique< ThisType >(checkerboard_function,
                                                  BaseType::create_matrix_function("diffusion_tensor", cfg),
                                                  BaseType::create_vector_function("force", cfg),
                                                  BaseType::create_vector_function("dirichlet", cfg),
                                                  BaseType::create_vector_function("neumann", cfg));
  } // ... create(...)

  Thermalblock(const std::shared_ptr< const CheckerboardFunctionType >& checkerboard_function,
               const std::shared_ptr< const DiffusionTensorType >& diffusion_tensor,
               const std::shared_ptr< const FunctionType >& force,
               const std::shared_ptr< const FunctionType >& dirichlet,
               const std::shared_ptr< const FunctionType >& neumann)
    : BaseType(checkerboard_function, diffusion_tensor, force, dirichlet, neumann)
  {}

  virtual std::string type() const override
  {
    return BaseType::BaseType::static_id() + ".thermalblock";
  }
}; // class Thermalblock< ..., 1 >


template< class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim = 1 >
class LocalThermalblock
{
  static_assert(AlwaysFalse< EntityImp >::value, "Not available for dimRange > 1!");
};


template< class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp >
class LocalThermalblock< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 >
  : public Default< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 >
{
  typedef Default< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 > BaseType;
  typedef LocalThermalblock< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 > ThisType;

  typedef Pymor::Functions::AffinelyDecomposableDefault< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 >
      AffinelyDecomposableDefaultFunctionType;
  typedef Stuff::Functions::Constant< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, domainDim, domainDim >
      ConstantMatrixFunctionType;
public:
  using typename BaseType::FunctionType;

  static std::string static_id()
  {
    return BaseType::BaseType::static_id() + ".localthermalblock";
  }

  static Stuff::Common::Configuration default_config(const std::string sub_name = "")
  {
    typedef Stuff::Functions::Constant< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 >
        ConstantFunctionType;
    Stuff::Common::Configuration config;
    Stuff::Common::Configuration constant_config = ConstantFunctionType::default_config();
    constant_config["type"] = ConstantFunctionType::static_id();
    constant_config["name"] = "force";
    constant_config["value"] = "1";
    config.add(constant_config, "force");
    constant_config["name"] = "dirichlet";
    constant_config["value"] = "0";
    config.add(constant_config, "dirichlet");
    constant_config["name"] = "neumann";
    config.add(constant_config, "neumann");
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
    const Stuff::Common::Configuration cfg = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    return Stuff::Common::make_unique< ThisType >(BaseType::create_vector_function("force", cfg),
                                                  BaseType::create_vector_function("dirichlet", cfg),
                                                  BaseType::create_vector_function("neumann", cfg));
  } // ... create(...)

  LocalThermalblock(const std::shared_ptr< const FunctionType >& force,
                    const std::shared_ptr< const FunctionType >& dirichlet,
                    const std::shared_ptr< const FunctionType >& neumann)
    : BaseType(create_diffusion_factor(),
               create_diffusion_tensor(),
               force,
               dirichlet,
               neumann)
  {}

  virtual std::string type() const override
  {
    return BaseType::BaseType::static_id() + ".localthermalblock";
  }

private:
  static std::shared_ptr< AffinelyDecomposableDefaultFunctionType > create_diffusion_factor()
  {
    typedef Stuff::Functions::DomainIndicator< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 > IndicatorFunctionType;
    const Pymor::ParameterType mu("diffusion", 3);

    auto ret = std::make_shared< AffinelyDecomposableDefaultFunctionType >("diffusion_factor");
    ret->register_component(new IndicatorFunctionType({{{{0.0, 0.0}, {0.5, 0.16}}, 1.0},
                                                       {{{0.0, 0.16}, {0.16, 0.33}}, 1.0},
                                                       {{{0.33, 0.16}, {0.5, 0.33}}, 1.0},
                                                       {{{0.0, 0.33}, {0.5, 1.0}}, 1.0}},
                                                      "left_block"),
                            new Pymor::ParameterFunctional(mu, "diffusion_factor[0]"));
    ret->register_component(new IndicatorFunctionType({{{{0.5, 0.0}, {1.0, 1.0}}, 1.0}},
                                                      "right_block"),
                            new Pymor::ParameterFunctional(mu, "diffusion_factor[1]"));
    ret->register_component(new IndicatorFunctionType({{{{0.16, 0.16}, {0.33, 0.33}}, 1.0}},
                                                      "small_block"),
                            new Pymor::ParameterFunctional(mu, "diffusion_factor[2]"));
    return ret;
  } // ... create_diffusion_factor()

  typedef Pymor::Functions::NonparametricDefault
      < EntityImp, DomainFieldImp, domainDim, RangeFieldImp, domainDim, domainDim > NonparametricFunctionType;

  static std::shared_ptr< NonparametricFunctionType > create_diffusion_tensor()
  {
    Stuff::Common::Configuration diffusion_tensor_config = ConstantMatrixFunctionType::default_config();
    diffusion_tensor_config["name"] = "diffusion_tensor";
    return std::make_shared< NonparametricFunctionType >(ConstantMatrixFunctionType::create(diffusion_tensor_config));
  } // ... create_diffusion_tensor()
}; // class LocalThermalblock< ..., 1 >


} // namespace Problems
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_PROBLEMS_THERMALBLOCK_HH
