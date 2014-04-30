// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_PROBLEMS_DEFAULT_HH
#define DUNE_HDD_LINEARELLIPTIC_PROBLEMS_DEFAULT_HH

#include "config.h"

#include <memory>

#include <dune/stuff/common/string.hh>
#include <dune/stuff/common/configtree.hh>
#include <dune/stuff/functions/constant.hh>

#include <dune/pymor/functions.hh>

#include "../../../linearelliptic/problems/interfaces.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Problems {


template< class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim >
class Default
  : public ProblemInterface< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim >
{
  typedef ProblemInterface< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim > BaseType;
public:
  typedef Default< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim > ThisType;

  using typename BaseType::EntityType;
  using typename BaseType::DomainFieldType;
  using BaseType::dimDomain;
  using typename BaseType::RangeFieldType;
  using BaseType::dimRange;

  using typename BaseType::DiffusionFactorType;
//  using typename BaseType::DiffusionTensorType;
  using typename BaseType::FunctionType;

  static const std::string static_id()
  {
    return BaseType::static_id() + ".default";
  }

  static Stuff::Common::ConfigTree default_config(const std::string sub_name = "")
  {
    Stuff::Common::ConfigTree config;
    // diffusion factor
    typedef Pymor::Function::Checkerboard< EntityType, DomainFieldType, dimDomain, RangeFieldType, 1 >
        CheckerBoardFunctionType;
    config.add(CheckerBoardFunctionType::defaultSettings(), "diffusion_factor");
    config["diffusion_factor.parameterName"] = "diffusion_factor";
    config["diffusion_factor.name"] = "diffusion_factor";
    config["diffusion_factor.type"] = CheckerBoardFunctionType::static_id();
//    // diffusion tensor
//    typedef Stuff::Function::Constant< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimDomain, dimDomain >
//        ConstantFunctionType;
//    config.add(ConstantFunctionType::defaultSettings(), "diffusion_tensor");
//    config["diffusion_tensor.name"] = "diffusion_tensor";
//    config["diffusion_tensor.type"] = ConstantFunctionType::static_id();
    // force
    typedef Stuff::Functions::Expression< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimDomain, dimDomain >
        ExpressionFunctionType;
    config["force.variable"] = "x";
    config["force.expression"] = "1.0";
    config["force.order"] = "0";
    config["force.name"] = "force";
    config["force.type"] = ExpressionFunctionType::static_id();
    // dirichlet values
    config["dirichlet.variable"] = "x";
    config["dirichlet.expression"] = "0.1*x[0]";
    config["dirichlet.order"] = "1";
    config["dirichlet.name"] = "dirichlet";
    config["dirichlet.type"] = ExpressionFunctionType::static_id();
    // neumann values
    config["neumann.variable"] = "x";
    config["neumann.expression"] = "1.0";
    config["neumann.order"] = "0";
    config["neumann.name"] = "neumann";
    config["neumann.type"] = ExpressionFunctionType::static_id();
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
    return std::unique_ptr< ThisType >(new ThisType(create_function< 1 >("diffusion_factor", cfg),
//                        create_function< dimDomain, dimDomain >("diffusion_tensor", cfg),
                        create_function< dimRange >("force", cfg),
                        create_function< dimRange >("dirichlet", cfg),
                        create_function< dimRange >("neumann", cfg)));
  } // ... create(...)

  Default(const std::shared_ptr< const DiffusionFactorType >& diff_fac,
//          const std::shared_ptr< const DiffusionTensorType >& diff_ten,
          const std::shared_ptr< const FunctionType >& forc,
          const std::shared_ptr< const FunctionType >& dir,
          const std::shared_ptr< const FunctionType >& neum)
    : diffusion_factor_(diff_fac)
//    , diffusion_tensor_(diff_ten)
    , force_(forc)
    , dirichlet_(dir)
    , neumann_(neum)
  {
    this->inherit_parameter_type(diffusion_factor_->parameter_type(), "diffusion_factor");
//    this->inherit_parameter_type(diffusion_tensor_->parameter_type(), "diffusion_tensor");
    this->inherit_parameter_type(force_->parameter_type(),     "force");
    this->inherit_parameter_type(dirichlet_->parameter_type(), "dirichlet");
    this->inherit_parameter_type(neumann_->parameter_type(),   "neumann");
  }

  Default(const ThisType& other) = delete;

  ThisType& operator=(const ThisType& other) = delete;

  virtual const DiffusionFactorType& diffusion_factor() const DS_OVERRIDE
  {
    return *diffusion_factor_;
  }

//  virtual const DiffusionTensorType& diffusion_tensor() const DS_OVERRIDE
//  {
//    return *diffusion_tensor_;
//  }

  virtual const FunctionType& force() const DS_OVERRIDE
  {
    return *force_;
  }

  virtual const FunctionType& dirichlet() const DS_OVERRIDE
  {
    return *dirichlet_;
  }

  virtual const FunctionType& neumann() const DS_OVERRIDE
  {
    return *neumann_;
  }

protected:
  template< int r, int rC = 1 >
      static std::shared_ptr< Pymor::AffinelyDecomposableFunctionInterface< EntityType, DomainFieldType, dimDomain, RangeFieldType, r, rC > >
  create_function(const std::string& id,
                  const Stuff::Common::ConfigTree& config)
  {
    typedef Stuff::Functions::Expression< EntityType, DomainFieldType, dimDomain, RangeFieldType, r, rC >
        ExpressionFunctionType;
    typedef Pymor::AffinelyDecomposableFunctions< EntityType, DomainFieldType, dimDomain, RangeFieldType, r, rC >
        FunctionsProvider;
    const Stuff::Common::ConfigTree cfg = config.sub(id);
    const std::string type = cfg.get("type", ExpressionFunctionType::static_id());
    return FunctionsProvider::create(type, cfg);
  } // ... create_function(...)

  std::shared_ptr< const DiffusionFactorType > diffusion_factor_;
//  std::shared_ptr< const DiffusionTensorType > diffusion_tensor_;
  std::shared_ptr< const FunctionType > force_;
  std::shared_ptr< const FunctionType > dirichlet_;
  std::shared_ptr< const FunctionType > neumann_;
}; // class Default

} // namespace Problems
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_PROBLEMS_DEFAULT_HH
