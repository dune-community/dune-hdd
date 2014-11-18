// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_PROBLEMS_DEFAULT_HH
#define DUNE_HDD_LINEARELLIPTIC_PROBLEMS_DEFAULT_HH

#include <memory>

#if HAVE_ALUGRID
# include <dune/stuff/common/disable_warnings.hh>
#   include <dune/grid/alugrid.hh>
# include <dune/stuff/common/reenable_warnings.hh>
#endif

#include <dune/stuff/common/string.hh>
#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/expression.hh>

#include <dune/pymor/functions.hh>
#include <dune/pymor/functions/default.hh>

#include "interfaces.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Problems {


template< class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim >
class Default
  : public ProblemInterface< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim >
{
protected:
  typedef ProblemInterface< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim > BaseType;
public:
  typedef Default< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim > ThisType;

  using typename BaseType::EntityType;
  using typename BaseType::DomainFieldType;
  using BaseType::dimDomain;
  using typename BaseType::RangeFieldType;
  using BaseType::dimRange;

  using typename BaseType::DiffusionFactorType;
  using typename BaseType::DiffusionTensorType;
  using typename BaseType::FunctionType;

  typedef typename DiffusionFactorType::NonparametricType NonparametricDiffusionFactorType;
  typedef typename DiffusionTensorType::NonparametricType NonparametricDiffusionTensorType;
  typedef typename FunctionType::NonparametricType NonparametricFunctionType;

private:
  typedef Pymor::Functions::NonparametricDefault
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, 1 >                    DiffusionFactorWrapperType;
  typedef Pymor::Functions::NonparametricDefault
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimDomain, dimDomain > DiffusiontensorWrapperType;
  typedef Pymor::Functions::NonparametricDefault
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >             FunctionWrapperType;

public:
  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::static_id() + ".default";
  }

  static Stuff::Common::Configuration default_config(const std::string sub_name = "")
  {
    Stuff::Common::Configuration config;
    // diffusion factor
    typedef Pymor::Functions::Checkerboard< EntityType, DomainFieldType, dimDomain, RangeFieldType, 1 >
        CheckerBoardFunctionType;
    config.add(CheckerBoardFunctionType::default_config(), "diffusion_factor");
    config["diffusion_factor.parameter_name"] = "diffusion_factor";
    config["diffusion_factor.name"] = "diffusion_factor";
    config["diffusion_factor.type"] = CheckerBoardFunctionType::static_id();
    // diffusion tensor
    typedef Stuff::Functions::Constant< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimDomain, dimDomain >
        ConstantFunctionType;
    config.add(ConstantFunctionType::default_config(), "diffusion_tensor");
    config["diffusion_tensor.name"] = "diffusion_tensor";
    config["diffusion_tensor.type"] = ConstantFunctionType::static_id();
    // force
    typedef Stuff::Functions::Expression< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >
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
      Stuff::Common::Configuration tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  static std::unique_ptr< ThisType > create(const Stuff::Common::Configuration config = default_config(),
                                            const std::string sub_name = static_id())
  {
    const Stuff::Common::Configuration cfg = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    return std::unique_ptr< ThisType >(new ThisType(create_scalar_function("diffusion_factor", cfg),
                                                    create_matrix_function("diffusion_tensor", cfg),
                                                    create_vector_function("force", cfg),
                                                    create_vector_function("dirichlet", cfg),
                                                    create_vector_function("neumann", cfg)));
  } // ... create(...)

  Default(const std::shared_ptr< const NonparametricDiffusionFactorType >& diff_fac,
          const std::shared_ptr< const NonparametricDiffusionTensorType >& diff_ten,
          const std::shared_ptr< const NonparametricFunctionType >& forc,
          const std::shared_ptr< const NonparametricFunctionType >& dir,
          const std::shared_ptr< const NonparametricFunctionType >& neum)
    : diffusion_factor_(std::make_shared< DiffusionFactorWrapperType >(diff_fac))
    , diffusion_tensor_(std::make_shared< DiffusiontensorWrapperType >(diff_ten))
    , force_(std::make_shared< FunctionWrapperType >(forc))
    , dirichlet_(std::make_shared< FunctionWrapperType >(dir))
    , neumann_(std::make_shared< FunctionWrapperType >(neum))
  {}

  Default(const std::shared_ptr< const DiffusionFactorType >& diff_fac,
          const std::shared_ptr< const DiffusionTensorType >& diff_ten,
          const std::shared_ptr< const FunctionType >& forc,
          const std::shared_ptr< const FunctionType >& dir,
          const std::shared_ptr< const FunctionType >& neum)
    : diffusion_factor_(diff_fac)
    , diffusion_tensor_(diff_ten)
    , force_(forc)
    , dirichlet_(dir)
    , neumann_(neum)
  {
    update_parameter_dependency();
  }

  Default(const ThisType& other)
    : diffusion_factor_(other.diffusion_factor_)
    , diffusion_tensor_(other.diffusion_tensor_)
    , force_(other.force_)
    , dirichlet_(other.dirichlet_)
    , neumann_(other.neumann_)
  {
    update_parameter_dependency();
  }

  ThisType& operator=(const ThisType& other) = delete;

  virtual std::string type() const override
  {
    return BaseType::static_id() + ".default";
  }

  virtual const std::shared_ptr< const DiffusionFactorType >& diffusion_factor() const override
  {
    return diffusion_factor_;
  }

  virtual const std::shared_ptr< const DiffusionTensorType >& diffusion_tensor() const override
  {
    return diffusion_tensor_;
  }

  virtual const std::shared_ptr< const FunctionType >& force() const override
  {
    return force_;
  }

  virtual const std::shared_ptr< const FunctionType >& dirichlet() const override
  {
    return dirichlet_;
  }

  virtual const std::shared_ptr< const FunctionType >& neumann() const override
  {
    return neumann_;
  }

protected:
      static std::shared_ptr< Pymor::AffinelyDecomposableFunctionInterface< EntityType, DomainFieldType, dimDomain, RangeFieldType, 1, 1 > >
  create_scalar_function(const std::string& id,
                  const Stuff::Common::Configuration& config)
  {
    typedef Stuff::Functions::Expression< EntityType, DomainFieldType, dimDomain, RangeFieldType, 1, 1 >
        ExpressionFunctionType;
    typedef Pymor::AffinelyDecomposableFunctionsProvider< EntityType, DomainFieldType, dimDomain, RangeFieldType, 1, 1 >
        FunctionsProvider;
    const Stuff::Common::Configuration cfg = config.sub(id);
    const std::string type = cfg.get("type", ExpressionFunctionType::static_id());
    return FunctionsProvider::create(type, cfg);
  } // ... create_scalar_function(...)

      static std::shared_ptr< Pymor::AffinelyDecomposableFunctionInterface< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 > >
  create_vector_function(const std::string& id,
                  const Stuff::Common::Configuration& config)
  {
    typedef Stuff::Functions::Expression< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 >
        ExpressionFunctionType;
    typedef Pymor::AffinelyDecomposableFunctionsProvider< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 >
        FunctionsProvider;
    const Stuff::Common::Configuration cfg = config.sub(id);
    const std::string type = cfg.get("type", ExpressionFunctionType::static_id());
    return FunctionsProvider::create(type, cfg);
  } // ... create_vector_function(...)

      static std::shared_ptr< Pymor::AffinelyDecomposableFunctionInterface< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimDomain, dimDomain > >
  create_matrix_function(const std::string& id, const Stuff::Common::Configuration& config)
  {
    typedef Stuff::Functions::Constant< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimDomain, dimDomain >
        ConstantFunctionType;
    typedef Pymor::AffinelyDecomposableFunctionsProvider< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimDomain, dimDomain >
        FunctionsProvider;
    const Stuff::Common::Configuration cfg = config.sub(id);
    const std::string type = cfg.get("type", ConstantFunctionType::static_id());
    return FunctionsProvider::create(type, cfg);
  } // ... create_matrix_function(...)

  void update_parameter_dependency()
  {
    this->inherit_parameter_type(diffusion_factor_->parameter_type(), "diffusion_factor");
    this->inherit_parameter_type(diffusion_tensor_->parameter_type(), "diffusion_tensor");
    this->inherit_parameter_type(force_->parameter_type(),     "force");
    this->inherit_parameter_type(dirichlet_->parameter_type(), "dirichlet");
    this->inherit_parameter_type(neumann_->parameter_type(),   "neumann");
  } // ... update_parameter_dependency()

  std::shared_ptr< const DiffusionFactorType > diffusion_factor_;
  std::shared_ptr< const DiffusionTensorType > diffusion_tensor_;
  std::shared_ptr< const FunctionType > force_;
  std::shared_ptr< const FunctionType > dirichlet_;
  std::shared_ptr< const FunctionType > neumann_;
}; // class Default


template< class P >
class ConvertToDefault
{
  typedef typename P::EntityType E;
  typedef typename P::DomainFieldType D;
  static const unsigned int d = P::dimDomain;
  typedef typename P::RangeFieldType R;
  static const unsigned int r = P::dimRange;
  static_assert(std::is_base_of< ProblemInterface< E, D, d, R, r >, P >::value,
                "P has to be derived from ProblemInterface!");
public:
  typedef Default< E, D, d, R, r > Type;
}; // class ConvertToDefault


} // namespace Problems
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_PROBLEMS_DEFAULT_HH
