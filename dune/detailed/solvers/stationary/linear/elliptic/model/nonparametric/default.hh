#ifndef DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_NONPARAMETRIC_DEFAULT_HH
#define DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_NONPARAMETRIC_DEFAULT_HH

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#include <dune/common/shared_ptr.hh>

#include <dune/stuff/function/nonparametric/expression.hh>
#include <dune/stuff/common/parameter/tree.hh>

#include "../interface.hh"

namespace Dune {
namespace Detailed {
namespace Solvers {
namespace Stationary {
namespace Linear {
namespace Elliptic {
namespace Model {


template< class DomainFieldImp, int domainDim,
          class RangeFieldImp, int rangeDim,
          class ParamFieldImp = double, int maxParamDim = 0 >
class NonparametricDefault;


template< class DomainFieldImp, int domainDim,
          class RangeFieldImp,
          class ParamFieldImp, int maxParamDim >
class NonparametricDefault< DomainFieldImp, domainDim, RangeFieldImp, 1, ParamFieldImp, maxParamDim >
  : public Interface< DomainFieldImp, domainDim, RangeFieldImp, 1, ParamFieldImp, maxParamDim >
{
public:
  typedef Interface<  DomainFieldImp, domainDim,
                      RangeFieldImp, 1,
                      ParamFieldImp, maxParamDim >            BaseType;
  typedef NonparametricDefault< DomainFieldImp, domainDim,
                                RangeFieldImp, 1,
                                ParamFieldImp, maxParamDim >  ThisType;

  typedef typename BaseType::DomainFieldType  DomainFieldType;
  static const int                            dimDomain = BaseType::dimDomain;

  typedef typename BaseType::RangeFieldType   RangeFieldType;
  static const int                            dimRange = BaseType::dimRange;

  typedef typename BaseType::ParamFieldType ParamFieldType;
  static const int                          maxDimParam = BaseType::maxDimParam;
  typedef typename BaseType::size_type      size_type;

  typedef typename BaseType::FunctionType FunctionType;

  static const std::string id()
  {
    return BaseType::id() + ".nonparametric.default";
  }

  NonparametricDefault(const Dune::shared_ptr< const FunctionType > _diffusion,
                       const Dune::shared_ptr< const FunctionType > _force,
                       const Dune::shared_ptr< const FunctionType > _dirichlet,
                       const Dune::shared_ptr< const FunctionType > _neumann)
    : diffusion_(_diffusion)
    , force_(_force)
    , dirichlet_(_dirichlet)
    , neumann_(_neumann)
  {}

  NonparametricDefault(const ThisType& _other)
    : diffusion_(_other.diffusion_)
    , force_(_other.force_)
    , dirichlet_(_other.dirichlet_)
    , neumann_(_other.neumann_)
  {}

  ThisType& operator=(const ThisType& _other)
  {
    if (this != &_other) {
      diffusion_ = _other.diffusion();
      force_ = _other.force();
      dirichlet_ = _other.dirichlet();
      neumann_ = _other.neumann();
    }
    return this;
  }

  static ThisType createFromParamTree(const Dune::ParameterTree& _paramTree, const std::string _subName = id())
  {
    // get correct paramTree
    Dune::Stuff::Common::ExtendedParameterTree paramTree;
    if (_paramTree.hasSub(_subName))
      paramTree = _paramTree.sub(_subName);
    else
      paramTree = _paramTree;
    // create functions
    typedef Dune::Stuff::Function::NonparametricExpression<   DomainFieldType, dimDomain,
                                                              RangeFieldType, dimRange > ExpressionFunctionType;
    // * diffusion
    Dune::shared_ptr< ExpressionFunctionType > _diffusion = Dune::make_shared< ExpressionFunctionType >(
        ExpressionFunctionType::createFromParamTree(paramTree.sub("diffusion")));
    // * force
    Dune::shared_ptr< ExpressionFunctionType > _force = Dune::make_shared< ExpressionFunctionType >(
        ExpressionFunctionType::createFromParamTree(paramTree.sub("force")));
    // * dirichlet
    Dune::shared_ptr< ExpressionFunctionType > _dirichlet = Dune::make_shared< ExpressionFunctionType >(
        ExpressionFunctionType::createFromParamTree(paramTree.sub("dirichlet")));
    // * neumann
    Dune::shared_ptr< ExpressionFunctionType > _neumann = Dune::make_shared< ExpressionFunctionType >(
        ExpressionFunctionType::createFromParamTree(paramTree.sub("neumann")));
    // create and return
    return ThisType(_diffusion, _force, _dirichlet, _neumann);
  } // static ThisType createFromParamTree(const Dune::ParameterTree& paramTree)

  virtual bool parametric() const
  {
    return false;
  }

  virtual Dune::shared_ptr< const FunctionType > diffusion() const
  {
    return diffusion_;
  }

  virtual Dune::shared_ptr< const FunctionType > force() const
  {
    return force_;
  }

  virtual Dune::shared_ptr< const FunctionType > dirichlet() const
  {
    return dirichlet_;
  }

  virtual Dune::shared_ptr< const FunctionType > neumann() const
  {
    return neumann_;
  }

private:
  Dune::shared_ptr< const FunctionType > diffusion_;
  Dune::shared_ptr< const FunctionType > force_;
  Dune::shared_ptr< const FunctionType > dirichlet_;
  Dune::shared_ptr< const FunctionType > neumann_;
}; // class Default

} // namespace Model
} // namespace Elliptic
} // namespace Linear
} // namespace Stationary
} // namespace Solvers
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_NONPARAMETRIC_DEFAULT_HH
