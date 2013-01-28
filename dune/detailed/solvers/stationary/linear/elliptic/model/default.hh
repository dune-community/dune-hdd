#ifndef DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_DEFAULT_HH
#define DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_DEFAULT_HH

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#include <dune/common/shared_ptr.hh>

#include <dune/stuff/function/expression.hh>
#include <dune/stuff/common/parameter/tree.hh>

#include "interface.hh"

namespace Dune {
namespace Detailed {
namespace Solvers {
namespace Stationary {
namespace Linear {
namespace Elliptic {
namespace Model {

template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim >
class Default
  : public Interface < DomainFieldImp, domainDim, RangeFieldImp, rangeDim >
{
public:
  typedef DomainFieldImp DomainFieldType;

  static const int dimDomain = domainDim;

  typedef RangeFieldImp RangeFieldType;

  static const int dimRange = rangeDim;

  typedef Interface< DomainFieldType, dimDomain, RangeFieldType, dimRange > BaseType;

  typedef Default< DomainFieldType, dimDomain, RangeFieldType, dimRange > ThisType;

  typedef typename BaseType::DiffusionType DiffusionType;

  typedef typename BaseType::ForceType ForceType;

  typedef typename BaseType::DirichletType DirichletType;

  typedef typename BaseType::NeumannType NeumannType;

  static const std::string id()
  {
    return BaseType::id() + ".default";
  }

  Default(const Dune::shared_ptr< const DiffusionType > _diffusion,
          const Dune::shared_ptr< const ForceType > _force,
          const Dune::shared_ptr< const DirichletType > _dirichlet,
          const Dune::shared_ptr< const NeumannType > _neumann,
          const int _diffusionOrder = -1,
          const int _forceOrder = -1,
          const int _dirichletOrder = -1,
          const int _neumannOrder = -1)
    : diffusion_(_diffusion)
    , force_(_force)
    , dirichlet_(_dirichlet)
    , neumann_(_neumann)
    , diffusionOrder_(_diffusionOrder)
    , forceOrder_(_forceOrder)
    , dirichletOrder_(_dirichletOrder)
    , neumannOrder_(_neumannOrder)
  {}

  Default(const ThisType& other)
    : diffusion_(other.diffusion_)
    , force_(other.force_)
    , dirichlet_(other.dirichlet_)
    , neumann_(other.neumann_)
    , diffusionOrder_(other.diffusionOrder_)
    , forceOrder_(other.forceOrder_)
    , dirichletOrder_(other.dirichletOrder_)
    , neumannOrder_(other.neumannOrder_)
  {}

  static ThisType createFromParamTree(const Dune::ParameterTree& paramTree, const std::string subName = id())
  {
    // get correct paramTree
    Dune::Stuff::Common::ExtendedParameterTree extendedParamTree;
    if (paramTree.hasSub(subName))
      extendedParamTree = paramTree.sub(subName);
    else
      extendedParamTree = paramTree;
    // build functions
    typedef Dune::Stuff::Function::Expression<  DomainFieldType, dimDomain,
                                                RangeFieldType, dimRange >
      ExpressionFunctionType;
    //   * diffusion
    if (!extendedParamTree.hasSub("diffusion"))
      DUNE_THROW(Dune::RangeError,
                 "\nError: sub 'diffusion' not found in the following Dune::ParameterTree:\n" << extendedParamTree.reportString("  "));
    //!TODO make_shared
    const Dune::shared_ptr< const ExpressionFunctionType > _diffusion = Dune::shared_ptr< ExpressionFunctionType >(new ExpressionFunctionType(
        ExpressionFunctionType::createFromParamTree(extendedParamTree.sub("diffusion"))));
    if (!extendedParamTree.sub("diffusion").hasKey("order"))
      std::cout << "Warning in " << id() << ": no key 'diffusion.order' given, defaulting to 0!" << std::endl;
    const int _diffusionOrder = extendedParamTree.sub("diffusion").get("order", 0);
    //   * force
    if (!extendedParamTree.hasSub("force"))
      DUNE_THROW(Dune::RangeError,
                 "\nError: sub 'force' not found in the following Dune::ParameterTree:\n" << extendedParamTree.reportString("  "));
    const Dune::shared_ptr< const ExpressionFunctionType > _force = Dune::shared_ptr< ExpressionFunctionType >(new ExpressionFunctionType(
        ExpressionFunctionType::createFromParamTree(extendedParamTree.sub("force"))));
    if (!extendedParamTree.sub("force").hasKey("order"))
      std::cout << "Warning in " << id() << ": no key 'force.order' given, defaulting to 0!" << std::endl;
    const int _forceOrder = extendedParamTree.sub("force").get("order", 0);
    //   * dirichlet
    if (!extendedParamTree.hasSub("dirichlet"))
      DUNE_THROW(Dune::RangeError,
                 "\nError: sub 'dirichlet' not found in the following Dune::ParameterTree:\n" << extendedParamTree.reportString("  "));
    const Dune::shared_ptr< const ExpressionFunctionType > _dirichlet = Dune::shared_ptr< ExpressionFunctionType >(new ExpressionFunctionType(
        ExpressionFunctionType::createFromParamTree(extendedParamTree.sub("dirichlet"))));
    if (!extendedParamTree.sub("dirichlet").hasKey("order"))
      std::cout << "Warning in " << id() << ": no key 'dirichlet.order' given, defaulting to 0!" << std::endl;
    const int _dirichletOrder = extendedParamTree.sub("dirichlet").get("order", 0);
    //   * neumann
    if (!extendedParamTree.hasSub("neumann"))
      DUNE_THROW(Dune::RangeError,
                 "\nError: sub 'neumann' not found in the following Dune::ParameterTree:\n" << extendedParamTree.reportString("  "));
    const Dune::shared_ptr< const ExpressionFunctionType > _neumann = Dune::shared_ptr< ExpressionFunctionType >(new ExpressionFunctionType(
        ExpressionFunctionType::createFromParamTree(extendedParamTree.sub("neumann"))));
    if (!extendedParamTree.sub("neumann").hasKey("order"))
      std::cout << "Warning in " << id() << ": no key 'neumann.order' given, defaulting to 0!" << std::endl;
    const int _neumannOrder = extendedParamTree.sub("neumann").get("order", 0);
    // create and return
    Default d(_diffusion, _force, _dirichlet, _neumann, _diffusionOrder, _forceOrder, _dirichletOrder, _neumannOrder);
    return d;
  } // static ThisType createFromParamTree(const Dune::ParameterTree& paramTree)

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      diffusion_ = other.diffusion();
      force_ = other.force();
      dirichlet_ = other.dirichlet();
      neumann_ = other.neumann();
      diffusionOrder_ = other.diffusionOrder();
      forceOrder_ = other.forceOrder();
      dirichletOrder_ = other.dirichletOrder();
      neumannOrder_ = other.neumannOrder();
    }
    return this;
  }

  virtual const Dune::shared_ptr< const DiffusionType > diffusion() const
  {
    return diffusion_;
  }

  virtual int diffusionOrder() const
  {
    return diffusionOrder_;
  }

  virtual const Dune::shared_ptr< const ForceType > force() const
  {
    return force_;
  }

  virtual int forceOrder() const
  {
    return forceOrder_;
  }

  virtual const Dune::shared_ptr< const DirichletType > dirichlet() const
  {
    return dirichlet_;
  }

  virtual int dirichletOrder() const
  {
    return dirichletOrder_;
  }

  virtual const Dune::shared_ptr< const NeumannType > neumann() const
  {
    return neumann_;
  }

  virtual int neumannOrder() const
  {
    return neumannOrder_;
  }

private:
  Dune::shared_ptr< const DiffusionType > diffusion_;
  Dune::shared_ptr< const ForceType > force_;
  Dune::shared_ptr< const DirichletType > dirichlet_;
  Dune::shared_ptr< const NeumannType > neumann_;
  int diffusionOrder_;
  int forceOrder_;
  int dirichletOrder_;
  int neumannOrder_;
}; // class Default

} // namespace Model
} // namespace Elliptic
} // namespace Linear
} // namespace Stationary
} // namespace Solvers
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_DEFAULT_HH
