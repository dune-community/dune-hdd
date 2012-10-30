#ifndef DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_DEFAULT_HH
#define DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_DEFAULT_HH

// dune-common
#include <dune/common/shared_ptr.hh>

// dune-stuff
#include <dune/stuff/function/expression.hh>
#include <dune/stuff/common/parameter/tree.hh>

// local
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

  Default(const Dune::Stuff::Common::ExtendedParameterTree paramTree)
    : diffusionOrder_(-1)
    , forceOrder_(-1)
    , dirichletOrder_(-1)
  {
    // check parametertree
    paramTree.assertSub("diffusion", id());
    paramTree.assertSub("force", id());
    paramTree.assertSub("dirichlet", id());
    paramTree.assertKey("diffusion.order", id());
    paramTree.assertKey("force.order", id());
    paramTree.assertKey("dirichlet.order", id());
    // build functions
    typedef Dune::Stuff::Function::Expression< DomainFieldType, dimDomain, RangeFieldType, dimRange > ExpressionFunctionType;
    diffusion_ = Dune::shared_ptr< DiffusionType >(new ExpressionFunctionType(paramTree.sub("diffusion")));
    force_ = Dune::shared_ptr< ForceType >(new ExpressionFunctionType(paramTree.sub("force")));
    dirichlet_ = Dune::shared_ptr< DirichletType >(new ExpressionFunctionType(paramTree.sub("dirichlet")));
    // set orders
    diffusionOrder_ = paramTree.sub("diffusion").get("order", -1);
    forceOrder_ = paramTree.sub("force").get("order", -1);
    dirichletOrder_ = paramTree.sub("dirichlet").get("order", -1);
  }

  Default(const ThisType& other)
    : diffusionOrder_(other.diffusionOrder_)
    , forceOrder_(other.forceOrder_)
    , dirichletOrder_(other.dirichletOrder_)
    , diffusion_(other.diffusion_)
    , force_(other.force_)
    , dirichlet_(other.dirichlet_)
  {}

  static const std::string id()
  {
    return "detailed.solvers.stationary.linear.elliptic.model.default";
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

private:
  ThisType& operator=(const ThisType&);

  int diffusionOrder_;
  int forceOrder_;
  int dirichletOrder_;
  Dune::shared_ptr< DiffusionType > diffusion_;
  Dune::shared_ptr< ForceType > force_;
  Dune::shared_ptr< DirichletType > dirichlet_;
}; // class Default

} // namespace Model

} // namespace Elliptic

} // namespace Linear

} // namespace Stationary

} // namespace Solvers

} // namespace Detailed

} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_DEFAULT_HH
