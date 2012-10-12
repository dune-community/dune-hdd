#ifndef DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_THERMALBLOCK_HH
#define DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_THERMALBLOCK_HH

// system
#include <vector>

// dune-common
#include <dune/common/shared_ptr.hh>

// dune-stuff
#include <dune/stuff/function/expression.hh>
#include <dune/stuff/function/interface.hh>
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
class Thermalblock
  : public Interface < DomainFieldImp, domainDim, RangeFieldImp, rangeDim >
{
public:
  typedef DomainFieldImp DomainFieldType;

  static const int dimDomain = domainDim;

  typedef RangeFieldImp RangeFieldType;

  static const int dimRange = rangeDim;

  typedef Interface< DomainFieldType, dimDomain, RangeFieldType, dimRange > BaseType;

  typedef Thermalblock< DomainFieldType, dimDomain, RangeFieldType, dimRange > ThisType;

  static const std::string id;

private:
  class Diffusion
    : public Dune::Stuff::Function::Interface< DomainFieldType, dimDomain, RangeFieldType, dimRange >
  {
  public:
    typedef Dune::FieldVector< DomainFieldType, dimDomain > DomainType;

    typedef Dune::FieldVector< RangeFieldType, dimRange > RangeType;

    Diffusion(const DomainType lowerLeft,
              const DomainType upperRight,
              const std::vector< unsigned int > numElements,
              const std::vector< RangeFieldType > components)
      : lowerLeft_(lowerLeft)
      , upperRight_(upperRight)
      , numElements_(numElements)
      , components_(components)
    {}

    Diffusion(const Diffusion& other)
      : lowerLeft_(other.lowerLeft_)
      , upperRight_(other.upperRight_)
      , numElements_(other.numElements_)
      , components_(other.components_)
    {}

    virtual void evaluate(const DomainType& x, RangeType& ret) const
    {

      assert(false && "Implement me!");

    } // virtual void evaluate(const DomainType& x, RangeType& ret) const

  private:
    const DomainType lowerLeft_;
    const DomainType upperRight_;
    const std::vector< unsigned int > numElements_;
    const std::vector< RangeFieldType > components_;
  }; // class Diffusion

public:
  typedef typename BaseType::DiffusionType DiffusionType;

  typedef typename BaseType::ForceType ForceType;

  typedef typename BaseType::DirichletType DirichletType;

  Thermalblock(const Dune::Stuff::Common::ExtendedParameterTree paramTree)
    : diffusionOrder_(-1)
    , forceOrder_(-1)
    , dirichletOrder_(-1)
  {
    // check parametertree
    paramTree.assertSub("diffusion", id);
    paramTree.assertSub("force", id);
    paramTree.assertSub("dirichlet", id);
    paramTree.assertKey("diffusion.order", id);
    paramTree.assertKey("force.order", id);
    paramTree.assertKey("dirichlet.order", id);
    // build functions
    typedef Dune::Stuff::Function::Expression< DomainFieldType, dimDomain, RangeFieldType, dimRange > ExpressionFunctionType;
    force_ = Dune::shared_ptr< ForceType >(new ExpressionFunctionType(paramTree.sub("force")));
    dirichlet_ = Dune::shared_ptr< DirichletType >(new ExpressionFunctionType(paramTree.sub("dirichlet")));
    // build the thermalblock diffusion
    Dune::FieldVector< DomainFieldType, dimDomain > lowerLeft(DomainFieldType(0));
    Dune::FieldVector< DomainFieldType, dimDomain > upperRight(DomainFieldType(0));
    std::vector< unsigned int > numElements(dimDomain, 0);
    std::vector< RangeFieldType > components;

    assert(false && "Implement me!");

    Dune::Stuff::Common::ExtendedParameterTree paramTreeDiffusion = paramTree.sub("diffusion");
    diffusion_ = Dune::shared_ptr< DiffusionType >(new Diffusion(lowerLeft, upperRight, numElements, components));
    // set orders
    diffusionOrder_ = paramTree.sub("diffusion").get("order", -1);
    forceOrder_ = paramTree.sub("force").get("order", -1);
    dirichletOrder_ = paramTree.sub("dirichlet").get("order", -1);
  }

  Thermalblock(const ThisType& other)
    : diffusionOrder_(other.diffusionOrder_)
    , forceOrder_(other.forceOrder_)
    , dirichletOrder_(other.dirichletOrder_)
    , diffusion_(other.diffusion_)
    , force_(other.force_)
    , dirichlet_(other.dirichlet_)
  {}

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
}; // class Thermalblock

template< class DomainFieldType, int dimDomain, class RangeFieldType, int dimRange >
const std::string Thermalblock< DomainFieldType, dimDomain, RangeFieldType, dimRange >::id = "detailed.solvers.stationary.linear.elliptic.model.thermalblock";

} // namespace Model

} // namespace Elliptic

} // namespace Linear

} // namespace Stationary

} // namespace Solvers

} // namespace Detailed

} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_THERMALBLOCK_HH
