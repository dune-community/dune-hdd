#ifndef DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_THERMALBLOCK_HH
#define DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_THERMALBLOCK_HH

// system
#include <vector>
#include <string>

// dune-common
#include <dune/common/shared_ptr.hh>

// dune-stuff
#include <dune/stuff/function/expression.hh>
#include <dune/stuff/function/interface.hh>
#include <dune/stuff/function/checkerboard.hh>
#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/print.hh>

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

  typedef typename BaseType::DiffusionType DiffusionType;

  typedef typename BaseType::ForceType ForceType;

  typedef typename BaseType::DirichletType DirichletType;

  typedef typename BaseType::NeumannType NeumannType;

  Thermalblock(const Dune::Stuff::Common::ExtendedParameterTree paramTree)
    : forceOrder_(-1)
    , dirichletOrder_(-1)
    , neumannOrder_(-1)
  {
    // check parametertree
    paramTree.assertSub("diffusion", id());
    paramTree.assertSub("force", id());
    paramTree.assertSub("dirichlet", id());
    paramTree.assertSub("neumann", id());
    paramTree.assertKey("diffusion.order", id());
    paramTree.assertKey("force.order", id());
    paramTree.assertKey("dirichlet.order", id());
    paramTree.assertKey("neumann.order", id());
    // build expression functions
    typedef typename Dune::Stuff::Function::Expression< DomainFieldType, dimDomain,
                                                        RangeFieldType, dimRange >
      ExpressionFunctionType;
    force_ = Dune::shared_ptr< ForceType >(new ExpressionFunctionType(paramTree.sub("force")));
    dirichlet_ = Dune::shared_ptr< DirichletType >(new ExpressionFunctionType(paramTree.sub("dirichlet")));
    neumann_ = Dune::shared_ptr< NeumannType >(new ExpressionFunctionType(paramTree.sub("neumann")));
    // build the thermalblock diffusion
    typedef typename Dune::Stuff::Function::Checkerboard< DomainFieldType, dimDomain,
                                                          RangeFieldType, dimRange >
      CheckerboardFunctionType;
    const Dune::Stuff::Common::ExtendedParameterTree subTree = paramTree.sub("diffusion");
    subTree.assertVector("lowerLeft");
    subTree.assertVector("upperRight");
    subTree.assertVector("numElements");
    subTree.assertVector("components");
    const std::vector< DomainFieldType > lowerLefts = subTree.getVector("lowerLeft", DomainFieldType(0));
    const std::vector< DomainFieldType > upperRights = subTree.getVector("upperRight", DomainFieldType(1));
    const std::vector< unsigned int > numElements = subTree.getVector("numElements", (unsigned int)(1));
    const std::vector< RangeFieldType > components = subTree.getVector("components", RangeFieldType(1));
    assert(int(lowerLefts.size()) >= dimDomain && "Given vector is too short!");
    assert(int(upperRights.size()) >= dimDomain && "Given vector is too short!");
    assert(int(numElements.size()) >= dimDomain && "Given vector is too short!");
    assert(int(components.size()) >= dimDomain && "Given vector is too short!");
    Dune::FieldVector< DomainFieldType, dimDomain > lowerLeft;
    Dune::FieldVector< DomainFieldType, dimDomain > upperRight;
    unsigned int numSubdomains = 1;
    for (int d = 0; d < dimDomain; ++d) {
      lowerLeft[d] = lowerLefts[d];
      upperRight[d] = upperRights[d];
      numSubdomains *= numElements[d];
    }
    diffusion_ = Dune::shared_ptr< DiffusionType >(new CheckerboardFunctionType(lowerLeft,
                                                                                upperRight,
                                                                                numElements,
                                                                                components));
    // set integration orders
    forceOrder_ = paramTree.sub("force").get("order", -1);
    dirichletOrder_ = paramTree.sub("dirichlet").get("order", -1);
    neumannOrder_ = paramTree.sub("neumann").get("order", -1);
  }

  Thermalblock(const ThisType& other)
    : forceOrder_(other.forceOrder_)
    , dirichletOrder_(other.dirichletOrder_)
    , neumannOrder_(other.nemannOrder_)
    , diffusion_(other.diffusion_)
    , force_(other.force_)
    , dirichlet_(other.dirichlet_)
    , neumann_(other.neumann_)
  {}

  static const std::string id()
  {
    return BaseType::id() + ".thermalblock";
  }

  virtual const Dune::shared_ptr< const DiffusionType > diffusion() const
  {
    return diffusion_;
  }

  virtual int diffusionOrder() const
  {
    return 0;
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
  ThisType& operator=(const ThisType&);

  int forceOrder_;
  int dirichletOrder_;
  int neumannOrder_;
  Dune::shared_ptr< DiffusionType > diffusion_;
  Dune::shared_ptr< ForceType > force_;
  Dune::shared_ptr< DirichletType > dirichlet_;
  Dune::shared_ptr< NeumannType > neumann_;
}; // class Thermalblock

} // namespace Model

} // namespace Elliptic

} // namespace Linear

} // namespace Stationary

} // namespace Solvers

} // namespace Detailed

} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_THERMALBLOCK_HH
