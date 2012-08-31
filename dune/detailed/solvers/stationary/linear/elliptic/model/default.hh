#ifndef DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_DEFAULT_HH
#define DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_DEFAULT_HH

// dune-common
#include <dune/common/parametertree.hh>
#include <dune/common/shared_ptr.hh>

// dune-stuff
#include <dune/stuff/function/expression.hh>
#include <dune/stuff/common/parameter/tree.hh>

namespace Dune
{

namespace Detailed {

namespace Solvers
{

namespace Stationary
{

namespace Linear
{

namespace Elliptic
{

namespace Model {

template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim >
class Default
{
public:
  typedef DomainFieldImp DomainFieldType;

  static const int dimDomain = domainDim;

  typedef RangeFieldImp RangeFieldType;

  static const int dimRange = rangeDim;

  typedef Default< DomainFieldType, dimDomain, RangeFieldType, dimRange > ThisType;

  static const std::string id;

  typedef Dune::Stuff::Function::Expression< DomainFieldType, dimDomain, RangeFieldType, dimRange > DiffusionType;

  typedef DiffusionType ForceType;

  typedef DiffusionType DirichletType;

  Default(const Dune::ParameterTree& paramTree)
    : diffusionOrder_(-1)
    , forceOrder_(-1)
    , dirichletOrder_(-1)
  {
    // check parametertree
    Dune::Stuff::Common::Parameter::Tree::assertSub(paramTree, "diffusion", id);
    Dune::Stuff::Common::Parameter::Tree::assertSub(paramTree, "force", id);
    Dune::Stuff::Common::Parameter::Tree::assertSub(paramTree, "dirichlet", id);
    // build functions
    diffusion_ = Dune::shared_ptr< DiffusionType >(new DiffusionType(paramTree.sub("diffusion")));
    force_ = Dune::shared_ptr< ForceType >(new ForceType(paramTree.sub("force")));
    dirichlet_ = Dune::shared_ptr< DirichletType >(new DirichletType(paramTree.sub("dirichlet")));
    // set orders
    Dune::Stuff::Common::Parameter::Tree::assertKey(paramTree, "diffusion.order", id);
    diffusionOrder_ = paramTree.sub("diffusion").get("order", -1);
    Dune::Stuff::Common::Parameter::Tree::assertKey(paramTree, "force.order", id);
    forceOrder_ = paramTree.sub("force").get("order", -1);
    Dune::Stuff::Common::Parameter::Tree::assertKey(paramTree, "dirichlet.order", id);
    dirichletOrder_ = paramTree.sub("dirichlet").get("order", -1);
  }

  const Dune::shared_ptr< DiffusionType > diffusion() const
  {
    return diffusion_;
  }

  int diffusionOrder() const
  {
    return diffusionOrder_;
  }

  const Dune::shared_ptr< ForceType > force() const
  {
    return force_;
  }

  int forceOrder() const
  {
    return forceOrder_;
  }

  const Dune::shared_ptr< DirichletType > dirichlet() const
  {
    return dirichlet_;
  }

  int dirichletOrder() const
  {
    return dirichletOrder_;
  }

private:
  int diffusionOrder_;
  int forceOrder_;
  int dirichletOrder_;
  Dune::shared_ptr< DiffusionType > diffusion_;
  Dune::shared_ptr< ForceType > force_;
  Dune::shared_ptr< DirichletType > dirichlet_;
}; // class Default

template< class DomainFieldType, int dimDomain, class RangeFieldType, int dimRange >
const std::string Default< DomainFieldType, dimDomain, RangeFieldType, dimRange >::id = "detailed.solvers.stationary.linear.elliptic.model.default";

} // namespace Model

} // namespace Elliptic

} // namespace Linear

} // namespace Stationary

} // namespace Solvers

} // namespace Detailed

} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_DEFAULT_HH
