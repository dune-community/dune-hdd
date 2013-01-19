#ifndef DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_SPE10_HH
#define DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_SPE10_HH

#include <dune/stuff/data/provider/spe10/model1.hh>
#include <dune/stuff/function/expression.hh>
#include <dune/stuff/common/parameter/tree.hh>

namespace Dune {

namespace Detailed {

namespace Solvers {

namespace Stationary {

namespace Linear {

namespace Elliptic {

namespace Spe10 {

template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim >
class One;

template< class DomainFieldImp, class RangeFieldImp >
class One< DomainFieldImp, 2, RangeFieldImp, 1 >
{
public:
  typedef DomainFieldImp DomainFieldType;

  static const int dimDomain = 2;

  typedef RangeFieldImp RangeFieldType;

  static const int dimRange = 1;

  typedef One< DomainFieldType, dimDomain, RangeFieldType, dimRange > ThisType;

  static const std::string id;

private:
  typedef Dune::Stuff::Data::Provider::Spe10::Model1::Permeability< DomainFieldType, dimDomain, RangeFieldType, dimRange > DiffusionCreator;

public:
  typedef typename DiffusionCreator::Function DiffusionType;

  typedef Dune::Stuff::Function::Expression< DomainFieldType, dimDomain, RangeFieldType, dimRange > ForceType;

  typedef ForceType DirichletType;

  One(const Dune::ParameterTree& paramTree)
    : forceOrder_(-1)
    , dirichletOrder_(-1)
  {
    // check parametertree
    Dune::Stuff::Common::Parameter::Tree::assertSub(paramTree, "diffusion", id);
    Dune::Stuff::Common::Parameter::Tree::assertSub(paramTree, "force", id);
    Dune::Stuff::Common::Parameter::Tree::assertSub(paramTree, "dirichlet", id);
    // build functions
    diffusion_ = DiffusionCreator::create(paramTree.sub("diffusion"), "  ", std::cout);
    force_ = Dune::shared_ptr< ForceType >(new ForceType(paramTree.sub("force")));
    dirichlet_ = Dune::shared_ptr< DirichletType >(new DirichletType(paramTree.sub("dirichlet")));
    // set orders
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
    return 0;
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
  int forceOrder_;
  int dirichletOrder_;
  Dune::shared_ptr< DiffusionType > diffusion_;
  Dune::shared_ptr< ForceType > force_;
  Dune::shared_ptr< DirichletType > dirichlet_;

}; // class One

template< class DomainFieldType, class RangeFieldType >
const std::string One< DomainFieldType, 2, RangeFieldType, 1 >::id = "detailed.solvers.stationary.linear.elliptic.model.spe10.one";

} // namespace Spe10

} // namespace Elliptic

} // namespace Stationary

} // namespace Stationary

} // namespace Solvers

} // namespace Detailed

} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_SPE10_HH
