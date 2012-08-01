#ifndef DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_HH
#define DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_HH

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

  Default(const Dune::ParameterTree& paramTree)
  {
    // check parametertree
    Dune::Stuff::Common::Parameter::Tree::assertSub(paramTree, "diffusion", id);
    Dune::Stuff::Common::Parameter::Tree::assertSub(paramTree, "force", id);
    // build functions
    diffusion_ = Dune::shared_ptr< DiffusionType >(new DiffusionType(paramTree.sub("diffusion")));
    force_ = Dune::shared_ptr< ForceType >(new ForceType(paramTree.sub("force")));
  }

  const Dune::shared_ptr< DiffusionType > diffusion() const
  {
    return diffusion_;
  }

  const Dune::shared_ptr< ForceType > force() const
  {
    return force_;
  }

private:
  Dune::shared_ptr< DiffusionType > diffusion_;
  Dune::shared_ptr< ForceType > force_;
}; // class Default

template< class DomainFieldType, int dimDomain, class RangeFieldType, int dimRange >
const std::string Default< DomainFieldType, dimDomain, RangeFieldType, dimRange >::id = "detailed.solvers.stationary.linear.elliptic.model.default";

template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim >
class TwoPhase
{
public:
  typedef DomainFieldImp DomainFieldType;

  static const int dimDomain = domainDim;

  typedef RangeFieldImp RangeFieldType;

  static const int dimRange = rangeDim;

  typedef TwoPhase< DomainFieldType, dimDomain, RangeFieldType, dimRange > ThisType;

  static const std::string id;

  typedef Dune::Stuff::Function::Expression< DomainFieldType, dimDomain, RangeFieldType, dimRange > ForceType;

  typedef ForceType TotalMobilityType;

  typedef ForceType PermeabilityType;

  TwoPhase(const Dune::ParameterTree& paramTree)
  {
    // check parametertree
    Dune::Stuff::Common::Parameter::Tree::assertSub(paramTree, "totalMobility", id);
    Dune::Stuff::Common::Parameter::Tree::assertSub(paramTree, "permeability", id);
    Dune::Stuff::Common::Parameter::Tree::assertSub(paramTree, "force", id);
    // build functions
    totalMobility_ = Dune::shared_ptr< TotalMobilityType >(new TotalMobilityType(paramTree.sub("totalMobility")));
    permeability_ = Dune::shared_ptr< PermeabilityType >(new PermeabilityType(paramTree.sub("permeability")));
    force_ = Dune::shared_ptr< ForceType >(new ForceType(paramTree.sub("force")));
  }

  const Dune::shared_ptr< TotalMobilityType > totalMobility() const
  {
    return totalMobility_;
  }

  const Dune::shared_ptr< PermeabilityType > permeability() const
  {
    return permeability_;
  }

  const Dune::shared_ptr< ForceType > force() const
  {
    return force_;
  }

private:
  Dune::shared_ptr< TotalMobilityType > totalMobility_;
  Dune::shared_ptr< PermeabilityType > permeability_;
  Dune::shared_ptr< ForceType > force_;
}; // class TwoPhase

template< class DomainFieldType, int dimDomain, class RangeFieldType, int dimRange >
const std::string TwoPhase< DomainFieldType, dimDomain, RangeFieldType, dimRange >::id = "detaileddsolvers.stationary.linear.elliptic.model.twophase";

} // namespace Model

} // namespace Elliptic

} // namespace Linear

} // namespace Stationary

} // namespace Solvers

} // namespace Detailed

} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_HH
