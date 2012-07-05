#ifndef DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_HH
#define DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_HH

// dune-common
#include <dune/common/parametertree.hh>
#include <dune/common/shared_ptr.hh>

// dune-helper-tools
#include <dune/helper-tools/function/expression.hh>
#include <dune/helper-tools/common/parametertree.hh>

namespace Dune
{

namespace DetailedSolvers
{

namespace Stationary
{

namespace Linear
{

namespace Elliptic
{

template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim >
class Model
{
public:
  typedef DomainFieldImp DomainFieldType;

  static const int dimDomain = domainDim;

  typedef RangeFieldImp RangeFieldType;

  static const int dimRange = rangeDim;

  typedef Model< DomainFieldType, dimDomain, RangeFieldType, dimRange > ThisType;

  static const std::string id;

  typedef Dune::HelperTools::Function::Expression< DomainFieldType, dimDomain, RangeFieldType, dimRange > DiffusionType;

  typedef DiffusionType ForceType;

  Model(const Dune::ParameterTree& paramTree)
  {
    // check parametertree
    Dune::HelperTools::Common::ParameterTree::assertSub(paramTree, "diffusion", id);
    Dune::HelperTools::Common::ParameterTree::assertSub(paramTree, "force", id);
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
}; // class Model


template< class DomainFieldType, int dimDomain, class RangeFieldType, int dimRange >
const std::string Model< DomainFieldType, dimDomain, RangeFieldType, dimRange >::id = "detailed-solvers.stationary.linear.elliptic.model";

} // namespace Elliptic

} // namespace Linear

} // namespace Stationary

} // namespace DetailedSolvers

} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_HH
