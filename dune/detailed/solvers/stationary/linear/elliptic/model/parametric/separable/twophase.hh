#ifndef DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_TWOPHASE_HH
#define DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_TWOPHASE_HH

#include <dune/stuff/function/product.hh>

#include "default.hh"

namespace Dune {
namespace Detailed {
namespace Solvers {
namespace Stationary {
namespace Linear {
namespace Elliptic {
namespace Model {


// forward, to allow for specialization
template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim >
class SeparableTwoPhase;


template< class DomainFieldImp, class RangeFieldImp >
class SeparableTwoPhase< DomainFieldImp, 2, RangeFieldImp, 1 >
  : public SeparableDefault< DomainFieldImp, 2, RangeFieldImp, 1 >
{
public:
  typedef SeparableTwoPhase< DomainFieldImp, 2, RangeFieldImp, 1 >  ThisType;
  typedef SeparableDefault < DomainFieldImp, 2, RangeFieldImp, 1 >  BaseType;

  typedef typename BaseType::DomainFieldType  DomainFieldType;
  static const int                            dimDomain = BaseType::dimDomain;

  typedef typename BaseType::RangeFieldType   RangeFieldType;
  static const int                            dimRange = BaseType::dimRange;

  typedef typename BaseType::ParamFieldType   ParamFieldType;
  static const int                            maxParamDim = BaseType::maxParamDim;
  typedef typename BaseType::ParamType        ParamType;

  typedef typename BaseType::FunctionType   FunctionType;
  typedef typename FunctionType::DomainType DomainType;
  typedef typename FunctionType::RangeType  RangeType;

private:
  typedef Dune::Stuff::Function::Product< DomainFieldType, dimDomain, RangeFieldType, dimRange,
                                          DomainFieldType, dimDomain, RangeFieldType, dimRange > ProductFunctionType;

public:
  static std::string id()
  {
    return BaseType::BaseType::id() + ".parametric.separable.twophase";
  }

  SeparableTwoPhase(const Dune::shared_ptr< const FunctionType > _totalMobility,
                    const Dune::shared_ptr< const FunctionType > _permeability,
                    const Dune::shared_ptr< const FunctionType > _force,
                    const Dune::shared_ptr< const FunctionType > _dirichlet,
                    const Dune::shared_ptr< const FunctionType > _neumann)
    : BaseType(Dune::make_shared< ProductFunctionType >(_totalMobility,
                                                        _permeability),
               _force,
               _dirichlet,
               _neumann)
  {}

  static Dune::ParameterTree createSampleDescription(const std::string subName = "")
  {
    Dune::ParameterTree description;
    description["totalMobility.type"] = "function.separable.default";
    description["totalMobility.name"] = "total mobility (two drops)";
    description["totalMobility.order"] = "2";
    description["totalMobility.component.0"] = "-1.0*exp(-1.0*(((x[0]-0.75)*(x[0]-0.75))/(2*0.1*0.1)))";
    description["totalMobility.component.1"] = "-1.0*exp(-1.0*(((x[0]-4.25)*(x[0]-4.25))/(2*0.1*0.1)))";
    description["totalMobility.component.2"] = "1.0";
    description["totalMobility.coefficient.0"] = "mu[0]";
    description["totalMobility.coefficient.1"] = "mu[1]";
    description["totalMobility.paramSize"] = "2";
    description["totalMobility.paramMin"] = "[0.0; 0.0]";
    description["totalMobility.paramMax"] = "[0.95; 0.95]";
    description["totalMobility.paramExplanation"] = "[drop_at_0.75; drop_at_4.25]";
    description["permeability.type"] = "function.spe10.model1";
    description["permeability.filename"] = "perm_case1.dat";
    description["permeability.name"] = "permeability (spe10.model1)";
    description["permeability.order"] = "0";
    description["permeability.lowerLeft"] = "[0.0; 0.0]";
    description["permeability.upperRight"] = "[5.0; 1.0]";
    description["force.type"] = "function.expression";
    description["force.name"] = "force (1)";
    description["force.order"] = "0";
    description["force.variable"] = "x";
    description["force.expression"] = "1.0";
    description["dirichlet.type"] = "function.expression";
    description["dirichlet.name"] = "dirichlet (0)";
    description["dirichlet.order"] = "0";
    description["dirichlet.variable"] = "x";
    description["dirichlet.expression"] = "0.0";
    description["neumann.type"] = "function.expression";
    description["neumann.name"] = "neumann (0)";
    description["neumann.order"] = "0";
    description["neumann.variable"] = "x";
    description["neumann.expression"] = "0.0";
    if (subName.empty())
      return description;
    else {
      Dune::Stuff::Common::ExtendedParameterTree extendedDescription;
      extendedDescription.add(description, subName);
      return extendedDescription;
    }
  } // ... createSampleDescription(...)

  static ThisType createFromDescription(const Dune::ParameterTree& _description, const std::string _subName = id())
  {
    // get correct paramTree
    Dune::Stuff::Common::ExtendedParameterTree description;
    if (_description.hasSub(_subName))
      description = _description.sub(_subName);
    else
      description = _description;
    return ThisType(BaseType::createFunction("totalMobility", description),
                    BaseType::createFunction("permeability", description),
                    BaseType::createFunction("force", description),
                    BaseType::createFunction("dirichlet", description),
                    BaseType::createFunction("neumann", description));
  } // ... createFromDescription(...)
}; // class SeparableTwoPhase< DomainFieldImp, 2, RangeFieldImp, 1 >


} // namespace Model
} // namespace Elliptic
} // namespace Linear
} // namespace Stationary
} // namespace Solvers
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_TWOPHASE_HH
