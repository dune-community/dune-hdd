#ifndef DUNE_DETAILED_SOLVERS_LINEARELLIPTIC_MODEL_AFFINEPARAMETRIC_TWOPHASE_HH
#define DUNE_DETAILED_SOLVERS_LINEARELLIPTIC_MODEL_AFFINEPARAMETRIC_TWOPHASE_HH

#include <memory>

#include <dune/stuff/function/product.hh>

#include "../default.hh"
#include "default.hh"

namespace Dune {
namespace DetailedSolvers {
namespace LinearElliptic {


// forward, to allow for specialization
template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, bool scalarDiffusion = true >
class ModelAffineParametricTwoPhase;


template< class DomainFieldImp, int domainDim, class RangeFieldImp >
class ModelAffineParametricTwoPhase< DomainFieldImp, domainDim, RangeFieldImp, 1, true >
  : public ModelAffineParametricDefault< DomainFieldImp, domainDim, RangeFieldImp, 1, true >
{
  typedef ModelAffineParametricDefault < DomainFieldImp, domainDim, RangeFieldImp, 1, true >  BaseType;
public:
  typedef ModelAffineParametricTwoPhase< DomainFieldImp, domainDim, RangeFieldImp, 1, true >  ThisType;

  typedef typename BaseType::DomainFieldType  DomainFieldType;
  static const int                            dimDomain = BaseType::dimDomain;

  typedef typename BaseType::RangeFieldType   RangeFieldType;
  static const int                            dimRange = BaseType::dimRange;

  typedef typename BaseType::ParamFieldType   ParamFieldType;
  static const int                            maxParamDim = BaseType::maxParamDim;
  typedef typename BaseType::ParamType        ParamType;

  typedef typename BaseType::ForceType      ForceType;
  typedef typename BaseType::DirichletType  DirichletType;
  typedef typename BaseType::NeumannType    NeumannType;

  typedef typename BaseType::DiffusionType                                                        TotalMobilityTypeType;
  typedef Dune::Stuff::FunctionInterface< DomainFieldType, dimDomain, RangeFieldType, dimRange >  PermeabilityType;
  typedef Dune::Stuff::FunctionProduct< DomainFieldType, dimDomain, RangeFieldType, dimRange,
                                        DomainFieldType, dimDomain, RangeFieldType, dimRange >    DiffusionType;

  static std::string id()
  {
    return BaseType::BaseType::id() + ".affineparametric.twophase";
  }

  ModelAffineParametricTwoPhase(const std::shared_ptr< const TotalMobilityTypeType > _totalMobility,
                                const std::shared_ptr< const PermeabilityType > _permeability,
                                const std::shared_ptr< const ForceType > _force,
                                const std::shared_ptr< const DirichletType > _dirichlet,
                                const std::shared_ptr< const NeumannType > _neumann)
    : BaseType(std::make_shared< DiffusionType >(_totalMobility, _permeability),
               _force,
               _dirichlet,
               _neumann)
  {}

  static Dune::ParameterTree createSampleDescription(const std::string subName = "")
  {
    Dune::ParameterTree description;
    description["totalMobility.type"] = "function.affineparametric.default";
    description["totalMobility.name"] = "total mobility (two drops)";
    description["totalMobility.order"] = "2";
    description["totalMobility.component.0"] = "-1.0*exp(-1.0*(((x[0]-0.75)*(x[0]-0.75))/(2*0.1*0.1)))";
    description["totalMobility.component.1"] = "-1.0*exp(-1.0*(((x[0]-4.25)*(x[0]-4.25))/(2*0.1*0.1)))";
    description["totalMobility.affinePart"] = "1.0";
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

  static ThisType* create(const Dune::ParameterTree& _settings, const std::string _subName = id())
  {
    // get correct paramTree
    Dune::Stuff::Common::ExtendedParameterTree settings;
    if (_settings.hasSub(_subName))
      settings = _settings.sub(_subName);
    else
      settings = _settings;
    return new ThisType(BaseType::createFunction("totalMobility", settings),
                        BaseType::createFunction("permeability", settings),
                        BaseType::createFunction("force", settings),
                        BaseType::createFunction("dirichlet", settings),
                        BaseType::createFunction("neumann", settings));
  } // ... create(...)
}; // class ModelAffineParametricTwoPhase


} // namespace LinearElliptic
} // namespace DetailedSolvers
} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_LINEARELLIPTIC_MODEL_AFFINEPARAMETRIC_TWOPHASE_HH
