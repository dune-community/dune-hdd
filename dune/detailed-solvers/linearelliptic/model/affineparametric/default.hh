#ifndef DUNE_DETAILED_SOLVERS_LINEARELLIPTIC_MODEL_AFFINEPARAMETRIC_DEFAULT_HH
#define DUNE_DETAILED_SOLVERS_LINEARELLIPTIC_MODEL_AFFINEPARAMETRIC_DEFAULT_HH

#include <vector>
#include <sstream>
#include <memory>

#include <dune/common/exceptions.hh>

#include <dune/stuff/common/string.hh>
#include <dune/stuff/common/color.hh>
#include <dune/stuff/common/parameter.hh>
#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/function.hh>
#include <dune/stuff/function/expression.hh>
#include <dune/stuff/function/affineparametric/coefficient.hh>

#include "../interface.hh"

namespace Dune {
namespace DetailedSolvers {
namespace LinearElliptic {


//// forwards, to allow for some friendlyness
//template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim >
//class ModelAffineParametricTwoPhase;
template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, bool scalarDiffusion >
class ModelAffineParametricThermalblock;



// forward to allow for specialization
template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, bool scalarDiffusion = true >
class ModelAffineParametricDefault;


template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim >
class ModelAffineParametricDefault< DomainFieldImp, domainDim, RangeFieldImp, rangeDim, true >
  : public ModelInterface< DomainFieldImp, domainDim, RangeFieldImp, rangeDim, true >
{
  typedef ModelInterface < DomainFieldImp, domainDim, RangeFieldImp, rangeDim, true >               BaseType;
public:
  typedef ModelAffineParametricDefault< DomainFieldImp, domainDim, RangeFieldImp, rangeDim, true >  ThisType;

  typedef typename BaseType::DomainFieldType  DomainFieldType;
  static const int                            dimDomain = BaseType::dimDomain;

  typedef typename BaseType::RangeFieldType   RangeFieldType;
  static const int                            dimRange = BaseType::dimRange;

  typedef Stuff::Common::Parameter::FieldType ParamFieldType;
  static const int                            maxParamDim = Stuff::Common::Parameter::maxDim;
  typedef Stuff::Common::Parameter::Type      ParamType;

  typedef typename BaseType::DiffusionType  DiffusionType;
  typedef typename BaseType::ForceType      ForceType;
  typedef typename BaseType::DirichletType  DirichletType;
  typedef typename BaseType::NeumannType    NeumannType;
//  typedef typename FunctionType::DomainType DomainType;
//  typedef typename FunctionType::RangeType  RangeType;

  static std::string id()
  {
    return BaseType::id() + ".affineparametric.default";
  }

  ModelAffineParametricDefault(const std::shared_ptr< const DiffusionType > _diffusion,
                               const std::shared_ptr< const ForceType > _force,
                               const std::shared_ptr< const DirichletType > _dirichlet,
                               const std::shared_ptr< const NeumannType > _neumann)
    : diffusion_(_diffusion)
    , force_(_force)
    , dirichlet_(_dirichlet)
    , neumann_(_neumann)
  {
    // do sanity checks
    if (!(diffusion_->parametric() || force_->parametric() || dirichlet_->parametric() || neumann_->parametric()))
      DUNE_THROW(Dune::RangeError,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " at least one of the given functions has to be parametric!");
    // compute paramSize
    paramSize_ = 0u;
    if (diffusion_->parametric())
      paramSize_ += diffusion_->paramSize();
    if (force_->parametric())
      paramSize_ += force_->paramSize();
    if (dirichlet_->parametric())
      paramSize_ += dirichlet_->paramSize();
    if (neumann_->parametric())
      paramSize_ += neumann_->paramSize();
    // compute parameter ranges and build param explanation
    ParamType paramMin(paramSize_);
    ParamType paramMax(paramSize_);
    size_t component = 0;
    std::vector< std::string > tmp;
    if (diffusion_->parametric()) {
      tmp = diffusion_->paramExplanation();
      assert(tmp.size() == diffusion_->paramSize());
      for (size_t qq = 0; qq < diffusion_->paramSize(); ++qq, ++component) {
        paramMin[component] = diffusion_->paramRange()[0][qq];
        paramMax[component] = diffusion_->paramRange()[1][qq];
        paramExplanation_.push_back("diffusion_" + Stuff::Common::toString(qq) + ": " + tmp[qq]);
      }
    }
    if (force_->parametric()) {
      tmp = force_->paramExplanation();
      assert(tmp.size() == force_->paramSize());
      for (size_t qq = 0; qq < force_->paramSize(); ++qq, ++component) {
        paramMin[component] = force_->paramRange()[0][qq];
        paramMax[component] = force_->paramRange()[1][qq];
        paramExplanation_.push_back("force_" + Stuff::Common::toString(qq) + ": " + tmp[qq]);
      }
    }
    if (dirichlet_->parametric()) {
      tmp = dirichlet_->paramExplanation();
      assert(tmp.size() == dirichlet_->paramSize());
      for (size_t qq = 0; qq < dirichlet_->paramSize(); ++qq, ++component) {
        paramMin[component] = dirichlet_->paramRange()[0][qq];
        paramMax[component] = dirichlet_->paramRange()[1][qq];
        paramExplanation_.push_back("dirichlet_" + Stuff::Common::toString(qq) + ": " + tmp[qq]);
      }
    }
    if (neumann_->parametric()) {
      tmp = neumann_->paramExplanation();
      assert(tmp.size() == neumann_->paramSize());
      for (size_t qq = 0; qq < neumann_->paramSize(); ++qq, ++component) {
        paramMin[component] = neumann_->paramRange()[0][qq];
        paramMax[component] = neumann_->paramRange()[1][qq];
        paramExplanation_.push_back("neumann_" + Stuff::Common::toString(qq) + ": " + tmp[qq]);
      }
    }
    paramRange_.push_back(paramMin);
    paramRange_.push_back(paramMax);
  } // Default(...)

  static Dune::ParameterTree defaultSettings(const std::string subName = "")
  {
    Dune::ParameterTree description;
    description["diffusion.type"] = "function.affineparametric.default";
    description["diffusion.name"] = "diffusion";
    description["diffusion.order"] = "0";
    description["diffusion.component.0"] = "x[0]";
    description["diffusion.component.1"] = "2*x[0]";
    description["diffusion.coefficient.0"] = "mu[0]";
    description["diffusion.coefficient.1"] = "mu[1]";
    description["diffusion.paramSize"] = "2";
    description["diffusion.paramMin"] = "[0.0; 0.0]";
    description["diffusion.paramMax"] = "[1.0; 1.0]";
    description["force.type"] = "function.affineparametric.default";
    description["force.name"] = "force";
    description["force.order"] = "0";
    description["force.component.0"] = "x[0]";
    description["force.component.1"] = "2*x[0]";
    description["force.coefficient.0"] = "mu[0]";
    description["force.coefficient.1"] = "mu[1]";
    description["force.paramSize"] = "2";
    description["force.paramMin"] = "[0.0; 0.0]";
    description["force.paramMax"] = "[1.0; 1.0]";
    description["dirichlet.type"] = "function.expression";
    description["dirichlet.name"] = "dirichlet";
    description["dirichlet.order"] = "0";
    description["dirichlet.variable"] = "x";
    description["dirichlet.expression"] = "[0.1*x[0]; 0.0]";
    description["neumann.type"] = "function.affineparametric.default";
    description["neumann.name"] = "neumann";
    description["neumann.order"] = "0";
    description["neumann.component.0"] = "x[0]";
    description["neumann.component.1"] = "2*x[0]";
    description["neumann.coefficient.0"] = "mu[0]";
    description["neumann.coefficient.1"] = "mu[1]";
    description["neumann.paramSize"] = "2";
    description["neumann.paramMin"] = "[0.0; 0.0]";
    description["neumann.paramMax"] = "[1.0; 1.0]";
    if (subName.empty())
      return description;
    else {
      Dune::Stuff::Common::ExtendedParameterTree extendedDescription;
      extendedDescription.add(description, subName);
      return extendedDescription;
    }
  } // defaultSettings

  static ThisType* create(const Dune::ParameterTree& _settiings, const std::string _subName = id())
  {
    // get correct paramTree
    Dune::Stuff::Common::ExtendedParameterTree settiings;
    if (_settiings.hasSub(_subName))
      settiings = _settiings.sub(_subName);
    else
      settiings = _settiings;
    return new ThisType(createFunction("diffusion", settiings),
                        createFunction("force", settiings),
                        createFunction("dirichlet", settiings),
                        createFunction("neumann", settiings));
  } // ... create(...)

  virtual size_t paramSize() const
  {
    return paramSize_;
  }

  virtual const std::vector< ParamType >& paramRange() const
  {
    return paramRange_;
  }

  virtual const std::vector< std::string >& paramExplanation() const
  {
    return paramExplanation_;
  }

  virtual ParamType mapParam(const ParamType& _globalMu, const std::string _id) const
  {
    ParamType localMu;
    // do some checks and initialize localMu tih correct size
    assert(_globalMu.size() == paramSize_);
    if (_id == "diffusion") {
      localMu = ParamType(diffusion()->paramSize());
      for (size_t pp = 0; pp < diffusion()->paramSize(); ++pp)
        localMu[pp] = _globalMu[pp];
    } else if (_id == "force") {
      localMu = ParamType(force()->paramSize());
      for (size_t pp = 0; pp < force()->paramSize(); ++pp)
        localMu[pp] = _globalMu[diffusion()->paramSize() + pp];
    } else if (_id == "dirichlet") {
      localMu = ParamType(dirichlet()->paramSize());
      for (size_t pp = 0; pp < dirichlet()->paramSize(); ++pp)
        localMu[pp] = _globalMu[diffusion()->paramSize() + force()->paramSize() + pp];
    } else if (_id == "neumann") {
      localMu = ParamType(neumann()->paramSize());
      for (size_t pp = 0; pp < neumann()->paramSize(); ++pp)
        localMu[pp] = _globalMu[diffusion()->paramSize() + force()->paramSize() + dirichlet()->paramSize() + pp];
    } else
      DUNE_THROW(Dune::RangeError,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " wrong '_id' given (has to be one of 'diffusion', 'force', 'dirichlet' or 'neumann', is '"
                 << _id << "')!");
    // the parameters' order is diffusion, force, dirichlet, neumann
    return localMu;
  } // ... mapParam(...)

  virtual std::shared_ptr< const DiffusionType > diffusion() const
  {
    return diffusion_;
  }

  virtual std::shared_ptr< const ForceType > force() const
  {
    return force_;
  }

  virtual std::shared_ptr< const DirichletType > dirichlet() const
  {
    return dirichlet_;
  }

  virtual std::shared_ptr< const NeumannType > neumann() const
  {
    return neumann_;
  }

private:
  typedef Dune::Stuff::GenericStationaryFunctionInterface< DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 >
    FunctionType;

  static std::shared_ptr< const FunctionType > createFunction(const std::string& _id,
                                                              const Dune::Stuff::Common::ExtendedParameterTree& _description)
  {
    const Dune::Stuff::Common::ExtendedParameterTree functionDescription = _description.sub(_id);
    std::string type;
    // we do this instead of using get() with 'function.expression' as a default in order to prevent get() from
    // throwing a warning
    if (functionDescription.hasKey("type"))
      type = functionDescription.get< std::string >("type");
    else
      type = "function.expression";
    return std::shared_ptr< const FunctionType >(
          Dune::Stuff::GenericStationaryFunctions<  DomainFieldType, dimDomain,
                                                    RangeFieldType, dimRange >::create(type,
                                                                                       functionDescription));
  } // ... createFunction(...)

//  friend class ModelAffineParametricTwoPhase< DomainFieldType, dimDomain, RangeFieldType, dimRange >;
  friend class ModelAffineParametricThermalblock< DomainFieldType, dimDomain, RangeFieldType, dimRange, true >;

  std::shared_ptr< const DiffusionType > diffusion_;
  std::shared_ptr< const ForceType > force_;
  std::shared_ptr< const DirichletType > dirichlet_;
  std::shared_ptr< const NeumannType > neumann_;
  size_t paramSize_;
  std::vector< ParamType > paramRange_;
  std::vector< std::string > paramExplanation_;
}; // class AffineParametricDefault


} // namespace LinearElliptic
} // namespace DetailedSolvers
} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_LINEARELLIPTIC_MODEL_AFFINEPARAMETRIC_DEFAULT_HH
