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


// forward, to allow for some friendlyness
template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim >
class ModelAffineParametricTwoPhase;


template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim >
class ModelAffineParametricDefault
  : public ModelInterface< DomainFieldImp, domainDim, RangeFieldImp, rangeDim >
{
public:
  typedef ModelAffineParametricDefault< DomainFieldImp, domainDim, RangeFieldImp, rangeDim >  ThisType;
  typedef ModelInterface < DomainFieldImp, domainDim, RangeFieldImp, rangeDim >               BaseType;

  typedef typename BaseType::DomainFieldType  DomainFieldType;
  static const int                            dimDomain = BaseType::dimDomain;

  typedef typename BaseType::RangeFieldType   RangeFieldType;
  static const int                            dimRange = BaseType::dimRange;

  typedef Stuff::Common::Parameter::FieldType ParamFieldType;
  static const int                            maxParamDim = Stuff::Common::Parameter::maxDim;
  typedef Stuff::Common::Parameter::Type      ParamType;

  typedef typename BaseType::FunctionType   FunctionType;
  typedef typename FunctionType::DomainType DomainType;
  typedef typename FunctionType::RangeType  RangeType;

  static std::string id()
  {
    return BaseType::id() + ".affineparametric.default";
  }

  ModelAffineParametricDefault(const std::shared_ptr< const FunctionType > _diffusion,
                               const std::shared_ptr< const FunctionType > _force,
                               const std::shared_ptr< const FunctionType > _dirichlet,
                               const std::shared_ptr< const FunctionType > _neumann)
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

  ModelAffineParametricDefault(const ThisType& other)
    : diffusion_(other.diffusion_)
    , force_(other.force_)
    , dirichlet_(other.dirichlet_)
    , neumann_(other.neumann_)
    , paramSize_(other.paramSize_)
    , paramRange_(other.paramRange_)
  {}

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      diffusion_ = other.diffusion();
      force_ = other.force();
      dirichlet_ = other.dirichlet();
      neumann_ = other.neumann();
      paramSize_ = other.paramSize();
      paramRange_ = other.paramRange();
    }
    return *this;
  } // ThisType& operator=(ThisType& other)

  static Dune::ParameterTree createSampleDescription(const std::string subName = "")
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
    description["dirichlet.type"] = "function.affineparametric.default";
    description["dirichlet.name"] = "dirichlet";
    description["dirichlet.order"] = "0";
    description["dirichlet.component.0"] = "x[0]";
    description["dirichlet.component.1"] = "2*x[0]";
    description["dirichlet.coefficient.0"] = "mu[0]";
    description["dirichlet.coefficient.1"] = "mu[1]";
    description["dirichlet.paramSize"] = "2";
    description["dirichlet.paramMin"] = "[0.0; 0.0]";
    description["dirichlet.paramMax"] = "[1.0; 1.0]";
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
  }

  static ThisType* create(const Dune::ParameterTree& _description, const std::string _subName = id())
  {
    // get correct paramTree
    Dune::Stuff::Common::ExtendedParameterTree description;
    if (_description.hasSub(_subName))
      description = _description.sub(_subName);
    else
      description = _description;
    return new ThisType(createFunction("diffusion", description),
                        createFunction("force", description),
                        createFunction("dirichlet", description),
                        createFunction("neumann", description));
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

  virtual std::shared_ptr< const FunctionType > diffusion() const
  {
    return diffusion_;
  }

  virtual std::shared_ptr< const FunctionType > force() const
  {
    return force_;
  }

  virtual std::shared_ptr< const FunctionType > dirichlet() const
  {
    return dirichlet_;
  }

  virtual std::shared_ptr< const FunctionType > neumann() const
  {
    return neumann_;
  }

private:
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
    return std::shared_ptr< const FunctionType >(Dune::Stuff::Functions<  DomainFieldType, dimDomain,
                                                                          RangeFieldType, dimRange >::create(type,
                                                                                                             functionDescription));
  } // ... createFunction(...)

  friend class ModelAffineParametricTwoPhase< DomainFieldType, dimDomain, RangeFieldType, dimRange >;

  std::shared_ptr< const FunctionType > diffusion_;
  std::shared_ptr< const FunctionType > force_;
  std::shared_ptr< const FunctionType > dirichlet_;
  std::shared_ptr< const FunctionType > neumann_;
  size_t paramSize_;
  std::vector< ParamType > paramRange_;
  std::vector< std::string > paramExplanation_;
}; // class AffineParametricDefault


} // namespace LinearElliptic
} // namespace DetailedSolvers
} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_LINEARELLIPTIC_MODEL_AFFINEPARAMETRIC_DEFAULT_HH
