// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_MODEL_DEFAULT_HH
#define DUNE_HDD_LINEARELLIPTIC_MODEL_DEFAULT_HH

#include <memory>

#include <dune/stuff/common/string.hh>
#include <dune/stuff/common/parameter/tree.hh>
#include <dune/pymor/functions.hh>

#include "interfaces.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Problem {


template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, bool scalarDiffusion = true >
class Default
  : public ProblemInterface< DomainFieldImp, domainDim, RangeFieldImp, rangeDim, scalarDiffusion >
{
  typedef ProblemInterface<  DomainFieldImp, domainDim, RangeFieldImp, rangeDim, scalarDiffusion >  BaseType;
public:
  typedef Default< DomainFieldImp, domainDim, RangeFieldImp, rangeDim, scalarDiffusion >            ThisType;

  typedef typename BaseType::DomainFieldType  DomainFieldType;
  static const int                            dimDomain = BaseType::dimDomain;

  typedef typename BaseType::RangeFieldType   RangeFieldType;
  static const int                            dimRange = BaseType::dimRange;

  typedef typename BaseType::DiffusionType  DiffusionType;
  typedef typename BaseType::ForceType      ForceType;
  typedef typename BaseType::DirichletType  DirichletType;
  typedef typename BaseType::NeumannType    NeumannType;

  Default(const std::shared_ptr< const DiffusionType > diff,
          const std::shared_ptr< const ForceType > forc,
          const std::shared_ptr< const DirichletType > dir,
          const std::shared_ptr< const NeumannType > neum)
    : diffusion_(diff)
    , force_(forc)
    , dirichlet_(dir)
    , neumann_(neum)
  {
    Pymor::Parametric::inherit_parameter_type(diffusion_->parameter_type(), "diffusion");
    Pymor::Parametric::inherit_parameter_type(force_->parameter_type(),     "force");
    Pymor::Parametric::inherit_parameter_type(dirichlet_->parameter_type(), "dirichlet");
    Pymor::Parametric::inherit_parameter_type(neumann_->parameter_type(),   "neumann");
  }

  static const std::string static_id()
  {
    return BaseType::static_id() + ".default";
  }

  virtual const std::string id()
  {
    return BaseType::static_id() + ".default";
  }

  static Dune::ParameterTree defaultSettings(const std::string subName = "")
  {
    Dune::Stuff::Common::ExtendedParameterTree settings;
    if (dimRange == 1) {
      typedef Pymor::Function::Checkerboard< DomainFieldType, dimDomain, RangeFieldType > CheckerBoardFunction;
      settings.add(CheckerBoardFunction::defaultSettings(),
                   "diffusion");
      settings["diffusion.type"] = CheckerBoardFunction::static_id();
      settings["diffusion.parameterName"] = "diffusion";
    } else {
      settings.add(Stuff::FunctionExpression< DomainFieldType, dimDomain, RangeFieldType, dimRange >::defaultSettings(),
                   "diffusion");
      settings["diffusion.expression"] = "1.0";
      settings["diffusion.order"] = "0";
    }
    settings["diffusion.name"] = "diffusion";
    settings.add(Stuff::FunctionExpression< DomainFieldType, dimDomain, RangeFieldType, dimRange >::defaultSettings(),
                 "force");
    settings["force.name"] = "force";
    settings["force.expression"] = "1.0";
    settings["force.order"] = "0";
    settings.add(Stuff::FunctionExpression< DomainFieldType, dimDomain, RangeFieldType, dimRange >::defaultSettings(),
                 "dirichlet");
    settings["dirichlet.name"] = "dirichlet";
    settings["dirichlet.expression"] = "0.0";
    settings["dirichlet.order"] = "1";
    settings.add(Stuff::FunctionExpression< DomainFieldType, dimDomain, RangeFieldType, dimRange >::defaultSettings(),
                 "neumann");
    settings["neumann.name"] = "neumann";
    settings["neumann.expression"] = "0.0";
    settings["neumann.order"] = "0";
    if (subName.empty())
      return settings;
    else {
      Dune::Stuff::Common::ExtendedParameterTree extendedSettings;
      extendedSettings.add(settings, subName);
      return extendedSettings;
    }
  } // ... defaultSettings(...)

  static ThisType* create(const Dune::ParameterTree settings = Dune::ParameterTree(),
                          const std::string subName = static_id())
  {
    // get correct description
    Dune::Stuff::Common::ExtendedParameterTree sett;
    if (settings.hasSub(subName))
      sett = settings.sub(subName);
    else
      sett = settings;
    return new ThisType(createVectorvaluedFunction< DomainFieldType, dimDomain, RangeFieldType, dimRange >("diffusion", sett),
                        createVectorvaluedFunction< DomainFieldType, dimDomain, RangeFieldType, dimRange >("force", sett),
                        createVectorvaluedFunction< DomainFieldType, dimDomain, RangeFieldType, dimRange >("dirichlet", sett),
                        createVectorvaluedFunction< DomainFieldType, dimDomain, RangeFieldType, dimRange >("neumann", sett));
  } // ... create(...)

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

protected:
  template< class D, int d, class R, int r >
  static std::shared_ptr< const Pymor::ParametricFunctionInterface< D, d, R, r, 1 > >
  createVectorvaluedFunction(const std::string& id,
                             const Dune::Stuff::Common::ExtendedParameterTree& settings)
  {

    const Dune::Stuff::Common::ExtendedParameterTree functionSettings = settings.sub(id);
    std::string type;
    // we do this instead of using get() with a default in order to prevent get() from throwing a warning
    if (functionSettings.hasKey("type"))
      type = functionSettings.get< std::string >("type");
    else
      type = Dune::Stuff::FunctionExpression< D, d, R, r >::static_id();
    return std::shared_ptr< const Pymor::ParametricFunctionInterface< D, d, R, r, 1 > >(
          Pymor::ParametricFunctions< D, d, R, r, 1 >::create(type, functionSettings));
  } // ... createVectorvaluedFunction(...)

  std::shared_ptr< const DiffusionType > diffusion_;
  std::shared_ptr< const ForceType > force_;
  std::shared_ptr< const DirichletType > dirichlet_;
  std::shared_ptr< const NeumannType > neumann_;
}; // class Default

} // namespace Problem
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_MODEL_DEFAULT_HH
