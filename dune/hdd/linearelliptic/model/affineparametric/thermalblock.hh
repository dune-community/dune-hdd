#ifndef DUNE_HDD_LINEARELLIPTIC_MODEL_AFFINEPARAMETRIC_THERMALBLOCK_HH
#define DUNE_HDD_LINEARELLIPTIC_MODEL_AFFINEPARAMETRIC_THERMALBLOCK_HH

#include <memory>
#include <vector>

#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/string.hh>
#include <dune/stuff/function/affineparametric/checkerboard.hh>
#include <dune/stuff/function.hh>

#include "../default.hh"
#include "default.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {


// forward, to allow for specialization
template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, bool scalarDiffusion = true >
class ModelAffineParametricThermalblock;


template< class DomainFieldImp, int domainDim, class RangeFieldImp >
class ModelAffineParametricThermalblock< DomainFieldImp, domainDim, RangeFieldImp, 1, true >
  : public ModelAffineParametricDefault< DomainFieldImp, domainDim, RangeFieldImp, 1, true >
{
  typedef ModelAffineParametricDefault< DomainFieldImp, domainDim, RangeFieldImp, 1, true >       BaseType;
public:
  typedef ModelAffineParametricThermalblock< DomainFieldImp, domainDim, RangeFieldImp, 1, true >  ThisType;

  typedef typename BaseType::DomainFieldType  DomainFieldType;
  static const int                            dimDomain = BaseType::dimDomain;
  typedef typename BaseType::RangeFieldType   RangeFieldType;
  static const int                            dimRange = BaseType::dimRange;
  typedef typename BaseType::ParamFieldType   ParamFieldType;
  static const int                            maxParamDim = BaseType::maxParamDim;

  typedef typename Dune::Stuff::AffineParametricFunctionCheckerboard< DomainFieldType, dimDomain,
                                                                      RangeFieldType, dimRange >
                                          CheckerboardFunctionType;
  typedef typename BaseType::ForceType      ForceType;
  typedef typename BaseType::DirichletType  DirichletType;
  typedef typename BaseType::NeumannType    NeumannType;

  static const std::string id()
  {
    return BaseType::BaseType::id() + ".affineparametric.thermalblock";
  }

  ModelAffineParametricThermalblock(const Dune::shared_ptr< const CheckerboardFunctionType > _diffusion,
                                    const Dune::shared_ptr< const ForceType > _force,
                                    const Dune::shared_ptr< const DirichletType > _dirichlet,
                                    const Dune::shared_ptr< const NeumannType > _neumann)
    : BaseType(_diffusion, _force, _dirichlet, _neumann)
  {}

  static Dune::ParameterTree defaultSettings(const std::string subName = "")
  {
    Dune::Stuff::Common::ExtendedParameterTree description;
    description.add(Dune::Stuff::AffineParametricFunctions< DomainFieldType, dimDomain,
                                                            RangeFieldType, dimRange >::defaultSettings("function.affineparametric.checkerboard"),
                    "diffusion");
    description["diffusion.name"] = "diffusion";
    description.add(Dune::Stuff::Functions< DomainFieldType, dimDomain,
                                            RangeFieldType, dimRange >::defaultSettings("function.expression"),
                    "force");
    description["force.name"] = "force";
    description["force.expression"] = "1.0";
    description["force.order"] = "0";
    description.add(Dune::Stuff::Functions< DomainFieldType, dimDomain,
                                            RangeFieldType, dimRange >::defaultSettings("function.expression"),
                    "neumann");
    description["neumann.name"] = "neumann";
    description["neumann.expression"] = "0.1";
    description["neumann.order"] = "0";
    description.add(Dune::Stuff::Functions< DomainFieldType, dimDomain,
                                            RangeFieldType, dimRange >::defaultSettings("function.expression"),
                    "dirichlet");
    description["dirichlet.name"] = "dirichlet";
    description["dirichlet.expression"] = "0.1*x[0]";
    description["dirichlet.order"] = "1";
    if (subName.empty())
      return description;
    else {
      Dune::Stuff::Common::ExtendedParameterTree extendedDescription;
      extendedDescription.add(description, subName);
      return extendedDescription;
    }
  } // ... defaultSettings(...)

  static ThisType* create(const Dune::ParameterTree& _settings, const std::string _subName = id())
  {
    // get correct settings
    typedef Stuff::Common::ExtendedParameterTree SettingsType;
    SettingsType settings;
    if (_settings.hasSub(_subName))
      settings = _settings.sub(_subName);
    else
      settings = _settings;
    // create the checkerboard diffusion
    SettingsType diffusionSettings = settings.sub("diffusion");
    diffusionSettings["diffusion.name"] = "checkerboard diffusion";
    const shared_ptr< const CheckerboardFunctionType >
        diffusion(CheckerboardFunctionType::create(diffusionSettings));
    // create the rest of the functions
    //   * therefore create a fake diffusion subSettings,
    settings.sub("diffusion") = CheckerboardFunctionType::defaultSettings();
    settings["diffusion.type"] = "function.affineparametric.checkerboard";
    settings["diffusion.name"] = "checkerboard diffusion";
    //   * create a default model
    const std::shared_ptr< const BaseType > baseModel(BaseType::create(settings));
    // create and return
    return new ThisType(diffusion,
                        baseModel->force(),
                        baseModel->dirichlet(),
                        baseModel->neumann());
  } // static ThisType create(...)
}; // class ModelAffineParametricThermalblock


} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_MODEL_AFFINEPARAMETRIC_THERMALBLOCK_HH
