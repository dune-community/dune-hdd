#ifndef DUNE_DETAILED_SOLVERS_LINEARELLIPTIC_MODEL_THERMALBLOCK_HH
#define DUNE_DETAILED_SOLVERS_LINEARELLIPTIC_MODEL_THERMALBLOCK_HH

#include <vector>
#include <string>
#include <memory>

#include <dune/stuff/function/checkerboard.hh>
#include <dune/stuff/common/parameter/tree.hh>

#include "default.hh"

namespace Dune {
namespace DetailedSolvers {
namespace LinearElliptic {


template< class DomainFieldImp, int domainDim, class RangeFieldImp >
class ModelThermalblock< DomainFieldImp, domainDim, RangeFieldImp, 1, true >
  : public ModelDefault< DomainFieldImp, domainDim, RangeFieldImp, 1, true >
{
  typedef ModelDefault< DomainFieldImp, domainDim, RangeFieldImp, 1, true >      BaseType;
public:
  typedef ModelThermalblock< DomainFieldImp, domainDim, RangeFieldImp, 1, true > ThisType;

  typedef typename BaseType::DomainFieldType  DomainFieldType;
  static const int                            dimDomain = BaseType::dimDomain;
  typedef typename BaseType::RangeFieldType   RangeFieldType;
  static const int                            dimRange = BaseType::dimRange;

  typedef typename BaseType::ForceType      ForceType;
  typedef typename BaseType::DirichletType  DirichletType;
  typedef typename BaseType::NeumannType    NeumannType;
  typedef typename Stuff::FunctionCheckerboard< DomainFieldType, dimDomain, RangeFieldType, dimRange >
                                            CheckerboardFunctionType;

  typedef typename CheckerboardFunctionType::DomainType   DomainType;


  static const std::string id()
  {
    return BaseType::BaseType::id() + ".thermalblock";
  }

  ModelThermalblock(const DomainType& _lowerLeft,
                    const DomainType& _upperRight,
                    const std::vector< size_t >& _numElements,
                    const std::vector< RangeFieldType >& _components,
                    const std::shared_ptr< const ForceType > _force,
                    const std::shared_ptr< const DirichletType > _dirichlet,
                    const std::shared_ptr< const NeumannType > _neumann)
    : BaseType(std::make_shared< CheckerboardFunctionType >(_lowerLeft,
                                                            _upperRight,
                                                            _numElements,
                                                            _components),
               _force,
               _dirichlet,
               _neumann)
  {}

  ModelThermalblock(const std::shared_ptr< const CheckerboardFunctionType > _diffusion,
                    const std::shared_ptr< const ForceType > _force,
                    const std::shared_ptr< const DirichletType > _dirichlet,
                    const std::shared_ptr< const NeumannType > _neumann)
    : BaseType(_diffusion, _force, _dirichlet, _neumann)
  {}

  static Dune::ParameterTree defaultSettings(const std::string subName = "")
  {
    Dune::Stuff::Common::ExtendedParameterTree description;
    description.add(Dune::Stuff::Functions< DomainFieldType, dimDomain,
                                            RangeFieldType, dimRange >::defaultSettings("function.checkerboard"),
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
  }

  static ThisType* create(const Dune::ParameterTree& _settings, const std::string _subName = id())
  {
    // get correct paramTree
    Dune::Stuff::Common::ExtendedParameterTree settings;
    if (_settings.hasSub(_subName))
      settings = _settings.sub(_subName);
    else
      settings = _settings;
    // create the correct diffusion
    return new ThisType(std::shared_ptr< CheckerboardFunctionType >(CheckerboardFunctionType::create(settings.sub("diffusion"))),
                        BaseType::template createFunction< ForceType >("force", settings),
                        BaseType::template createFunction< DirichletType >("dirichlet", settings),
                        BaseType::template createFunction< NeumannType >("neumann", settings));
  } // static ThisType createFromParamTree(const Dune::ParameterTree& paramTree, const std::string subName = id())

private:

}; // class ModelThermalblock


} // namespace LinearElliptic
} // namespace DetailedSolvers
} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_LINEARELLIPTIC_MODEL_THERMALBLOCK_HH
