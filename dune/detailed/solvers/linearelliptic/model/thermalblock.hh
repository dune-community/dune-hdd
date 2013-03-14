#ifndef DUNE_DETAILED_SOLVERS_LINEARELLIPTIC_MODEL_THERMALBLOCK_HH
#define DUNE_DETAILED_SOLVERS_LINEARELLIPTIC_MODEL_THERMALBLOCK_HH

#include <vector>
#include <string>
#include <memory>

#include <dune/stuff/function/checkerboard.hh>
#include <dune/stuff/common/parameter/tree.hh>

#include "default.hh"

namespace Dune {
namespace Detailed {
namespace Solvers {
namespace LinearElliptic {


template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim >
class ModelThermalblock;


template< class DomainFieldImp, int domainDim, class RangeFieldImp >
class ModelThermalblock< DomainFieldImp, domainDim, RangeFieldImp, 1 >
  : public ModelDefault< DomainFieldImp, domainDim, RangeFieldImp, 1 >
{
public:
  typedef ModelThermalblock< DomainFieldImp, domainDim, RangeFieldImp, 1 > ThisType;
  typedef ModelDefault< DomainFieldImp, domainDim, RangeFieldImp, 1 >      BaseType;

  typedef typename BaseType::DomainFieldType  DomainFieldType;
  static const int                            dimDomain = BaseType::dimDomain;
  typedef typename BaseType::RangeFieldType   RangeFieldType;
  static const int                            dimRange = BaseType::dimRange;

  typedef typename BaseType::FunctionType   FunctionType;
  typedef typename Stuff::FunctionCheckerboard< DomainFieldType, dimDomain, RangeFieldType, dimRange >
                                            CheckerboardFunctionType;
  typedef typename FunctionType::DomainType DomainType;

  static const std::string id()
  {
    return BaseType::BaseType::id() + ".thermalblock";
  }

  ModelThermalblock(const DomainType& _lowerLeft,
                    const DomainType& _upperRight,
                    const std::vector< size_t >& _numElements,
                    const std::vector< RangeFieldType >& _components,
                    const std::shared_ptr< const FunctionType > _force,
                    const std::shared_ptr< const FunctionType > _dirichlet,
                    const std::shared_ptr< const FunctionType > _neumann)
    : BaseType(std::make_shared< CheckerboardFunctionType >(_lowerLeft,
                                                            _upperRight,
                                                            _numElements,
                                                            _components),
               _force,
               _dirichlet,
               _neumann)
  {}

  ModelThermalblock(const std::shared_ptr< const CheckerboardFunctionType > _diffusion,
                    const std::shared_ptr< const FunctionType > _force,
                    const std::shared_ptr< const FunctionType > _dirichlet,
                    const std::shared_ptr< const FunctionType > _neumann)
    : BaseType(_diffusion, _force, _dirichlet, _neumann)
  {}

  ModelThermalblock(const ThisType& _other)
    :BaseType(_other.diffusion(),
              _other.force(),
              _other.dirichlet(),
              _other.neumann())
  {}

  ThisType& operator=(const ThisType& _other)
  {
    if (this != &_other) {
      BaseType::operator=(_other);
    }
    return this;
  } // ThisType& operator=(const ThisType& other)

  static Dune::ParameterTree createSampleDescription(const std::string subName = "")
  {
    Dune::Stuff::Common::ExtendedParameterTree description;
    description.add(Dune::Stuff::Functions< DomainFieldType, dimDomain,
                                            RangeFieldType, dimRange >::createSampleDescription("function.checkerboard"),
                    "diffusion");
    description["diffusion.name"] = "diffusion";
    description.add(Dune::Stuff::Functions< DomainFieldType, dimDomain,
                                            RangeFieldType, dimRange >::createSampleDescription("function.expression"),
                    "force");
    description["force.name"] = "force";
    description["force.expression"] = "1.0";
    description["force.order"] = "0";
    description.add(Dune::Stuff::Functions< DomainFieldType, dimDomain,
                                            RangeFieldType, dimRange >::createSampleDescription("function.expression"),
                    "neumann");
    description["neumann.name"] = "neumann";
    description["neumann.expression"] = "0.1";
    description["neumann.order"] = "0";
    description.add(Dune::Stuff::Functions< DomainFieldType, dimDomain,
                                            RangeFieldType, dimRange >::createSampleDescription("function.expression"),
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

  static ThisType* create(const Dune::ParameterTree& _description, const std::string _subName = id())
  {
    // get correct paramTree
    Dune::Stuff::Common::ExtendedParameterTree description;
    if (_description.hasSub(_subName))
      description = _description.sub(_subName);
    else
      description = _description;
    // create the correct diffusion
    const Dune::Stuff::Common::ExtendedParameterTree& diffusionDescription = description.sub("diffusion");
    const std::shared_ptr< const CheckerboardFunctionType >
        diffusion(CheckerboardFunctionType::create(diffusionDescription));
    return new ThisType(diffusion,
                        createFunction("force", description),
                        createFunction("dirichlet", description),
                        createFunction("neumann", description));
  } // static ThisType createFromParamTree(const Dune::ParameterTree& paramTree, const std::string subName = id())

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
  }
}; // class ModelThermalblock


} // namespace LinearElliptic
} // namespace Solvers
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_LINEARELLIPTIC_MODEL_THERMALBLOCK_HH
