#ifndef DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_NONPARAMETRIC_THERMALBLOCK_HH
#define DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_NONPARAMETRIC_THERMALBLOCK_HH

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#include <vector>
#include <string>

#include <dune/common/shared_ptr.hh>

#include <dune/stuff/function/checkerboard.hh>
#include <dune/stuff/common/parameter/tree.hh>

#include "default.hh"

namespace Dune {
namespace Detailed {
namespace Solvers {
namespace Stationary {
namespace Linear {
namespace Elliptic {
namespace Model {


template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim >
class Thermalblock;


template< class DomainFieldImp, int domainDim, class RangeFieldImp >
class Thermalblock< DomainFieldImp, domainDim, RangeFieldImp, 1 >
  : public Default< DomainFieldImp, domainDim, RangeFieldImp, 1 >
{
public:
  typedef Thermalblock< DomainFieldImp, domainDim, RangeFieldImp, 1 > ThisType;
  typedef Default< DomainFieldImp, domainDim, RangeFieldImp, 1 >      BaseType;

  typedef typename BaseType::DomainFieldType  DomainFieldType;
  static const int                            dimDomain = BaseType::dimDomain;
  typedef typename BaseType::RangeFieldType   RangeFieldType;
  static const int                            dimRange = BaseType::dimRange;

  typedef typename BaseType::FunctionType   FunctionType;
  typedef typename Stuff::Function::Checkerboard< DomainFieldType, dimDomain, RangeFieldType, dimRange >
                                            CheckerboardFunctionType;
  typedef typename FunctionType::DomainType DomainType;

  static const std::string id()
  {
    return BaseType::BaseType::id() + ".thermalblock";
  }

  Thermalblock(const DomainType& _lowerLeft,
               const DomainType& _upperRight,
               const std::vector< size_t >& _numElements,
               const std::vector< RangeFieldType >& _components,
               const Dune::shared_ptr< const FunctionType > _force,
               const Dune::shared_ptr< const FunctionType > _dirichlet,
               const Dune::shared_ptr< const FunctionType > _neumann)
    : BaseType(Dune::shared_ptr< const CheckerboardFunctionType >(new CheckerboardFunctionType(_lowerLeft,
                                                                                               _upperRight,
                                                                                               _numElements,
                                                                                               _components)),
               _force,
               _dirichlet,
               _neumann)
  {}

  Thermalblock(const Dune::shared_ptr< const CheckerboardFunctionType > _diffusion,
                            const Dune::shared_ptr< const FunctionType > _force,
                            const Dune::shared_ptr< const FunctionType > _dirichlet,
                            const Dune::shared_ptr< const FunctionType > _neumann)
    : BaseType(_diffusion, _force, _dirichlet, _neumann)
  {}

  Thermalblock(const ThisType& _other)
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
    description.add(Dune::Stuff::Function::createSampleDescription< DomainFieldType, dimDomain,
                                                                    RangeFieldType, dimRange >("function.checkerboard"),
                    "diffusion");
    description["diffusion.name"] = "diffusion";
    description.add(Dune::Stuff::Function::createSampleDescription< DomainFieldType, dimDomain,
                                                                    RangeFieldType, dimRange >("function.expression"),
                    "force");
    description["force.name"] = "force";
    description["force.expression"] = "1.0";
    description["force.order"] = "0";
    description.add(Dune::Stuff::Function::createSampleDescription< DomainFieldType, dimDomain,
                                                                    RangeFieldType, dimRange >("function.expression"),
                    "neumann");
    description["neumann.name"] = "neumann";
    description["neumann.expression"] = "0.1";
    description["neumann.order"] = "0";
    description.add(Dune::Stuff::Function::createSampleDescription< DomainFieldType, dimDomain,
                                                                    RangeFieldType, dimRange >("function.expression"),
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

  static ThisType createFromDescription(const Dune::ParameterTree& _description, const std::string _subName = id())
  {
    // get correct paramTree
    Dune::Stuff::Common::ExtendedParameterTree description;
    if (_description.hasSub(_subName))
      description = _description.sub(_subName);
    else
      description = _description;
    // create the correct diffusion
    const Dune::Stuff::Common::ExtendedParameterTree& diffusionDescription = description.sub("diffusion");
    const Dune::shared_ptr< const CheckerboardFunctionType >
        diffusion(new CheckerboardFunctionType(CheckerboardFunctionType::createFromDescription(diffusionDescription)));
    return ThisType(diffusion,
                    createFunction("force", description),
                    createFunction("dirichlet", description),
                    createFunction("neumann", description));
  } // static ThisType createFromParamTree(const Dune::ParameterTree& paramTree, const std::string subName = id())

private:
  static Dune::shared_ptr< const FunctionType > createFunction(const std::string& _id,
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
    return Dune::shared_ptr< const FunctionType >(Dune::Stuff::Function::create< DomainFieldType, dimDomain, RangeFieldType, dimRange >(type,
                                                                                                 functionDescription));
  }

}; // class Thermalblock


} // namespace Model
} // namespace Elliptic
} // namespace Linear
} // namespace Stationary
} // namespace Solvers
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_NONPARAMETRIC_THERMALBLOCK_HH
