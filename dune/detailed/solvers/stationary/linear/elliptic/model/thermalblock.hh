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
    // create the corrent diffusion
    const Dune::Stuff::Common::ExtendedParameterTree& diffusionDescription = description.sub("diffusion");
    std::vector< DomainFieldType > lowerLefts = diffusionDescription.getVector("lowerLeft", DomainFieldType(0), dimDomain);
    std::vector< DomainFieldType > upperRights = diffusionDescription.getVector("upperRight", DomainFieldType(1), dimDomain);
    std::vector< size_t > numElements = diffusionDescription.getVector("numElements", size_t(1), dimDomain);
    std::vector< RangeFieldType > components = diffusionDescription.getVector("components", RangeFieldType(1), dimDomain);
    assert(int(lowerLefts.size()) >= dimDomain && "Given vector is too short!");
    assert(int(upperRights.size()) >= dimDomain && "Given vector is too short!");
    assert(int(numElements.size()) >= dimDomain && "Given vector is too short!");
    assert(int(components.size()) >= dimDomain && "Given vector is too short!");
    Dune::FieldVector< DomainFieldType, dimDomain > lowerLeft;
    Dune::FieldVector< DomainFieldType, dimDomain > upperRight;
    size_t numSubdomains = 1;
    for (int d = 0; d < dimDomain; ++d) {
      lowerLeft[d] = lowerLefts[d];
      upperRight[d] = upperRights[d];
      numSubdomains *= numElements[d];
    }
    Dune::shared_ptr< CheckerboardFunctionType > _diffusion(new CheckerboardFunctionType(lowerLeft,
                                                                                         upperRight,
                                                                                         numElements,
                                                                                         components));
    // use method from base to create the rest, therefore create a copy of the paramtree to fake the diffusion
    Dune::Stuff::Common::ExtendedParameterTree fakeDescription(description);
    fakeDescription["diffusion.variable"] = "x";
    fakeDescription["diffusion.expression"] = "[1.0; 1.0; 1.0]";
    fakeDescription["diffusion.order"] = "0";
    fakeDescription["diffusion.name"] = "diffusion";
    const BaseType base = BaseType::createFromDescription(fakeDescription);
    // create and return
    return ThisType(_diffusion, base.force(), base.dirichlet(), base.neumann());
  } // static ThisType createFromParamTree(const Dune::ParameterTree& paramTree, const std::string subName = id())
}; // class Thermalblock


} // namespace Model
} // namespace Elliptic
} // namespace Linear
} // namespace Stationary
} // namespace Solvers
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_NONPARAMETRIC_THERMALBLOCK_HH
