#ifndef DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_THERMALBLOCK_HH
#define DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_THERMALBLOCK_HH

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#include <vector>
#include <string>

#include <dune/common/shared_ptr.hh>

#include <dune/stuff/function/nonparametric/checkerboard.hh>

#include "default.hh"

namespace Dune {
namespace Detailed {
namespace Solvers {
namespace Stationary {
namespace Linear {
namespace Elliptic {
namespace Model {

template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim >
class NonparametricThermalblock
  : public NonparametricDefault< DomainFieldImp, domainDim, RangeFieldImp, rangeDim >
{
public:
  typedef NonparametricDefault< DomainFieldImp, domainDim, RangeFieldImp, rangeDim >      BaseType;
  typedef NonparametricThermalblock< DomainFieldImp, domainDim, RangeFieldImp, rangeDim > ThisType;

  typedef typename BaseType::DomainFieldType  DomainFieldType;
  static const int                            dimDomain = BaseType::dimDomain;
  typedef typename BaseType::RangeFieldType   RangeFieldType;
  static const int                            dimRange = BaseType::dimRange;

  typedef typename BaseType::FunctionType   FunctionType;
  typedef typename Dune::Stuff::Function::NonparametricCheckerboard<  DomainFieldType, dimDomain,
                                                                      RangeFieldType, dimRange >
                                            CheckerboardFunctionType;
  typedef typename FunctionType::DomainType DomainType;

  static const std::string id()
  {
    return BaseType::BaseType::id() + ".nonparametric.thermalblock";
  }

  NonparametricThermalblock(const DomainType& _lowerLeft,
                            const DomainType& _upperRight,
                            const std::vector< unsigned int > _numElements,
                            const std::vector< RangeFieldType > _components,
                            const Dune::shared_ptr< const FunctionType > _force,
                            const Dune::shared_ptr< const FunctionType > _dirichlet,
                            const Dune::shared_ptr< const FunctionType > _neumann)
    : BaseType(Dune::shared_ptr< const CheckerboardFunctionType >(new CheckerboardFunctionType(_lowerLeft,
                                                                                               _upperRight,
                                                                                               _numElements,
                                                                                               _components)),
               _force, _dirichlet, _neumann)
  {}

  NonparametricThermalblock(const Dune::shared_ptr< const CheckerboardFunctionType > _diffusion,
                            const Dune::shared_ptr< const FunctionType > _force,
                            const Dune::shared_ptr< const FunctionType > _dirichlet,
                            const Dune::shared_ptr< const FunctionType > _neumann)
    : BaseType(_diffusion, _force, _dirichlet, _neumann)
  {}

  NonparametricThermalblock(const ThisType& _other)
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

  static ThisType createFromParamTree(const Dune::ParameterTree& _paramTree, const std::string _subName = id())
  {
    // get correct paramTree
    Dune::Stuff::Common::ExtendedParameterTree paramTree;
    if (_paramTree.hasSub(_subName))
      paramTree = _paramTree.sub(_subName);
    else
      paramTree = _paramTree;
    // create the corrent diffusion
    const Dune::Stuff::Common::ExtendedParameterTree& diffusionTree = paramTree.sub("diffusion");
    std::vector< DomainFieldType > lowerLefts = diffusionTree.getVector("lowerLeft", DomainFieldType(0), dimDomain);
    std::vector< DomainFieldType > upperRights = diffusionTree.getVector("upperRight", DomainFieldType(1), dimDomain);
    std::vector< unsigned int > numElements = diffusionTree.getVector("numElements", 1u, dimDomain);
    std::vector< RangeFieldType > components = diffusionTree.getVector("components", RangeFieldType(1), dimDomain);
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
    Dune::Stuff::Common::ExtendedParameterTree fakeParamTree(paramTree);
    fakeParamTree["diffusion.variable"] = "x";
    fakeParamTree["diffusion.expression"] = "[1.0; 1.0; 1.0]";
    fakeParamTree["diffusion.order"] = "0";
    fakeParamTree["diffusion.name"] = "diffusion";
    const BaseType base = BaseType::createFromParamTree(fakeParamTree);
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

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_THERMALBLOCK_HH
