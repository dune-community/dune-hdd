#ifndef DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_THERMALBLOCK_HH
#define DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_THERMALBLOCK_HH

#include <vector>
#include <string>

#include <dune/common/shared_ptr.hh>

#include <dune/stuff/function/expression.hh>
#include <dune/stuff/function/interface.hh>
#include <dune/stuff/function/checkerboard.hh>
#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/print.hh>

#include "default.hh"

namespace Dune {
namespace Detailed {
namespace Solvers {
namespace Stationary {
namespace Linear {
namespace Elliptic {
namespace Model {

template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim >
class Thermalblock
  : public Default < DomainFieldImp, domainDim, RangeFieldImp, rangeDim >
{
public:
  typedef DomainFieldImp DomainFieldType;

  static const int dimDomain = domainDim;

  typedef RangeFieldImp RangeFieldType;

  static const int dimRange = rangeDim;

  typedef Default< DomainFieldType, dimDomain, RangeFieldType, dimRange > BaseType;

  typedef Thermalblock< DomainFieldType, dimDomain, RangeFieldType, dimRange > ThisType;

  typedef typename Dune::Stuff::Function::Checkerboard< DomainFieldType, dimDomain,
                                                        RangeFieldType, dimRange >
    CheckerboardFunctionType;

  typedef typename BaseType::DiffusionType DiffusionType;

  typedef typename DiffusionType::DomainType DomainType;

  typedef typename BaseType::ForceType ForceType;

  typedef typename BaseType::DirichletType DirichletType;

  typedef typename BaseType::NeumannType NeumannType;

  static const std::string id()
  {
    return BaseType::BaseType::id() + ".thermalblock";
  }

  Thermalblock(const DomainType& _lowerLeft,
               const DomainType& _upperRight,
               const std::vector< unsigned int > _numElements,
               const std::vector< RangeFieldType > _components,
               const Dune::shared_ptr< const ForceType > _force,
               const Dune::shared_ptr< const DirichletType > _dirichlet,
               const Dune::shared_ptr< const NeumannType > _neumann,
               const int _diffusionOrder = 0,
               const int _forceOrder = -1,
               const int _dirichletOrder = -1,
               const int _neumannOrder = -1)
    : BaseType(Dune::shared_ptr< const CheckerboardFunctionType >(new CheckerboardFunctionType(_lowerLeft,
                                                                                               _upperRight,
                                                                                               _numElements,
                                                                                               _components)),
               _force, _dirichlet, _neumann, _diffusionOrder, _forceOrder, _dirichletOrder, _neumannOrder)
  {}

  Thermalblock(const Dune::shared_ptr< const CheckerboardFunctionType > _diffusion,
               const Dune::shared_ptr< const ForceType > _force,
               const Dune::shared_ptr< const DirichletType > _dirichlet,
               const Dune::shared_ptr< const NeumannType > _neumann,
               const int _diffusionOrder = 0,
               const int _forceOrder = -1,
               const int _dirichletOrder = -1,
               const int _neumannOrder = -1)
    : BaseType(_diffusion, _force, _dirichlet, _neumann, _diffusionOrder, _forceOrder, _dirichletOrder, _neumannOrder)
  {}

  Thermalblock(const ThisType& other)
    :BaseType(other.diffusion(),
              other.force(),
              other.dirichlet(),
              other.neumann(),
              other.diffusionOrder(),
              other.forceOrder(),
              other.dirichletOrder(),
              other.neumannOrder())
  {}

  static ThisType createFromParamTree(const Dune::ParameterTree& paramTree, const std::string subName = id())
  {
    // get correct paramTree
    Dune::Stuff::Common::ExtendedParameterTree extendedParamTree;
    if (paramTree.hasSub(subName))
      extendedParamTree = paramTree.sub(subName);
    else
      extendedParamTree = paramTree;
    // create the corrent diffusion
    if (!extendedParamTree.hasSub("diffusion"))
      DUNE_THROW(Dune::RangeError,
                 "\nError: sub 'diffusion' not found in the following Dune::ParameterTree\n:" << extendedParamTree.reportString("  "));
    const Dune::Stuff::Common::ExtendedParameterTree& diffusionTree = extendedParamTree.sub("diffusion");
    std::vector< DomainFieldType > lowerLefts;
    std::vector< DomainFieldType > upperRights;
    std::vector< unsigned int > numElements;
    std::vector< RangeFieldType > components;
    if (!diffusionTree.hasKey("lowerLeft"))
      std::cout << "Warning in " << id() << ": neither vector nor key 'lowerLeft' given, defaulting to 0.0!" << std::endl;
    if (!diffusionTree.hasKey("upperRight"))
      std::cout << "Warning in " << id() << ": neither vector nor key 'upperRight' given, defaulting to 1.0!" << std::endl;
    if (!diffusionTree.hasKey("numElements"))
      std::cout << "Warning in " << id() << ": neither vector nor key 'numElements' given, defaulting to 1!" << std::endl;
    if (!diffusionTree.hasKey("components"))
      std::cout << "Warning in " << id() << ": neither vector nor key 'components' given, defaulting to 1!" << std::endl;
    if (!diffusionTree.hasKey("order"))
      std::cout << "Warning in " << id() << ": no key 'order' given, defaulting to 0!" << std::endl;
    const unsigned int _diffusionOrder = diffusionTree.get("order", 0u);
    lowerLefts = diffusionTree.getVector("lowerLeft", DomainFieldType(0), dimDomain);
    upperRights = diffusionTree.getVector("upperRight", DomainFieldType(1), dimDomain);
    numElements = diffusionTree.getVector("numElements", 1u, dimDomain);
    components = diffusionTree.getVector("components", RangeFieldType(1), dimDomain);
    assert(int(lowerLefts.size()) >= dimDomain && "Given vector is too short!");
    assert(int(upperRights.size()) >= dimDomain && "Given vector is too short!");
    assert(int(numElements.size()) >= dimDomain && "Given vector is too short!");
    assert(int(components.size()) >= dimDomain && "Given vector is too short!");
    Dune::FieldVector< DomainFieldType, dimDomain > lowerLeft;
    Dune::FieldVector< DomainFieldType, dimDomain > upperRight;
    unsigned int numSubdomains = 1;
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
    Dune::Stuff::Common::ExtendedParameterTree fakeParamTree(extendedParamTree);
    fakeParamTree["diffusion.variable"] = "x";
    fakeParamTree["diffusion.expression"] = "[1.0; 1.0; 1.0]";
    fakeParamTree["diffusion.order"] = "0";
    const BaseType base = BaseType::createFromParamTree(fakeParamTree);
    // create and return
    return ThisType(_diffusion, base.force(), base.dirichlet(), base.neumann(),
                    _diffusionOrder, base.forceOrder(), base.dirichletOrder(), base.neumannOrder());
  } // static ThisType createFromParamTree(const Dune::ParameterTree& paramTree, const std::string subName = id())

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      BaseType::operator=(other);
    }
    return this;
  } // ThisType& operator=(const ThisType& other)

private:
  Dune::shared_ptr< const DiffusionType > diffusion_;
  Dune::shared_ptr< const ForceType > force_;
  Dune::shared_ptr< const DirichletType > dirichlet_;
  Dune::shared_ptr< const NeumannType > neumann_;
  int forceOrder_;
  int dirichletOrder_;
  int neumannOrder_;
}; // class Thermalblock

} // namespace Model

} // namespace Elliptic

} // namespace Linear

} // namespace Stationary

} // namespace Solvers

} // namespace Detailed

} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_THERMALBLOCK_HH
