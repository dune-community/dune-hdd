#ifndef DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_INTERFACE_HH
#define DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_INTERFACE_HH

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#include <dune/common/shared_ptr.hh>

#include <dune/stuff/function/interface.hh>

namespace Dune {
namespace Detailed {
namespace Solvers {
namespace Stationary {
namespace Linear {
namespace Elliptic {
namespace Model {

template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim >
class Interface
{
public:
  typedef DomainFieldImp DomainFieldType;

  static const int dimDomain = domainDim;

  typedef RangeFieldImp RangeFieldType;

  static const int dimRange = rangeDim;

  typedef Dune::Stuff::Function::Interface< DomainFieldType, dimDomain, RangeFieldType, dimRange > DiffusionType;

  typedef Dune::Stuff::Function::Interface< DomainFieldType, dimDomain, RangeFieldType, dimRange > ForceType;

  typedef Dune::Stuff::Function::Interface< DomainFieldType, dimDomain, RangeFieldType, dimRange > DirichletType;

  typedef Dune::Stuff::Function::Interface< DomainFieldType, dimDomain, RangeFieldType, dimRange > NeumannType;

  static const std::string id()
  {
    return "detailed.solvers.stationary.linear.elliptic.model";
  }

  virtual const Dune::shared_ptr< const DiffusionType > diffusion() const = 0;

  virtual int diffusionOrder() const = 0;

  virtual const Dune::shared_ptr< const ForceType > force() const = 0;

  virtual int forceOrder() const = 0;

  virtual const Dune::shared_ptr< const DirichletType > dirichlet() const = 0;

  virtual int dirichletOrder() const = 0;

  virtual const Dune::shared_ptr< const NeumannType > neumann() const = 0;

  virtual int neumannOrder() const = 0;
}; // class Interface

} // namespace Model
} // namespace Elliptic
} // namespace Linear
} // namespace Stationary
} // namespace Solvers
} // namespace Detailed
} // namespace Dune

#include "default.hh"
#include "thermalblock.hh"

namespace Dune {
namespace Detailed {
namespace Solvers {
namespace Stationary {
namespace Linear {
namespace Elliptic {
namespace Model {

template< class DomainFieldType, int dimDomain, class RangeFieldType, int dimRange >
Interface< DomainFieldType, dimDomain, RangeFieldType, dimRange >* create(const std::string type, const Dune::ParameterTree paramTree = Dune::ParameterTree())
{
  // choose model
  if (type == "detailed.solvers.stationary.linear.elliptic.model.default") {
    typedef Dune::Detailed::Solvers
        ::Stationary
        ::Linear
        ::Elliptic
        ::Model
        ::Default< DomainFieldType, dimDomain, RangeFieldType, dimRange >
      DefaultModelType;
    return new DefaultModelType(DefaultModelType::createFromParamTree(paramTree));
  } else if (type == "detailed.solvers.stationary.linear.elliptic.model.thermalblock") {
    typedef Dune::Detailed::Solvers
        ::Stationary
        ::Linear
        ::Elliptic
        ::Model
        ::Thermalblock< DomainFieldType, dimDomain, RangeFieldType, dimRange >
      ThermalblockModelType;
    return new ThermalblockModelType(ThermalblockModelType::createFromParamTree(paramTree));
  } else
    DUNE_THROW(Dune::RangeError,
               "\nError: unknown model '" << type << "'given!");
} // Interface* create(const std::string type, const Dune::ParameterTree paramTree = Dune::ParameterTree())

} // namespace Model
} // namespace Elliptic
} // namespace Linear
} // namespace Stationary
} // namespace Solvers
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_INTERFACE_HH
