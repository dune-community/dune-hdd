#ifndef DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_HH
#define DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_HH

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#include "model/interface.hh"
#include "model/default.hh"
#include "model/thermalblock.hh"

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

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_HH
