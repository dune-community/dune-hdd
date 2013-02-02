#ifndef DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_HH
#define DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_HH

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parametertree.hh>

#include <dune/stuff/common/color.hh>

#include "model/interface.hh"
#include "model/nonparametric/default.hh"
#include "model/nonparametric/thermalblock.hh"

namespace Dune {
namespace Detailed {
namespace Solvers {
namespace Stationary {
namespace Linear {
namespace Elliptic {
namespace Model {


template< class DomainFieldType, int dimDomain,
          class RangeFieldType, int dimRange,
          class ParamFieldType = double, int maxParams = 0 >
Dune::shared_ptr< Interface<  DomainFieldType, dimDomain,
                              RangeFieldType, dimRange,
                              ParamFieldType, maxParams > >
  create(const std::string type, const Dune::ParameterTree paramTree = Dune::ParameterTree())
{
  // choose model
  if (type == "model.stationary.linear.elliptic.nonparametric.default") {
    typedef Dune::Detailed::Solvers
        ::Stationary
        ::Linear
        ::Elliptic
        ::Model
        ::NonparametricDefault< DomainFieldType, dimDomain, RangeFieldType, dimRange, ParamFieldType, maxParams >
      ModelType;
    return Dune::make_shared< ModelType >(ModelType::createFromParamTree(paramTree));
  } else if (type == "model.stationary.linear.elliptic.nonparametric.thermalblock") {
    typedef Dune::Detailed::Solvers
        ::Stationary
        ::Linear
        ::Elliptic
        ::Model
        ::NonparametricThermalblock< DomainFieldType, dimDomain, RangeFieldType, dimRange >
      ThermalblockModelType;
    return Dune::make_shared< ThermalblockModelType >(ThermalblockModelType::createFromParamTree(paramTree));
  } else
    DUNE_THROW(Dune::RangeError,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " unknown model '" << type << "' requested!");
} // ... create(...)


} // namespace Model
} // namespace Elliptic
} // namespace Linear
} // namespace Stationary
} // namespace Solvers
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_HH
