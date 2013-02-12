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
#include "model/default.hh"
#include "model/thermalblock.hh"
#include "model/parametric/separable/default.hh"
#include "model/parametric/separable/thermalblock.hh"

namespace Dune {
namespace Detailed {
namespace Solvers {
namespace Stationary {
namespace Linear {
namespace Elliptic {
namespace Model {


std::vector< std::string > types()
{
  std::vector< std::string > ret;
  ret.push_back("model.stationary.linear.elliptic.default");
  ret.push_back("model.stationary.linear.elliptic.thermalblock");
  ret.push_back("model.stationary.linear.elliptic.parametric.separable.default");
  ret.push_back("model.stationary.linear.elliptic.parametric.separable.thermalblock");
  return ret;
}


template< class DomainFieldType, int dimDomain, class RangeFieldType, int dimRange >
Dune::shared_ptr< Interface< DomainFieldType, dimDomain, RangeFieldType, dimRange > >
  create(const std::string type, const Dune::ParameterTree description = Dune::ParameterTree())
{
  // choose model
  if (type == "model.stationary.linear.elliptic.default") {
    typedef Stationary::Linear::Elliptic::Model::Default< DomainFieldType, dimDomain, RangeFieldType, dimRange >
      ModelType;
    return Dune::make_shared< ModelType >(ModelType::createFromDescription(description));
  } else if (type == "model.stationary.linear.elliptic.thermalblock") {
    typedef Stationary::Linear::Elliptic::Model::Thermalblock<  DomainFieldType, dimDomain, RangeFieldType, dimRange >
      ModelType;
    return Dune::make_shared< ModelType >(ModelType::createFromDescription(description));
  } else if (type == "model.stationary.linear.elliptic.parametric.separable.default") {
    typedef Stationary::Linear::Elliptic::Model::SeparableDefault<  DomainFieldType, dimDomain,
                                                                    RangeFieldType, dimRange >
      ModelType;
    return Dune::make_shared< ModelType >(ModelType::createFromDescription(description));
  } else if (type == "model.stationary.linear.elliptic.parametric.separable.thermalblock") {
    typedef Stationary::Linear::Elliptic::Model::SeparableThermalblock< DomainFieldType, dimDomain,
                                                                        RangeFieldType, dimRange >
      ModelType;
    return Dune::make_shared< ModelType >(ModelType::createFromDescription(description));
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
