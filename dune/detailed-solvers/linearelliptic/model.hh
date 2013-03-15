#ifndef DUNE_DETAILED_SOLVERS_LINEARELLIPTIC_MODEL_HH
#define DUNE_DETAILED_SOLVERS_LINEARELLIPTIC_MODEL_HH

#include <dune/common/exceptions.hh>
#include <dune/common/parametertree.hh>

#include <dune/stuff/common/color.hh>

#include "model/interface.hh"
#include "model/default.hh"
#include "model/thermalblock.hh"
#include "model/affineparametric/default.hh"
#include "model/affineparametric/twophase.hh"
#include "model/affineparametric/thermalblock.hh"

namespace Dune {
namespace DetailedSolvers {
namespace LinearElliptic {


template< class DomainFieldType, int dimDomain, class RangeFieldType, int dimRange >
class Models
{
public:
  static std::vector< std::string > available()
  {
    return {
        "model.linearelliptic.default"
        , "model.linearelliptic.thermalblock"
        , "model.linearelliptic.affineparametric.default"
        , "model.linearelliptic.affineparametric.twophase"
        , "model.linearelliptic.affineparametric.thermalblock"
    };
  } // ... available()

  static Dune::ParameterTree createSampleDescription(const std::string type, const std::string subname = "")
  {
    if (type == "model.linearelliptic.default")
      return ModelDefault< DomainFieldType, dimDomain, RangeFieldType, dimRange >::createSampleDescription(subname);
    else if (type == "model.linearelliptic.thermalblock")
      return ModelThermalblock< DomainFieldType, dimDomain, RangeFieldType, dimRange >::createSampleDescription(subname);
    else if (type == "model.linearelliptic.affineparametric.default")
      return ModelAffineParametricDefault<  DomainFieldType, dimDomain,
                                            RangeFieldType, dimRange >::createSampleDescription(subname);
    else if (type == "model.linearelliptic.affineparametric.twophase")
      return ModelAffineParametricTwoPhase< DomainFieldType, dimDomain,
                                            RangeFieldType, dimRange >::createSampleDescription(subname);
    else if (type == "model.linearelliptic.affineparametric.thermalblock")
      return ModelAffineParametricThermalblock< DomainFieldType, dimDomain,
                                                RangeFieldType, dimRange >::createSampleDescription(subname);
    else
      DUNE_THROW(Dune::RangeError,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " unknown model '" << type << "' requested!");
  } // ... createSampleDescription(...)

  static ModelInterface< DomainFieldType, dimDomain, RangeFieldType, dimRange >*
  create(const std::string type = available()[0],
           const Dune::ParameterTree description = Dune::ParameterTree())
  {
    if (type == "model.linearelliptic.default")
      return ModelDefault< DomainFieldType, dimDomain, RangeFieldType, dimRange >::create(description);
    else if (type == "model.linearelliptic.thermalblock")
      return ModelThermalblock< DomainFieldType, dimDomain, RangeFieldType, dimRange >::create(description);
    else if (type == "model.linearelliptic.affineparametric.default")
      return ModelAffineParametricDefault<  DomainFieldType, dimDomain,
                                            RangeFieldType, dimRange >::create(description);
    else if (type == "model.linearelliptic.affineparametric.twophase")
      return ModelAffineParametricTwoPhase< DomainFieldType, dimDomain,
                                            RangeFieldType, dimRange >::create(description);
    else if (type == "model.linearelliptic.affineparametric.thermalblock")
      return ModelAffineParametricThermalblock< DomainFieldType, dimDomain,
                                                RangeFieldType, dimRange >::create(description);
    else
      DUNE_THROW(Dune::RangeError,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " unknown model '" << type << "' requested!");
  } // ... create(...)
}; // class Models


} // namespace LinearElliptic
} // namespace DetailedSolvers
} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_LINEARELLIPTIC_MODEL_HH
