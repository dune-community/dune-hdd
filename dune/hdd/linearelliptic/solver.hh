#ifndef DUNE_HDD_LINEARELLIPTIC_SOLVER_HH
#define DUNE_HDD_LINEARELLIPTIC_SOLVER_HH

#include <vector>
#include <string>

#include <dune/common/parametertree.hh>
#include <dune/common/exceptions.hh>

#include <dune/stuff/common/color.hh>

#include "solver/interface.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {


// forward, to allow for specialization
template< class DomainFieldType, int dimDomain, class RangeFieldType, int dimRange >
class Solvers
{
public:
  Solvers() = delete;
};


template< class DomainFieldType, int dimDomain, class RangeFieldType >
class Solvers< DomainFieldType, dimDomain, RangeFieldType, 1 >
{
public:
  static std::vector< std::string > available()
  {
    return {
        "solver.linearelliptic.cg.dd"
    };
  }

  static Dune::ParameterTree createSampleDescription(const std::string type)
  {
//    if (type == "solver.linearelliptic.cg.dd")
//      return SolverContinuousGalerkinDD< DomainFieldType, dimDomain, RangeFieldType, dimRange >::createSampleDescription();
//    else
      DUNE_THROW(Dune::RangeError,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " unknown solver '" << type << "' requested!");
  } // ... createSampleDescription(...)

  static SolverInterface< DomainFieldType, dimDomain, RangeFieldType, 1 >*
    create(const std::string type,
           const Dune::ParameterTree description = Dune::ParameterTree())
  {
//    if (type == "solver.linearelliptic.cg.dd")
//      return SolverContinuousGalerkinDD< DomainFieldType, dimDomain, RangeFieldType, 1 >::create(description);
//    else
//      DUNE_THROW(Dune::RangeError,
//                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " unknown solver '" << type << "' requested!");
  } // ... create(...)
}; // class Solvers< DomainFieldType, dimDomain, RangeFieldType, 1 >


} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_SOLVER_HH
