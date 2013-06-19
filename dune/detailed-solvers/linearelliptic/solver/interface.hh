#ifndef DUNE_DETAILED_SOLVERS_LINEARELLIPTIC_SOLVER_INTERFACE_HH
#define DUNE_DETAILED_SOLVERS_LINEARELLIPTIC_SOLVER_INTERFACE_HH

#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/common/bartonnackmanifcheck.hh>

#include <dune/grid/part/interface.hh>

#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/common/color.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/parameter/tree.hh>

#include "../model/interface.hh"

namespace Dune {
namespace DetailedSolvers {
namespace LinearElliptic {


template< class Traits >
class SolverInterface
{
public:
  typedef SolverInterface< Traits >     ThisType;
  typedef typename Traits::derived_type derived_type;

  typedef Dune::grid::Part::Interface< typename Traits::GridPartTraits > GridPartType;

  static const int polOrder = Traits::polOrder;

  typedef typename GridPartType::ctype    DomainFieldType;
  static const int                        dimDomain = GridPartType::dimension;
  typedef typename Traits::RangeFieldType RangeFieldType;
  static const int                        dimRange = Traits::dimRange;

  typedef Dune::Stuff::GridboundaryInterface< typename GridPartType::GridViewType > BoundaryInfoType;
  typedef ModelInterface< DomainFieldType, dimDomain, RangeFieldType, dimRange >            ModelType;

  typedef typename Traits::VectorType VectorType;
  typedef Dune::Stuff::Common::ExtendedParameterTree SettingsType;

  static const std::string id()
  {
    return "solver.linearelliptic";
  }

  std::shared_ptr< const GridPartType > gridPart() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().gridPart());
    return asImp().gridPart();
  }

  std::shared_ptr< const BoundaryInfoType > boundaryInfo() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().boundaryInfo());
    return asImp().boundaryInfo();
  }

  std::shared_ptr< const ModelType > model() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().model());
    return asImp().model();
  }

  std::shared_ptr< VectorType > createVector() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().createVector());
    return asImp().createVector();
  }

  void visualize(const std::shared_ptr< const VectorType > vector,
                 const std::string filename,
                 const std::string name,
                 std::ostream& out = Dune::Stuff::Common::Logger().devnull(),
                 const std::string prefix = "") const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().visualize(vector, filename, name, out, prefix));
    asImp().visualize(vector, filename, name, out, prefix);
  }

  void init(std::ostream& out = Dune::Stuff::Common::Logger().devnull(),
            const std::string prefix = "")
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().init(out, prefix));
    asImp().init(out, prefix);
  }

  bool initialized() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().initialized());
    return asImp().initialized();
  }

  void solve(std::shared_ptr< VectorType > solution,
             const std::string linearSolverType,
             const double linearSolverPrecision,
             const size_t linearSolverMaxIterations,
             std::ostream& out = Dune::Stuff::Common::Logger().devnull(),
             const std::string prefix = "") const
  {
    if (!initialized())
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " call init() before calling solve()!");
    CHECK_INTERFACE_IMPLEMENTATION(asImp().solve(solution,
                                                 linearSolverType,
                                                 linearSolverPrecision,
                                                 linearSolverMaxIterations,
                                                 out,
                                                 prefix));
    asImp().solve(solution, linearSolverType, linearSolverPrecision, linearSolverMaxIterations, out, prefix);
  }

  derived_type& asImp()
  {
    return static_cast< derived_type& >(*this);
  }

  const derived_type& asImp() const
  {
    return static_cast< const derived_type& >(*this);
  }
}; // class SolverInterface< DomainFieldType, dimDomain, RangeFieldType, 1 >


template< class Traits >
class SolverParametricInterface
{
public:
  typedef SolverParametricInterface< Traits > ThisType;
  typedef typename Traits::derived_type       derived_type;

  typedef Dune::grid::Part::Interface< typename Traits::GridPartTraits > GridPartType;

  static const int polOrder = Traits::polOrder;

  typedef typename GridPartType::ctype    DomainFieldType;
  static const int                        dimDomain = GridPartType::dimension;
  typedef typename Traits::RangeFieldType RangeFieldType;
  static const int                        dimRange = Traits::dimRange;

  typedef Dune::Stuff::GridboundaryInterface< typename GridPartType::GridViewType > BoundaryInfoType;
  typedef ModelInterface< DomainFieldType, dimDomain, RangeFieldType, dimRange >            ModelType;

  typedef typename ModelType::ParamType ParamType;

  typedef typename Traits::VectorType VectorType;

  static const std::string id()
  {
    return "solver.linearelliptic";
  }

  std::shared_ptr< const GridPartType > gridPart() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().gridPart());
    return asImp().gridPart();
  }

  std::shared_ptr< const BoundaryInfoType > boundaryInfo() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().boundaryInfo());
    return asImp().boundaryInfo();
  }

  std::shared_ptr< const ModelType > model() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().model());
    return asImp().model();
  }

  std::shared_ptr< VectorType > createVector() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().createVector());
    return asImp().createVector();
  }

  void visualize(const std::shared_ptr< const VectorType > vector,
                 const std::string filename,
                 const std::string name,
                 std::ostream& out = Dune::Stuff::Common::Logger().devnull(),
                 const std::string prefix = "") const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().visualize(vector, filename, name, out, prefix));
    asImp().visualize(vector, filename, name, out, prefix);
  }

  void init(std::ostream& out = Dune::Stuff::Common::Logger().devnull(),
            const std::string prefix = "")
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().init(out, prefix));
    asImp().init(out, prefix);
  }

  bool initialized() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().initialized());
    return asImp().initialized();
  }

  void solve(std::shared_ptr< VectorType > solution,
             const ParamType& mu,
             const std::string linearSolverType,
             const double linearSolverPrecision,
             const size_t linearSolverMaxIterations,
             std::ostream& out = Dune::Stuff::Common::Logger().devnull(),
             const std::string prefix = "") const
  {
    if (!initialized())
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " call init() before calling solve()!");
    CHECK_INTERFACE_IMPLEMENTATION(asImp().solve(solution,
                                                 mu,
                                                 linearSolverType,
                                                 linearSolverPrecision,
                                                 linearSolverMaxIterations,
                                                 out,
                                                 prefix));
    asImp().solve(solution, mu, linearSolverType, linearSolverPrecision, linearSolverMaxIterations, out, prefix);
  }

  derived_type& asImp()
  {
    return static_cast< derived_type& >(*this);
  }

  const derived_type& asImp() const
  {
    return static_cast< const derived_type& >(*this);
  }
}; // class SolverInterface< DomainFieldType, dimDomain, RangeFieldType, 1 >


} // namespace LinearElliptic
} // namespace DetailedSolvers
} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_LINEARELLIPTIC_SOLVER_INTERFACE_HH
