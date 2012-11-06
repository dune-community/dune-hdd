#ifndef DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_FINITEVOLUME_DUNE_PDELAB_HH
#define DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_FINITEVOLUME_DUNE_PDELAB_HH

// system
#include <sstream>

// dune-common
#include <dune/common/exceptions.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/timer.hh>

// dune-stuff
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/grid/boundaryinfo.hh>

namespace Dune {

namespace Detailed {

namespace Solvers {

namespace Stationary {

namespace Linear {

namespace Elliptic {

namespace FiniteVolume {

/**
 *  \todo Add method to create a discrete function.
 *  \todo Change id -> static id()
 */
template< class ModelImp, class GridViewImp, class BoundaryInfoImp, int polynomialOrder >
class DunePdelab
{
public:
  typedef ModelImp ModelType;

  typedef GridViewImp GridViewType;

  typedef BoundaryInfoImp BoundaryInfoType;

  static const int polOrder = polynomialOrder;

  typedef DunePdelab< ModelType, GridViewType, BoundaryInfoType, polOrder > ThisType;

  DunePdelab(const Dune::shared_ptr< const ModelType > model,
                              const Dune::shared_ptr< const GridViewType > gridView,
                              const Dune::shared_ptr< const BoundaryInfoType > boundaryInfo)
    : model_(model)
    , gridView_(gridView)
    , boundaryInfo_(boundaryInfo)
    , initialized_(false)
  {}

  static const std::string id()
  {
    return "detailed.solvers.stationary.linear.elliptic.finitevolume.dune-pdelab";
  }

  const Dune::shared_ptr< const ModelType > model() const
  {
    return model_;
  }

  const Dune::shared_ptr< const GridViewType > gridView() const
  {
    return gridView_;
  }

  const Dune::shared_ptr< const BoundaryInfoType > boundaryInfo() const
  {
    return boundaryInfo_;
  }

  void init(const std::string prefix = "", std::ostream& out = Dune::Stuff::Common::Logger().debug())
  {
    if (!initialized_) {
      // timer
      Dune::Timer timer;

      // spaces
      out << prefix << "initializing spaces... " << std::flush;
      timer.reset();

      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      // done
      initialized_ = true;
    } // if !(initialized_)
  } // void init()

//  Dune::shared_ptr< VectorBackendType > createVector() const
//  {
//    assert(initialized_ && "A vector can only be created after init() has been called!");
//  }

//  Dune::shared_ptr< DiscreteFunctionType > createDiscreteFunction(const std::string name = "discrete_function") const
//  {
//    assert(initialized_ && "Please call init() beafore calling createDiscreteFunction()!");
//  } // ... createDiscreteFunction(...)

//  Dune::shared_ptr< DiscreteFunctionType > createDiscreteFunction(VectorBackendType vector,
//                                                                  const std::string name = "discrete_function") const
//  {
//    assert(initialized_ && "Please call init() beafore calling createDiscreteFunction()!");
//  } // ... createDiscreteFunction(...)

//  void solve(VectorBackendType& solution,
//             const Dune::ParameterTree paramTree = Dune::ParameterTree(),
//             const std::string prefix = "",
//             std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
//  {
//    // preparations
//    assert(initialized_ && "The system can only be solved after init() has been called! ");
//    assert(solution.size() == ansatzSpace_->map().size() && "Given vector has wrong size!");

//    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
//  } // void solve()

//  void visualize(VectorBackendType& vector,
//                 const std::string filename = "solution",
//                 const std::string name = "solution",
//                 const std::string prefix = "",
//                 std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
//  {
//    // preparations
//    assert(initialized_ && "A vector can only be visualized after init() has been called! ");

//    Dune::Timer timer;
//    out << prefix << "writing '" << name << "' to '" << filename;
//    if (dimDomain == 1)
//      out << ".vtp";
//    else
//      out << ".vtu";
//    out << "'... " << std::flush;

//    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
//  }

//  void visualize(Dune::shared_ptr< DiscreteFunctionType > discreteFunction,
//                 const std::string filename = "solution",
//                 const std::string prefix = "",
//                 std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
//  {
//    // preparations
//    assert(initialized_ && "A discrete function can only be visualized after init() has been called! ");
//    Dune::Timer timer;
//    out << prefix << "writing '" << discreteFunction->name() << "' to '" << filename;
//    if (dimDomain == 1)
//      out << ".vtp";
//    else
//      out << ".vtu";
//    out << "'... " << std::flush;
//    typedef Dune::VTKWriter< typename AnsatzSpaceType::GridViewType > VTKWriterType;
//    VTKWriterType vtkWriter(ansatzSpace_->gridView());
//    vtkWriter.addVertexData(discreteFunction);
//    vtkWriter.write(filename);
//    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
//  }

//  const Dune::shared_ptr< const MatrixBackendType > systemMatrix() const
//  {
//    assert(initialized_);
//    return matrix_;
//  }

//  Dune::shared_ptr< MatrixBackendType > systemMatrix()
//  {
//    assert(initialized_);
//    return matrix_;
//  }

//  const Dune::shared_ptr< const VectorBackendType > rightHandSide() const
//  {
//    assert(initialized_);
//    return rhs_;
//  }

//  Dune::shared_ptr< VectorBackendType > rightHandSide()
//  {
//    assert(initialized_);
//    return rhs_;
//  }

//  const Dune::shared_ptr< const PatternType > pattern() const
//  {
//    assert(initialized_);
//    return pattern_;
//  }

private:
  const Dune::shared_ptr< const ModelType > model_;
  const Dune::shared_ptr< const GridViewType > gridView_;
  const Dune::shared_ptr< const BoundaryInfoType > boundaryInfo_;
  bool initialized_;
//  Dune::shared_ptr< const PatternType > pattern_;
//  Dune::shared_ptr< MatrixBackendType > matrix_;
//  Dune::shared_ptr< VectorBackendType > rhs_;
}; // class DuneDetailedDiscretizations

} // namespace ContinuousGalerkin

} // namespace Elliptic

} // namespace Linear

} // namespace Stationary

} // namespace Solvers

} // namespace Detailed

} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_FINITEVOLUME_DUNE_PDELAB_HH
