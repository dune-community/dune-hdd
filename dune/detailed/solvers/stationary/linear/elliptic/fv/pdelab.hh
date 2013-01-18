#ifndef DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_FINITEVOLUME_DUNE_PDELAB_HH
#define DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_FINITEVOLUME_DUNE_PDELAB_HH

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#include <sstream>
#include <math.h>
#include <iostream>
#include <vector>
#include <map>
#include <string>

#include <dune/common/shared_ptr.hh>
#include <dune/common/fvector.hh>
#include <dune/common/static_assert.hh>
#include <dune/common/timer.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/io/file/vtk/vtuwriter.hh>

#include <dune/pdelab/finiteelementmap/p0fem.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/backend/eigenvectorbackend.hh>
#include <dune/pdelab/backend/eigensolverbackend.hh>
#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/flags.hh>

#include <dune/grid/part/interface.hh>

#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/la/container/eigen.hh>

#include "../model/interface.hh"

namespace Dune {
namespace Detailed {
namespace Solvers {
namespace Stationary {
namespace Linear {
namespace Elliptic {
namespace FV {

template< class GridPartImp, int polynomialOrder, class RangeFieldImp, int dimensionRange >
class Pdelab;

//template< class ModelImp, class GridViewImp, class BoundaryInfoImp, int polynomialOrder >
template< class GridPartImp, int polynomialOrder, class RangeFieldImp >
class Pdelab< GridPartImp, polynomialOrder, RangeFieldImp, 1 >
{
public:
  typedef Pdelab< GridPartImp, polynomialOrder, RangeFieldImp, 1 > ThisType;

  typedef Dune::grid::Part::Interface< typename GridPartImp::Traits > GridPartType;

private:
  typedef typename GridPartType::GridViewType GridViewType;

public:
  static const int polOrder = polynomialOrder;

  typedef typename GridPartType::ctype DomainFieldType;

  static const int dimDomain = GridPartType::dimension;

  typedef RangeFieldImp RangeFieldType;

  static const int dimRange = 1;

  typedef Dune::Stuff::Grid::BoundaryInfo::Interface< typename GridPartType::GridViewType > BoundaryInfoType;

  typedef typename Dune::Detailed::Solvers
      ::Stationary
      ::Linear
      ::Elliptic
      ::Model::Interface< DomainFieldType, dimDomain, RangeFieldType, dimRange > ModelType;

private:
  typedef Dune::PDELab::P0LocalFiniteElementMap< DomainFieldType, RangeFieldType, dimDomain > FiniteElementMapType;

  typedef Dune::PDELab::NoConstraints ConstraintsType;

  typedef Dune::PDELab::EigenVectorBackend PdelabVectorBackendType;

public:
  typedef Dune::PDELab::GridFunctionSpace< GridViewType, FiniteElementMapType, ConstraintsType, PdelabVectorBackendType >
    AnsatzSpaceType;

  typedef AnsatzSpaceType TestSpaceType;

private:
  /** a local operator for solving the equation
   *
   *  -\nabla \cdot (a*\nabla u) + b*u = f   in \Omega
   *                                 u = g   on \Gamma_D\subseteq\partial\Omega
   *               -a*\nabla u \cdot n = j   on \Gamma_N = \partial\Omega\setminus\Gamma_D
   *
   * with cell-centered finite volumes on axiparallel, structured grids
   *
   * \tparam B a function indicating the type of boundary condition
   * \tparam G a function for the values of the Dirichlet boundary condition
   */
  class EllipticLocalOperator :    // implement jacobian evaluation in base classes
    public Dune::PDELab::NumericalJacobianApplyVolume< EllipticLocalOperator >,
    public Dune::PDELab::NumericalJacobianVolume< EllipticLocalOperator >,
    public Dune::PDELab::NumericalJacobianApplySkeleton< EllipticLocalOperator >,
    public Dune::PDELab::NumericalJacobianSkeleton< EllipticLocalOperator >,
    public Dune::PDELab::NumericalJacobianApplyBoundary< EllipticLocalOperator >,
    public Dune::PDELab::NumericalJacobianBoundary< EllipticLocalOperator >,
    public Dune::PDELab::FullSkeletonPattern,                     // matrix entries skeleton
    public Dune::PDELab::FullVolumePattern,
    public Dune::PDELab::LocalOperatorDefaultFlags
  {
  public:
    // pattern assembly flags
    enum { doPatternVolume = true };
    enum { doPatternSkeleton = true };

    // residual assembly flags
    enum { doAlphaVolume  = true };
    enum { doAlphaSkeleton  = true };                             // assemble skeleton term
    enum { doAlphaBoundary  = true };

    // constructor
    EllipticLocalOperator(const ModelType& model, const BoundaryInfoType& boundaryInfo)
        : localModel_(model), localBoundaryInfo_(boundaryInfo)
    {}

    // volume integral depending on test and ansatz functions
    template<typename EntityType, typename TrialFunctionSpaceType, typename TrialVector,
                typename TestFunctionSpaceType, typename ResidualType>
    void alpha_volume (const EntityType& entity, const TrialFunctionSpaceType& trialFunction,
                       const TrialVector& trialVector, const TestFunctionSpaceType& testFunction,
                       ResidualType& residual) const
    {
      // domain and range field type
      typedef typename TrialFunctionSpaceType::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::DomainFieldType DomainFieldType;
      typedef  Dune::FieldVector<DomainFieldType,dimDomain> DomainType;
      typedef  Dune::FieldVector<RangeFieldType,dimRange> RangeType;

      // evaluate reaction term
      const DomainType center = entity.geometry().center();
//      RangeType b = 0.0;
      const RangeType f = localModel_.force()->evaluate(center);
      residual.accumulate(trialFunction, 0, (/*b*trialVector(trialFunction, 0)*/ - f)*entity.geometry().volume());
    } // void alpha_volume()

    // skeleton integral depending on test and ansatz functions
    // each face is only visited ONCE!
    template<typename IntersectionType, typename TrialFunctionSpaceType, typename TrialVector,
                typename TestFunctionSpaceType, typename ResidualType>
    void alpha_skeleton (const IntersectionType& intersection,
                         const TrialFunctionSpaceType& trialFunction_s, const TrialVector& trialVector_s,
                         const TestFunctionSpaceType& testFunction_s, const TrialFunctionSpaceType& trialFunction_n,
                         const TrialVector& trialVector_n, const TestFunctionSpaceType& testFunction_n,
                         ResidualType& residual_s, ResidualType& residual_n) const
    {
      // domain and range field type
      typedef typename TrialFunctionSpaceType::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::DomainFieldType DomainFieldType;
      typedef  Dune::FieldVector<RangeFieldType,dimRange> RangeType;

      // distance between cell centers in global coordinates
      Dune::FieldVector<DomainFieldType,dimDomain> inside_global = intersection.inside()->geometry().center();
      Dune::FieldVector<DomainFieldType,dimDomain> outside_global = intersection.outside()->geometry().center();
      outside_global -= inside_global;
      RangeType distance = outside_global.two_norm();

      // face geometry
      RangeType face_volume = intersection.geometry().volume();

      // diffusive flux for both sides
      RangeType a;
      localModel_.diffusion()->evaluate(inside_global,a);
      residual_s.accumulate(trialFunction_s,0,-a*(trialVector_n(trialFunction_n,0)-trialVector_s(trialFunction_s,0))
                                                                                                *face_volume/distance);
      residual_n.accumulate(trialFunction_n,0, a*(trialVector_n(trialFunction_n,0)-trialVector_s(trialFunction_s,0))
                                                                                                *face_volume/distance);
    } // void alpha_skeleton()

    // skeleton integral depending on test and ansatz functions
    // Here Dirichlet and Neumann boundary conditions are evaluated
    template<typename IntersectionType, typename TrialFunctionSpaceType, typename TrialVector,
                typename TestFunctionSpaceType, typename ResidualType>
    void alpha_boundary (const IntersectionType& intersection,
                         const TrialFunctionSpaceType& trialFunction_s, const TrialVector& trialVector_s,
                         const TestFunctionSpaceType& /*testFunction_s*/, ResidualType& residual_s) const
    {
      // domain and range field type
      typedef typename TrialFunctionSpaceType::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::DomainFieldType DomainFieldType;
      typedef  Dune::FieldVector<RangeFieldType,dimRange> RangeType;

      // face geometry
      Dune::FieldVector<DomainFieldType,dimDomain> face_center = intersection.geometry().center();
      RangeType face_volume = intersection.geometry().volume();

      // do things depending on boundary condition type (call intersection() to get original intersection)
      if (localBoundaryInfo_.dirichlet(intersection.intersection())) {
          const RangeType g = localModel_.dirichlet()->evaluate(face_center);
          const RangeType a = localModel_.diffusion()->evaluate(face_center);
          Dune::FieldVector<DomainFieldType,dimDomain> inside_global = intersection.inside()->geometry().center();
          inside_global -= face_center;
          RangeType distance = inside_global.two_norm();
          residual_s.accumulate(trialFunction_s, 0, -a*(g-trialVector_s(trialFunction_s,0))*face_volume/distance);
          return;
      } else if (localBoundaryInfo_.neumann(intersection.intersection())) {
          const RangeType neumannValue = localModel_.neumann()->evaluate(face_center);
          residual_s.accumulate(trialFunction_s, 0, -neumannValue * face_volume);
          return;
      }
    } // void alpha_boundary()

  private:
    const ModelType& localModel_;
    const BoundaryInfoType& localBoundaryInfo_;
  }; // class EllipticLocalOperator

  // Make grid operator (left hand side)
  typedef EllipticLocalOperator LocalOperatorType;

  typedef PdelabVectorBackendType::MatrixBackend PdelabMatrixBackendType;

  typedef Dune::PDELab::EmptyTransformation ConstraintsContainerType;

  typedef Dune::PDELab::GridOperator< AnsatzSpaceType, TestSpaceType, LocalOperatorType,
                                      PdelabMatrixBackendType, DomainFieldType, RangeFieldType, RangeFieldType,
                                      ConstraintsContainerType, ConstraintsContainerType > GridOperatorType;

  typedef typename GridOperatorType::Traits::Domain PdelabVectorType;

  typedef typename GridOperatorType::template MatrixContainer< RangeFieldType >::Type PdelabMatrixType;

  typedef typename GridOperatorType::Traits::TrialGridFunctionSpace TrialDiscreteFunctionSpaceType;

  typedef Dune::PDELab::EigenBackend_BiCGSTAB_Diagonal LinearSolverType;

//  typedef typename GridOperatorType::Traits::Domain VectorBackendType;

public:
  typedef Dune::Stuff::LA::Container::Eigen::RowMajorSparseMatrix< RangeFieldType > MatrixType;

  typedef Dune::Stuff::LA::Container::Eigen::DenseVector< RangeFieldType > VectorType;

private:
  typedef Dune::PDELab::DiscreteGridFunction< AnsatzSpaceType, PdelabVectorType > PdelabDiscreteAnsatzFunctionType;

  typedef Dune::PDELab::DiscreteGridFunction< TestSpaceType, PdelabVectorType > PdelabDiscreteTestFunctionType;

public:
  typedef Dune::PDELab::VTKGridFunctionAdapter< PdelabDiscreteAnsatzFunctionType > DiscreteAnsatzFunctionType;

  typedef Dune::PDELab::VTKGridFunctionAdapter< PdelabDiscreteTestFunctionType > DiscreteTestFunctionType;

  Pdelab(const Dune::shared_ptr< const GridPartType > _gridPart,
         const Dune::shared_ptr< const BoundaryInfoType > _boundaryInfo,
         const Dune::shared_ptr< const ModelType > _model)
    : gridPart_(_gridPart)
    , gridView_(new GridViewType(gridPart_->gridView()))
    , boundaryInfo_(_boundaryInfo)
    , model_(_model)
    , initialized_(false)
  {
    // sanity checks
    if (model_->diffusionOrder() < 0)
      DUNE_THROW(Dune::RangeError,
                 "\nERROR: negative integration order given in model.diffusionOrder()!");
    if (model_->forceOrder() < 0)
      DUNE_THROW(Dune::RangeError,
                 "\nERROR: negative integration order given in model.forceOrder()!");
    if (model_->dirichletOrder() < 0)
      DUNE_THROW(Dune::RangeError,
                 "\nERROR: negative integration order given in model.dirichletOrder()!");
    if (model_->neumannOrder() < 0)
      DUNE_THROW(Dune::RangeError,
                 "\nERROR: negative integration order given in model.neumannOrder()!");
  }

  static const std::string id()
  {
    return "detailed.solvers.stationary.linear.elliptic.fv.pdelab";
  }

  const Dune::shared_ptr< const ModelType > model() const
  {
    return model_;
  }

  const Dune::shared_ptr< const GridPartType > gridPart() const
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
      Dune::Timer timer;

      out << prefix << "initializing function space... " << std::flush;
      timer.reset();
      finiteElementMap_ = Dune::shared_ptr< FiniteElementMapType >(new FiniteElementMapType(
          Dune::GeometryType(gridView_->template begin< 0 >()->geometry().type())));
      ansatzSpace_ = Dune::shared_ptr< AnsatzSpaceType >(new AnsatzSpaceType(*gridView_, *finiteElementMap_));
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      out << prefix << "initializing operator... " << std::flush;
      timer.reset();
      localOperator_ = Dune::shared_ptr< LocalOperatorType >(new LocalOperatorType(*model_, *boundaryInfo_));
      gridOperator_ = Dune::shared_ptr< GridOperatorType >(new GridOperatorType(*ansatzSpace_,
                                                                                *ansatzSpace_,
                                                                                *localOperator_));
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      out << prefix << "initializing matrix and vectors... " << std::flush;
      timer.reset();
      matrix_ = Dune::shared_ptr< PdelabMatrixType >(new PdelabMatrixType(*gridOperator_));
      *matrix_ = 0.0;
      PdelabVectorType initialGuess(*ansatzSpace_, 0.0);
      PdelabVectorType residual(*ansatzSpace_, 0.0);
      rhs_ = Dune::shared_ptr< PdelabVectorType >(new PdelabVectorType(*ansatzSpace_, 0.0));
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      out << prefix << "assembling system... " << std::flush;
      timer.reset();
      // given a start vector initialGuess
      // assemble jacobian system matrix_ = \nabla R(initialGuess) (with residualoperator R)
      gridOperator_->jacobian(initialGuess, *matrix_);
      // compute the residual R(initialGuess)
      gridOperator_->residual(initialGuess, residual);
      *rhs_ -= residual;
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      // done
      initialized_ = true;
    } // if !(initialized_)
  } // void init()

  const AnsatzSpaceType& ansatzSpace() const
  {
    assert(initialized_);
    return *ansatzSpace_;
  }

  const TestSpaceType& testSpace() const
  {
    assert(initialized_);
    return *ansatzSpace_;
  }

  Dune::shared_ptr< VectorType > createAnsatzVector() const
  {
    assert(initialized_ && "Please call init() before calling createAnsatzVector()!");
    PdelabVectorType pdelabVector(*ansatzSpace_, 0.0);
    Dune::shared_ptr< VectorType > vector(new VectorType(pdelabVector));
    return vector;
  } // Dune::shared_ptr< VectorType > createAnsatzVector() const

  Dune::shared_ptr< VectorType > createTestVector() const
  {
    return createAnsatzVector();
  } // Dune::shared_ptr< VectorType > createTestVector() const

  Dune::shared_ptr< DiscreteAnsatzFunctionType > createAnsatzFunction(const std::string name = "ansatzFunction") const
  {
    assert(initialized_ && "Please call init() before calling createAnsatzFunction()!");
    Dune::shared_ptr< PdelabVectorType > pdelabVector(new PdelabVectorType(*ansatzSpace_, 0.0));
    Dune::shared_ptr< PdelabDiscreteAnsatzFunctionType > pdelabAnsatzFunction(
          new PdelabDiscreteAnsatzFunctionType(ansatzSpace_, pdelabVector));
    Dune::shared_ptr< DiscreteAnsatzFunctionType > ansatzFunction(new DiscreteAnsatzFunctionType(pdelabAnsatzFunction,
                                                                                                 name));
    return ansatzFunction;
  } // Dune::shared_ptr< DiscreteAnsatzFunctionType > createAnsatzFunction(...) const

  Dune::shared_ptr< DiscreteAnsatzFunctionType > createAnsatzFunction(Dune::shared_ptr< VectorType > vector,
                                                                      const std::string name = "ansatzFunction") const
  {
    assert(initialized_ && "Please call init() before calling createAnsatzFunction()!");
    Dune::shared_ptr< PdelabVectorType > pdelabVector(new PdelabVectorType(*ansatzSpace_, 0.0));
    pdelabVector->base() = vector->base();
    Dune::shared_ptr< PdelabDiscreteAnsatzFunctionType > pdelabAnsatzFunction(
          new PdelabDiscreteAnsatzFunctionType(ansatzSpace_, pdelabVector));
    Dune::shared_ptr< DiscreteAnsatzFunctionType > ansatzFunction(new DiscreteAnsatzFunctionType(pdelabAnsatzFunction,
                                                                                                 name));
    return ansatzFunction;
  } // Dune::shared_ptr< DiscreteAnsatzFunctionType > createAnsatzFunction(...) const

  Dune::shared_ptr< DiscreteTestFunctionType > createTestFunction(const std::string name = "testFunction") const
  {
    assert(initialized_ && "Please call init() before calling createTestFunction()!");
    Dune::shared_ptr< PdelabVectorType > pdelabVector(new PdelabVectorType(*ansatzSpace_, 0.0));
    Dune::shared_ptr< PdelabDiscreteTestFunctionType > pdelabTestFunction(
          new PdelabDiscreteTestFunctionType(ansatzSpace_, pdelabVector));
    Dune::shared_ptr< DiscreteTestFunctionType > testFunction(new DiscreteTestFunctionType(pdelabTestFunction,
                                                                                                 name));
    return testFunction;
  } // Dune::shared_ptr< DiscreteTestFunctionType > createTestFunction(...) const

  Dune::shared_ptr< DiscreteTestFunctionType > createTestFunction(Dune::shared_ptr< VectorType > vector,
                                                                      const std::string name = "testFunction") const
  {
    assert(initialized_ && "Please call init() before calling createTestFunction()!");
    Dune::shared_ptr< PdelabVectorType > pdelabVector(new PdelabVectorType(*ansatzSpace_, 0.0));
    pdelabVector->base() = vector->base();
    Dune::shared_ptr< PdelabDiscreteTestFunctionType > pdelabTestFunction(
          new PdelabDiscreteTestFunctionType(ansatzSpace_, pdelabVector));
    Dune::shared_ptr< DiscreteTestFunctionType > testFunction(new DiscreteTestFunctionType(pdelabTestFunction,
                                                                                                 name));
    return testFunction;
  } // Dune::shared_ptr< DiscreteTestFunctionType > createTestFunction(...) const

  void solve(VectorType& solutionVector,
             /*const std::string linearSolverType = "eigen.iterative.bicgstab.diagonal",*/
             const unsigned int linearSolverMaxIter = 5000,
             const double linearSolverPrecision = 1e-12,
             const std::string prefix = "",
             std::ostream& out = Dune::Stuff::Common::Logger().debug())
  {
    assert(initialized_ && "Please call init() before calling solve()!");
    out << prefix << "solving linear system (of size " << matrix_->rows()
        << "x" << matrix_->cols() << ")... " << std::flush;
    Dune::Timer timer;
    LinearSolverType linearSolver(linearSolverMaxIter);
    typename PdelabVectorType::ElementType defect = linearSolver.norm(*rhs_);
    typename PdelabVectorType::ElementType mindefect = 1e-99;
    PdelabVectorType update(gridOperator_->trialGridFunctionSpace(), 0.0);
    typename PdelabVectorType::ElementType red = std::min(linearSolverPrecision, defect/(mindefect));
    // solve matrix_ * update = residual = -*rhs_ for update
    linearSolver.apply(*matrix_, update, *rhs_, red);
    // and store the solution u = u0 + update = update (u0 = 0)
    solutionVector.base() = update.base();
    out << prefix << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void solve(...)

  void visualizeAnsatzVector(VectorType& vector,
                             const std::string filename = id() + ".ansatzVector",
                             const std::string name = id() + ".ansatzVector",
                             const std::string prefix = "",
                             std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
  {
    // preparations
    assert(initialized_ && "Please call init() before calling visualize()");
    Dune::Timer timer;
    out << prefix << "writing '" << name << "'" << std::endl;
    out << prefix << "     to '" << filename;
    if (dimDomain == 1)
      out << ".vtp";
    else
      out << ".vtu";
    out << "'... " << std::flush;
    Dune::shared_ptr< VectorType > tmpVectorPtr(new VectorType(vector));
    const Dune::shared_ptr< const DiscreteAnsatzFunctionType > discreteFunction
        = createAnsatzFunction(tmpVectorPtr, name);
    visualizeFunction(discreteFunction, filename, "", Dune::Stuff::Common::Logger().devnull());
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void visualizeAnsatzVector(...)

  void visualizeTestVector(VectorType& vector,
                             const std::string filename = id() + ".testVector",
                             const std::string name = id() + ".testVector",
                             const std::string prefix = "",
                             std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
  {
    // preparations
    assert(initialized_ && "Please call init() before calling visualize()");
    Dune::Timer timer;
    out << prefix << "writing '" << name << "'" << std::endl;
    out << prefix << "     to '" << filename;
    if (dimDomain == 1)
      out << ".vtp";
    else
      out << ".vtu";
    out << "'... " << std::flush;
    Dune::shared_ptr< VectorType > tmpVectorPtr(new VectorType(vector));
    const Dune::shared_ptr< const DiscreteTestFunctionType > discreteFunction
        = createTestFunction(tmpVectorPtr, name);
    visualizeFunction(discreteFunction, filename, "", Dune::Stuff::Common::Logger().devnull());
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void visualizeAnsatzVector(...)

  void visualizeFunction(const Dune::shared_ptr< const DiscreteAnsatzFunctionType > discreteFunction,
                         const std::string filename = id() + ".discreteAnsatzFunction",
                         const std::string prefix = "",
                         std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
  {
    // preparations
    assert(initialized_ && "Please call init() before calling visualize()");
    Dune::Timer timer;
    out << prefix << "writing '" << discreteFunction->name() << "' to '" << filename;
    if (dimDomain == 1)
      out << ".vtp";
    else
      out << ".vtu";
    out << "'... " << std::flush;
    typedef Dune::VTKWriter< GridViewType > VTKWriterType;
    VTKWriterType vtkWriter(*gridView_, Dune::VTK::conforming);
    vtkWriter.addCellData(discreteFunction);
    vtkWriter.write(filename, Dune::VTK::appendedraw);
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void visualizeFunction(...)

//  void visualizeFunction(const Dune::shared_ptr< const DiscreteTestFunctionType > discreteFunction,
//                         const std::string filename = id() + ".discreteTestFunction",
//                         const std::string prefix = "",
//                         std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
//  {
//    // preparations
//    assert(initialized_ && "Please call init() before calling visualize()");
//    Dune::Timer timer;
//    out << prefix << "writing '" << discreteFunction->name() << "' to '" << filename;
//    if (dimDomain == 1)
//      out << ".vtp";
//    else
//      out << ".vtu";
//    out << "'... " << std::flush;
//    typedef Dune::VTKWriter< GridViewType > VTKWriterType;
//    VTKWriterType vtkWriter(*gridView_, Dune::VTK::conforming);
//    vtkWriter.addCellData(discreteFunction);
//    vtkWriter.write(filename, Dune::VTK::appendedraw);
//    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
//  } // void visualizeFunction(...)

private:
  const Dune::shared_ptr< const GridPartType > gridPart_;
  const Dune::shared_ptr< const GridViewType > gridView_;
  const Dune::shared_ptr< const BoundaryInfoType > boundaryInfo_;
  const Dune::shared_ptr< const ModelType > model_;
  bool initialized_;
  Dune::shared_ptr< FiniteElementMapType > finiteElementMap_;
  Dune::shared_ptr< AnsatzSpaceType > ansatzSpace_;
  Dune::shared_ptr< LocalOperatorType> localOperator_;
  Dune::shared_ptr< GridOperatorType> gridOperator_;
  Dune::shared_ptr< PdelabMatrixType > matrix_;
  Dune::shared_ptr< PdelabVectorType > rhs_;
}; // class Pdelab

} // namespace FV
} // namespace Elliptic
} // namespace Linear
} // namespace Stationary
} // namespace Solvers
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_FINITEVOLUME_DUNE_PDELAB_HH
