#ifndef DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_FINITEVOLUME_DUNE_PDELAB_HH
#define DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_FINITEVOLUME_DUNE_PDELAB_HH

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

// system
#include <sstream>
#include <math.h>
#include <iostream>
#include <vector>
#include <map>
#include <string>

// dune-common
#include <dune/common/exceptions.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/timer.hh>
#include <dune/common/mpihelper.hh>
#include <dune/common/fvector.hh>
#include <dune/common/static_assert.hh>

// dune-stuff
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/grid/boundaryinfo.hh>

// dune-geometry
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>

// dune-grid
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/yaspgrid.hh>
#if HAVE_ALBERTA
#include <dune/grid/albertagrid.hh>
#include <dune/grid/albertagrid/dgfparser.hh>
#endif
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif
#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#endif

// dune-pdelab
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/common/geometrywrapper.hh>
#include <dune/pdelab/finiteelementmap/p0fem.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/constraints/constraints.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/backend/eigenvectorbackend.hh>
#include <dune/pdelab/backend/eigenmatrixbackend.hh>
#include <dune/pdelab/backend/eigensolverbackend.hh>
#include <dune/pdelab/backend/backendselector.hh>
#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/flags.hh>

// dune-istl
//#include <dune/istl/bvector.hh>
//#include <dune/istl/operators.hh>
//#include <dune/istl/solvers.hh>
//#include <dune/istl/preconditioners.hh>
//#include <dune/istl/io.hh>
//#include <dune/istl/superlu.hh>


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

  static const int dimDomain = ModelType::dimDomain;

  //TODO: which types should be private?

  // Choose domain and range field type
  typedef typename GridViewType::Grid::ctype CoordinateType;
  typedef double RangeFieldType;
  static const int dim = GridViewType::dimension;

  // Make grid function space
  typedef Dune::PDELab::P0LocalFiniteElementMap<CoordinateType,RangeFieldType,dim> FiniteElementMap;
  typedef Dune::PDELab::NoConstraints ConstraintsType;
  typedef Dune::PDELab::EigenVectorBackend VectorBackendType;
  typedef Dune::PDELab::GridFunctionSpace<GridViewType,FiniteElementMap,ConstraintsType,VectorBackendType> GridFunctionSpace;

  /** a local operator for solving the equation
   *
   *   - \Delta u + a*u = f   in \Omega
   *                  u = g   on \Gamma_D\subseteq\partial\Omega
   *  -\nabla u \cdot n = j   on \Gamma_N = \partial\Omega\setminus\Gamma_D
   *
   * with cell-centered finite volumes on axiparallel, structured grids
   *
   * \tparam B a function indicating the type of boundary condition
   * \tparam G a function for the values of the Dirichlet boundary condition
   */
  class EllipticLocalOperator :    // implement jacobian evaluation in base classes
    public Dune::PDELab::NumericalJacobianApplyVolume<EllipticLocalOperator>,
    public Dune::PDELab::NumericalJacobianVolume<EllipticLocalOperator>,
    public Dune::PDELab::NumericalJacobianApplySkeleton<EllipticLocalOperator>,
    public Dune::PDELab::NumericalJacobianSkeleton<EllipticLocalOperator>,
    public Dune::PDELab::NumericalJacobianApplyBoundary<EllipticLocalOperator>,
    public Dune::PDELab::NumericalJacobianBoundary<EllipticLocalOperator>,
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

    // volume integral depending on test and ansatz functions
    template<typename EntityType, typename TrialFunctionSpaceType, typename TrialVector, typename TestFunctionSpaceType, typename ResidualType>
    void alpha_volume (const EntityType& entity, const TrialFunctionSpaceType& trialFunction, const TrialVector& trialVector,
                       const TestFunctionSpaceType& testFunction, ResidualType& residual) const
    {
      // domain and range field type
      typedef typename TrialFunctionSpaceType::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::DomainFieldType DomainFieldType;
      typedef typename TrialFunctionSpaceType::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::RangeFieldType RangeType;

      // evaluate reaction term
      RangeType a = 0.0;
      RangeType f = 0.0;

      residual.accumulate(trialFunction,0,(a*trialVector(trialFunction,0)-f)*entity.geometry().volume());
    } // void alpha_volume()

    // skeleton integral depending on test and ansatz functions
    // each face is only visited ONCE!
    template<typename IntersectionType, typename TrialFunctionSpaceType, typename TrialVector, typename TestFunctionSpaceType, typename ResidualType>
    void alpha_skeleton (const IntersectionType& intersection,
                         const TrialFunctionSpaceType& trialFunction_s, const TrialVector& trialVector_s,
                         const TestFunctionSpaceType& testFunction_s, const TrialFunctionSpaceType& trialFunction_n,
                         const TrialVector& trialVector_n, const TestFunctionSpaceType& testFunction_n,
                         ResidualType& residual_s, ResidualType& residual_n) const
    {
      // domain and range field type
      typedef typename TrialFunctionSpaceType::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::DomainFieldType DomainFieldType;
      typedef typename TrialFunctionSpaceType::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::RangeFieldType RangeType;

      // distance between cell centers in global coordinates
      Dune::FieldVector<DomainFieldType,dim> inside_global = intersection.inside()->geometry().center();
      Dune::FieldVector<DomainFieldType,dim> outside_global = intersection.outside()->geometry().center();
      inside_global -= outside_global;
      RangeType distance = inside_global.two_norm();

      // face geometry
      RangeType face_volume = intersection.geometry().volume();

      // diffusive flux for both sides
      residual_s.accumulate(trialFunction_s,0,-(trialVector_n(trialFunction_n,0)-trialVector_s(trialFunction_s,0))*face_volume/distance);
      residual_n.accumulate(trialFunction_n,0,(trialVector_n(trialFunction_n,0)-trialVector_s(trialFunction_s,0))*face_volume/distance);
    } // void alpha_skeleton()

    // skeleton integral depending on test and ansatz functions
    // Here Dirichlet and Neumann boundary conditions are evaluated
    template<typename IntersectionType, typename TrialFunctionSpaceType, typename TrialVector, typename TestFunctionSpaceType, typename ResidualType>
    void alpha_boundary (const IntersectionType& intersection,
                         const TrialFunctionSpaceType& trialFunction_s, const TrialVector& trialVector_s,
                         const TestFunctionSpaceType& testFunction_s, ResidualType& residual_s) const
    {
      // domain and range field type
      typedef typename TrialFunctionSpaceType::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::DomainFieldType DomainFieldType;
      typedef typename TrialFunctionSpaceType::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::RangeFieldType RangeType;

      // face geometry
      Dune::FieldVector<DomainFieldType,dim> face_center = intersection.geometry().center();
      RangeType face_volume = intersection.geometry().volume();

      // evaluate boundary condition type
      int boundaryType;
      if (face_center[0]>1.0-1e-6)
        boundaryType = 0; // Neumann
      else
        boundaryType = 1; // Dirichlet

      // do things depending on boundary condition type
      if (boundaryType==0) // Neumann boundary
        {
          RangeType j; if (face_center[1]<0.5) j = 1.0; else j = -1.0;
          residual_s.accumulate(trialFunction_s,0,j*face_volume);
          return;
        }

      if (boundaryType==1) // Dirichlet boundary
        {
          RangeType g;
          if (face_center[0]<1E-6 && face_center[1]>0.25 && face_center[1]<0.5)
            g = 1.0;
          else
            g = 0.0;
          Dune::FieldVector<DomainFieldType,dim> inside_global = intersection.inside()->geometry().center();
          inside_global -= face_center;
          RangeType distance = inside_global.two_norm();
          residual_s.accumulate(trialFunction_s,0,-(g-trialVector_s(trialFunction_s,0))*face_volume/distance);
          return;
        }
    } // void alpha_boundary()

  }; // class EllipticLocalOperator


  // Make grid operator (left hand side)
  typedef EllipticLocalOperator LocalOperatorType;
  typedef VectorBackendType::MatrixBackend MatrixBackendType;
  typedef Dune::PDELab::EmptyTransformation ConstraintsContainer;
  typedef Dune::PDELab::GridOperator<GridFunctionSpace,GridFunctionSpace,LocalOperatorType,MatrixBackendType,RangeFieldType,RangeFieldType,RangeFieldType,ConstraintsContainer,ConstraintsContainer> GridOperator;

  typedef typename GridOperator::template MatrixContainer<RangeFieldType>::Type MatrixContainer;
  typedef typename GridOperator::Traits::TrialGridFunctionSpace TrialGridFunctionSpace;
  //TODO: replace DiscreteTrialFunction and/or TrialVectorType with VectorBackendType?
  typedef typename Dune::PDELab::BackendVectorSelector<TrialGridFunctionSpace,RangeFieldType>::Type DiscreteTrialFunction;

  typedef Dune::PDELab::EigenBackend_BiCGSTAB_Diagonal LinearSolverType;

  typedef typename GridOperator::Traits::Domain TrialVectorType; // Dof-vector for a discrete function

  typedef Dune::PDELab::DiscreteGridFunction<GridFunctionSpace,TrialVectorType> DiscreteGridFunction;

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

      // function space
      out << prefix << "initializing function space... " << std::flush;
      timer.reset();
      finiteElementMap_ = Dune::shared_ptr< FiniteElementMap >(new FiniteElementMap(Dune::GeometryType(gridView_->template begin< 0 >()->geometry().type())));
      gridFunctionSpace_ = Dune::shared_ptr< GridFunctionSpace >(new GridFunctionSpace(*gridView_, *finiteElementMap_));
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      // grid operator (left hand side)
      out << prefix << "initializing operator... " << std::flush;
      timer.reset();
      LocalOperatorType localOperator;
      gridOperator_ = Dune::shared_ptr< GridOperator >(new GridOperator(*gridFunctionSpace_,*gridFunctionSpace_,localOperator));
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      // Newton method (5 steps):
      // (1) given a start vector u0 (here: *solution_ )
      solution_ = Dune::shared_ptr<TrialVectorType> (new TrialVectorType(*gridFunctionSpace_,0.0));

      // (2) assemble jacobian system A = \nabla R(u0) (here: matrix_ = \nabla R(*solution_) )
      out << prefix << "assembling matrix... " << std::flush;
      timer.reset();
      matrix_ = Dune::shared_ptr< MatrixContainer >(new MatrixContainer(*gridOperator_));
      *matrix_ = 0.0;
      gridOperator_->jacobian(*solution_,*matrix_);
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      // (3) compute the residual r = R(u0)
      out << prefix << "assembling residual... " << std::flush;
      timer.reset();
      residual_ = Dune::shared_ptr< DiscreteTrialFunction >(new DiscreteTrialFunction(gridOperator_->testGridFunctionSpace(),0.0));
      gridOperator_->residual(*solution_,*residual_);
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
//    assert(initialized_ && "Please call init() before calling createDiscreteFunction()!");
//  } // ... createDiscreteFunction(...)

//  Dune::shared_ptr< DiscreteFunctionType > createDiscreteFunction(VectorBackendType vector,
//                                                                  const std::string name = "discrete_function") const
//  {
//    assert(initialized_ && "Please call init() before calling createDiscreteFunction()!");
//  } // ... createDiscreteFunction(...)


    //TODO: adapt function to interface
  void solve(//VectorBackendType& solution,
             const Dune::ParameterTree paramTree = Dune::ParameterTree(),
             const std::string prefix = "",
             std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
  {
    // preparations
    assert(initialized_ && "The system can only be solved after init() has been called! ");
    // timer
    Dune::Timer timer;

    // Assemble solver
    LinearSolverType linearSolver(100);

    typename TrialVectorType::ElementType defect = linearSolver.norm(*residual_);
    // die beiden folgenden Werte wurden im pdelab-Beispiel als Argumente übergeben
    typename TrialVectorType::ElementType mindefect = 1e-99;
    typename TrialVectorType::ElementType reduction = 1e-10;

    // compute correction
    timer.reset();
    TrialVectorType update(gridOperator_->trialGridFunctionSpace(),0.0);
    typename TrialVectorType::ElementType red = std::min(reduction,defect/(mindefect));

    // (4) solve A*update = residual_ for update
    linearSolver.apply(*matrix_,update,*residual_,red); // solver makes right hand side consistent

    // (5) update u = u0 - update
    *solution_ -= update;

    //TODO: *rhs_ and *residual_ must be of the same type
    //*rhs_ = -*residual_;
    out << prefix << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void solve()


  void visualize(//VectorBackendType& vector,
                 const std::string filename = "solution",
                 const std::string name = "solution",
                 const std::string prefix = "",
                 std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
  {
    // preparations
    assert(initialized_ && "A vector can only be visualized after init() has been called! ");
    // timer
    Dune::Timer timer;

    out << prefix << "writing '" << name << "' to '" << filename;
    if (dimDomain == 1)
      out << ".vtp";
    else
      out << ".vtu";
    out << "'... " << std::flush;

    DiscreteGridFunction discreteGridFunction(*gridFunctionSpace_,*solution_);
    Dune::VTKWriter<GridViewType> vtkwriter(*gridView_,Dune::VTK::conforming);
    //TODO: fix
    vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<DiscreteGridFunction>(discreteGridFunction,name));
    vtkwriter.write(filename,Dune::VTK::appendedraw);

    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
  }


#if 0
  //TODO: replace DiscreteFunctionType with DiscreteGridFunction?
  void visualize(Dune::shared_ptr< DiscreteFunctionType > discreteFunction,
                 const std::string filename = "solution",
                 const std::string prefix = "",
                 std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
  {
    // preparations
    assert(initialized_ && "A discrete function can only be visualized after init() has been called! ");
    Dune::Timer timer;
    out << prefix << "writing '" << discreteFunction->name() << "' to '" << filename;
    if (dimDomain == 1)
      out << ".vtp";
    else
      out << ".vtu";
    out << "'... " << std::flush;
    typedef Dune::VTKWriter< typename AnsatzSpaceType::GridViewType > VTKWriterType;
    VTKWriterType vtkWriter(ansatzSpace_->gridView());
    vtkWriter.addVertexData(discreteFunction);
    vtkWriter.write(filename);
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;

  typedef Dune::PDELab::DiscreteGridFunction<GridFunctionSpace,TrialVectorType> DiscreteGridFunction;
  DiscreteGridFunction discreteGridFunction(gridFunctionSpace_,solution);
  Dune::VTKWriter<GridViewType> vtkwriter(gridView_,Dune::VTK::conforming);
  vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<DiscreteGridFunction>(discreteGridFunction,"solution"));
  vtkwriter.write("dune-pdelab",Dune::VTK::appendedraw);
  }

  const Dune::shared_ptr< const MatrixBackendType > systemMatrix() const
  {
    assert(initialized_);
    return matrix_;
  }

  Dune::shared_ptr< MatrixBackendType > systemMatrix()
  {
    assert(initialized_);
    return matrix_;
  }

  const Dune::shared_ptr< const VectorBackendType > rightHandSide() const
  {
    assert(initialized_);
    return rhs_;
  }

  Dune::shared_ptr< VectorBackendType > rightHandSide()
  {
    assert(initialized_);
    return rhs_;
  }

  const Dune::shared_ptr< const PatternType > pattern() const
  {
    assert(initialized_);
    return pattern_;
  }

#endif // visualize, ...


private:
  const Dune::shared_ptr< const ModelType > model_;
  const Dune::shared_ptr< const GridViewType > gridView_;
  const Dune::shared_ptr< const BoundaryInfoType > boundaryInfo_;
  bool initialized_;
//  Dune::shared_ptr< const PatternType > pattern_;
  Dune::shared_ptr< MatrixContainer> matrix_;
  Dune::shared_ptr< VectorBackendType > rhs_;
  Dune::shared_ptr< TrialVectorType > solution_;
  Dune::shared_ptr< DiscreteTrialFunction > residual_;
  Dune::shared_ptr< FiniteElementMap > finiteElementMap_;
  Dune::shared_ptr< GridFunctionSpace > gridFunctionSpace_;
  Dune::shared_ptr< GridOperator> gridOperator_;

}; // class DunePdelab

} // namespace FiniteVolume

} // namespace Elliptic

} // namespace Linear

} // namespace Stationary

} // namespace Solvers

} // namespace Detailed

} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_FINITEVOLUME_DUNE_PDELAB_HH
