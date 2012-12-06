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
private:
  typedef ModelImp ModelType;
  typedef GridViewImp GridViewType;
  typedef BoundaryInfoImp BoundaryInfoType;
  static const int polOrder = polynomialOrder;
  typedef std::map<unsigned int, std::set<unsigned int> > PatternType;

  typedef DunePdelab< ModelType, GridViewType, BoundaryInfoType, polOrder > ThisType;

  static const int dimDomain = ModelType::dimDomain;
  static const int dimRange = ModelType::dimRange;

  // Choose domain and range field type
  typedef typename GridViewType::Grid::ctype CoordinateType;
  typedef double RangeFieldType;

  // Make grid function space
  typedef Dune::PDELab::P0LocalFiniteElementMap<CoordinateType,RangeFieldType,dimDomain> FiniteElementMapType;
  typedef Dune::PDELab::NoConstraints ConstraintsType;
  typedef Dune::PDELab::EigenVectorBackend PdelabVectorBackendType;
  typedef Dune::PDELab::GridFunctionSpace<GridViewType,FiniteElementMapType,ConstraintsType,PdelabVectorBackendType> GridFunctionSpaceType;
public:
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


    // constructor
    EllipticLocalOperator(const ModelType& model, const BoundaryInfoType& boundaryInfo)
        : localModel_(model), localBoundaryInfo_(boundaryInfo)
    {}

    // volume integral depending on test and ansatz functions
    template<typename EntityType, typename TrialFunctionSpaceType, typename TrialVector, typename TestFunctionSpaceType, typename ResidualType>
    void alpha_volume (const EntityType& entity, const TrialFunctionSpaceType& trialFunction, const TrialVector& trialVector,
                       const TestFunctionSpaceType& testFunction, ResidualType& residual) const
    {
      // domain and range field type
      typedef typename TrialFunctionSpaceType::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::DomainFieldType DomainFieldType;
      typedef  Dune::FieldVector<DomainFieldType,dimDomain> DomainType;
      typedef  Dune::FieldVector<RangeFieldType,dimRange> RangeType;

      // evaluate reaction term
      DomainType center = entity.geometry().center();
      RangeType b = 0.0;
      RangeType f = 1.0;
      localModel_.force()->evaluate(center,f);

      residual.accumulate(trialFunction,0,(b*trialVector(trialFunction,0)-f)*entity.geometry().volume());
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
      residual_s.accumulate(trialFunction_s,0,-a*(trialVector_n(trialFunction_n,0)-trialVector_s(trialFunction_s,0))*face_volume/distance);
      residual_n.accumulate(trialFunction_n,0, a*(trialVector_n(trialFunction_n,0)-trialVector_s(trialFunction_s,0))*face_volume/distance);
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
      typedef  Dune::FieldVector<RangeFieldType,dimRange> RangeType;

      // face geometry
      Dune::FieldVector<DomainFieldType,dimDomain> face_center = intersection.geometry().center();
      RangeType face_volume = intersection.geometry().volume();

      // do things depending on boundary condition type
      if (localBoundaryInfo_.dirichlet(intersection))
      {
          RangeType g = 0.0;
          localModel_.dirichlet()->evaluate(face_center,g);
          RangeType a = 1.0;
          localModel_.diffusion()->evaluate(face_center,a);
          Dune::FieldVector<DomainFieldType,dimDomain> inside_global = intersection.inside()->geometry().center();
          inside_global -= face_center;
          RangeType distance = inside_global.two_norm();
          residual_s.accumulate(trialFunction_s,0,-a*(g-trialVector_s(trialFunction_s,0))*face_volume/distance);
          return;
      }
      else if (localBoundaryInfo_.neumann(intersection))
      {
          RangeType j; if (face_center[1]<0.5) j = 1.0; else j = -1.0;
          residual_s.accumulate(trialFunction_s,0,j*face_volume);
          return;
      }


/*
      // evaluate boundary condition type
      int boundaryType;
//      if (face_center[0]>1.0-1e-6)
//        boundaryType = 0; // Neumann
//      else
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
          RangeType g = 0.0;
          localModel_.dirichlet()->evaluate(face_center,g);
          RangeType a = 1.0;
          localModel_.diffusion()->evaluate(face_center,a);
          Dune::FieldVector<DomainFieldType,dimDomain> inside_global = intersection.inside()->geometry().center();
          inside_global -= face_center;
          RangeType distance = inside_global.two_norm();
          residual_s.accumulate(trialFunction_s,0,-a*(g-trialVector_s(trialFunction_s,0))*face_volume/distance);
          return;
        }
        */
    } // void alpha_boundary()

  private:
    const ModelType& localModel_;
    const BoundaryInfoType& localBoundaryInfo_;

  }; // class EllipticLocalOperator

private:
  // Make grid operator (left hand side)
  typedef EllipticLocalOperator LocalOperatorType;
  typedef PdelabVectorBackendType::MatrixBackend PdelabMatrixBackendType;
  typedef Dune::PDELab::EmptyTransformation ConstraintsContainerType;
  typedef Dune::PDELab::GridOperator<GridFunctionSpaceType,GridFunctionSpaceType,LocalOperatorType,PdelabMatrixBackendType,
                                     RangeFieldType,RangeFieldType,RangeFieldType,ConstraintsContainerType,
                                     ConstraintsContainerType> GridOperatorType;

  typedef typename GridOperatorType::template MatrixContainer<RangeFieldType>::Type MatrixBackendType;

  typedef typename GridOperatorType::Traits::TrialGridFunctionSpace TrialGridFunctionSpaceType;

  typedef Dune::PDELab::EigenBackend_BiCGSTAB_Diagonal LinearSolverType;
public:
  typedef typename GridOperatorType::Traits::Domain VectorBackendType;
  typedef Dune::PDELab::DiscreteGridFunction<GridFunctionSpaceType,VectorBackendType> DiscreteGridFunctionType;


public:
  DunePdelab(const Dune::shared_ptr< const ModelType > model,
                              const Dune::shared_ptr< const GridViewType > gridView,
                              const Dune::shared_ptr< const BoundaryInfoType > boundaryInfo)
    : model_(model)
    , gridView_(gridView)
    , boundaryInfo_(boundaryInfo)
    , initialized_(false)
  {
//      std::string a = *rhs_;

      // ##VectorBackendType##
//       Dune::PDELab::EigenVectorBackend::VectorContainer
//             <
//               Dune::PDELab::GridFunctionSpace
//                <Dune::GridView
//                  <Dune::DefaultLeafGridViewTraits
//                    <const Dune::SGrid<2, 2>, (Dune::PartitionIteratorType)4u>
//                  >, Dune::PDELab::P0LocalFiniteElementMap<double, double, 2>, Dune::PDELab::NoConstraints,
//                        Dune::PDELab::EigenVectorBackend, Dune::PDELab::GridFunctionGeneralMapper
//                >,
//                  double
//             >

      // ##DiscreteGridFunctionType## :
//          Dune::Detailed::Solvers::Stationary::Linear::Elliptic::FiniteVolume::DunePdelab
//               <Dune::Detailed::Solvers::Stationary::Linear::Elliptic::Model::Interface<double, 2, double, 1>
//                  , Dune::GridView
//                    <Dune::DefaultLeafGridViewTraits
//                        <const Dune::SGrid<2, 2>, (Dune::PartitionIteratorType)4u>
//                    >,
//                 Dune::Stuff::Grid::BoundaryInfo::IdBased, 1
//                >::DiscreteGridFunctionType
//          {aka Dune::PDELab::DiscreteGridFunction
//                   <

//                     Dune::PDELab::GridFunctionSpace
//                      <Dune::GridView
//                        <Dune::DefaultLeafGridViewTraits
//                          <const Dune::SGrid<2, 2>, (Dune::PartitionIteratorType)4u>
//                        >, Dune::PDELab::P0LocalFiniteElementMap<double, double, 2>, Dune::PDELab::NoConstraints,
//                            Dune::PDELab::EigenVectorBackend, Dune::PDELab::GridFunctionGeneralMapper
//                      >,

//                          Dune::PDELab::EigenVectorBackend::VectorContainer
//                            <

//                              Dune::PDELab::GridFunctionSpace
//                                <Dune::GridView
//                                  <Dune::DefaultLeafGridViewTraits
//                                     <const Dune::SGrid<2, 2>, (Dune::PartitionIteratorType)4u>
//                                  >, Dune::PDELab::P0LocalFiniteElementMap<double, double, 2>, Dune::PDELab::NoConstraints,
//                                      Dune::PDELab::EigenVectorBackend, Dune::PDELab::GridFunctionGeneralMapper
//                                >,

//                          double
//                      >
//                    >}

         //  kann man den MatrixBackendType auch aus dem PdelabMatrixBackendType bekommen?
         // und was ist mit PdelabVectorBackendType und VectorBackendType ?
         // Pdelab... =..., nur ohne den Container?
  }

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
      finiteElementMap_ = Dune::shared_ptr< FiniteElementMapType >
              (new FiniteElementMapType(Dune::GeometryType(gridView_->template begin< 0 >()->geometry().type())));
      gridFunctionSpace_ = Dune::shared_ptr< GridFunctionSpaceType >(new GridFunctionSpaceType(*gridView_, *finiteElementMap_));
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      // grid operator (left hand side)
      out << prefix << "initializing operator... " << std::flush;
      timer.reset();
      LocalOperatorType localOperator(*model_, *boundaryInfo_);
      gridOperator_ = Dune::shared_ptr< GridOperatorType >(new GridOperatorType(*gridFunctionSpace_,*gridFunctionSpace_,localOperator));
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      // Newton method (5 steps):
      // (1) given a start vector u0=0
      VectorBackendType u0(*gridFunctionSpace_,0.0);

      // (2) assemble jacobian system matrix_ = \nabla R(u0) (with residualoperator R)
      out << prefix << "assembling matrix... " << std::flush;
      timer.reset();
      matrix_ = Dune::shared_ptr< MatrixBackendType >(new MatrixBackendType(*gridOperator_));
      *matrix_ = 0.0;
      gridOperator_->jacobian(u0,*matrix_);
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      // (3) compute the residual R(u0)
      out << prefix << "assembling residual... " << std::flush;
      timer.reset();
      VectorBackendType residual(*gridFunctionSpace_,0.0);
      gridOperator_->residual(u0,residual);
      rhs_ = Dune::shared_ptr< VectorBackendType >(new VectorBackendType(gridOperator_->testGridFunctionSpace(),0.0));
      *rhs_ -= residual;
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;
      // done
      initialized_ = true;
    } // if !(initialized_)
  } // void init()

  Dune::shared_ptr< VectorBackendType > createVector() const
  {
    assert(initialized_ && "A vector can only be created after init() has been called!");
    Dune::shared_ptr<VectorBackendType> vector = Dune::shared_ptr<VectorBackendType> (new VectorBackendType(*gridFunctionSpace_,0.0));
    return vector;
  }

  Dune::shared_ptr< DiscreteGridFunctionType > createDiscreteFunction(const std::string name = "discrete_function") const
  {
    assert(initialized_ && "Please call init() before calling createDiscreteFunction()!");
    VectorBackendType vector(*gridFunctionSpace_,0.0);
    Dune::shared_ptr< DiscreteGridFunctionType > discreteGridFunction = Dune::shared_ptr<DiscreteGridFunctionType> (new DiscreteGridFunctionType(*gridFunctionSpace_,vector));
    return discreteGridFunction;
    // wofür brauche ich hier den string name? unten will ich zwar discFunc->name() aufrufen, aber die Klasse
    // DiscreteGridFunctionType hat gar kein member named 'name'
    // dito für die zweite Methode createDiscreteFunction()

  } // ... createDiscreteFunction()

  Dune::shared_ptr< DiscreteGridFunctionType > createDiscreteFunction(VectorBackendType& vector,
                                                                  const std::string name = "discrete_function") const
  {
    assert(initialized_ && "Please call init() before calling createDiscreteFunction()!");
    Dune::shared_ptr< DiscreteGridFunctionType > discreteGridFunction = Dune::shared_ptr<DiscreteGridFunctionType> (new DiscreteGridFunctionType(*gridFunctionSpace_,vector));
    return discreteGridFunction;
  } // ... createDiscreteFunction()


  void solve(VectorBackendType& solution,
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

    typename VectorBackendType::ElementType defect = linearSolver.norm(*rhs_);
    typename VectorBackendType::ElementType mindefect = 1e-99;
    typename VectorBackendType::ElementType reduction = 1e-10;

    // compute correction
    timer.reset();
    VectorBackendType update(gridOperator_->trialGridFunctionSpace(),0.0);
    typename VectorBackendType::ElementType red = std::min(reduction,defect/(mindefect));

    // (4) solve matrix_*update = residual = -*rhs_ for update
    linearSolver.apply(*matrix_,update,*rhs_,red); // solver makes right hand side consistent

    // (5) and store the solution u = u0 + update = update
    solution = update;
    out << prefix << "done (took " << timer.elapsed() << " sec)" << std::endl;

  } // void solve()


  void visualize(VectorBackendType& vector,
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

    DiscreteGridFunctionType discreteGridFunction(*gridFunctionSpace_,vector);
    Dune::VTKWriter<GridViewType> vtkwriter(*gridView_,Dune::VTK::conforming);
    vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<DiscreteGridFunctionType>(discreteGridFunction,name));
    vtkwriter.write(filename,Dune::VTK::appendedraw);

    out << "done (took " << timer.elapsed() << " sec)" << std::endl;

    //alternativ Aufrufen der anderen visualize-Methode
//    Dune::shared_ptr< DiscreteGridFunctionType > discreteGridFunction = Dune::shared_ptr<DiscreteGridFunctionType> (new DiscreteGridFunctionType(*gridFunctionSpace_,vector));
//    this->visualize(discreteGridFunction, filename, prefix, out);

  } // void visualize()


  void visualize(Dune::shared_ptr< DiscreteGridFunctionType > discreteGridFunction,
                 const std::string filename = "solution",
                 const std::string prefix = "",
                 std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
  {
    // preparations
    assert(initialized_ && "A discrete function can only be visualized after init() has been called! ");
    Dune::Timer timer;

    out << prefix << "writing '" << "name" /*discreteGridFunction->name()*/ << "' to '" << filename;    //es gibt kein discreteGridFunction->name() !?!
    if (dimDomain == 1)
      out << ".vtp";
    else
      out << ".vtu";
    out << "'... " << std::flush;

    Dune::VTKWriter<GridViewType> vtkwriter(*gridView_,Dune::VTK::conforming);
    vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<DiscreteGridFunctionType>(*discreteGridFunction,"solution")); /*discreteGridFunction->name()*/
    vtkwriter.write(filename,Dune::VTK::appendedraw);

    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void visualize()



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
    assert(false && "Implement me!");
    return pattern_;
  }

  Dune::shared_ptr< FiniteElementMapType >finiteElementMap() const
  {
      assert(initialized_);
      return finiteElementMap_;
  }

  Dune::shared_ptr< GridFunctionSpaceType > gridFunctionSpace() const
  {
      assert(initialized_);
      return gridFunctionSpace_;
  }

  Dune::shared_ptr< GridOperatorType > gridOperator() const
  {
      assert(initialized_);
      return gridOperator_;
  }


private:
  const Dune::shared_ptr< const ModelType > model_;
  const Dune::shared_ptr< const GridViewType > gridView_;
  const Dune::shared_ptr< const BoundaryInfoType > boundaryInfo_;
  bool initialized_;
  Dune::shared_ptr< const PatternType > pattern_;
  Dune::shared_ptr< MatrixBackendType> matrix_;
  Dune::shared_ptr< VectorBackendType > rhs_;
  Dune::shared_ptr< FiniteElementMapType > finiteElementMap_;
  Dune::shared_ptr< GridFunctionSpaceType > gridFunctionSpace_;
  Dune::shared_ptr< GridOperatorType> gridOperator_;

}; // class DunePdelab

} // namespace FiniteVolume

} // namespace Elliptic

} // namespace Linear

} // namespace Stationary

} // namespace Solvers

} // namespace Detailed

} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_FINITEVOLUME_DUNE_PDELAB_HH
