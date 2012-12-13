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
#include <dune/common/shared_ptr.hh>
#include <dune/common/fvector.hh>
#include <dune/common/static_assert.hh>

// dune-geometry
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>

// dune-grid
#include <dune/grid/io/file/vtk/vtuwriter.hh>

// dune-pdelab
#include <dune/pdelab/finiteelementmap/p0fem.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/backend/eigenvectorbackend.hh>
#include <dune/pdelab/backend/eigensolverbackend.hh>
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

template< class ModelImp, class GridViewImp, class BoundaryInfoImp, int polynomialOrder >
class DunePdelab
{
public:
  typedef ModelImp ModelType;
  typedef GridViewImp GridViewType;
  typedef BoundaryInfoImp BoundaryInfoType;
  static const int polOrder = polynomialOrder;
  typedef std::map<unsigned int, std::set<unsigned int> > PatternType;

  typedef DunePdelab< ModelType, GridViewType, BoundaryInfoType, polOrder > ThisType;

  static const int dimDomain = ModelType::dimDomain;
  static const int dimRange = ModelType::dimRange;
private:
  // Choose domain and range field type
  typedef typename GridViewType::Grid::ctype CoordinateType;
  typedef double RangeFieldType;

  // Make grid function space
  typedef Dune::PDELab::P0LocalFiniteElementMap<CoordinateType,RangeFieldType,dimDomain> FiniteElementMapType;
  typedef Dune::PDELab::NoConstraints ConstraintsType;
  typedef Dune::PDELab::EigenVectorBackend PdelabVectorBackendType;
  typedef Dune::PDELab::GridFunctionSpace<GridViewType,FiniteElementMapType,ConstraintsType,PdelabVectorBackendType>
                                                                                                  DiscreteFunctionSpaceType;
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
      DomainType center = entity.geometry().center();
      RangeType b = 0.0;
      RangeType f = 1.0;
      localModel_.force()->evaluate(center,f);

      residual.accumulate(trialFunction,0,(b*trialVector(trialFunction,0)-f)*entity.geometry().volume());
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
          residual_s.accumulate(trialFunction_s,0,face_volume);
          return;
      }
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
  typedef Dune::PDELab::GridOperator<DiscreteFunctionSpaceType,DiscreteFunctionSpaceType,LocalOperatorType,
                                     PdelabMatrixBackendType,RangeFieldType,RangeFieldType,RangeFieldType,
                                     ConstraintsContainerType,ConstraintsContainerType> GridOperatorType;
  typedef typename GridOperatorType::template MatrixContainer<RangeFieldType>::Type MatrixBackendType;
  typedef typename GridOperatorType::Traits::TrialGridFunctionSpace TrialDiscreteFunctionSpaceType;
  typedef Dune::PDELab::EigenBackend_BiCGSTAB_Diagonal LinearSolverType;
public:
  typedef typename GridOperatorType::Traits::Domain VectorBackendType;
  typedef Dune::PDELab::DiscreteGridFunction<DiscreteFunctionSpaceType,VectorBackendType> DiscreteFunctionType;


public:
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
      finiteElementMap_ = Dune::shared_ptr< FiniteElementMapType >
              (new FiniteElementMapType(Dune::GeometryType(gridView_->template begin< 0 >()->geometry().type())));
      gridFunctionSpace_ = Dune::shared_ptr< DiscreteFunctionSpaceType >(new DiscreteFunctionSpaceType
                                                                        (*gridView_, *finiteElementMap_));
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      // grid operator (left hand side)
      out << prefix << "initializing operator... " << std::flush;
      timer.reset();
      LocalOperatorType localOperator(*model_, *boundaryInfo_);
      gridOperator_ = Dune::shared_ptr< GridOperatorType >(new GridOperatorType(*gridFunctionSpace_,
                                                                                *gridFunctionSpace_,localOperator));
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
    Dune::shared_ptr<VectorBackendType> vector = Dune::shared_ptr<VectorBackendType>
                                                    (new VectorBackendType(*gridFunctionSpace_,0.0));
    return vector;
  }

  Dune::shared_ptr< DiscreteFunctionType > createDiscreteFunction() const
  {
    assert(initialized_ && "Please call init() before calling createDiscreteFunction()!");
    VectorBackendType vector(*gridFunctionSpace_,0.0);
    Dune::shared_ptr< DiscreteFunctionType > discreteGridFunction = Dune::shared_ptr<DiscreteFunctionType>
                                (new DiscreteFunctionType(*gridFunctionSpace_,vector));
    return discreteGridFunction;
  } // ... createDiscreteFunction()

  Dune::shared_ptr< DiscreteFunctionType > createDiscreteFunction(VectorBackendType& vector) const
  {
    assert(initialized_ && "Please call init() before calling createDiscreteFunction()!");
    Dune::shared_ptr< DiscreteFunctionType > discreteGridFunction = Dune::shared_ptr<DiscreteFunctionType>
                                (new DiscreteFunctionType(*gridFunctionSpace_,vector));
    return discreteGridFunction;
  } // ... createDiscreteFunction()


  void solve(VectorBackendType& solution,
             const Dune::ParameterTree paramTree = Dune::ParameterTree(),
             const unsigned int maxIter = 500,
             const double precision = 1e-12,
             const std::string prefix = "",
             std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
  {
    // preparations
    assert(initialized_ && "The system can only be solved after init() has been called! ");
    // timer
    Dune::Timer timer;

    // Assemble solver
    LinearSolverType linearSolver(maxIter);

    typename VectorBackendType::ElementType defect = linearSolver.norm(*rhs_);
    typename VectorBackendType::ElementType mindefect = 1e-99;

    // compute correction
    timer.reset();
    VectorBackendType update(gridOperator_->trialGridFunctionSpace(),0.0);
    typename VectorBackendType::ElementType red = std::min(precision,defect/(mindefect));

    // (4) solve matrix_ * update = residual = -*rhs_ for update
    linearSolver.apply(*matrix_,update,*rhs_,red); // solver makes right hand side consistent

    // (5) and store the solution u = u0 + update = update
    solution = update;
    out << prefix << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void solve()

  void visualize(VectorBackendType& vector,
                 const std::string filename = id() + "solution",
                 const std::string name = "solution",
                 const std::string prefix = "",
                 std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
  {
    Dune::shared_ptr< DiscreteFunctionType > discreteGridFunction = Dune::shared_ptr<DiscreteFunctionType>
                                (new DiscreteFunctionType(*gridFunctionSpace_,vector));
    this->visualize(discreteGridFunction, filename, name, prefix, out);
  } // void visualize()


  void visualize(Dune::shared_ptr< DiscreteFunctionType > discreteGridFunction,
                 const std::string filename = id() + "solution",
                 const std::string name = "solution",
                 const std::string prefix = "",
                 std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
  {
    // preparations
    assert(initialized_ && "A discrete function can only be visualized after init() has been called! ");
    Dune::Timer timer;

    out << prefix << "writing '" << name << "' to '" << filename;
    if (dimDomain == 1)
      out << ".vtp";
    else
      out << ".vtu";
    out << "'... " << std::flush;

    Dune::VTKWriter<GridViewType> vtkwriter(*gridView_,Dune::VTK::conforming);
    vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<DiscreteFunctionType>(*discreteGridFunction,name));
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

  Dune::shared_ptr< DiscreteFunctionSpaceType > discreteFunctionSpace() const
  {
    assert(initialized_);
    return gridFunctionSpace_;
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
  Dune::shared_ptr< DiscreteFunctionSpaceType > gridFunctionSpace_;
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
