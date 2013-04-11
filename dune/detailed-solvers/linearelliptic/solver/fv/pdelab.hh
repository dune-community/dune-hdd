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
#include <dune/stuff/la/container/affineparametric.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/function/expression.hh>
#include <dune/stuff/la/solver.hh>

#include "../../model/interface.hh"
#include "../interface.hh"

namespace Dune {
namespace DetailedSolvers {
namespace LinearElliptic {


// forward of the solver, to be used in the traits and allow for specialization
template< class GridPartImp, class RangeFieldImp, int rangeDim >
class SolverFiniteVolumePdelab
{
public:
  SolverFiniteVolumePdelab() = delete;
};

/**
 *  \brief  Traits for SolverFiniteVolumePdelab
 */
template< class GridPartImp, class RangeFieldImp, int rangeDim >
class SolverFiniteVolumePdelabTraits
{
public:
  typedef SolverFiniteVolumePdelab< GridPartImp, RangeFieldImp, rangeDim > derived_type;
  typedef typename GridPartImp::Traits  GridPartTraits;
  static const int                      polOrder = 0;
  typedef RangeFieldImp                 RangeFieldType;
  static const int                      dimRange = rangeDim;
  typedef Dune::Stuff::LA::Container::EigenDenseVector< RangeFieldType > VectorType;
}; // class SolverFiniteVolumePdelabTraits

/**
 *  \brief  Solver of linear elliptic pdes using a finite volume discretization provided by dune-pdelab
 */
template< class GridPartImp, class RangeFieldImp >
class SolverFiniteVolumePdelab< GridPartImp, RangeFieldImp, 1 >
  : public SolverInterface< SolverFiniteVolumePdelabTraits< GridPartImp, RangeFieldImp, 1 > >
  , public SolverParametricInterface< SolverFiniteVolumePdelabTraits< GridPartImp, RangeFieldImp, 1 > >
{
public:
  typedef SolverFiniteVolumePdelab< GridPartImp, RangeFieldImp, 1 >       ThisType;
  typedef SolverFiniteVolumePdelabTraits< GridPartImp, RangeFieldImp, 1 > Traits;
  typedef SolverInterface< Traits >                                       BaseType;
  typedef SolverParametricInterface< Traits >                             ParametricBaseType;
  typedef Dune::grid::Part::Interface< typename Traits::GridPartTraits >  GridPartType;

private:
  typedef typename GridPartType::GridViewType GridViewType;

public:
  typedef typename GridPartType::ctype    DomainFieldType;
  static const int                        dimDomain = GridPartType::dimension;
  typedef typename Traits::RangeFieldType RangeFieldType;
  static const int                        dimRange = Traits::dimRange;
  static const int                        polOrder = Traits::polOrder;

  typedef Dune::Stuff::GridboundaryInterface< GridViewType >                      BoundaryInfoType;
  typedef ModelInterface< DomainFieldType, dimDomain, RangeFieldType, dimRange >  ModelType;
  typedef typename ModelType::ParamType                                           ParamType;
  typedef Dune::Stuff::LA::Container::EigenRowMajorSparseMatrix< RangeFieldType > MatrixType;
  typedef typename Traits::VectorType                                             VectorType;

private:
  typedef Dune::Stuff::LA::Container::AffineParametric< MatrixType > AffineParametricMatrixType;
  typedef Dune::Stuff::LA::Container::AffineParametric< VectorType > AffineParametricVectorType;

  typedef Dune::PDELab::P0LocalFiniteElementMap< DomainFieldType, RangeFieldType, dimDomain > FiniteElementMapType;
  typedef Dune::PDELab::NoConstraints                                                         ConstraintsType;
  typedef Dune::PDELab::EigenVectorBackend                                                    PdelabVectorBackendType;

public:
  typedef Dune::PDELab::GridFunctionSpace< GridViewType, FiniteElementMapType, ConstraintsType, PdelabVectorBackendType >
                                                                                                  AnsatzSpaceType;

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

    // volume integral depending on test and trial/ansatz functions
    template<typename EntityType, typename TrialFunctionSpaceType, typename TrialVector,
                typename TestFunctionSpaceType, typename ResidualType>
    void alpha_volume (const EntityType& entity, const TrialFunctionSpaceType& trialFunction,
                       const TrialVector& /*trialVector*/, const TestFunctionSpaceType& /*testFunction*/,
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

    // skeleton integral depending on test and trial/ansatz functions
    // each face is only visited ONCE!
    template<typename IntersectionType, typename TrialFunctionSpaceType, typename TrialVector,
             typename TestFunctionSpaceType, typename ResidualType>
    void alpha_skeleton (const IntersectionType& intersection,
                         const TrialFunctionSpaceType& trialFunction_s, const TrialVector& trialVector_s,
                         const TestFunctionSpaceType& /*testFunction_s*/,
                         const TrialFunctionSpaceType& trialFunction_n, const TrialVector& trialVector_n,
                         const TestFunctionSpaceType& /*testFunction_n*/,
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
      RangeType a = localModel_.diffusion()->evaluate(inside_global);
      residual_s.accumulate(trialFunction_s,0,-a*(trialVector_n(trialFunction_n,0)-trialVector_s(trialFunction_s,0))
                                                                                                *face_volume/distance);
      residual_n.accumulate(trialFunction_n,0, a*(trialVector_n(trialFunction_n,0)-trialVector_s(trialFunction_s,0))
                                                                                                *face_volume/distance);
    } // void alpha_skeleton()

    // skeleton integral depending on test and trial/ansatz functions
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
  typedef EllipticLocalOperator                  LocalOperatorType;
  typedef PdelabVectorBackendType::MatrixBackend PdelabMatrixBackendType;
  typedef Dune::PDELab::EmptyTransformation      ConstraintsContainerType;

  typedef Dune::PDELab::GridOperator< AnsatzSpaceType, AnsatzSpaceType, LocalOperatorType,
                                      PdelabMatrixBackendType, DomainFieldType, RangeFieldType, RangeFieldType,
                                      ConstraintsContainerType, ConstraintsContainerType > GridOperatorType;

  typedef typename GridOperatorType::Traits::Domain                                   PdelabVectorType;
  typedef typename GridOperatorType::template MatrixContainer< RangeFieldType >::Type PdelabMatrixType;

  typedef Dune::PDELab::EigenBackend_BiCGSTAB_Diagonal LinearSolverType;

  typedef Dune::PDELab::DiscreteGridFunction< AnsatzSpaceType, PdelabVectorType > DiscreteFunctionType;

public:
  SolverFiniteVolumePdelab(const Dune::shared_ptr< const GridPartType > _gridPart,
                           const Dune::shared_ptr< const BoundaryInfoType > _boundaryInfo,
                           const Dune::shared_ptr< const ModelType > _model)
    : gridPart_(_gridPart)
    , gridView_(new GridViewType(gridPart_->gridView()))
    , boundaryInfo_(_boundaryInfo)
    , model_(_model)
    , initialized_(false)
  {
    // sanity checks
    std::stringstream msg;
    unsigned int throw_up = 0;
    // * integration orders
    if (model_->diffusion()->order() < 0) {
      msg << "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
          << " negative integration order given for the diffusion!";
      ++throw_up;
    }
    if (model_->force()->order() < 0) {
      msg << "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
          << " negative integration order given for the force!";
      ++throw_up;
    }
    if (model_->neumann()->order() < 0) {
      msg << "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
          << " negative integration order given for the neumann values!";
      ++throw_up;
    }
    // * parametrization
    if (model_->parametric() && !model_->affineparametric()) {
      msg << "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
          << " only implemented for nonparametric or affineparametric models!";
      ++throw_up;
    }
    if (throw_up)
      DUNE_THROW(Dune::InvalidStateException, msg.str());
    // function spaces
    finiteElementMap_ = Dune::shared_ptr< FiniteElementMapType >(new FiniteElementMapType(
        Dune::GeometryType(gridView_->template begin< 0 >()->geometry().type())));
    ansatzSpace_ = Dune::shared_ptr< AnsatzSpaceType >(new AnsatzSpaceType(*gridView_, *finiteElementMap_));
  }

  static const std::string id()
  {
    return BaseType::id() + ".fv.pdelab";
  }

  bool initialized() const
  {
    return initialized_;
  }

  const Dune::shared_ptr< const ModelType > model() const
  {
    return model_;
  }

  const Dune::shared_ptr< const GridPartType > gridPart() const
  {
    return gridPart_;
  }

  const Dune::shared_ptr< const GridViewType > gridView() const
  {
    return gridView_;
  }

  const Dune::shared_ptr< const BoundaryInfoType > boundaryInfo() const
  {
    return boundaryInfo_;
  }

  void init(std::ostream& out = Dune::Stuff::Common::Logger().debug(), const std::string prefix = "")
  {
    if (!initialized_) {
      Dune::Timer timer;
      // check
      if(model_->dirichlet()->parametric() || model_->force()->parametric() || model_->neumann()->parametric())
        DUNE_THROW(Dune::NotImplemented,
                   "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                   << " only implemented for nonparametric force, dirichlet and neumann functions!");
      typedef Dune::DetailedSolvers::LinearElliptic::
                ModelDefault<DomainFieldType, dimDomain, RangeFieldType, dimRange> DefaultModelType;
      typedef Dune::Stuff::FunctionExpression<DomainFieldType, dimDomain, RangeFieldType, dimRange> zeroFunctionType;
      auto zero = Dune::shared_ptr< zeroFunctionType >(new zeroFunctionType("x", "0.0", 0));
      // start vector initialGuess
      PdelabVectorType initialGuess(*ansatzSpace_, 0.0);

      out << prefix << "initializing matrices and vectors:" << std::endl;
      // * diffusion matrices and dirichlet vectors
      std::shared_ptr< AffineParametricMatrixType > diffusionMatrix;
      std::shared_ptr< AffineParametricVectorType > dirichletVector;
      // if the diffusion is not parametric
      if (!model_->diffusion()->parametric()) {
        out << prefix << "  1 diffusion matrix and 1 dirichlet vector...    " << std::flush;
        timer.reset();
        std::shared_ptr< MatrixType > diffusionMatrices(std::make_shared<MatrixType>());
        std::shared_ptr< VectorType > dirichletVectors(std::make_shared<VectorType>());
        // force and neumann are always zero
        DefaultModelType model(model_->diffusion(), zero, model_->dirichlet(), zero);
        auto localOperator = Dune::shared_ptr< LocalOperatorType >(new LocalOperatorType(model, *boundaryInfo_));
        auto gridOperator  = Dune::shared_ptr< GridOperatorType >(new GridOperatorType(*ansatzSpace_,
                                                                                       *ansatzSpace_,
                                                                                       *localOperator));
        auto matrix = Dune::shared_ptr< PdelabMatrixType >(new PdelabMatrixType(*gridOperator));
        *matrix = 0.0;
        PdelabVectorType residual(*ansatzSpace_, 0.0);
        auto rhs = Dune::shared_ptr< PdelabVectorType >(new PdelabVectorType(*ansatzSpace_, 0.0));

        // assemble jacobian system matrix = \nabla R(initialGuess) (with residualoperator R)
        gridOperator->jacobian(initialGuess, *matrix);
        // compute residual R(initialGuess)
        gridOperator->residual(initialGuess, residual);
        *rhs -= residual;
        // store matrix and rhs
        diffusionMatrices->backend() = matrix->base();
        dirichletVectors->backend() = rhs->base();
        diffusionMatrix = std::make_shared< AffineParametricMatrixType >(diffusionMatrices);
        dirichletVector = std::make_shared< AffineParametricVectorType >(dirichletVectors);
        out << "done (took " << timer.elapsed() << " sec)" << std::endl;
      } else {
        // we are separable (see constructor), loop over all components
        out << prefix << "  " << model_->diffusion()->numComponents()
                      << " diffusion matrices and dirichlet vectors...   " << std::flush;
        timer.reset();
        std::vector< std::shared_ptr< MatrixType > > diffusionMatrices;
        std::vector< std::shared_ptr< VectorType > > dirichletVectors;
        for (size_t qq = 0; qq < model_->diffusion()->numComponents(); ++qq) {
          // one model for each component (force and neumann are always zero)
          DefaultModelType model_q(model_->diffusion()->components()[qq], zero, model_->dirichlet(), zero);
          auto localOperator = Dune::shared_ptr< LocalOperatorType >(new LocalOperatorType(model_q, *boundaryInfo_));
          auto gridOperator  = Dune::shared_ptr< GridOperatorType >(new GridOperatorType(*ansatzSpace_,
                                                                                    *ansatzSpace_,
                                                                                    *localOperator));
          auto matrix = Dune::shared_ptr< PdelabMatrixType >(new PdelabMatrixType(*gridOperator));
          *matrix = 0.0;
          PdelabVectorType residual(*ansatzSpace_, 0.0);
          auto rhs = Dune::shared_ptr< PdelabVectorType >(new PdelabVectorType(*ansatzSpace_, 0.0));

          // assemble jacobian system matrix = \nabla R(initialGuess) (with residualoperator R)
          gridOperator->jacobian(initialGuess, *matrix);
          // compute residual R(initialGuess)
          gridOperator->residual(initialGuess, residual);
          *rhs -= residual;
          // store matrix and rhs
          diffusionMatrices.push_back(std::make_shared<MatrixType>());
          diffusionMatrices[qq]->backend() = matrix->base();
          dirichletVectors.push_back(std::make_shared<VectorType>());
          dirichletVectors[qq]->backend() = rhs->base();
          } // loop over all components
          diffusionMatrix = std::make_shared< AffineParametricMatrixType >(model_->diffusion()->paramSize(),
                                                                           diffusionMatrices,
                                                                           model_->diffusion()->coefficients());
          dirichletVector = std::make_shared< AffineParametricVectorType >(model_->diffusion()->paramSize(),
                                                                           dirichletVectors,
                                                                           model_->diffusion()->coefficients());
          out << "done (took " << timer.elapsed() << " sec)" << std::endl;
      } // if (!model_->diffusion()->parametric())

      // * force vector
      out << prefix << "  1 force vector...    " << std::flush;
      timer.reset();
      std::shared_ptr< AffineParametricVectorType > forceVector;
      std::shared_ptr< VectorType > forceVectors(std::make_shared<VectorType>());
      // diffusion, dirichlet and neumann are always zero
      DefaultModelType model_f(zero, model_->force(), zero, zero);
      auto localOperator_f = Dune::shared_ptr< LocalOperatorType >(new LocalOperatorType(model_f, *boundaryInfo_));
      auto gridOperator_f  = Dune::shared_ptr< GridOperatorType >(new GridOperatorType(*ansatzSpace_,
                                                                                       *ansatzSpace_,
                                                                                       *localOperator_f));
      PdelabVectorType residual(*ansatzSpace_, 0.0);
      auto rhs = Dune::shared_ptr< PdelabVectorType >(new PdelabVectorType(*ansatzSpace_, 0.0));
      // compute residual R(initialGuess)
      gridOperator_f->residual(initialGuess, residual);
      *rhs -= residual;
      forceVectors->backend() = rhs->base();
      forceVector = std::make_shared< AffineParametricVectorType >(forceVectors);
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      // * neumann vector
      out << prefix << "  1 neumann vector...    " << std::flush;
      timer.reset();
      std::shared_ptr< AffineParametricVectorType > neumannVector;
      std::shared_ptr< VectorType > neumannVectors(std::make_shared<VectorType>());
      // diffusion, force and dirichlet are always zero
      DefaultModelType model_gN(zero, zero, zero, model_->neumann());
      auto localOperator_gN = Dune::shared_ptr< LocalOperatorType >(new LocalOperatorType(model_gN, *boundaryInfo_));
      auto gridOperator_gN  = Dune::shared_ptr< GridOperatorType >(new GridOperatorType(*ansatzSpace_,
                                                                                        *ansatzSpace_,
                                                                                        *localOperator_gN));
      residual=0.0;
      *rhs = 0.0;
      // compute residual R(initialGuess)
      gridOperator_gN->residual(initialGuess, residual);
      *rhs -= residual;
      neumannVectors->backend() = rhs->base();
      neumannVector = std::make_shared< AffineParametricVectorType >(neumannVectors);
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      vectors_.insert(std::make_pair("dirichlet", dirichletVector));
      matrices_.insert(std::make_pair("diffusion", diffusionMatrix));
      vectors_.insert(std::make_pair("force", forceVector));
      vectors_.insert(std::make_pair("neumann", neumannVector));

      // done
      initialized_ = true;
    } // if !(initialized_)
  } // void init()

  Dune::shared_ptr< VectorType > createVector() const
  {
    assert(initialized_ && "Please call init() before calling createVector()!");
    //!TODO don't use naked refs
    PdelabVectorType pdelabVector(*ansatzSpace_, 0.0);
    Dune::shared_ptr< VectorType > vector(new VectorType(pdelabVector));
    return vector;
  } // Dune::shared_ptr< VectorType > createVector() const

private:
  void generic_solve(std::shared_ptr< VectorType >& solutionVector,
                     const ParamType& mu,
                     const std::string& linearSolverType,
                     const size_t& linearSolverMaxIter,
                     const double linearSolverPrecision,
                     const std::string prefix,
                     std::ostream& out) const
  {
    assert(initialized_ && "Please call init() before calling solve()!");
    // first of all, get the corect parameters (the model returns empty ones for nonparametric functions)
    const ParamType muDiffusion = model_->mapParam(mu, "diffusion");
    const ParamType muForce = model_->mapParam(mu, "force");
    const ParamType muDirichlet = model_->mapParam(mu, "dirichlet");
    const ParamType muNeumann = model_->mapParam(mu, "neumann");
    assert(matrices_.find("diffusion") != matrices_.end());
    assert(vectors_.find("force") != vectors_.end());
    assert(vectors_.find("neumann") != vectors_.end());
    assert(vectors_.find("dirichlet") != vectors_.end());
    const AffineParametricMatrixType& diffusionMatrix = *(matrices_.find("diffusion")->second);
    const AffineParametricVectorType& forceVector = *(vectors_.find("force")->second);
    const AffineParametricVectorType& neumannVector = *(vectors_.find("neumann")->second);
    const AffineParametricVectorType& dirichletVector = *(vectors_.find("dirichlet")->second);
    Dune::Timer timer;
    out << prefix << "computing system matrix...   " << std::flush;
    std::shared_ptr< MatrixType > systemMatrix = diffusionMatrix.fix(muDiffusion);
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
    out << prefix << "computing right hand side... " << std::flush;
    timer.reset();
    VectorType rightHandSide;
    rightHandSide.backend() = forceVector.fix(muForce)->backend() + neumannVector.fix(muNeumann)->backend()
                                  + dirichletVector.fix(muDiffusion)->backend();
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;

    out << prefix << "solving linear system (of size " << systemMatrix->rows()
        << "x" << systemMatrix->cols() << ")... " << std::endl;
    out << prefix << "  using '" << linearSolverType << "'... " << std::flush;
    timer.reset();
    typedef typename Dune::Stuff::LA::Solver::Interface< MatrixType, VectorType > SolverType;
    const std::shared_ptr< const SolverType >
                                   solver(Dune::Stuff::LA::Solver::create< MatrixType, VectorType >(linearSolverType));
    const unsigned int failure = solver->apply(*systemMatrix,
                                               rightHandSide,
                                               *solutionVector,
                                               linearSolverMaxIter,
                                               linearSolverPrecision);
    if (failure)
      DUNE_THROW(Dune::MathError,
                 "\n"
                 << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " linear solver '" << linearSolverType << "' reported error code " << failure << "!\n"
                 << "  1: did not converge\n"
                 << "  2: had numerical issues\n"
                 << "  3: dude, I have no idea");
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void generic_solve(...)

public:
  void solve(std::shared_ptr< VectorType > solutionVector,
             const std::string linearSolverType = "bicgstab.ilut",
             const double linearSolverPrecision = 1e-12,
             const size_t linearSolverMaxIter = 5000,
             std::ostream& out = Dune::Stuff::Common::Logger().devnull(),
             const std::string prefix = "") const
  {
    if (!initialized_)
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " call init() before calling solve()!");
    // check, that we are really in the nonparametric setting!
    if (model_->parametric())
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " nonparametric solve() called for a parametric model!");
    generic_solve(solutionVector,
                  ParamType(),
                  linearSolverType, linearSolverMaxIter, linearSolverPrecision,
                  prefix, out);
  } // ... solve(...)

  void solve(std::shared_ptr< VectorType > solutionVector,
             const ParamType& mu,
             const std::string linearSolverType = "bicgstab.ilut",
             const double linearSolverPrecision = 1e-12,
             const size_t linearSolverMaxIter = 5000,
             std::ostream& out = Dune::Stuff::Common::Logger().devnull(),
             const std::string prefix = "") const
  {
    if (!initialized_)
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " call init() before calling solve()!");
    // check, that we are really in the parametric setting!
    if (!model_->parametric())
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " parametric solve() called for a nonparametric model!");
    generic_solve(solutionVector,
                  mu,
                  linearSolverType, linearSolverMaxIter, linearSolverPrecision,
                  prefix, out);
  } // ... solve(..., mu, ...)

  void visualize(const std::shared_ptr< const VectorType >& vector,
                 const std::string filename = id() + ".vector",
                 const std::string name = id() + ".vector",
                 std::ostream& out = Dune::Stuff::Common::Logger().debug(),
                 const std::string prefix = "") const
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
    PdelabVectorType pdelabVector(*ansatzSpace_, 0.0);
    pdelabVector.base() = vector->backend();
    DiscreteFunctionType discreteFunction(*ansatzSpace_, pdelabVector);

    typedef Dune::VTKWriter< GridViewType > VTKWriterType;
    VTKWriterType vtkWriter(*gridView_, Dune::VTK::conforming);
    vtkWriter.addCellData(std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DiscreteFunctionType > >
                                                                                          (discreteFunction, name));
    vtkWriter.write(filename, Dune::VTK::appendedraw);
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void visualize(...)

private:
  const Dune::shared_ptr< const GridPartType > gridPart_;
  const Dune::shared_ptr< const GridViewType > gridView_;
  const Dune::shared_ptr< const BoundaryInfoType > boundaryInfo_;
  const Dune::shared_ptr< const ModelType > model_;
  bool initialized_;
  Dune::shared_ptr< FiniteElementMapType > finiteElementMap_;
  Dune::shared_ptr< AnsatzSpaceType > ansatzSpace_;
  std::map< const std::string, std::shared_ptr< AffineParametricMatrixType > > matrices_;
  std::map< const std::string, std::shared_ptr< AffineParametricVectorType > > vectors_;
}; // class SolverFiniteVolumePdelab


} // namespace LinearElliptic
} // namespace DetailedSolvers
} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_FINITEVOLUME_DUNE_PDELAB_HH
