
#ifndef DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_BOUNDARY_DUNE_DETAILED_DISCRETIZATIONS_HH
#define DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_BOUNDARY_DUNE_DETAILED_DISCRETIZATIONS_HH

// dune-common
#include <dune/common/shared_ptr.hh>
#include <dune/common/timer.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/exceptions.hh>

// dune-fem
#include <dune/fem/space/common/functionspace.hh>

// dune-detailed-discretizations
#include <dune/detailed/discretizations/evaluation/local/binary/ipdgfluxes.hh>
#include <dune/detailed/discretizations/discreteoperator/local/codim1/boundaryintegral.hh>
#include <dune/detailed/discretizations/assembler/local/codim1/matrix.hh>
#include <dune/detailed/discretizations/assembler/multiscale/boundary.hh>

// dune-stuff
#include <dune/stuff/common/parameter/tree.hh>

namespace Dune {

namespace Detailed {

namespace Solvers {

namespace Stationary {

namespace Linear {

namespace Elliptic {

namespace Boundary {

template< class BoundaryGridPartImp, class BoundaryInfoImp, class LocalSolverImp, class ContainerFactoryImp >
class DuneDetailedDiscretizations
{
public:
  typedef BoundaryGridPartImp BoundaryGridPartType;

  typedef BoundaryInfoImp BoundaryInfoType;

  typedef LocalSolverImp LocalSolverType;

  typedef ContainerFactoryImp ContainerFactory;

  typedef DuneDetailedDiscretizations< BoundaryGridPartType, BoundaryInfoType, LocalSolverType, ContainerFactory > ThisType;

  static const std::string id;

  typedef typename ContainerFactory::SparseMatrixType MatrixBackendType;

  typedef typename ContainerFactory::DenseVectorType VectorBackendType;

private:
  typedef typename LocalSolverType::ModelType ModelType;

  typedef typename ModelType::DomainFieldType DomainFieldType;

  static const int dimDomain = ModelType::dimDomain;

  typedef typename ModelType::RangeFieldType RangeFieldType;

  static const int dimRange = ModelType::dimRange;

public:
  typedef typename LocalSolverType::AnsatzSpaceType::PatternType PatternType;

  DuneDetailedDiscretizations(const Dune::shared_ptr< const BoundaryGridPartType > boundaryGridPart,
                              const Dune::shared_ptr< const BoundaryInfoType > boundaryInfo,
                              const Dune::shared_ptr< const LocalSolverType > localSolver,
                              const Dune::Stuff::Common::ExtendedParameterTree paramTree)
    : boundaryGridPart_(boundaryGridPart)
    , boundaryInfo_(boundaryInfo)
    , localSolver_(localSolver)
    , model_(localSolver_->model())
    , initialized_(false)
  {
    // get penalty factor
    std::string key = "discretization.penaltyFactor";
    paramTree.assertKey(key, id);
    penaltyFactor_ = paramTree.get(key, RangeFieldType(-1.0));
    if (!penaltyFactor_ > 0) {
      std::stringstream msg;
      msg << "Error in " << id << ":"
          << "wrong '" << key << "' given (should be posititve double, is '" << penaltyFactor_ << "') in the following Dune::ParameterTree:" << std::endl;
      paramTree.report(msg);
      DUNE_THROW(Dune::InvalidStateException, msg.str());
    }
    // test orders
    assert(model_->diffusionOrder() >= 0 && "Please provide a nonnegative order for the diffusion!");
    assert(model_->forceOrder() >= 0 && "Please provide a nonnegative order for the force!");
  }

  void init(const std::string prefix = "", std::ostream& out = Dune::Stuff::Common::Logger().debug())
  {
    if (!initialized_) {
      // prepare
      Dune::Timer timer;
      out << prefix << "initializing matrix container (of size "
          << localSolver_->ansatzSpace().map().size() << "x" << localSolver_->testSpace().map().size()
          << ")... " << std::flush;
      timer.reset();
      // compute pattern
      pattern_ = localSolver_->ansatzSpace().computeLocalPattern(*boundaryGridPart_,
                                                                 localSolver_->testSpace());
      // create matrix
      matrix_ = Dune::shared_ptr< MatrixBackendType >(new MatrixBackendType(
          ContainerFactory::createSparseMatrix(localSolver_->ansatzSpace().map().size(),
                                               localSolver_->testSpace().map().size(),
                                               *pattern_)));
      // create vector
      rhs_ = Dune::shared_ptr< VectorBackendType >(new VectorBackendType(
          ContainerFactory::createDenseVector(localSolver_->testSpace().map().size())));
      out<< "done (took " << timer.elapsed() << " sek)" << std::endl;
      out << prefix << "assembling matrix... " << std::flush;
      // evaluation
      typedef Dune::FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, dimRange > FunctionSpaceType;
      typedef typename ModelType::DiffusionType DiffusionType;
      typedef Dune::Detailed::Discretizations::Evaluation::Local::Binary::IPDGfluxes::Dirichlet< FunctionSpaceType, DiffusionType > IPDGfluxType;
      const IPDGfluxType ipdgFlux(model_->diffusion(), model_->diffusionOrder(), penaltyFactor_);
      // operator
      typedef Dune::Detailed::Discretizations::DiscreteOperator::Local::Codim1::BoundaryIntegral< IPDGfluxType > DirichletOperatorType;
      const DirichletOperatorType dirichletOperator(ipdgFlux);
      // local matrix assembler
      typedef Dune::Detailed::Discretizations::Assembler::Local::Codim1::Boundary< DirichletOperatorType > DirichletMatrixAssemblerType;
      const DirichletMatrixAssemblerType dirichletMatrixAssembler(dirichletOperator);
      // evaluation

      // functional

      // local vector assemblerh

      // assemble boundary matrix
      typedef Dune::Detailed::Discretizations::Assembler::Multiscale::Boundary  <
          BoundaryGridPartType,
          BoundaryInfoType,
          typename LocalSolverType::AnsatzSpaceType,
          typename LocalSolverType::TestSpaceType >
        BoundaryAssemblerType;
      const BoundaryAssemblerType boundaryAssembler(*boundaryGridPart_,
                                                    boundaryInfo_,
                                                    localSolver_->ansatzSpace(),
                                                    localSolver_->testSpace());
      boundaryAssembler.assembleMatrices(dirichletMatrixAssembler,
                                         *matrix_);
      out<< "done (took " << timer.elapsed() << " sek)" << std::endl;
      initialized_ = true;
    } // if (!initialized_)
  } // void init()

  const Dune::shared_ptr< const PatternType > pattern() const
  {
    return pattern_;
  }

  const Dune::shared_ptr< const MatrixBackendType > systemMatrix() const
  {
    return matrix_;
  }

  Dune::shared_ptr< MatrixBackendType > systemMatrix()
  {
    return matrix_;
  }

  const Dune::shared_ptr< const VectorBackendType > rightHandSide() const
  {
    return rhs_;
  }

  Dune::shared_ptr< VectorBackendType > rightHandSide()
  {
    return rhs_;
  }

private:
  const Dune::shared_ptr< const BoundaryGridPartType > boundaryGridPart_;
  const Dune::shared_ptr< const BoundaryInfoType > boundaryInfo_;
  const Dune::shared_ptr< const LocalSolverType > localSolver_;
  const Dune::shared_ptr< const ModelType > model_;
  bool initialized_;
  RangeFieldType penaltyFactor_;
  Dune::shared_ptr< const PatternType > pattern_;
  Dune::shared_ptr< MatrixBackendType > matrix_;
  Dune::shared_ptr< VectorBackendType > rhs_;
}; // DuneDetailedDiscretizations

template< class BoundaryGridPartType, class BoundaryInfoType, class LocalSolverType, class ContainerFactoryType >
const std::string DuneDetailedDiscretizations< BoundaryGridPartType, BoundaryInfoType, LocalSolverType, ContainerFactoryType >::id = "detailed.solvers.stationary.linear.elliptic.boundary";

} // namespace Boundary

} // namespace Elliptic

} // namespace Linear

} // namespace Solvers

} // namespace Stationary

} // namespace Detailed

} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_BOUNDARY_DUNE_DETAILED_DISCRETIZATIONS_HH
