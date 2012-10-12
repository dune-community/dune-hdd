
#ifndef DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_COUPLING_PRIMAL_DUNE_DETAILED_DISCRETIZATIONS_HH
#define DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_COUPLING_PRIMAL_DUNE_DETAILED_DISCRETIZATIONS_HH

// dune-common
#include <dune/common/shared_ptr.hh>
#include <dune/common/timer.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/exceptions.hh>

// dune-fem
#include <dune/fem/space/common/functionspace.hh>

// dune-detailed-discretizations
#include <dune/detailed/discretizations/discreteoperator/local/codim1/innerintegral.hh>
#include <dune/detailed/discretizations/evaluation/local/quaternary/ipdgfluxes.hh>
#include <dune/detailed/discretizations/assembler/local/codim1/matrix.hh>
#include <dune/detailed/discretizations/assembler/multiscale/coupling.hh>

// dune-stuff
#include <dune/stuff/common/parameter/tree.hh>

namespace Dune {

namespace Detailed {

namespace Solvers {

namespace Stationary {

namespace Linear {

namespace Elliptic {

namespace Coupling {

namespace Primal {

template< class CouplingGridPartImp, class InnerSolverImp, class OuterSolverImp, class ContainerFactoryImp >
class DuneDetailedDiscretizations
{
public:
  typedef CouplingGridPartImp CouplingGridPartType;

  typedef InnerSolverImp InnerSolverType;

  typedef OuterSolverImp OuterSolverType;

  typedef ContainerFactoryImp ContainerFactory;

  typedef DuneDetailedDiscretizations< CouplingGridPartType, InnerSolverType, OuterSolverType, ContainerFactory > ThisType;

  static const std::string id;

  typedef typename ContainerFactory::SparseMatrixType MatrixBackendType;

private:
  typedef typename InnerSolverType::ModelType ModelType;

  typedef typename ModelType::DomainFieldType DomainFieldType;

  static const int dimDomain = ModelType::dimDomain;

  typedef typename ModelType::RangeFieldType RangeFieldType;

  static const int dimRange = ModelType::dimRange;

public:
  typedef typename InnerSolverType::AnsatzSpaceType::PatternType PatternType;

  DuneDetailedDiscretizations(const Dune::shared_ptr< const CouplingGridPartType > innerOuterCouplingGridPart,
                              const Dune::shared_ptr< const CouplingGridPartType > outerInnerCouplingGridPart,
                              const Dune::shared_ptr< const InnerSolverType > innerSolver,
                              const Dune::shared_ptr< const OuterSolverType > outerSolver,
                              const Dune::Stuff::Common::ExtendedParameterTree paramTree)
    : innerOuterCouplingGridPart_(innerOuterCouplingGridPart)
    , outerInnerCouplingGridPart_(outerInnerCouplingGridPart)
    , innerSolver_(innerSolver)
    , outerSolver_(outerSolver)
    , model_(innerSolver_->model())
    , initialized_(false)
  {
    // check compatibility between the models
    assert(innerSolver_->model()->id == outerSolver_->model()->id && "The solvers models have to be compatible");
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
  }

  void init(const std::string prefix = "", std::ostream& out = Dune::Stuff::Common::Logger().debug())
  {
    if (!initialized_) {
      // prepare
      Dune::Timer timer;
      out << prefix << "initializing 4 matrix containers (of sizes ~ "
          << std::max(innerSolver_->ansatzSpace().map().size(), outerSolver_->ansatzSpace().map().size())
          << "x"
          << std::max(innerSolver_->testSpace().map().size(), outerSolver_->testSpace().map().size())
          << ")... " << std::flush;
      timer.reset();
      // computer patterns
      innerInnerPattern_ = innerSolver_->ansatzSpace().computeLocalPattern(*innerOuterCouplingGridPart_,
                                                                           innerSolver_->testSpace());
      innerOuterPattern_ = innerSolver_->ansatzSpace().computeCouplingPattern(*innerOuterCouplingGridPart_,
                                                                              outerSolver_->testSpace());
      outerInnerPattern_ = outerSolver_->ansatzSpace().computeCouplingPattern(*outerInnerCouplingGridPart_,
                                                                              innerSolver_->testSpace());
      outerOuterPattern_ = outerSolver_->ansatzSpace().computeLocalPattern(*outerInnerCouplingGridPart_,
                                                                           outerSolver_->testSpace());
      // create matrices
      innerInnerMatrix_ = Dune::shared_ptr< MatrixBackendType >(new MatrixBackendType(
          ContainerFactory::createSparseMatrix(innerSolver_->ansatzSpace().map().size(),
                                               innerSolver_->testSpace().map().size(),
                                               *innerInnerPattern_)));
      innerOuterMatrix_ = Dune::shared_ptr< MatrixBackendType >(new MatrixBackendType(
          ContainerFactory::createSparseMatrix(innerSolver_->ansatzSpace().map().size(),
                                               outerSolver_->testSpace().map().size(),
                                               *innerOuterPattern_)));
      outerInnerMatrix_ = Dune::shared_ptr< MatrixBackendType >(new MatrixBackendType(
          ContainerFactory::createSparseMatrix(outerSolver_->ansatzSpace().map().size(),
                                               innerSolver_->testSpace().map().size(),
                                               *outerInnerPattern_)));
      outerOuterMatrix_ = Dune::shared_ptr< MatrixBackendType >(new MatrixBackendType(
          ContainerFactory::createSparseMatrix(outerSolver_->ansatzSpace().map().size(),
                                               outerSolver_->testSpace().map().size(),
                                               *outerOuterPattern_)));
      initialized_ = true;
      out<< "done (took " << timer.elapsed() << " sek)" << std::endl;
      out << prefix << "assembling 4 matrices... " << std::flush;
      // evaluation
      typedef Dune::FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, dimRange > FunctionSpaceType;
      typedef typename ModelType::DiffusionType DiffusionType;
      typedef Dune::Detailed::Discretizations::Evaluation::Local::Quaternary::IPDGfluxes::Inner< FunctionSpaceType, DiffusionType > IPDGfluxType;
      const IPDGfluxType ipdgFlux(model_->diffusion(), model_->diffusionOrder(), penaltyFactor_);
      // operator
      typedef Dune::Detailed::Discretizations::DiscreteOperator::Local::Codim1::InnerIntegral< IPDGfluxType > IPDGoperatorType;
      const IPDGoperatorType ipdgOperator(ipdgFlux);
      // local assembler
      typedef Dune::Detailed::Discretizations::Assembler::Local::Codim1::Inner< IPDGoperatorType > CouplingMatrixAssemblerType;
      const CouplingMatrixAssemblerType couplingMatrixAssembler(ipdgOperator);
      // assemble coupling matrix
      typedef Dune::Detailed::Discretizations::Assembler::Multiscale::Coupling::Primal<
          CouplingGridPartType,
          typename InnerSolverType::AnsatzSpaceType,
          typename InnerSolverType::TestSpaceType,
          typename OuterSolverType::AnsatzSpaceType,
          typename OuterSolverType::TestSpaceType >
        CouplingAssemblerType;
      const CouplingAssemblerType couplingAssembler(*innerOuterCouplingGridPart_,
                                                    innerSolver_->ansatzSpace(),
                                                    innerSolver_->testSpace(),
                                                    outerSolver_->ansatzSpace(),
                                                    outerSolver_->testSpace());
      couplingAssembler.assembleMatrices(couplingMatrixAssembler,
                                         *innerInnerMatrix_,
                                         *innerOuterMatrix_,
                                         *outerInnerMatrix_,
                                         *outerOuterMatrix_);
      out<< "done (took " << timer.elapsed() << " sek)" << std::endl;
    } // if (!initialized_)
  } // void init()

  const Dune::shared_ptr< const PatternType > innerInnerPattern() const
  {
    return innerInnerPattern_;
  }

  const Dune::shared_ptr< const PatternType > innerOuterPattern() const
  {
    return innerOuterPattern_;
  }

  const Dune::shared_ptr< const PatternType > outerInnerPattern() const
  {
    return outerInnerPattern_;
  }

  const Dune::shared_ptr< const PatternType > outerOuterPattern() const
  {
    return outerOuterPattern_;
  }

  const Dune::shared_ptr< const MatrixBackendType > innerInnerMatrix() const
  {
    return innerInnerMatrix_;
  }

  Dune::shared_ptr< MatrixBackendType > innerInnerMatrix()
  {
    return innerInnerMatrix_;
  }

  const Dune::shared_ptr< const MatrixBackendType > innerOuterMatrix() const
  {
    return innerOuterMatrix_;
  }

  Dune::shared_ptr< MatrixBackendType > innerOuterMatrix()
  {
    return innerOuterMatrix_;
  }

  const Dune::shared_ptr< const MatrixBackendType > outerInnerMatrix() const
  {
    return outerInnerMatrix_;
  }

  Dune::shared_ptr< MatrixBackendType > outerInnerMatrix()
  {
    return outerInnerMatrix_;
  }

  const Dune::shared_ptr< const MatrixBackendType > outerOuterMatrix() const
  {
    return outerOuterMatrix_;
  }

  Dune::shared_ptr< MatrixBackendType > outerOuterMatrix()
  {
    return outerOuterMatrix_;
  }

private:
  const Dune::shared_ptr< const CouplingGridPartType > innerOuterCouplingGridPart_;
  const Dune::shared_ptr< const CouplingGridPartType > outerInnerCouplingGridPart_;
  const Dune::shared_ptr< const InnerSolverType > innerSolver_;
  const Dune::shared_ptr< const OuterSolverType > outerSolver_;
  const Dune::shared_ptr< const ModelType > model_;
  bool initialized_;
  RangeFieldType penaltyFactor_;
  Dune::shared_ptr< const PatternType > innerInnerPattern_;
  Dune::shared_ptr< const PatternType > innerOuterPattern_;
  Dune::shared_ptr< const PatternType > outerInnerPattern_;
  Dune::shared_ptr< const PatternType > outerOuterPattern_;
  Dune::shared_ptr< MatrixBackendType > innerInnerMatrix_;
  Dune::shared_ptr< MatrixBackendType > innerOuterMatrix_;
  Dune::shared_ptr< MatrixBackendType > outerInnerMatrix_;
  Dune::shared_ptr< MatrixBackendType > outerOuterMatrix_;
}; // DuneDetailedDiscretizations

template< class CouplingGridPartType, class InnerSolverType, class OuterSolverType, class ContainerFactoryType >
const std::string DuneDetailedDiscretizations< CouplingGridPartType, InnerSolverType, OuterSolverType, ContainerFactoryType >::id = "detailed.solvers.stationary.linear.elliptic.coupling.primal";

} // namespace Primal

} // namespace Coupling

} // namespace Elliptic

} // namespace Linear

} // namespace Solvers

} // namespace Stationary

} // namespace Detailed

} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_COUPLING_PRIMAL_DUNE_DETAILED_DISCRETIZATIONS_HH
