#ifndef DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MULTISCALE_SEMICONTINUOUSGALERKIN_DUNE_DETAILED_DISCRETIZATIONS_HH
#define DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MULTISCALE_SEMICONTINUOUSGALERKIN_DUNE_DETAILED_DISCRETIZATIONS_HH

// system
#include <vector>
#include <sstream>

// dune-common
#include <dune/common/exceptions.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/timer.hh>

// dune-fem
#include <dune/fem/space/common/functionspace.hh>

// dune-grid-multiscale
#include <dune/grid/multiscale/default.hh>

// dune-detailed-discretizations
#include <dune/detailed/discretizations/la/factory/eigen.hh>

// dune-detailed-solvers
#include <dune/detailed/solvers/stationary/linear/elliptic/continuousgalerkin/dune-detailed-discretizations.hh>

// dune-stuff
#include <dune/stuff/common/logging.hh>

namespace Dune {

namespace Detailed {

namespace Solvers {

namespace Stationary {

namespace Linear {

namespace Elliptic {

namespace Multiscale {

namespace SemicontinuousGalerkin {

template< class ModelImp, class MsGridImp, int polynomialOrder >
class DuneDetailedDiscretizations
{
public:
  typedef ModelImp ModelType;

  typedef MsGridImp MsGridType;

  typedef typename MsGridType::GlobalGridPartType GlobalGridPartType;

  typedef GlobalGridPartType GridPartType;

  static const int polOrder = polynomialOrder;

  typedef DuneDetailedDiscretizations< ModelType, MsGridType, polOrder > ThisType;

  static const std::string id;

private:
  typedef typename ModelType::DomainFieldType DomainFieldType;

  static const int dimDomain = ModelType::dimDomain;

  typedef typename ModelType::RangeFieldType RangeFieldType;

  static const int dimRange = ModelType::dimRange;

  typedef Dune::Detailed::Discretizations::LA::Factory::Eigen< RangeFieldType > ContainerFactory;

  typedef typename ContainerFactory::SparseMatrixType MatrixBackendType;

  typedef typename ContainerFactory::DenseVectorType VectorBackendType;

  typedef Dune::FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, dimRange > FunctionSpaceType;

  typedef typename MsGridType::LocalGridPartType LocalGridPartType;

  typedef typename Dune::Detailed::Solvers::Stationary::Linear::Elliptic::ContinuousGalerkin::DuneDetailedDiscretizations< ModelType, LocalGridPartType, polOrder > LocalSolverType;

public:
  typedef typename MatrixBackendType::StorageType MatrixType;

  typedef typename VectorBackendType::StorageType VectorType;

  DuneDetailedDiscretizations(const ModelType& model, const MsGridType& msGrid)
    : model_(model),
      msGrid_(msGrid)
  {}

  void init(Dune::ParameterTree paramTree = Dune::ParameterTree())
  {
    // logging
    const std::string prefix = paramTree.get("prefix", "");
    Dune::Stuff::Common::LogStream& debug = Dune::Stuff::Common::Logger().debug();

    // timer
    Dune::Timer timer;

    // create local solvers for each subdomain
    const unsigned int subdomains = msGrid_.size();
    debug << std::endl << prefix << id << ".init:" << std::endl;
    debug << prefix << "initializing " << subdomains << " local continuous galerkin solvers" << std::flush;
    for (unsigned int subdomain = 0; subdomain < subdomains; ++subdomain) {
      debug << "." << std::flush;
      localGridParts_.push_back(msGrid_.localGridPart(subdomain));
      localSolvers_.push_back(Dune::shared_ptr< LocalSolverType >(new LocalSolverType(model_, *(localGridParts_[subdomain]))));
      localSolvers_[subdomain]->init();
    }
    debug << " done (took " << timer.elapsed() << " sek)" << std::endl;

//    // function spaces
//    if (verbose) {
//      std::cout << prefix << "setting up function spaces... " << std::flush;
//      timer.reset();
//    }
//    discreteH1_ = Dune::shared_ptr< DiscreteH1Type >(new DiscreteH1Type(gridPart_));
//    ansatzSpace_ = Dune::shared_ptr< AnsatzSpaceType >(new AnsatzSpaceType(*discreteH1_));
//    testSpace_ = Dune::shared_ptr< TestSpaceType >(new TestSpaceType(*discreteH1_));
//    if (verbose)
//      std::cout << "done (took " << timer.elapsed() << " sec)" << std::endl;

//    // left hand side (operator)
//    if (verbose) {
//      std::cout << prefix << "setting up operator and functional... " << std::flush;
//      timer.reset();
//    }
//    typedef typename ModelType::DiffusionType DiffusionType;
//    typedef Dune::Detailed::Discretizations::Evaluation::Local::Binary::Elliptic< FunctionSpaceType, DiffusionType > EllipticEvaluationType;
//    const EllipticEvaluationType ellipticEvaluation(model_.diffusion());
//    typedef Dune::Detailed::Discretizations::DiscreteOperator::Local::Codim0::Integral< EllipticEvaluationType > EllipticOperatorType;
//    const EllipticOperatorType ellipticOperator(ellipticEvaluation);
//    typedef Dune::Detailed::Discretizations::Evaluation::Local::Quaternary::IPDGfluxes::Inner< FunctionSpaceType, DiffusionType > InnerIPDGEvaluationType;
//    const InnerIPDGEvaluationType innerIPDGEvaluation(model_.diffusion());
//    typedef Dune::Detailed::Discretizations::DiscreteOperator::Local::Codim1::InnerIntegral< InnerIPDGEvaluationType > InnerIPDGOperatorType;
//    const InnerIPDGOperatorType innerIPDGOperator(innerIPDGEvaluation);

//    // right hand side (functional)
//    typedef typename ModelType::ForceType ForceType;
//    typedef Dune::Detailed::Discretizations::Evaluation::Local::Unary::Scale< FunctionSpaceType, ForceType > ProductEvaluationType;
//    const ProductEvaluationType productEvaluation(model_.force());
//    typedef Dune::Detailed::Discretizations::DiscreteFunctional::Local::Codim0::Integral< ProductEvaluationType > L2FunctionalType;
//    const L2FunctionalType l2Functional(productEvaluation);
//    if (verbose)
//      std::cout << "done (took " << timer.elapsed() << " sec)" << std::endl;

//    // system matrix and right hand side
//    if (verbose) {
//      std::cout << prefix << "setting up matrix and vector container... " << std::flush;
//      timer.reset();
//    }
//    matrix_ = Dune::shared_ptr< MatrixBackendType >(new MatrixBackendType(ContainerFactory::createSparseMatrix(*ansatzSpace_, *testSpace_)));
//    rhs_ = Dune::shared_ptr< VectorBackendType >(new VectorBackendType(ContainerFactory::createDenseVector(*testSpace_)));
//    if (verbose)
//      std::cout << "done (took " << timer.elapsed() << " sec)" << std::endl;

//    // assemble system
//    if (verbose) {
//      std::cout << prefix << "assembling system... " << std::flush;
//      timer.reset();
//    }
//    typedef Dune::Detailed::Discretizations::Assembler::Local::Codim0::Matrix< EllipticOperatorType > LocalCodim0MatrixAssemblerType;
//    const LocalCodim0MatrixAssemblerType localCodim0MatrixAssembler(ellipticOperator);
//    typedef Dune::Detailed::Discretizations::Assembler::Local::Codim1::Matrix< InnerIPDGOperatorType > LocalCodim1MatrixAssemblerType;
//    const LocalCodim1MatrixAssemblerType localCodim1MatrixAssembler(innerIPDGOperator);
//    typedef Dune::Detailed::Discretizations::Assembler::Local::Combined::Matrix< LocalCodim0MatrixAssemblerType, LocalCodim1MatrixAssemblerType > LocalMatrixAssemblerType;
//    const LocalMatrixAssemblerType localMatrixAssembler(localCodim0MatrixAssembler, localCodim1MatrixAssembler);
//    typedef Dune::Detailed::Discretizations::Assembler::Local::Codim0::Vector< L2FunctionalType > LocalVectorAssemblerType;
//    const LocalVectorAssemblerType localVectorAssembler(l2Functional);
//    typedef Dune::Detailed::Discretizations::Assembler::System::Constrained< AnsatzSpaceType, TestSpaceType > SystemAssemblerType;
//    const SystemAssemblerType systemAssembler(*ansatzSpace_, *testSpace_);
//    systemAssembler.assembleSystem(localMatrixAssembler, *matrix_, localVectorAssembler, *rhs_);
//    if (verbose)
//      std::cout << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void init(Dune::ParameterTree paramTree = Dune::ParameterTree())

//  void solve(Dune::shared_ptr< VectorType >& solution, Dune::ParameterTree paramTree = Dune::ParameterTree()) const
//  {
//    // preparations
//    const bool verbose = paramTree.get("verbose", false);
//    const std::string prefix = paramTree.get("prefix", "");
//    const std::string type = paramTree.get("type", "eigen.cg.diagonal.upper");
//    const unsigned int maxIter = paramTree.get("maxIter", 5000);
//    const double precision = paramTree.get("precision", 1e-12);
//    Dune::Timer timer;
//    VectorBackendType vector(solution);
//    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::BicgstabIlut BicgstabIlutSolver;
//    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::BicgstabDiagonal BicgstabDiagonalSolver;
//    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::CgDiagonalUpper CgDiagonalUpperSolver;
//    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::CgDiagonalLower CgDiagonalLowerSolver;
//    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::SimplicialcholeskyUpper SimplicialcholeskyUpperSolver;
//    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::SimplicialcholeskyLower SimplicialcholeskyLowerSolver;
//    if (verbose)
//      std::cout << prefix << "solving linear system of size " << matrix_->rows() << "x" << matrix_->cols() << std::endl
//                << prefix << "using " << type << "... " << std::flush;
//    if (type == "eigen.bicgstab.incompletelut"){
//      BicgstabIlutSolver::apply(*matrix_, vector, *rhs_, maxIter, precision);
//    } else if (type == "eigen.bicgstab.diagonal"){
//      BicgstabDiagonalSolver::apply(*matrix_, vector, *rhs_, maxIter, precision);
//    } else if (type == "eigen.cg.diagonal.upper"){
//      CgDiagonalUpperSolver::apply(*matrix_, vector, *rhs_, maxIter, precision);
//    } else if (type == "eigen.cg.diagonal.lower"){
//      CgDiagonalLowerSolver::apply(*matrix_, vector, *rhs_, maxIter, precision);
//    } else if (type == "eigen.simplicialcholesky.upper"){
//      SimplicialcholeskyUpperSolver::apply(*matrix_, vector, *rhs_, maxIter, precision);
//    } else if (type == "eigen.simplicialcholesky.lower"){
//      SimplicialcholeskyLowerSolver::apply(*matrix_, vector, *rhs_, maxIter, precision);
//    } else {
//      std::stringstream msg;
//      msg << "Error";
//      if (id != "") {
//        msg << " in " << id;
//      }
//      msg << ": solver type '" << type << "not supported!" << std::endl;
//      DUNE_THROW(Dune::InvalidStateException, msg.str());
//    }
//    if (verbose)
//      std::cout << "done (took " << timer.elapsed() << " sec)" << std::endl;
//  } // void solve(Dune::shared_ptr< VectorType >& solution, Dune::ParameterTree paramTree = Dune::ParameterTree()) const

//  void visualize(const Dune::shared_ptr< VectorType >& vector, Dune::ParameterTree paramTree = Dune::ParameterTree()) const
//  {
//    // preparations
//    const bool verbose = paramTree.get("verbose", false);
//    const std::string prefix = paramTree.get("prefix", "");
//    const std::string name = paramTree.get("name", id);
//    const std::string filename = paramTree.get("filename", "visualization");
//    Dune::Timer timer;
//    if (verbose) {
//        std::cout << prefix << "writing '" << name << "' to '" << filename;
//      if (dimDomain == 1)
//        std::cout << ".vtp";
//      else
//        std::cout << ".vtu";
//      std::cout << "'... " << std::flush;
//    }
//    const VectorBackendType vectorBackend(vector);
//    typedef Dune::Detailed::Discretizations::DiscreteFunction::Default< AnsatzSpaceType, VectorBackendType > DiscreteFunctionType;
//    Dune::shared_ptr< DiscreteFunctionType > discreteFunction(new DiscreteFunctionType(*ansatzSpace_, vectorBackend, name));
//    typedef Dune::VTKWriter< typename AnsatzSpaceType::GridViewType > VTKWriterType;
//    VTKWriterType vtkWriter(ansatzSpace_->gridView());
//    vtkWriter.addVertexData(discreteFunction);
//    vtkWriter.write(filename);
//    if (verbose)
//      std::cout << "done (took " << timer.elapsed() << " sec)" << std::endl;
//  }

//  Dune::shared_ptr< VectorType > createVector() const
//  {
//    VectorBackendType tmp = ContainerFactory::createDenseVector(*testSpace_);
//    return tmp.storage();
//  }

private:
  const ModelType& model_;
  const MsGridType& msGrid_;
  std::vector< Dune::shared_ptr< const LocalGridPartType > > localGridParts_;
  std::vector< Dune::shared_ptr< LocalSolverType > > localSolvers_;
}; // class DuneDetailedDiscretizations

template< class ModelType, class GridPartType, int polOrder >
const std::string DuneDetailedDiscretizations< ModelType, GridPartType, polOrder >::id = "detailed.solvers.stationary.linear.elliptic.multiscale.semicontinuousgalerkin";

} // namespace SemicontinuousGalerkin

} // namespace Multiscale

} // namespace Elliptic

} // namespace Linear

} // namespace Stationary

} // namespace Solvers

} // namespace Detailed

} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MULTISCALE_SEMICONTINUOUSGALERKIN_DUNE_DETAILED_DISCRETIZATIONS_HH
