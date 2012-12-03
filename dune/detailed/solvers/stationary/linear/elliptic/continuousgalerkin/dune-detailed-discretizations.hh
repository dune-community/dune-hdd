#ifndef DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_CONTINUOUSGALERKIN_DUNE_DETAILED_DISCRETIZATIONS_HH
#define DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_CONTINUOUSGALERKIN_DUNE_DETAILED_DISCRETIZATIONS_HH

#include <sstream>

#include <dune/common/exceptions.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/timer.hh>

#include <dune/fem/space/common/functionspace.hh>

#include <dune/detailed/discretizations/discretefunctionspace/continuous/lagrange.hh>
#include <dune/detailed/discretizations/discretefunctionspace/sub/linear.hh>
#include <dune/detailed/discretizations/discretefunction/default.hh>
#include <dune/detailed/discretizations/discretefunctionspace/sub/affine.hh>
#include <dune/detailed/discretizations/evaluation/local/binary/elliptic.hh>
#include <dune/detailed/discretizations/discreteoperator/local/codim0/integral.hh>
#include <dune/detailed/discretizations/evaluation/local/unary/scale.hh>
#include <dune/detailed/discretizations/discretefunctional/local/codim0/integral.hh>
#include <dune/detailed/discretizations/la/factory/eigen.hh>
#include <dune/detailed/discretizations/assembler/local/codim0/matrix.hh>
#include <dune/detailed/discretizations/assembler/local/codim0/vector.hh>
#include <dune/detailed/discretizations/assembler/system/constrained.hh>
#include <dune/detailed/discretizations/la/backend/solver/eigen.hh>
#include <dune/detailed/discretizations/discretefunction/default.hh>

#include <dune/stuff/common/logging.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/discretefunction/projection/dirichlet.hh>

namespace Dune {
namespace Detailed {
namespace Solvers {
namespace Stationary {
namespace Linear {
namespace Elliptic {
namespace ContinuousGalerkin {

/**
 *  \todo Add method to create a discrete function.
 *  \todo Change id -> static id()
 */
template< class ModelImp, class GridPartImp, class BoundaryInfoImp, int polynomialOrder >
class DuneDetailedDiscretizations
{
public:
  typedef ModelImp ModelType;

  typedef GridPartImp GridPartType;

  typedef BoundaryInfoImp BoundaryInfoType;

  static const int polOrder = polynomialOrder;

  typedef DuneDetailedDiscretizations< ModelType, GridPartType, BoundaryInfoType, polOrder > ThisType;

private:
  typedef typename ModelType::DomainFieldType DomainFieldType;

  static const int dimDomain = ModelType::dimDomain;

  typedef typename ModelType::RangeFieldType RangeFieldType;

  static const int dimRange = ModelType::dimRange;

  typedef Dune::Detailed::Discretizations::LA::Factory::Eigen< RangeFieldType > ContainerFactory;

  typedef typename ContainerFactory::SparseMatrixType MatrixBackendType;

public:
  typedef typename ContainerFactory::DenseVectorType VectorBackendType;

private:
  typedef Dune::FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, dimRange > FunctionSpaceType;

  typedef Dune::Detailed::Discretizations::DiscreteFunctionSpace::Continuous::Lagrange< FunctionSpaceType, GridPartType, polOrder > DiscreteH1Type;

public:
  typedef Dune::Detailed::Discretizations::DiscreteFunctionSpace::Sub::Linear::Dirichlet< DiscreteH1Type, BoundaryInfoType > TestSpaceType;

  typedef Dune::Detailed::Discretizations::DiscreteFunction::Default< TestSpaceType, VectorBackendType > DiscreteFunctionType;

  typedef TestSpaceType AnsatzSpaceType;

  typedef typename AnsatzSpaceType::PatternType PatternType;

  static const std::string id()
  {
    return "detailed.solvers.stationary.linear.elliptic.continuousgalerkin";
  }

  DuneDetailedDiscretizations(const Dune::shared_ptr< const ModelType > model,
                              const Dune::shared_ptr< const GridPartType > gridPart,
                              const Dune::shared_ptr< const BoundaryInfoType > boundaryInfo)
    : model_(model)
    , gridPart_(gridPart)
    , boundaryInfo_(boundaryInfo)
    , initialized_(false)
  {
  }

  const Dune::shared_ptr< const ModelType > model() const
  {
    return model_;
  }

  const Dune::shared_ptr< const GridPartType > gridPart() const
  {
    return gridPart_;
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

      // function spaces
      out << prefix << "initializing function spaces... " << std::flush;
      timer.reset();
      discreteH1_ = Dune::shared_ptr< const DiscreteH1Type >(new DiscreteH1Type(*gridPart_));
      testSpace_ = Dune::shared_ptr< TestSpaceType >(new TestSpaceType(*discreteH1_, boundaryInfo_));
//      Dune::shared_ptr< DiscreteFunctionType > dirichlet(new DiscreteFunctionType(*testSpace_, "dirichlet"));
//      Dune::Stuff::DiscreteFunction::Projection::Dirichlet::project(testSpace_->boundaryInfo(),
//                                                                    *(model_->dirichlet()),
//                                                                    *dirichlet);
      ansatzSpace_ = testSpace_;//Dune::shared_ptr< AnsatzSpaceType >(new AnsatzSpaceType(*testSpace_, dirichlet));
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;
      // left hand side (operator)
      out << prefix << "initializing operator and functional... " << std::flush;
      timer.reset();
      typedef typename ModelType::DiffusionType DiffusionType;
      typedef Dune::Detailed::Discretizations::Evaluation::Local::Binary::Elliptic< FunctionSpaceType, DiffusionType > EllipticEvaluationType;
      const EllipticEvaluationType ellipticEvaluation(model_->diffusion(), model_->diffusionOrder());
      typedef Dune::Detailed::Discretizations::DiscreteOperator::Local::Codim0::Integral< EllipticEvaluationType > EllipticOperatorType;
      const EllipticOperatorType ellipticOperator(ellipticEvaluation);
      // right hand side (functional)
      typedef typename ModelType::ForceType ForceType;
      typedef Dune::Detailed::Discretizations::Evaluation::Local::Unary::Scale< FunctionSpaceType, ForceType > ProductEvaluationType;
      const ProductEvaluationType productEvaluation(model_->force(), model_->forceOrder());
      typedef Dune::Detailed::Discretizations::DiscreteFunctional::Local::Codim0::Integral< ProductEvaluationType > L2FunctionalType;
      const L2FunctionalType l2Functional(productEvaluation);
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;
      // system matrix and right hand side
      out << prefix << "initializing matrix and vector container (of size "
          << ansatzSpace_->map().size() << "x" << testSpace_->map().size() << ")... " << std::flush;
      timer.reset();
      pattern_ = ansatzSpace_->computePattern(*testSpace_);
      matrix_ = Dune::shared_ptr< MatrixBackendType >(new MatrixBackendType(ContainerFactory::createSparseMatrix(ansatzSpace_->map().size(), testSpace_->map().size(), *pattern_)));
      rhs_ = Dune::shared_ptr< VectorBackendType >(new VectorBackendType(ContainerFactory::createDenseVector(*testSpace_)));
      affineShift_ = Dune::shared_ptr< VectorBackendType >(new VectorBackendType(ContainerFactory::createDenseVector(*ansatzSpace_)));
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;
      // assemble system
      out << prefix << "assembling system... " << std::flush;
      timer.reset();
      typedef Dune::Detailed::Discretizations::Assembler::Local::Codim0::Matrix< EllipticOperatorType > LocalMatrixAssemblerType;
      const LocalMatrixAssemblerType localmatrixAssembler(ellipticOperator);
      typedef Dune::Detailed::Discretizations::Assembler::Local::Codim0::Vector< L2FunctionalType > LocalVectorAssemblerType;
      const LocalVectorAssemblerType localVectorAssembler(l2Functional);
      typedef Dune::Detailed::Discretizations::Assembler::System::Constrained< AnsatzSpaceType, TestSpaceType > SystemAssemblerType;
      const SystemAssemblerType systemAssembler(*ansatzSpace_, *testSpace_);
      systemAssembler.assemble(localmatrixAssembler, *matrix_, localVectorAssembler, *rhs_);
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;
      // done
      initialized_ = true;
    } // if !(initialized_)
  } // void init()

  Dune::shared_ptr< VectorBackendType > createVector() const
  {
    assert(initialized_ && "A vector can only be created after init() has been called!");
    Dune::shared_ptr< VectorBackendType > ret(new VectorBackendType(ContainerFactory::createDenseVector(*ansatzSpace_)));
    return ret;
  }

  Dune::shared_ptr< DiscreteFunctionType > createDiscreteFunction(const std::string name = "discrete_function") const
  {
    assert(initialized_ && "Please call init() beafore calling createDiscreteFunction()!");
    return createDiscreteFunction(ContainerFactory::createDenseVector(*ansatzSpace_), name);
  } // ... createDiscreteFunction(...)

  Dune::shared_ptr< DiscreteFunctionType > createDiscreteFunction(VectorBackendType vector,
                                                                  const std::string name = "discrete_function") const
  {
    assert(initialized_ && "Please call init() beafore calling createDiscreteFunction()!");
    assert(vector.size() == ansatzSpace_->map().size());
    Dune::shared_ptr< DiscreteFunctionType > discreteFunction(new DiscreteFunctionType(*ansatzSpace_,
                                                                                       vector,
                                                                                       name));
    return discreteFunction;
  } // ... createDiscreteFunction(...)

  void solve(VectorBackendType& solution,
             const std::string type = "eigen.bicgstab.incompletelut",
             const unsigned int maxIter = 5000,
             const double precision = 1e-12,
             const std::string prefix = "",
             std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
  {
    // preparations
    assert(initialized_ && "The system can only be solved after init() has been called! ");
    assert(solution.size() == ansatzSpace_->map().size() && "Given vector has wrong size!");
    Dune::Timer timer;
    // solve
    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::BicgstabIlut BicgstabIlutSolver;
    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::BicgstabDiagonal BicgstabDiagonalSolver;
    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::CgDiagonalUpper CgDiagonalUpperSolver;
    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::CgDiagonalLower CgDiagonalLowerSolver;
    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::SimplicialcholeskyUpper SimplicialcholeskyUpperSolver;
    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::SimplicialcholeskyLower SimplicialcholeskyLowerSolver;
    out << prefix << "solving linear system of size " << matrix_->rows() << "x" << matrix_->cols() << std::endl
        << prefix << "  using " << type << "... " << std::flush;
    if (type == "eigen.bicgstab.incompletelut"){
      BicgstabIlutSolver::apply(*matrix_, solution, *rhs_, maxIter, precision);
    } else if (type == "eigen.bicgstab.diagonal"){
      BicgstabDiagonalSolver::apply(*matrix_, solution, *rhs_, maxIter, precision);
    } else if (type == "eigen.cg.diagonal.upper"){
      CgDiagonalUpperSolver::apply(*matrix_, solution, *rhs_, maxIter, precision);
    } else if (type == "eigen.cg.diagonal.lower"){
      CgDiagonalLowerSolver::apply(*matrix_, solution, *rhs_, maxIter, precision);
    } else if (type == "eigen.simplicialcholesky.upper"){
      SimplicialcholeskyUpperSolver::apply(*matrix_, solution, *rhs_, maxIter, precision);
    } else if (type == "eigen.simplicialcholesky.lower"){
      SimplicialcholeskyLowerSolver::apply(*matrix_, solution, *rhs_, maxIter, precision);
    } else {
      DUNE_THROW(Dune::RangeError,
                 "\nError: wrong 'type' given (" << type << ")!");
    }
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void solve()

  void visualize(VectorBackendType& vector,
                 const std::string filename = "solution",
                 const std::string name = "solution",
                 const std::string prefix = "",
                 std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
  {
    // preparations
    assert(initialized_ && "A vector can only be visualized after init() has been called! ");
    assert(vector.size() == ansatzSpace_->map().size() && "Given vector has wrong size!");
    Dune::Timer timer;
    out << prefix << "writing '" << name << "' to '" << filename;
    if (dimDomain == 1)
      out << ".vtp";
    else
      out << ".vtu";
    out << "'... " << std::flush;
    Dune::shared_ptr< DiscreteFunctionType > discreteFunction(new DiscreteFunctionType(*ansatzSpace_, vector, name));
    visualize(discreteFunction, filename, "", Dune::Stuff::Common::Logger().devnull());
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
  }

  void visualize(Dune::shared_ptr< DiscreteFunctionType > discreteFunction,
                 const std::string filename = "solution",
                 const std::string prefix = "",
                 std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
  {
    // preparations
    assert(initialized_ && "A discrete function can only be visualized after init() has been called! ");
    assert(discreteFunction->size() == ansatzSpace_->map().size() && "Given vector has wrong size!");
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
  }

  const AnsatzSpaceType& ansatzSpace() const
  {
    assert(initialized_);
    return *ansatzSpace_;
  }

  const TestSpaceType& testSpace() const
  {
    assert(initialized_);
    return *testSpace_;
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

private:
  const Dune::shared_ptr< const ModelType > model_;
  const Dune::shared_ptr< const GridPartType > gridPart_;
  const Dune::shared_ptr< const BoundaryInfoType > boundaryInfo_;
  bool initialized_;
  Dune::shared_ptr< const DiscreteH1Type > discreteH1_;
  Dune::shared_ptr< const TestSpaceType > testSpace_;
  Dune::shared_ptr< const AnsatzSpaceType > ansatzSpace_;
  Dune::shared_ptr< const PatternType > pattern_;
  Dune::shared_ptr< MatrixBackendType > matrix_;
  Dune::shared_ptr< VectorBackendType > rhs_;
  Dune::shared_ptr< VectorBackendType > affineShift_;
}; // class DuneDetailedDiscretizations

} // namespace ContinuousGalerkin
} // namespace Elliptic
} // namespace Linear
} // namespace Stationary
} // namespace Solvers
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_CONTINUOUSGALERKIN_DUNE_DETAILED_DISCRETIZATIONS_HH
