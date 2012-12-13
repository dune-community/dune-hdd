#ifndef DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_CG_DETAILED_DISCRETIZATIONS_HH
#define DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_CG_DETAILED_DISCRETIZATIONS_HH

#include <sstream>

#include <dune/common/exceptions.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/timer.hh>

#include <dune/stuff/common/logging.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/discretefunction/projection/dirichlet.hh>

#include <dune/grid/part/interface.hh>

#include <dune/detailed/discretizations/discretefunctionspace/continuous/lagrange.hh>
#include <dune/detailed/discretizations/discretefunctionspace/sub/linear.hh>
#include <dune/detailed/discretizations/la/container/factory/eigen.hh>
#include <dune/detailed/discretizations/discretefunction/default.hh>
#include <dune/detailed/discretizations/discretefunctionspace/sub/affine.hh>
#include <dune/detailed/discretizations/evaluation/local/binary/elliptic.hh>
#include <dune/detailed/discretizations/discreteoperator/local/codim0/integral.hh>
#include <dune/detailed/discretizations/evaluation/local/unary/scale.hh>
#include <dune/detailed/discretizations/discretefunctional/local/codim0/integral.hh>
#include <dune/detailed/discretizations/discretefunctional/local/codim1/integral.hh>
#include <dune/detailed/discretizations/assembler/local/codim0/matrix.hh>
#include <dune/detailed/discretizations/assembler/local/codim0/vector.hh>
#include <dune/detailed/discretizations/assembler/local/codim1/vector.hh>
#include <dune/detailed/discretizations/assembler/system.hh>

#include "../model/interface.hh"

namespace Dune {
namespace Detailed {
namespace Solvers {
namespace Stationary {
namespace Linear {
namespace Elliptic {
namespace CG {


template< class GridPartImp, int polynomialOrder, class RangeFieldImp, int dimensionRange >
class DetailedDiscretizations;

template< class GridPartImp, int polynomialOrder, class RangeFieldImp >
class DetailedDiscretizations< GridPartImp, polynomialOrder, RangeFieldImp, 1 >
{
public:
  typedef DetailedDiscretizations< GridPartImp, polynomialOrder, RangeFieldImp, 1 > ThisType;

  typedef Dune::grid::Part::Interface< typename GridPartImp::Traits > GridPartType;

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

  typedef Dune::FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, dimRange > FunctionSpaceType;

private:
  typedef Dune::Detailed::Discretizations
      ::DiscreteFunctionSpace
      ::Continuous
      ::Lagrange< FunctionSpaceType, GridPartType, polOrder >
    LagrangeSpaceType;

public:
  typedef Dune::Detailed::Discretizations
      ::DiscreteFunctionSpace
      ::Sub
      ::Linear
      ::Dirichlet< LagrangeSpaceType >
    TestSpaceType;

private:
  typedef typename Dune::Detailed::Discretizations::LA::Container::Factory::Eigen< RangeFieldType > ContainerFactory;

public:
  typedef typename ContainerFactory::DenseVectorType VectorType;

  typedef typename ContainerFactory::SparseMatrixType MatrixType;

  typedef typename TestSpaceType::PatternType PatternType;

  typedef Dune::Detailed::Discretizations
      ::DiscreteFunction
      ::Default< LagrangeSpaceType, VectorType >
    DiscreteFunctionType;

  typedef typename Dune::Detailed::Discretizations
      ::DiscreteFunctionSpace
      ::Sub
      ::Affine
      ::Dirichlet< TestSpaceType, VectorType >
    AnsatzSpaceType;

  static const std::string id()
  {
    return "detailed.solvers.stationary.linear.elliptic.continuousgalerkin";
  }

  DetailedDiscretizations(const Dune::shared_ptr< const GridPartType > _gridPart,
                          const Dune::shared_ptr< const BoundaryInfoType > _boundaryInfo,
                          const Dune::shared_ptr< const ModelType > _model)
    : gridPart_(_gridPart)
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
  } // DetailedDiscretizations

  const Dune::shared_ptr< const GridPartType > gridPart() const
  {
    return gridPart_;
  }

  const Dune::shared_ptr< const BoundaryInfoType > boundaryInfo() const
  {
    return boundaryInfo_;
  }

  const Dune::shared_ptr< const ModelType > model() const
  {
    return model_;
  }

  void init(const std::string prefix = "",
            std::ostream& out = Dune::Stuff::Common::Logger().debug())
  {
    if (!initialized_) {
      Dune::Timer timer;

      out << prefix << "initializing discrete function spaces... " << std::flush;
      lagrangeSpace_ = Dune::shared_ptr< const LagrangeSpaceType >(new LagrangeSpaceType(*gridPart_));
      testSpace_ = Dune::shared_ptr< const TestSpaceType >(new TestSpaceType(*lagrangeSpace_, boundaryInfo_));
      Dune::shared_ptr< DiscreteFunctionType > discreteDirichlet(new DiscreteFunctionType(*lagrangeSpace_,
                                                                                          "dirichlet"));
      Dune::Stuff::DiscreteFunction::Projection::Dirichlet::project(*boundaryInfo_,
                                                                    *(model_->dirichlet()),
                                                                    *discreteDirichlet);
      ansatzSpace_ = Dune::shared_ptr< const AnsatzSpaceType >(new AnsatzSpaceType(*testSpace_, discreteDirichlet));
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      out << prefix << "initializing operator and functionals... " << std::flush;
      timer.reset();
      // * left hand side
      //   * elliptic operator
      typedef Dune::Detailed::Discretizations
          ::Evaluation
          ::Local
          ::Binary
          ::Elliptic< FunctionSpaceType, typename ModelType::DiffusionType >
        EllipticEvaluationType;
      const EllipticEvaluationType ellipticEvaluation(model_->diffusion(), model_->diffusionOrder());
      typedef Dune::Detailed::Discretizations
          ::DiscreteOperator
          ::Local
          ::Codim0
          ::Integral< EllipticEvaluationType >
        EllipticOperatorType;
      const EllipticOperatorType ellipticOperator(ellipticEvaluation);
      // * right hand side
      //   * L2 force functional
      typedef Dune::Detailed::Discretizations
          ::Evaluation
          ::Local
          ::Unary
          ::Scale< FunctionSpaceType, typename ModelType::ForceType >
        ForceEvaluationType;
      const ForceEvaluationType forceEvaluation(model_->force(), model_->forceOrder());
      typedef Dune::Detailed::Discretizations
          ::DiscreteFunctional
          ::Local
          ::Codim0
          ::Integral< ForceEvaluationType >
        L2VolumeFunctionalType;
      const L2VolumeFunctionalType forceFunctional(forceEvaluation);
      //   * L2 neumann functional
      typedef Dune::Detailed::Discretizations
          ::Evaluation
          ::Local
          ::Unary
          ::Scale< FunctionSpaceType, typename ModelType::NeumannType >
        NeumannEvaluationType;
      const NeumannEvaluationType neumannEvaluation(model_->neumann(), model_->neumannOrder());
      typedef typename Dune::Detailed::Discretizations
          ::DiscreteFunctional
          ::Local
          ::Codim1
          ::Integral
          ::Boundary< NeumannEvaluationType >
        L2BoundaryFunctionalType;
      const L2BoundaryFunctionalType neumannFunctional(neumannEvaluation);
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      out << prefix << "initializing matrix (of size " << testSpace_->map().size() << "x" << ansatzSpace_->map().size()
           << ") and vectors... " << std::flush;
      timer.reset();
      // * create the left hand side matrix
      //   * therefore create the pattern
      Dune::shared_ptr< const PatternType > diffusionPattern = testSpace_->computePattern(*ansatzSpace_);
      patterns_.insert(std::pair< const std::string, Dune::shared_ptr< const PatternType > >(
                         "diffusion",
                         diffusionPattern));
      //   * and the matrix
      Dune::shared_ptr< MatrixType > diffusionMatrix = Dune::shared_ptr< MatrixType >(
            new MatrixType(testSpace_->map().size(), ansatzSpace_->map().size(), *diffusionPattern));
      matrices_.insert(std::pair< const std::string, Dune::shared_ptr< MatrixType > >(
                         "diffusion",
                         diffusionMatrix));
      // create the right hand side vectors
      Dune::shared_ptr< VectorType > forceVector = ContainerFactory::createDenseVector(*testSpace_);
      vectors_.insert(std::pair< const std::string, Dune::shared_ptr< VectorType > >(
                        "force",
                        forceVector));
      Dune::shared_ptr< VectorType > neumannVector = ContainerFactory::createDenseVector(*testSpace_);
      vectors_.insert(std::pair< const std::string, Dune::shared_ptr< VectorType > >(
                        "neumann",
                        neumannVector));
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      out << prefix << "assembing system... " << std::flush;
      timer.reset();
      // * local matrix assembler
      typedef Dune::Detailed::Discretizations::Assembler::Local::Codim0::Matrix< EllipticOperatorType >
          LocalMatrixAssemblerType;
      const Dune::shared_ptr< const LocalMatrixAssemblerType > localMatrixAssembler(
            new LocalMatrixAssemblerType(ellipticOperator));
      // * local vector assemblers
      //   * force vector
      typedef Dune::Detailed::Discretizations::Assembler::Local::Codim0::Vector< L2VolumeFunctionalType >
          LocalVolumeVectorAssemblerType;
      const Dune::shared_ptr< const LocalVolumeVectorAssemblerType > localforceVectorAssembler(
            new LocalVolumeVectorAssemblerType(forceFunctional));
      //   * neumann vector
      typedef Dune::Detailed::Discretizations::Assembler::Local::Codim1::Vector::Neumann< L2BoundaryFunctionalType,
                                                                                          BoundaryInfoType >
          LocalNeumannVectorAssemblerType;
      const Dune::shared_ptr< const LocalNeumannVectorAssemblerType > localNeumannVectorAssembler(
            new LocalNeumannVectorAssemblerType(neumannFunctional, boundaryInfo_));
      typedef Dune::Detailed::Discretizations::Assembler::System< TestSpaceType, AnsatzSpaceType > SystemAssemblerType;
      // * system assembler
      SystemAssemblerType systemAssembler(*testSpace_, *ansatzSpace_);
      systemAssembler.addLocalMatrixAssembler(localMatrixAssembler, diffusionMatrix);
      systemAssembler.addLocalVectorAssembler(localforceVectorAssembler, forceVector);
      systemAssembler.addLocalVectorAssembler(localNeumannVectorAssembler, neumannVector);
      systemAssembler.assemble();

      // done
      initialized_ = true;
    } // if !(initialized_)
  } // void init()

//  Dune::shared_ptr< VectorBackendType > createVector() const
//  {
//    assert(initialized_ && "A vector can only be created after init() has been called!");
//    Dune::shared_ptr< VectorBackendType > ret(new VectorBackendType(ContainerFactory::createDenseVector(*ansatzSpace_)));
//    return ret;
//  }

//  Dune::shared_ptr< DiscreteFunctionType > createDiscreteFunction(const std::string name = "discrete_function") const
//  {
//    assert(initialized_ && "Please call init() beafore calling createDiscreteFunction()!");
//    return createDiscreteFunction(ContainerFactory::createDenseVector(*ansatzSpace_), name);
//  } // ... createDiscreteFunction(...)

//  Dune::shared_ptr< DiscreteFunctionType > createDiscreteFunction(VectorBackendType vector,
//                                                                  const std::string name = "discrete_function") const
//  {
//    assert(initialized_ && "Please call init() beafore calling createDiscreteFunction()!");
//    assert(vector.size() == ansatzSpace_->map().size());
//    Dune::shared_ptr< DiscreteFunctionType > discreteFunction(new DiscreteFunctionType(*ansatzSpace_,
//                                                                                       vector,
//                                                                                       name));
//    return discreteFunction;
//  } // ... createDiscreteFunction(...)

//  void solve(VectorBackendType& solution,
//             const std::string type = "eigen.bicgstab.incompletelut",
//             const unsigned int maxIter = 5000,
//             const double precision = 1e-12,
//             const std::string prefix = "",
//             std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
//  {
//    // preparations
//    assert(initialized_ && "The system can only be solved after init() has been called! ");
//    assert(solution.size() == ansatzSpace_->map().size() && "Given vector has wrong size!");
//    Dune::Timer timer;
//    // solve
//    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::BicgstabIlut BicgstabIlutSolver;
//    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::BicgstabDiagonal BicgstabDiagonalSolver;
//    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::CgDiagonalUpper CgDiagonalUpperSolver;
//    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::CgDiagonalLower CgDiagonalLowerSolver;
//    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::SimplicialcholeskyUpper SimplicialcholeskyUpperSolver;
//    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::SimplicialcholeskyLower SimplicialcholeskyLowerSolver;
//    out << prefix << "solving linear system of size " << matrix_->rows() << "x" << matrix_->cols() << std::endl
//        << prefix << "  using " << type << "... " << std::flush;
//    if (type == "eigen.bicgstab.incompletelut"){
//      BicgstabIlutSolver::apply(*matrix_, solution, *rhs_, maxIter, precision);
//    } else if (type == "eigen.bicgstab.diagonal"){
//      BicgstabDiagonalSolver::apply(*matrix_, solution, *rhs_, maxIter, precision);
//    } else if (type == "eigen.cg.diagonal.upper"){
//      CgDiagonalUpperSolver::apply(*matrix_, solution, *rhs_, maxIter, precision);
//    } else if (type == "eigen.cg.diagonal.lower"){
//      CgDiagonalLowerSolver::apply(*matrix_, solution, *rhs_, maxIter, precision);
//    } else if (type == "eigen.simplicialcholesky.upper"){
//      SimplicialcholeskyUpperSolver::apply(*matrix_, solution, *rhs_, maxIter, precision);
//    } else if (type == "eigen.simplicialcholesky.lower"){
//      SimplicialcholeskyLowerSolver::apply(*matrix_, solution, *rhs_, maxIter, precision);
//    } else {
//      DUNE_THROW(Dune::RangeError,
//                 "\nError: wrong 'type' given (" << type << ")!");
//    }
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
//    assert(vector.size() == ansatzSpace_->map().size() && "Given vector has wrong size!");
//    Dune::Timer timer;
//    out << prefix << "writing '" << name << "' to '" << filename;
//    if (dimDomain == 1)
//      out << ".vtp";
//    else
//      out << ".vtu";
//    out << "'... " << std::flush;
//    Dune::shared_ptr< DiscreteFunctionType > discreteFunction(new DiscreteFunctionType(*ansatzSpace_, vector, name));
//    visualize(discreteFunction, filename, "", Dune::Stuff::Common::Logger().devnull());
//    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
//  }

//  void visualize(Dune::shared_ptr< DiscreteFunctionType > discreteFunction,
//                 const std::string filename = "solution",
//                 const std::string prefix = "",
//                 std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
//  {
//    // preparations
//    assert(initialized_ && "A discrete function can only be visualized after init() has been called! ");
//    assert(discreteFunction->size() == ansatzSpace_->map().size() && "Given vector has wrong size!");
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

//  const AnsatzSpaceType& ansatzSpace() const
//  {
//    assert(initialized_);
//    return *ansatzSpace_;
//  }

//  const TestSpaceType& testSpace() const
//  {
//    assert(initialized_);
//    return *testSpace_;
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
  const Dune::shared_ptr< const GridPartType > gridPart_;
  const Dune::shared_ptr< const BoundaryInfoType > boundaryInfo_;
  const Dune::shared_ptr< const ModelType > model_;
  bool initialized_;
  Dune::shared_ptr< const LagrangeSpaceType > lagrangeSpace_;
  Dune::shared_ptr< const TestSpaceType > testSpace_;
  Dune::shared_ptr< const AnsatzSpaceType > ansatzSpace_;
  std::map< const std::string, Dune::shared_ptr< const PatternType > > patterns_;
  std::map< const std::string, Dune::shared_ptr< MatrixType > > matrices_;
  std::map< const std::string, Dune::shared_ptr< VectorType > > vectors_;
}; // class DetailedDiscretizations

} // namespace CG
} // namespace Elliptic
} // namespace Linear
} // namespace Stationary
} // namespace Solvers
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_CG_DETAILED_DISCRETIZATIONS_HH
