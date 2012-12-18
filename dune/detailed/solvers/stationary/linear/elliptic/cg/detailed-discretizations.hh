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
#include <dune/stuff/la/solver/eigen.hh>

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

#if HAVE_DUNE_RB
  #include <dune/rb/model/stationary/linear/elliptic/interface.hh>
#endif // HAVE_DUNE_RB

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

  typedef typename ContainerFactory::RowMajorSparseMatrixType MatrixType;

  typedef typename TestSpaceType::PatternType PatternType;

  typedef Dune::Detailed::Discretizations
      ::DiscreteFunction
      ::Default< TestSpaceType, VectorType >
    DiscreteTestFunctionType;

private:
  typedef Dune::Detailed::Discretizations
      ::DiscreteFunction
      ::Default< LagrangeSpaceType, VectorType >
    DiscreteFunctionType;

public:
  typedef typename Dune::Detailed::Discretizations
      ::DiscreteFunctionSpace
      ::Sub
      ::Affine
      ::Dirichlet< TestSpaceType, VectorType >
    AnsatzSpaceType;

  typedef Dune::Detailed::Discretizations
      ::DiscreteFunction
      ::Default< AnsatzSpaceType, VectorType >
    DiscreteAnsatzFunctionType;

  static const std::string id()
  {
    return "detailed.solvers.stationary.linear.elliptic.cg.detailed_discretizations";
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
    return *testSpace_;
  }

  Dune::shared_ptr< VectorType > createAnsatzVector() const
  {
    assert(initialized_ && "Please call init() before calling createAnsatzVector()!");
    Dune::shared_ptr< VectorType > vector(ContainerFactory::createDenseVector(*ansatzSpace_));
    return vector;
  } // Dune::shared_ptr< VectorType > createAnsatzVector() const

  Dune::shared_ptr< VectorType > createTestVector() const
  {
    assert(initialized_ && "Please call init() before calling createAnsatzVector()!");
    Dune::shared_ptr< VectorType > vector(ContainerFactory::createDenseVector(*testSpace_));
    return vector;
  } // Dune::shared_ptr< VectorType > createAnsatzVector() const

  Dune::shared_ptr< DiscreteAnsatzFunctionType > createAnsatzFunction(const std::string name = "ansatzFunction") const
  {
    assert(initialized_ && "Please call init() before calling createAnsatzFunction()!");
    Dune::shared_ptr< DiscreteAnsatzFunctionType > ansatzFunction(new DiscreteAnsatzFunctionType(
                                                                    *ansatzSpace_,
                                                                    name));
    return ansatzFunction;
  } // Dune::shared_ptr< DiscreteAnsatzFunctionType > createAnsatzFunction(...) const

  Dune::shared_ptr< DiscreteAnsatzFunctionType > createAnsatzFunction(Dune::shared_ptr< VectorType > vector,
                                                                      const std::string name = "ansatzFunction") const
  {
    assert(initialized_ && "Please call init() before calling createAnsatzFunction()!");
    Dune::shared_ptr< DiscreteAnsatzFunctionType > ansatzFunction(new DiscreteAnsatzFunctionType(
                                                                    *ansatzSpace_,
                                                                    vector,
                                                                    name));
    return ansatzFunction;
  } // Dune::shared_ptr< DiscreteAnsatzFunctionType > createAnsatzFunction(...) const

  Dune::shared_ptr< DiscreteTestFunctionType > createTestFunction(const std::string name = "testFunction") const
  {
    assert(initialized_ && "Please call init() before calling createTestFunction()!");
    Dune::shared_ptr< DiscreteTestFunctionType > testFunction(new DiscreteTestFunctionType(
                                                                *testSpace_,
                                                                name));
    return testFunction;
  } // Dune::shared_ptr< DiscreteAnsatzFunctionType > createTestFunction(...) const

  Dune::shared_ptr< DiscreteTestFunctionType > createTestFunction(Dune::shared_ptr< VectorType > vector,
                                                                  const std::string name = "testFunction") const
  {
    assert(initialized_ && "Please call init() before calling createTestFunction()!");
    Dune::shared_ptr< DiscreteTestFunctionType > testFunction(new DiscreteTestFunctionType(
                                                                *testSpace_,
                                                                vector,
                                                                name));
    return testFunction;
  } // Dune::shared_ptr< DiscreteTestFunctionType > createTestFunction(...) const

  void solve(VectorType& solutionVector,
             const std::string linearSolverType = "eigen.iterative.bicgstab.diagonal",
             const unsigned int linearSolverMaxIter = 5000,
             const double linearSolverPrecision = 1e-12,
             const std::string prefix = "",
             std::ostream& out = Dune::Stuff::Common::Logger().debug())
  {
    assert(initialized_ && "Please call init() before calling solve()!");
    out << prefix << "applying constraints... " << std::flush;
    Dune::Timer timer;
    // compute system matrix and right hand side
    // * therefore get matrix and vectors
    MatrixType& systemMatrix = *(matrices_["diffusion"]);
    const VectorType& forceVector = *(vectors_["force"]);
    const VectorType& neumannVector = *(vectors_["neumann"]);
    const VectorType& dirichletVector = *(ansatzSpace_->affineShift()->vector());
    // * compute right hand side
    VectorType rhsVector(testSpace_->map().size());
    rhsVector.base() = forceVector.base()
        + neumannVector.base()
        - systemMatrix.base() * dirichletVector.base();
    // * create a system assembler
    typedef Dune::Detailed::Discretizations::Assembler::System< TestSpaceType, AnsatzSpaceType > SystemAssemblerType;
    SystemAssemblerType systemAssembler(*testSpace_, *ansatzSpace_);
    // * and apply the constraints
    systemAssembler.applyConstraints(systemMatrix, rhsVector);
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;

    out << prefix << "solving linear system (of size " << systemMatrix.rows()
        << "x" << systemMatrix.cols() << ")" << std::endl;
    out << prefix << "  using '" << linearSolverType << "'... " << std::flush;
    timer.reset();
    typedef typename Dune::Stuff::LA::Solver::Eigen::Interface< MatrixType, VectorType > SolverType;
    SolverType* solver = Dune::Stuff::LA::Solver::Eigen::create< MatrixType, VectorType >(linearSolverType);
    solver->init(systemMatrix);
    const bool success = solver->apply(rhsVector,
                                       solutionVector,
                                       linearSolverMaxIter,
                                       linearSolverPrecision);
    if (!success)
      DUNE_THROW(Dune::MathError,
                 "\nERROR: linear solver '" << linearSolverType << "' reported a problem!");
    if (solutionVector.size() != ansatzSpace_->map().size())
      DUNE_THROW(Dune::MathError,
                 "\nERROR: linear solver '" << linearSolverType << "' produced a solution of wrong size (is "
                 << solutionVector.size() << ", should be " << ansatzSpace_->map().size() << ")!");
    solutionVector.base() += dirichletVector.base();
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void solve(...)

  void visualizeAnsatzVector(VectorType& vector,
                             const std::string filename = id() + ".ansatzVector",
                             const std::string name = id() + "ansatzVector",
                             const std::string prefix = "",
                             std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
  {
    // preparations
    assert(initialized_ && "Please call init() before calling visualize()");
    assert(vector.size() == ansatzSpace_->map().size() && "Given vector has wrong size!");
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
                           const std::string name = id() + "testVector",
                           const std::string prefix = "",
                           std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
  {
    // preparations
    assert(initialized_ && "Please call init() before calling visualize()");
    assert(vector.size() == ansatzSpace_->map().size() && "Given vector has wrong size!");
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
    typedef Dune::VTKWriter< typename DiscreteAnsatzFunctionType::DiscreteFunctionSpaceType::GridViewType > VTKWriterType;
    VTKWriterType vtkWriter(discreteFunction->space().gridView());
    vtkWriter.addVertexData(discreteFunction);
    vtkWriter.write(filename);
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void visualizeFunction(...)

  void visualizeFunction(const Dune::shared_ptr< const DiscreteTestFunctionType > discreteFunction,
                         const std::string filename = id() + ".discreteTestFunction",
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
    typedef Dune::VTKWriter< typename DiscreteAnsatzFunctionType::DiscreteFunctionSpaceType::GridViewType > VTKWriterType;
    VTKWriterType vtkWriter(discreteFunction->space().gridView());
    vtkWriter.addVertexData(discreteFunction);
    vtkWriter.write(filename);
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void visualizeFunction(...)

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


#if HAVE_DUNE_RB
namespace Parametric {

template< class GridPartImp, int polynomialOrder,
          class RangeFieldImp, int dimensionRange,
          class ParamFieldImp, int maxNumParams >
class DetailedDiscretizations;

template< class GridPartImp, int polynomialOrder, class RangeFieldImp, class ParamFieldImp, int maxNumParams >
class DetailedDiscretizations< GridPartImp, polynomialOrder, RangeFieldImp, 1, ParamFieldImp, maxNumParams >
{
public:
  typedef DetailedDiscretizations< GridPartImp, polynomialOrder, RangeFieldImp, 1, ParamFieldImp, maxNumParams > ThisType;

  typedef Dune::grid::Part::Interface< typename GridPartImp::Traits > GridPartType;

  static const int polOrder = polynomialOrder;

  typedef typename GridPartType::ctype DomainFieldType;

  static const int dimDomain = GridPartType::dimension;

  typedef RangeFieldImp RangeFieldType;

  static const int dimRange = 1;

  typedef ParamFieldImp ParamFieldType;

  static const int maxParams = maxNumParams;

  typedef Dune::RB
      ::Model
      ::Stationary
      ::Linear
      ::Elliptic::Interface< DomainFieldType, dimDomain, RangeFieldType, dimRange, ParamFieldType, maxParams >
    ModelType;

  typedef typename ModelType::ParamType ParamType;

  typedef Dune::Stuff::Grid::BoundaryInfo::Interface< typename GridPartType::GridViewType > BoundaryInfoType;

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
      ::Default< TestSpaceType, VectorType >
    DiscreteTestFunctionType;

private:
  typedef Dune::Detailed::Discretizations
      ::DiscreteFunction
      ::Default< LagrangeSpaceType, VectorType >
    DiscreteFunctionType;

public:
  typedef typename Dune::Detailed::Discretizations
      ::DiscreteFunctionSpace
      ::Sub
      ::Affine
      ::Dirichlet< TestSpaceType, VectorType >
    AnsatzSpaceType;

  typedef Dune::Detailed::Discretizations
      ::DiscreteFunction
      ::Default< AnsatzSpaceType, VectorType >
    DiscreteAnsatzFunctionType;

  static const std::string id()
  {
    return "detailed.solvers.stationary.linear.elliptic.cg.parametric.detailed_discretizations";
  }

  DetailedDiscretizations(const shared_ptr< const ModelType > _model,
                          const shared_ptr< const GridPartType > _gridPart,
                          const shared_ptr< const BoundaryInfoType > _boundaryInfo)
    : model_(_model)
    , gridPart_(_gridPart)
    , boundaryInfo_(_boundaryInfo)
    , initialized_(false)
  {
    // allow only separable parametric or nonparametric data functions
    unsigned int throw_up = 0u;
    std::stringstream msg;
    msg << "\nERROR: only separable or nonparametric data functions allowed!" << std::endl;
    if (model_->diffusion()->parametric() && !model_->diffusion()->separable()) {
      ++throw_up;
      msg << "       - But 'model.diffusion()' is not!" << std::endl;
    }
    if (model_->force()->parametric() && !model_->force()->separable()) {
      ++throw_up;
      msg << "       - But 'model.force()' is not!" << std::endl;
    }
    if (model_->dirichlet()->parametric() && !model_->dirichlet()) {
      ++throw_up;
      msg << "       - But 'model.dirichlet()' is not!" << std::endl;
    }
    if (model_->neumann()->parametric() && !model_->neumann()) {
      ++throw_up;
      msg << "       - But 'model.neumann()' is not!" << std::endl;
    }
    if (throw_up)
      DUNE_THROW(Dune::InvalidStateException, msg.str());
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
  } // DuneDetailedSolvers(...)

  const shared_ptr< const ModelType > model() const
  {
    return model_;
  }

  const shared_ptr< const GridPartType > gridPart() const
  {
    return gridPart_;
  }

  const shared_ptr< const BoundaryInfoType > boundaryInfo() const
  {
    return boundaryInfo_;
  }

  void init(const std::string prefix = "",
            std::ostream& out = Dune::Stuff::Common::Logger().debug())
  {
    if (!initialized_) {
      Dune::Timer timer;

      if (model_->dirichlet()->parametric())
        DUNE_THROW(Dune::InvalidStateException,
                   "\nERROR: not implemented for parametric dirichlet values!");

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

      out << prefix << "initializing operators and functionals:" << std::endl;
      timer.reset();
      // * left hand side
      //   * elliptic operator
      typedef Dune::Detailed::Discretizations
          ::Evaluation
          ::Local
          ::Binary
          ::Elliptic< FunctionSpaceType, typename ModelType::DiffusionType::ComponentType >
        EllipticEvaluationType;
      typedef Dune::Detailed::Discretizations
          ::DiscreteOperator
          ::Local
          ::Codim0
          ::Integral< EllipticEvaluationType >
        EllipticOperatorType;
      std::vector< Dune::shared_ptr< const typename ModelType::DiffusionType::ComponentType > > diffusions;
      std::vector< const EllipticEvaluationType* > ellipticEvaluations;
      std::vector< const EllipticOperatorType* > ellipticOperators;
      if (model_->diffusion()->separable()) {
        out << prefix << "  " << model_->diffusion()->numComponents() << " diffusion operators... "
            << std::flush;
        for (unsigned int qq = 0; qq < model_->diffusion()->numComponents(); ++qq) {
          diffusions.push_back(model_->diffusion()->components()[qq]);
        }
      } else {
        out << prefix << "  1 diffusion operator... " << std::flush;
        diffusions.push_back(model_->diffusion());
      }
      for (unsigned int qq = 0; qq < diffusions.size(); ++qq) {
        ellipticEvaluations.push_back(new EllipticEvaluationType(diffusions[qq],
                                                                 model_->diffusionOrder()));
        ellipticOperators.push_back(new EllipticOperatorType(*(ellipticEvaluations[qq])));
      }
      out << "done" << std::endl;
      // * right hand side
      //   * L2 force functional
      typedef Dune::Detailed::Discretizations
          ::Evaluation
          ::Local
          ::Unary
          ::Scale< FunctionSpaceType, typename ModelType::ForceType::ComponentType >
        ForceEvaluationType;
      typedef Dune::Detailed::Discretizations
          ::DiscreteFunctional
          ::Local
          ::Codim0
          ::Integral< ForceEvaluationType >
        L2ForceFunctionalType;
      std::vector< Dune::shared_ptr< const typename ModelType::ForceType::ComponentType > > forces;
      std::vector< const ForceEvaluationType* > forceEvaluations;
      std::vector< const L2ForceFunctionalType* > forceFunctionals;
      if (model_->force()->separable()) {
        out << prefix << "  " << model_->force()->numComponents() << " force functionals... "
            << std::flush;
        for (unsigned int qq = 0; qq < model_->force()->numComponents(); ++qq) {
          forces.push_back(model_->force()->components()[qq]);
        }
      } else {
        out << prefix << "  1 force functional... " << std::flush;
        forces.push_back(model_->force());
      }
      for (unsigned int qq = 0; qq < forces.size(); ++qq) {
        forceEvaluations.push_back(new ForceEvaluationType(forces[qq],
                                                           model_->forceOrder()));
        forceFunctionals.push_back(new L2ForceFunctionalType(*(forceEvaluations[qq])));
      }
      out << "done" << std::endl;
      //   * L2 neumann functional
      typedef Dune::Detailed::Discretizations
          ::Evaluation
          ::Local
          ::Unary
          ::Scale< FunctionSpaceType, typename ModelType::NeumannType::ComponentType >
        NeumannEvaluationType;
      typedef typename Dune::Detailed::Discretizations
          ::DiscreteFunctional
          ::Local
          ::Codim1
          ::Integral
          ::Boundary< NeumannEvaluationType >
        L2NeumannFunctionalType;
      std::vector< Dune::shared_ptr< const typename ModelType::NeumannType::ComponentType > > neumanns;
      std::vector< const NeumannEvaluationType* > neumannEvaluations;
      std::vector< const L2NeumannFunctionalType* > neumannFunctionals;
      if (model_->neumann()->separable()) {
        out << prefix << "  " << model_->neumann()->numComponents() << " neumann functionals... "
            << std::flush;
        for (unsigned int qq = 0; qq < model_->neumann()->numComponents(); ++qq) {
          neumanns.push_back(model_->neumann()->components()[qq]);
        }
      } else {
        out << prefix << "  1 neumann functional... " << std::flush;
        neumanns.push_back(model_->neumann());
      }
      for (unsigned int qq = 0; qq < neumanns.size(); ++qq) {
        neumannEvaluations.push_back(new NeumannEvaluationType(neumanns[qq],
                                                               model_->neumannOrder()));
        neumannFunctionals.push_back(new L2NeumannFunctionalType(*(neumannEvaluations[qq])));
      }
      out << "done" << std::endl;

      out << prefix << "initializing matrices and vectors:" << std::endl;
      timer.reset();
      // * create left hand side matrices
      //   * therefore create the pattern
      Dune::shared_ptr< const PatternType > diffusionPattern = testSpace_->computePattern(*ansatzSpace_);
      //   * and the matrices
      if (ellipticOperators.size() == 1)
        out << prefix << "  1 diffusion matrix... " << std::flush;
      else
        out << prefix << "  " << ellipticOperators.size() << " diffusion matrices... " << std::flush;
      for (unsigned int qq = 0; qq < ellipticOperators.size(); ++ qq) {
        patterns_.insert(std::pair< const std::string, Dune::shared_ptr< const PatternType > >(
                           "diffusion_" + Dune::Stuff::Common::toString(qq),
                           diffusionPattern));
        Dune::shared_ptr< MatrixType > diffusionMatrix = Dune::shared_ptr< MatrixType >(
              new MatrixType(testSpace_->map().size(), ansatzSpace_->map().size(), *diffusionPattern));
        matrices_.insert(std::pair< const std::string, Dune::shared_ptr< MatrixType > >(
                           "diffusion_" + Dune::Stuff::Common::toString(qq),
                           diffusionMatrix));
      }
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;
      // create the right hand side vectors
      //   * for the force
      if (forceFunctionals.size() == 1)
        out << prefix << "  1 force vector... " << std::flush;
      else
        out << prefix << "  " << forceFunctionals.size() << " force vectors... " << std::flush;
      timer.reset();
      for (unsigned int qq = 0; qq < forceFunctionals.size(); ++qq)
        vectors_.insert(std::pair< const std::string, Dune::shared_ptr< VectorType > >(
                          "force_" + Dune::Stuff::Common::toString(qq),
                          Dune::shared_ptr< VectorType >(ContainerFactory::createDenseVector(*testSpace_))));
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;
      //   * for the neumann
      if (neumannFunctionals.size() == 1)
        out << prefix << "  1 neumann vector... " << std::flush;
      else
        out << prefix << "  " << neumannFunctionals.size() << " neumann vectors... " << std::flush;
      timer.reset();
      for (unsigned int qq = 0; qq < neumannFunctionals.size(); ++qq)
        vectors_.insert(std::pair< const std::string, Dune::shared_ptr< VectorType > >(
                          "neumann_" + Dune::Stuff::Common::toString(qq),
                          Dune::shared_ptr< VectorType >(ContainerFactory::createDenseVector(*testSpace_))));
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      out << prefix << "assembing system... " << std::flush;
      timer.reset();
      typedef Dune::Detailed::Discretizations::Assembler::System< TestSpaceType, AnsatzSpaceType > SystemAssemblerType;
      SystemAssemblerType systemAssembler(*testSpace_, *ansatzSpace_);
      // * local matrix assembler
      typedef Dune::Detailed::Discretizations::Assembler::Local::Codim0::Matrix< EllipticOperatorType >
          LocalMatrixAssemblerType;
      for (unsigned int qq = 0; qq < ellipticOperators.size(); ++qq) {
        const Dune::shared_ptr< const LocalMatrixAssemblerType > localMatrixAssembler(
              new LocalMatrixAssemblerType(*(ellipticOperators[qq])));
        Dune::shared_ptr< MatrixType > diffusionMatrix
            = matrices_.find("diffusion_" + Dune::Stuff::Common::toString(qq))->second;
        systemAssembler.addLocalMatrixAssembler(localMatrixAssembler, diffusionMatrix);
      }
      // * local vector assemblers
      //   * force vector
      typedef Dune::Detailed::Discretizations::Assembler::Local::Codim0::Vector< L2ForceFunctionalType >
          LocalForceVectorAssemblerType;
      for (unsigned int qq = 0; qq < forceFunctionals.size(); ++qq) {
        const Dune::shared_ptr< const LocalForceVectorAssemblerType > localForceVectorAssembler(
              new LocalForceVectorAssemblerType(*(forceFunctionals[qq])));
        Dune::shared_ptr< VectorType > forceVector
            = vectors_.find("force_" + Dune::Stuff::Common::toString(qq))->second;
        systemAssembler.addLocalVectorAssembler(localForceVectorAssembler, forceVector);
      }
      //   * neumann vector
      typedef Dune::Detailed::Discretizations::Assembler::Local::Codim1::Vector::Neumann< L2NeumannFunctionalType,
                                                                                          BoundaryInfoType >
          LocalNeumannVectorAssemblerType;
      for (unsigned int qq = 0; qq < neumannFunctionals.size(); ++qq) {
        const Dune::shared_ptr< const LocalNeumannVectorAssemblerType > localNeumannVectorAssembler(
              new LocalNeumannVectorAssemblerType(*(neumannFunctionals[qq]), boundaryInfo_));
        Dune::shared_ptr< VectorType > neumannVector
            = vectors_.find("neumann_" + Dune::Stuff::Common::toString(qq))->second;
        systemAssembler.addLocalVectorAssembler(localNeumannVectorAssembler, neumannVector);
      }
      // * system assembler
      systemAssembler.assemble();
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      initialized_ = true;
    }
  } // void init(...)

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

  Dune::shared_ptr< VectorType > createAnsatzVector() const
  {
    assert(initialized_ && "Please call init() before calling createAnsatzVector()!");
    Dune::shared_ptr< VectorType > vector(ContainerFactory::createDenseVector(*ansatzSpace_));
    return vector;
  } // Dune::shared_ptr< VectorType > createAnsatzVector() const

  Dune::shared_ptr< VectorType > createTestVector() const
  {
    assert(initialized_ && "Please call init() before calling createAnsatzVector()!");
    Dune::shared_ptr< VectorType > vector(ContainerFactory::createDenseVector(*testSpace_));
    return vector;
  } // Dune::shared_ptr< VectorType > createAnsatzVector() const

  Dune::shared_ptr< DiscreteAnsatzFunctionType > createAnsatzFunction(const std::string name = "ansatzFunction") const
  {
    assert(initialized_ && "Please call init() before calling createAnsatzFunction()!");
    Dune::shared_ptr< DiscreteAnsatzFunctionType > ansatzFunction(new DiscreteAnsatzFunctionType(
                                                                    *ansatzSpace_,
                                                                    name));
    return ansatzFunction;
  } // Dune::shared_ptr< DiscreteAnsatzFunctionType > createAnsatzFunction(...) const

  Dune::shared_ptr< DiscreteAnsatzFunctionType > createAnsatzFunction(Dune::shared_ptr< VectorType > vector,
                                                                      const std::string name = "ansatzFunction") const
  {
    assert(initialized_ && "Please call init() before calling createAnsatzFunction()!");
    Dune::shared_ptr< DiscreteAnsatzFunctionType > ansatzFunction(new DiscreteAnsatzFunctionType(
                                                                    *ansatzSpace_,
                                                                    vector,
                                                                    name));
    return ansatzFunction;
  } // Dune::shared_ptr< DiscreteAnsatzFunctionType > createAnsatzFunction(...) const

  Dune::shared_ptr< DiscreteTestFunctionType > createTestFunction(const std::string name = "testFunction") const
  {
    assert(initialized_ && "Please call init() before calling createTestFunction()!");
    Dune::shared_ptr< DiscreteTestFunctionType > testFunction(new DiscreteTestFunctionType(
                                                                *testSpace_,
                                                                name));
    return testFunction;
  } // Dune::shared_ptr< DiscreteAnsatzFunctionType > createTestFunction(...) const

  Dune::shared_ptr< DiscreteTestFunctionType > createTestFunction(Dune::shared_ptr< VectorType > vector,
                                                                  const std::string name = "testFunction") const
  {
    assert(initialized_ && "Please call init() before calling createTestFunction()!");
    Dune::shared_ptr< DiscreteTestFunctionType > testFunction(new DiscreteTestFunctionType(
                                                                *testSpace_,
                                                                vector,
                                                                name));
    return testFunction;
  } // Dune::shared_ptr< DiscreteTestFunctionType > createTestFunction(...) const

  void solve(const ParamType& mu,
             VectorType& solutionVector,
             const std::string linearSolverType = "eigen.iterative.bicgstab.diagonal",
             const unsigned int linearSolverMaxIter = 5000,
             const double linearSolverPrecision = 1e-12,
             const std::string prefix = "",
             std::ostream& out = Dune::Stuff::Common::Logger().debug())
  {
    assert(initialized_ && "Please call init() before calling solve()!");
    out << prefix << "computing system matrix and right hand side... " << std::flush;
    Dune::Timer timer;
    // compute the system matrix
    Dune::shared_ptr< MatrixType > systemMatrix;
    if (model_->diffusion()->separable()) {
      const ParamType muDiffusion = model_->getDiffusionParam(mu);
      systemMatrix = fixMatrix(*(model_->diffusion()), "diffusion", muDiffusion);
    } else {
      systemMatrix = Dune::shared_ptr< MatrixType >(new MatrixType(testSpace_->map().size(),
                                                                   ansatzSpace_->map().size()));
      systemMatrix->base() = matrices_.find("diffusion_0")->second->base();
    }
    // compute the right hand side
    Dune::shared_ptr< VectorType > rhsVector = ContainerFactory::createDenseVector(*testSpace_);
    // * add up force
    if (model_->force()->separable()) {
      const ParamType muForce = model_->getForceParam(mu);
      fixAndAddVector(*(model_->force()), "force", muForce, *rhsVector);
    } else
      rhsVector->base() = vectors_.find("force_0")->second->base();
    // * add up neumann
    if (model_->neumann()->separable()) {
      const ParamType muNeumann = model_->getNeumannParam(mu);
      fixAndAddVector(*(model_->neumann()), "neumann", muNeumann, *rhsVector);
    } else
      rhsVector->base() += vectors_.find("neumann_0")->second->base();
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;

    out << prefix << "applying constraints... " << std::flush;
    typedef Dune::Detailed::Discretizations::Assembler::System< TestSpaceType, AnsatzSpaceType > SystemAssemblerType;
    SystemAssemblerType systemAssembler(*testSpace_, *ansatzSpace_);
    systemAssembler.applyConstraints(*systemMatrix, *rhsVector);
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;

    out << prefix << "solving linear system (of size " << systemMatrix->rows()
        << "x" << systemMatrix->cols() << ")" << std::endl;
    out << prefix << "  using '" << linearSolverType << "'... " << std::flush;
    timer.reset();
    typedef typename Dune::Stuff::LA::Solver::Eigen::Interface< MatrixType > SolverType;
    SolverType* solver = Dune::Stuff::LA::Solver::Eigen::create< MatrixType >(linearSolverType);
    solver->init(*systemMatrix);
    const bool success = solver->apply(*rhsVector,
                                       solutionVector,
                                       linearSolverMaxIter,
                                       linearSolverPrecision);
    if (!success)
      DUNE_THROW(Dune::MathError,
                 "\nERROR: linear solver '" << linearSolverType << "' reported a problem!");
    if (solutionVector.size() != ansatzSpace_->map().size())
      DUNE_THROW(Dune::MathError,
                 "\nERROR: linear solver '" << linearSolverType << "' produced a solution of wrong size (is "
                 << solutionVector.size() << ", should be " << ansatzSpace_->map().size() << ")!");
//    solutionVector.base() += dirichletVector.base();
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void solve(...)

  void visualizeAnsatzVector(VectorType& vector,
                             const std::string filename = id() + ".ansatzVector",
                             const std::string name = id() + "ansatzVector",
                             const std::string prefix = "",
                             std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
  {
    // preparations
    assert(initialized_ && "Please call init() before calling visualize()");
    assert(vector.size() == ansatzSpace_->map().size() && "Given vector has wrong size!");
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
                           const std::string name = id() + "testVector",
                           const std::string prefix = "",
                           std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
  {
    // preparations
    assert(initialized_ && "Please call init() before calling visualize()");
    assert(vector.size() == ansatzSpace_->map().size() && "Given vector has wrong size!");
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
    typedef Dune::VTKWriter< typename DiscreteAnsatzFunctionType::DiscreteFunctionSpaceType::GridViewType > VTKWriterType;
    VTKWriterType vtkWriter(discreteFunction->space().gridView());
    vtkWriter.addVertexData(discreteFunction);
    vtkWriter.write(filename);
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void visualizeFunction(...)

  void visualizeFunction(const Dune::shared_ptr< const DiscreteTestFunctionType > discreteFunction,
                         const std::string filename = id() + ".discreteTestFunction",
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
    typedef Dune::VTKWriter< typename DiscreteAnsatzFunctionType::DiscreteFunctionSpaceType::GridViewType > VTKWriterType;
    VTKWriterType vtkWriter(discreteFunction->space().gridView());
    vtkWriter.addVertexData(discreteFunction);
    vtkWriter.write(filename);
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void visualizeFunction(...)

private:
  template< class FunctionType >
  Dune::shared_ptr< MatrixType > fixMatrix(const FunctionType& function,
                                           const std::string& name,
                                           const ParamType& mu)
  {
    Dune::shared_ptr< MatrixType > ret(new MatrixType(testSpace_->map().size(),
                                                      ansatzSpace_->map().size(),
                                                      *(patterns_.begin()->second)));
    // summ up all components which have a coefficient
    for (unsigned int qq = 0; qq < function.numCoefficients(); ++qq) {
      const Dune::shared_ptr< const MatrixType >
          componentMatrix = matrices_.find(name + "_" + Dune::Stuff::Common::toString(qq))->second;
      const ParamType coefficient = function.coefficients()[qq]->evaluate(mu);
      assert(coefficient.size() == 1);
      ret->base() += componentMatrix->base() * coefficient[0];
    }
    // and add the nonparametric contribution
    if (function.numComponents() > function.numCoefficients()) {
      const Dune::shared_ptr< const MatrixType >
          componentMatrix = matrices_.find(name + "_" + Dune::Stuff::Common::toString(function.numComponents()))->second;
      ret->base() += componentMatrix->base();
    }
    return ret;
  } // Dune::shared_ptr< MatrixType > completeMatrix(...)

  template< class FunctionType >
  void fixAndAddVector(const FunctionType& function,
                       const std::string& name,
                       const ParamType& mu,
                       VectorType& vector) const
  {
    for (unsigned int qq = 0; qq < function.numCoefficients(); ++qq) {
      const VectorType& componentVector = *(vectors_.find(name + "_" + Dune::Stuff::Common::toString(qq))->second);
      const ParamType coefficient = function.coefficients()[qq]->evaluate(mu);
      assert(coefficient.size() == 1);
      vector += componentVector * coefficient[0];
    }
    if (function.numComponents() > function.numCoefficients()) {
      const VectorType&
          componentVector= *(vectors_.find(name + "_" + Dune::Stuff::Common::toString(function.numComponents()))->second);
      vector += componentVector;
    }
  } // void fixAndAddVector(...)

  const shared_ptr< const ModelType > model_;
  const shared_ptr< const GridPartType > gridPart_;
  const shared_ptr< const BoundaryInfoType > boundaryInfo_;
  bool initialized_;
  Dune::shared_ptr< const LagrangeSpaceType > lagrangeSpace_;
  Dune::shared_ptr< const TestSpaceType > testSpace_;
  Dune::shared_ptr< const AnsatzSpaceType > ansatzSpace_;
  std::map< const std::string, Dune::shared_ptr< const PatternType > > patterns_;
  std::map< const std::string, Dune::shared_ptr< MatrixType > > matrices_;
  std::map< const std::string, Dune::shared_ptr< VectorType > > vectors_;
}; // class DetailedDiscretizations


} // namespace Parametric
#endif // HAVE_DUNE_RB

} // namespace CG
} // namespace Elliptic
} // namespace Linear
} // namespace Stationary
} // namespace Solvers
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_CG_DETAILED_DISCRETIZATIONS_HH
