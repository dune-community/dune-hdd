#ifndef DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_CG_DETAILED_DISCRETIZATIONS_HH
#define DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_CG_DETAILED_DISCRETIZATIONS_HH

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#include <sstream>

#include <dune/common/exceptions.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/timer.hh>

#include <dune/stuff/common/logging.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/discretefunction/projection/dirichlet.hh>
#include <dune/stuff/la/solver.hh>

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
class DetailedDiscretizations;//!TODO this meant to be empty imp base for specialization??

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

private:
  typedef Dune::FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, dimRange > FunctionSpaceType;

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
      DiscreteFunctionType discreteDirichlet(*lagrangeSpace_, "dirichlet");
      Dune::Stuff::DiscreteFunction::Projection::Dirichlet::project(*boundaryInfo_,
                                                                    *(model_->dirichlet()),
                                                                    discreteDirichlet);
      ansatzSpace_ = Dune::shared_ptr< const AnsatzSpaceType >(new AnsatzSpaceType(*testSpace_, discreteDirichlet.createConst()));
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
      auto diffusionMatrix = Dune::make_shared< MatrixType >(testSpace_->map().size(), ansatzSpace_->map().size(), *diffusionPattern);
      matrices_.insert(std::make_pair( "diffusion", diffusionMatrix));
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
    //!TODO don't use naked refs
    return *ansatzSpace_;
  }

  const TestSpaceType& testSpace() const
  {
    assert(initialized_);
    //!TODO don't use naked refs
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

  void solve(Dune::shared_ptr< VectorType > solutionVector,
             const std::string linearSolverType = "eigen.iterative.bicgstab.diagonal",
             const unsigned int linearSolverMaxIter = 5000,
             const double linearSolverPrecision = 1e-12,
             const std::string prefix = "",
             std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
  {
    assert(initialized_ && "Please call init() before calling solve()!");
    out << prefix << "applying constraints... " << std::flush;
    Dune::Timer timer;
    // compute system matrix and right hand side
    // * therefore get matrix and vectors
    assert(matrices_.find("diffusion") != matrices_.end());
    assert(vectors_.find("force") != vectors_.end());
    assert(vectors_.find("neumann") != vectors_.end());
    const MatrixType& diffusionMatrix = *(matrices_.find("diffusion")->second);
    const VectorType& forceVector = *(vectors_.find("force")->second);
    const VectorType& neumannVector = *(vectors_.find("neumann")->second);
    const VectorType& dirichletVector = *(ansatzSpace_->affineShift()->vector());
    // * compute the system matrix
    MatrixType systemMatrix = diffusionMatrix;
    // * compute right hand side
    VectorType rhsVector(testSpace_->map().size());
    rhsVector.backend() = forceVector.backend()
        + neumannVector.backend()
        - systemMatrix.backend() * dirichletVector.backend();
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
    typedef typename Dune::Stuff::LA::Solver::Interface< MatrixType, VectorType > SolverType;
    Dune::shared_ptr< SolverType > solver = Dune::Stuff::LA::Solver::create< MatrixType, VectorType >(linearSolverType);
    const bool success = solver->apply(systemMatrix,
                                       rhsVector,
                                       *solutionVector,
                                       linearSolverMaxIter,
                                       linearSolverPrecision);
    if (!success)
      DUNE_THROW(Dune::MathError,
                 "\nERROR: linear solver '" << linearSolverType << "' reported a problem!");
    if (solutionVector->size() != ansatzSpace_->map().size())
      DUNE_THROW(Dune::MathError,
                 "\nERROR: linear solver '" << linearSolverType << "' produced a solution of wrong size (is "
                 << solutionVector->size() << ", should be " << ansatzSpace_->map().size() << ")!");
    solutionVector->backend() += dirichletVector.backend();
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void solve(...)

  void visualizeAnsatzVector(Dune::shared_ptr< VectorType > vector,
                             const std::string filename = id() + ".ansatzVector",
                             const std::string name = id() + ".ansatzVector",
                             const std::string prefix = "",
                             std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
  {
    // preparations
    assert(initialized_ && "Please call init() before calling visualize()");
    assert(vector->size() == ansatzSpace_->map().size() && "Given vector has wrong size!");
    Dune::Timer timer;
    out << prefix << "writing '" << name << "'" << std::endl;
    out << prefix << "     to '" << filename;
    if (dimDomain == 1)
      out << ".vtp";
    else
      out << ".vtu";
    out << "'... " << std::flush;
    const Dune::shared_ptr< const DiscreteAnsatzFunctionType > discreteFunction
        = createAnsatzFunction(vector, name);
    visualizeFunction(discreteFunction, filename, "", Dune::Stuff::Common::Logger().devnull());
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void visualizeAnsatzVector(...)

  void visualizeTestVector(Dune::shared_ptr< VectorType > vector,
                           const std::string filename = id() + ".testVector",
                           const std::string name = id() + ".testVector",
                           const std::string prefix = "",
                           std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
  {
    // preparations
    assert(initialized_ && "Please call init() before calling visualize()");
    assert(vector->size() == ansatzSpace_->map().size() && "Given vector has wrong size!");
    Dune::Timer timer;
    out << prefix << "writing '" << name << "'" << std::endl;
    out << prefix << "     to '" << filename;
    if (dimDomain == 1)
      out << ".vtp";
    else
      out << ".vtu";
    out << "'... " << std::flush;
    const Dune::shared_ptr< const DiscreteTestFunctionType > discreteFunction
        = createTestFunction(vector, name);
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

  typedef typename ContainerFactory::RowMajorSparseMatrixType MatrixType;

  typedef typename Dune::Stuff::LA::Container::Separable< MatrixType, ParamFieldType, maxParams > SeparableMatrixType;

  typedef typename Dune::Stuff::LA::Container::Separable< VectorType, ParamFieldType, maxParams > SeparableVectorType;

private:
  typedef typename TestSpaceType::PatternType PatternType;

public:
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
    assert(model_->diffusion()->separable()
           && "Please wrap nonparametric and nonseparable functions using Dune::RB::Function::Parametric::Separable::Wrapper!");
    assert(model_->force()->separable()
           && "Please wrap nonparametric and nonseparable functions using Dune::RB::Function::Parametric::Separable::Wrapper!");
    assert(model_->dirichlet()->separable()
           && "Please wrap nonparametric and nonseparable functions using Dune::RB::Function::Parametric::Separable::Wrapper!");
    assert(model_->neumann()->separable()
           && "Please wrap nonparametric and nonseparable functions using Dune::RB::Function::Parametric::Separable::Wrapper!");
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

  //!TODO shared_ptr namespace
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
      ansatzSpace_ = Dune::shared_ptr< const AnsatzSpaceType >(new AnsatzSpaceType(*testSpace_, discreteDirichlet->createConst()));
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      out << prefix << "initializing operators and functionals:" << std::endl;
      timer.reset();
      // * left hand side
      //   * elliptic operator
      out << prefix << "  " << model_->diffusion()->numComponents() << " diffusion operator";
      if (model_->diffusion()->numComponents() > 1)
        out << "s";
      out << "... " << std::flush;
      typedef Dune::Detailed::Discretizations
          ::Evaluation
          ::Local
          ::Binary
          ::Elliptic< FunctionSpaceType, typename ModelType::DiffusionType::ComponentType >
        DiffusionEvaluationType;
      typedef Dune::Detailed::Discretizations
          ::DiscreteOperator
          ::Local
          ::Codim0
          ::Integral< DiffusionEvaluationType >
        DiffusionOperatorType;
      std::vector< const DiffusionEvaluationType* > diffusionEvaluationPtrs;
      std::vector< Dune::shared_ptr< const DiffusionOperatorType > > diffusionOperatorPtrs;
      for (unsigned int qq = 0; qq < model_->diffusion()->numComponents(); ++qq) {
        diffusionEvaluationPtrs.push_back(new DiffusionEvaluationType(model_->diffusion()->components()[qq],
                                                                     model_->diffusionOrder()));
        diffusionOperatorPtrs.push_back(Dune::shared_ptr< DiffusionOperatorType >(new DiffusionOperatorType(
                                                                                    *(diffusionEvaluationPtrs[qq]))));
      }
      typedef Dune::Stuff
          ::Common
          ::SeparableContainer< DiffusionOperatorType, typename ModelType::DiffusionType::CoefficientType >
        SeparableDiffusionOperatorType;
      const SeparableDiffusionOperatorType separableDiffusionOperator(model_->diffusion()->parametric() ? model_->diffusion()->paramSize() : 0,
                                                                      diffusionOperatorPtrs,
                                                                      model_->diffusion()->coefficients());
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;
      // * right hand side
      //   * L2 force functional
      out << prefix << "  " << model_->force()->numComponents() << " force functional";
      if (model_->force()->numComponents() > 1)
        out << "s";
      out << "... " << std::flush;
      timer.reset();
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
        ForceFunctionalType;
      std::vector< const ForceEvaluationType* > forceEvaluationPtrs;//!TODO leaks
      std::vector< Dune::shared_ptr< const ForceFunctionalType > > forceFunctionalPtrs;
      for (unsigned int qq = 0; qq < model_->force()->numComponents(); ++qq) {
        forceEvaluationPtrs.push_back(new ForceEvaluationType(model_->force()->components()[qq],
                                                              model_->forceOrder()));
        forceFunctionalPtrs.push_back(Dune::shared_ptr< ForceFunctionalType >(new ForceFunctionalType(
                                                                                  //!TODO don't use naked refs
                                                                                *(forceEvaluationPtrs[qq]))));
      }
      typedef Dune::Stuff
          ::Common
          ::SeparableContainer< ForceFunctionalType, typename ModelType::ForceType::CoefficientType >
        SeparableForceFunctionalType;
      const SeparableForceFunctionalType separableForceFunctional(model_->force()->parametric() ? model_->force()->paramSize() : 0,
                                                                  forceFunctionalPtrs,
                                                                  model_->force()->coefficients());
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;
      //   * L2 neumann functional
      out << prefix << "  " << model_->neumann()->numComponents() << " neumann functional";
      if (model_->neumann()->numComponents() > 1)
        out << "s";
      out << "... " << std::flush;
      timer.reset();
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
        NeumannFunctionalType;
      std::vector< const NeumannEvaluationType* > neumannEvaluationPtrs;//!TODO leaks
      std::vector< Dune::shared_ptr< const NeumannFunctionalType > > neumannFunctionalPtrs;
      for (unsigned int qq = 0; qq < model_->neumann()->numComponents(); ++qq) {
        neumannEvaluationPtrs.push_back(new NeumannEvaluationType(model_->neumann()->components()[qq],
                                                                  model_->neumannOrder()));
        neumannFunctionalPtrs.push_back(Dune::shared_ptr< NeumannFunctionalType >(new NeumannFunctionalType(
                                                                                      //!TODO don't use naked refs
                                                                                    *(neumannEvaluationPtrs[qq]))));
      }
      typedef Dune::Stuff
          ::Common
          ::SeparableContainer< NeumannFunctionalType, typename ModelType::NeumannType::CoefficientType >
        SeparableNeumannFunctionalType;
      const SeparableNeumannFunctionalType separableNeumannFunctional(model_->neumann()->parametric() ? model_->neumann()->paramSize() : 0,
                                                                      neumannFunctionalPtrs,
                                                                      model_->neumann()->coefficients());
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      out << prefix << "initializing matrices and vectors:" << std::endl;
      timer.reset();
      // * create left hand side matrices
      //   * therefore create the pattern
      Dune::shared_ptr< const PatternType > diffusionPattern = testSpace_->computePattern(*ansatzSpace_);
      patterns_.insert(std::pair< const std::string, Dune::shared_ptr< const PatternType > >("diffusion",
                                                                                             diffusionPattern));
      //   * and the matrices
      out << prefix << "  " << separableDiffusionOperator.numComponents() << " diffusion matri";
      if (separableDiffusionOperator.numComponents() > 1)
        out << "ces";
      else
        out << "x";
      out << "... " << std::flush;
      std::vector< Dune::shared_ptr< MatrixType > > diffusionMatrixPtrs;
      for (unsigned int qq = 0; qq < separableDiffusionOperator.numComponents(); ++ qq)
        diffusionMatrixPtrs.push_back(Dune::shared_ptr< MatrixType >(new MatrixType(
            testSpace_->map().size(),
            ansatzSpace_->map().size(),
            *diffusionPattern)));
      Dune::shared_ptr< SeparableMatrixType > diffusionMatrix(new SeparableMatrixType(separableDiffusionOperator.paramSize(),
                                                                                      diffusionMatrixPtrs,
                                                                                      separableDiffusionOperator.coefficients()));
      matrices_.insert(std::pair< const std::string, Dune::shared_ptr< SeparableMatrixType > >("diffusion",
                                                                                               diffusionMatrix));
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;
      // * create the right hand side vectors
      //   * for the force
      out << prefix << "  " << separableForceFunctional.numComponents() << " force vector";
      if (separableForceFunctional.numComponents() > 1)
        out << "s";
      out << "... " << std::flush;
      timer.reset();
      std::vector< Dune::shared_ptr< VectorType > > forceVectorPtrs;
      for (unsigned int qq = 0; qq < separableForceFunctional.numComponents(); ++qq)
        forceVectorPtrs.push_back(Dune::shared_ptr< VectorType >(ContainerFactory::createDenseVector(*testSpace_)));
      Dune::shared_ptr< SeparableVectorType > forceVector(new SeparableVectorType(separableForceFunctional.paramSize(),
                                                                                  forceVectorPtrs,
                                                                                  separableForceFunctional.coefficients()));
      vectors_.insert(std::pair< const std::string, Dune::shared_ptr< SeparableVectorType > >("force", forceVector));
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;
      //   * for the neumann
      out << prefix << "  " << separableNeumannFunctional.numComponents() << " neumann vector";
      if (separableNeumannFunctional.numComponents() > 1)
        out << "s";
      out << "... " << std::flush;
      timer.reset();
      std::vector< Dune::shared_ptr< VectorType > > neumannVectorPtrs;
      for (unsigned int qq = 0; qq < separableNeumannFunctional.numComponents(); ++qq)
        neumannVectorPtrs.push_back(Dune::shared_ptr< VectorType >(ContainerFactory::createDenseVector(*testSpace_)));
      Dune::shared_ptr< SeparableVectorType > neumannVector(new SeparableVectorType(separableNeumannFunctional.paramSize(),
                                                                                    neumannVectorPtrs,
                                                                                    separableNeumannFunctional.coefficients()));
      vectors_.insert(std::pair< const std::string, Dune::shared_ptr< SeparableVectorType > >("neumann", neumannVector));
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      out << prefix << "assembing system... " << std::flush;
      timer.reset();
      typedef Dune::Detailed::Discretizations::Assembler::System< TestSpaceType, AnsatzSpaceType > SystemAssemblerType;
      SystemAssemblerType systemAssembler(*testSpace_, *ansatzSpace_);
      // * local matrix assemblers
      typedef Dune::Detailed::Discretizations::Assembler::Local::Codim0::Matrix< DiffusionOperatorType >
          LocalDiffusionMatrixAssemblerType;
      for (unsigned int qq = 0; qq < separableDiffusionOperator.numComponents(); ++qq) {
        const Dune::shared_ptr< const LocalDiffusionMatrixAssemblerType > localDiffusionMatrixAssembler(
              new LocalDiffusionMatrixAssemblerType(*(separableDiffusionOperator.components()[qq])));
        systemAssembler.addLocalMatrixAssembler(localDiffusionMatrixAssembler, diffusionMatrix->components()[qq]);
      }
      // * local vector assemblers
      //   * force vector
      typedef Dune::Detailed::Discretizations::Assembler::Local::Codim0::Vector< ForceFunctionalType >
          LocalForceVectorAssemblerType;
      for (unsigned int qq = 0; qq < separableForceFunctional.numComponents(); ++qq) {
        const Dune::shared_ptr< const LocalForceVectorAssemblerType > localForceVectorAssembler(
              new LocalForceVectorAssemblerType(*(separableForceFunctional.components()[qq])));
        systemAssembler.addLocalVectorAssembler(localForceVectorAssembler, forceVector->components()[qq]);
      }
      //   * neumann vector
      typedef Dune::Detailed::Discretizations::Assembler::Local::Codim1::Vector::Neumann< NeumannFunctionalType,
                                                                                          BoundaryInfoType >
          LocalNeumannVectorAssemblerType;
      for (unsigned int qq = 0; qq < separableNeumannFunctional.numComponents(); ++qq) {
        const Dune::shared_ptr< const LocalNeumannVectorAssemblerType > localNeumannVectorAssembler(
              new LocalNeumannVectorAssemblerType(*(separableNeumannFunctional.components()[qq]), boundaryInfo_));
        systemAssembler.addLocalVectorAssembler(localNeumannVectorAssembler, neumannVector->components()[qq]);
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
             Dune::shared_ptr< VectorType > solutionVector,
             const std::string linearSolverType = "eigen.iterative.bicgstab.diagonal",
             const unsigned int linearSolverMaxIter = 5000,
             const double linearSolverPrecision = 1e-12,
             const std::string prefix = "",
             std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
  {
    assert(initialized_ && "Please call init() before calling solve()!");
    out << prefix << "computing system matrix and right hand side... " << std::flush;
    Dune::Timer timer;
    // compute the system matrix
    Dune::shared_ptr< MatrixType > systemMatrix;
    assert(matrices_.find("diffusion") != matrices_.end());
    const Dune::shared_ptr< const SeparableMatrixType > diffusionMatrix =  matrices_.find("diffusion")->second;
    if (model_->diffusion()->parametric())
      systemMatrix = diffusionMatrix->fix(model_->getDiffusionParam(mu));
    else
      systemMatrix = diffusionMatrix->fix();
    // compute the right hand side
    Dune::shared_ptr< VectorType > rhsVector = ContainerFactory::createDenseVector(*testSpace_);
    // * add up force
    assert(vectors_.find("force") != vectors_.end());
    Dune::shared_ptr< SeparableVectorType > forceVector = vectors_.find("force")->second;
    if (model_->force()->parametric())
      rhsVector->backend() += forceVector->fix(model_->getForceParam(mu))->backend();
    else
      rhsVector->backend() += forceVector->fix()->backend();
    // * add up neumann
    assert(vectors_.find("neumann") != vectors_.end());
    Dune::shared_ptr< SeparableVectorType > neumannVector = vectors_.find("neumann")->second;
    if (model_->neumann()->parametric())
      rhsVector->backend() += neumannVector->fix(model_->getNeumannParam(mu))->backend();
    else
      rhsVector->backend() += neumannVector->fix()->backend();
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
    typedef typename Dune::Stuff::LA::Solver::Interface< MatrixType, VectorType > SolverType;
    Dune::shared_ptr< SolverType > solver = Dune::Stuff::LA::Solver::create< MatrixType, VectorType >(linearSolverType);
    const bool success = solver->apply(*systemMatrix,
                                       *rhsVector,
                                       *solutionVector,
                                       linearSolverMaxIter,
                                       linearSolverPrecision);
    if (!success)
      DUNE_THROW(Dune::MathError,
                 "\nERROR: linear solver '" << linearSolverType << "' reported a problem!");
    if (solutionVector->size() != ansatzSpace_->map().size())
      DUNE_THROW(Dune::MathError,
                 "\nERROR: linear solver '" << linearSolverType << "' produced a solution of wrong size (is "
                 << solutionVector->size() << ", should be " << ansatzSpace_->map().size() << ")!");
    solutionVector->backend() += ansatzSpace_->affineShift()->vector()->backend();
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void solve(...)

  void visualizeAnsatzVector(Dune::shared_ptr< VectorType > vector,
                             const std::string filename = id() + ".ansatzVector",
                             const std::string name = id() + "ansatzVector",
                             const std::string prefix = "",
                             std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
  {
    // preparations
    assert(initialized_ && "Please call init() before calling visualize()");
    assert(vector->size() == ansatzSpace_->map().size() && "Given vector has wrong size!");
    Dune::Timer timer;
    out << prefix << "writing '" << name << "'" << std::endl;
    out << prefix << "     to '" << filename;
    if (dimDomain == 1)
      out << ".vtp";
    else
      out << ".vtu";
    out << "'... " << std::flush;
    const Dune::shared_ptr< const DiscreteAnsatzFunctionType > discreteFunction
        = createAnsatzFunction(vector, name);
    visualizeFunction(discreteFunction, filename, "", Dune::Stuff::Common::Logger().devnull());
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void visualizeAnsatzVector(...)

  void visualizeTestVector(Dune::shared_ptr< VectorType > vector,
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
    const Dune::shared_ptr< const DiscreteTestFunctionType > discreteFunction
        = createTestFunction(vector, name);
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

  Dune::shared_ptr< MatrixType > getSystemMatrix(const ParamType& mu) const
  {
    // compute the system matrix
    Dune::shared_ptr< MatrixType > systemMatrix;
    assert(matrices_.find("diffusion") != matrices_.end());
    const Dune::shared_ptr< const SeparableMatrixType > diffusionMatrix =  matrices_.find("diffusion")->second;
    if (model_->diffusion()->parametric())
      systemMatrix = diffusionMatrix->fix(model_->getDiffusionParam(mu));
    else
      systemMatrix = diffusionMatrix->fix();
    // apply constraints
    typedef Dune::Detailed::Discretizations::Assembler::System< TestSpaceType, AnsatzSpaceType > SystemAssemblerType;
    SystemAssemblerType systemAssembler(*testSpace_, *ansatzSpace_);
    systemAssembler.applyMatrixConstraints(*systemMatrix);
    return systemMatrix;
  }

private:
  const shared_ptr< const ModelType > model_;
  const shared_ptr< const GridPartType > gridPart_;
  const shared_ptr< const BoundaryInfoType > boundaryInfo_;
  bool initialized_;
  Dune::shared_ptr< const LagrangeSpaceType > lagrangeSpace_;
  Dune::shared_ptr< const TestSpaceType > testSpace_;
  Dune::shared_ptr< const AnsatzSpaceType > ansatzSpace_;
  std::map< const std::string, Dune::shared_ptr< const PatternType > > patterns_;
  std::map< const std::string, Dune::shared_ptr< SeparableMatrixType > > matrices_;
  std::map< const std::string, Dune::shared_ptr< SeparableVectorType > > vectors_;
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
