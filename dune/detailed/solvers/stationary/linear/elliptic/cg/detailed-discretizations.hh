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
#include <dune/stuff/la/container/separable.hh>
//#include <dune/stuff/common/separable-container.hh>
#include <dune/stuff/common/color.hh>

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
  static const int                                                                  polOrder = polynomialOrder;

  typedef Dune::grid::Part::Interface< typename GridPartImp::Traits > GridPartType;
  typedef typename GridPartType::ctype                                DomainFieldType;
  static const int                                                    dimDomain = GridPartType::dimension;
  typedef RangeFieldImp                                               RangeFieldType;
  static const int                                                    dimRange = 1;

  typedef Dune::Stuff::Grid::BoundaryInfo::Interface< typename GridPartType::GridViewType > BoundaryInfoType;
  typedef typename Model::Interface< DomainFieldType, dimDomain, RangeFieldType, dimRange > ModelType;
  typedef typename ModelType::ParamFieldType                                                ParamFieldType;
  static const int maxDimParam = ModelType::maxDimParam;

private:
  typedef Dune::FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, dimRange >               FunctionSpaceType;
  typedef Dune::Detailed::Discretizations::DiscreteFunctionSpace::Continuous::Lagrange< FunctionSpaceType,
                                                                                        GridPartType,
                                                                                        polOrder >  LagrangeSpaceType;
public:
  typedef Dune::Detailed::Discretizations::DiscreteFunctionSpace::Sub::Linear::Dirichlet< LagrangeSpaceType >
                                                                                                    TestSpaceType;
  typedef TestSpaceType                                                                             AnsatzSpaceType;

private:
  typedef typename Dune::Detailed::Discretizations::LA::Container::Factory::Eigen< RangeFieldType > ContainerFactory;
public:
  typedef typename ContainerFactory::RowMajorSparseMatrixType                                       MatrixType;
  typedef Dune::Stuff::LA::Container::Separable< MatrixType, ParamFieldType, maxDimParam >          SeparableMatrixType;
  typedef typename ContainerFactory::DenseVectorType                                                VectorType;
  typedef Dune::Stuff::LA::Container::Separable< VectorType, ParamFieldType, maxDimParam >          SeparableVectorType;
  typedef typename TestSpaceType::PatternType                                                       PatternType;

  typedef Dune::Detailed::Discretizations::DiscreteFunction::Default< TestSpaceType, VectorType >
                                                                                          DiscreteTestFunctionType;
  typedef Dune::Detailed::Discretizations::DiscreteFunction::Default< AnsatzSpaceType, VectorType >
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
    , systemComputed_(false)
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
    if (!model_->parametric()) {
      if (model_->diffusion()->parametric()) {
        msg << "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
            << " nonparametric model has a parametric diffusion!";
        ++throw_up;
      }
      if (model_->force()->parametric()) {
        msg << "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
            << " nonparametric model has a parametric force!";
        ++throw_up;
      }
      if (model_->dirichlet()->parametric()) {
        msg << "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
            << " nonparametric model has parametric dirichlet values!";
        ++throw_up;
      }
      if (model_->neumann()->parametric()) {
        msg << "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
            << " nonparametric model has parametric neumann values!";
        ++throw_up;
      }
    } else {
      if (model_->diffusion()->parametric() && !model_->diffusion()->separable()) {
        msg << "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
            << " only implemented for nonparametric or separable-parametric diffusion!";
        ++throw_up;
      }
      if (model_->force()->parametric() && !model_->force()->separable()) {
        msg << "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
            << " only implemented for nonparametric or separable-parametric force!";
        ++throw_up;
      }
      if (model_->dirichlet()->parametric() && !model_->dirichlet()->separable()) {
        msg << "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
            << " only implemented for nonparametric or separable-parametric dirichlet values!";
        ++throw_up;
      }
      if (model_->neumann()->parametric() && !model_->neumann()->separable()) {
        msg << "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
            << " only implemented for nonparametric or separable-parametric neumann!";
        ++throw_up;
      }
      if (throw_up)
        DUNE_THROW(Dune::InvalidStateException, msg.str());
    }
  } // DetailedDiscretizations

  Dune::shared_ptr< const GridPartType > gridPart() const
  {
    return gridPart_;
  }

  Dune::shared_ptr< const BoundaryInfoType > boundaryInfo() const
  {
    return boundaryInfo_;
  }

  Dune::shared_ptr< const ModelType > model() const
  {
    return model_;
  }

  void init(const std::string prefix = "",
            std::ostream& out = Dune::Stuff::Common::Logger().debug())
  {
    if (!initialized_) {
      Dune::Timer timer;

      out << prefix << "initializing discrete function spaces... " << std::flush;
      lagrangeSpace_ = Dune::make_shared< const LagrangeSpaceType >(*gridPart_);
      testSpace_ = Dune::make_shared< const TestSpaceType >(*lagrangeSpace_, boundaryInfo_);
      ansatzSpace_ = Dune::make_shared< const AnsatzSpaceType >(*lagrangeSpace_, boundaryInfo_);
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      out << prefix << "projecting dirichlet boundary values... " << std::flush;
      Dune::shared_ptr< SeparableVectorType > dirichletVector;
      // if the dirichlet values are not parametric
      if (!model_->diffusion()->parametric()) {
        // project them
        DiscreteAnsatzFunctionType dirichlet(*lagrangeSpace_);
        Dune::Stuff::DiscreteFunction::Projection::Dirichlet::project(*boundaryInfo_,
                                                                      *(model_->dirichlet()),
                                                                      dirichlet);
        dirichletVector = Dune::make_shared< SeparableVectorType >(dirichlet.vector());
      } else {
        // we can assume they are separable (see constructor)
        // so we project each component
        std::vector< Dune::shared_ptr< VectorType > > dirichletComponents;
        for (std::size_t qq = 0; qq < model_->diffusion()->numComponents(); ++qq) {
          DiscreteAnsatzFunctionType dirichlet(*lagrangeSpace_); // <- this is supposed to be here and not before for!
          Dune::Stuff::DiscreteFunction::Projection::Dirichlet::project(*boundaryInfo_,
                                                                        *(model_->dirichlet()->components()[qq]),
                                                                        dirichlet);
          dirichletComponents.push_back(dirichlet.vector());
        }
        dirichletVector = Dune::make_shared< SeparableVectorType >(model_->dirichlet()->paramSize(),
                                                                   dirichletComponents,
                                                                   model_->dirichlet()->coefficients());

      } // if the dirichlet values are not parametric
      vectors_.insert(std::make_pair("dirichlet", dirichletVector));
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      out << prefix << "initializing operators and functionals:" << std::endl;
      // * left hand side
      //   * diffusion operators
      typedef Dune::Detailed::Discretizations::Evaluation::Local::Binary::Elliptic< FunctionSpaceType,
                                                                                    typename ModelType::FunctionType >
        EllipticEvaluationType;
      std::vector< EllipticEvaluationType* > diffusionEvaluations;
      typedef Dune::Detailed::Discretizations::DiscreteOperator::Local::Codim0::Integral< EllipticEvaluationType >
        EllipticOperatorType;
      std::vector< EllipticOperatorType* > diffusionOperators;
      timer.reset();
      if (!model_->diffusion()->parametric()) {
        out << prefix << "  1 diffusion operator...   " << std::flush;
        diffusionEvaluations.push_back(new EllipticEvaluationType(model_->diffusion(), model_->diffusion()->order()));
        diffusionOperators.push_back(new EllipticOperatorType(*(diffusionEvaluations[0])));
      } else {
        assert(false && "Implement me for parametric models!");
      }
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;
      // * right hand side
      //   * L2 force functional
      typedef Dune::Detailed::Discretizations::Evaluation::Local::Unary::Scale< FunctionSpaceType,
                                                                                typename ModelType::FunctionType >
        ProductEvaluationType;
      std::vector< ProductEvaluationType* > forceEvaluations;
      typedef Dune::Detailed::Discretizations::DiscreteFunctional::Local::Codim0::Integral< ProductEvaluationType >
        L2VolumeFunctionalType;
      std::vector< L2VolumeFunctionalType* > forceFunctionals;
      timer.reset();
      if (!model_->force()->parametric())
      {
        out << prefix << "  1 force     functional... " << std::flush;
        forceEvaluations.push_back(new ProductEvaluationType(model_->force(), model_->force()->order()));
        forceFunctionals.push_back(new L2VolumeFunctionalType(*(forceEvaluations[0])));
      } else {
        assert(false && "Implement me for parametric models!");
      }
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;
      timer.reset();
      //   * L2 neumann functional
      std::vector< ProductEvaluationType* > neumannEvaluations;
      std::vector< L2VolumeFunctionalType* > neumannFunctionals;
      if (!model_->neumann()->parametric())
      {
        out << prefix << "  1 neumann   functional... " << std::flush;
        neumannEvaluations.push_back(new ProductEvaluationType(model_->neumann(), model_->neumann()->order()));
        neumannFunctionals.push_back(new L2VolumeFunctionalType(*(neumannEvaluations[0])));
      } else {
        assert(false && "Implement me for parametric models!");
      }
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      out << prefix << "initializing matrices and vectors:" << std::endl;
      timer.reset();
      // create the left hand side matrices for the diffusion
      // * therefore create the pattern
      Dune::shared_ptr< const PatternType > diffusionPattern = testSpace_->computePattern(*ansatzSpace_);
      patterns_.insert(std::make_pair("diffusion", diffusionPattern));
      // * and the matrices
      Dune::shared_ptr< SeparableMatrixType > diffusionMatrices;
      if (!model_->diffusion()->parametric()) {
        out << prefix << "  1 diffusion matrix (of size "
            << testSpace_->map().size() << "x" << ansatzSpace_->map().size() << ")... " << std::flush;
        diffusionMatrices
            = Dune::make_shared < SeparableMatrixType >(ContainerFactory::createRowMajorSparseMatrix(*testSpace_,
                                                                                                     *ansatzSpace_,
                                                                                                     *diffusionPattern));
      } else {
        assert(false && "Implement me for parametric models!");
      }
      matrices_.insert(std::make_pair("diffusion", diffusionMatrices));
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;
      // create the right hand side vectors
      // * for the force
      Dune::shared_ptr< SeparableVectorType > forceVectors;
      if (!model_->force()->parametric()) {
        out << prefix << "  1 force     vector (of size " << testSpace_->map().size() << ")... "
            << " " << Dune::Stuff::Common::whitespaceify(ansatzSpace_->map().size()) << std::flush;
        forceVectors
            = Dune::make_shared < SeparableVectorType >(ContainerFactory::createDenseVector(*testSpace_));
      } else {
        assert(false && "Implement me for parametric models!");
      }
      vectors_.insert(std::make_pair("force", forceVectors));
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;
      // * for the neumann values
      Dune::shared_ptr< SeparableVectorType > neumannVectors;
      if (!model_->neumann()->parametric()) {
        out << prefix << "  1 neumann   vector (of size " << testSpace_->map().size() << ")... "
            << " " << Dune::Stuff::Common::whitespaceify(ansatzSpace_->map().size()) << std::flush;
        neumannVectors
            = Dune::make_shared < SeparableVectorType >(ContainerFactory::createDenseVector(*testSpace_));
      } else {
        assert(false && "Implement me for parametric models!");
      }
      vectors_.insert(std::make_pair("neumann", neumannVectors));
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      out << prefix << "initialiting assemblers:" << std::endl;
      timer.reset();
      // * local matrix assembler
      typedef Dune::Detailed::Discretizations::Assembler::Local::Codim0::Matrix< EllipticOperatorType >
          LocalDiffusionMatrixAssemblerType;
      std::vector< Dune::shared_ptr< const LocalDiffusionMatrixAssemblerType > > localDiffusionMatrixAssemblers;
      if (!model_->diffusion()->parametric()) {
        out << prefix << "  1 diffusion matrix assembler... " << std::flush;
        localDiffusionMatrixAssemblers.push_back(
              Dune::make_shared< LocalDiffusionMatrixAssemblerType >(*(diffusionOperators[0])));
      } else {
        assert(false && "Implement me for parametric models!");
      }
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;
      timer.reset();
      // * local vector assemblers
      //   * force vector
      typedef Dune::Detailed::Discretizations::Assembler::Local::Codim0::Vector< L2VolumeFunctionalType >
          LocalFunctionalVectorAssemblerType;
      std::vector< Dune::shared_ptr< const LocalFunctionalVectorAssemblerType > > localForceVectorAssemblers;
      if (!model_->force()->parametric()) {
        out << prefix << "  1 force     vector assembler... " << std::flush;
        localForceVectorAssemblers.push_back(
              Dune::make_shared< LocalFunctionalVectorAssemblerType >(*(forceFunctionals[0])));
      } else {
        assert(false && "Implement me for parametric models!");
      }
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;
      timer.reset();
      //   * neumann vector
      std::vector< Dune::shared_ptr< const LocalFunctionalVectorAssemblerType > > localNeumannVectorAssemblers;
      if (!model_->neumann()->parametric()) {
        out << prefix << "  1 neumann   vector assembler... " << std::flush;
        localNeumannVectorAssemblers.push_back(
              Dune::make_shared< LocalFunctionalVectorAssemblerType >(*(neumannFunctionals[0])));
      } else {
        assert(false && "Implement me for parametric models!");
      }
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;
      // * system assembler
      out << prefix << "assembling system... " << std::flush;
      typedef Dune::Detailed::Discretizations::Assembler::System< TestSpaceType, AnsatzSpaceType > SystemAssemblerType;
      SystemAssemblerType systemAssembler(*testSpace_, *ansatzSpace_);
      for (std::size_t qq = 0; qq < localDiffusionMatrixAssemblers.size(); ++qq)
        systemAssembler.addLocalMatrixAssembler(localDiffusionMatrixAssemblers[qq],
                                                diffusionMatrices->components()[qq]);
      for (std::size_t qq = 0; qq < localForceVectorAssemblers.size(); ++qq)
        systemAssembler.addLocalVectorAssembler(localForceVectorAssemblers[qq],
                                                forceVectors->components()[qq]);
      for (std::size_t qq = 0; qq < localNeumannVectorAssemblers.size(); ++qq)
        systemAssembler.addLocalVectorAssembler(localNeumannVectorAssemblers[qq],
                                                neumannVectors->components()[qq]);
      systemAssembler.assemble();
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      // clean up
      for (std::size_t ii = 0; ii < neumannFunctionals.size(); ++ii)
        delete neumannFunctionals[ii];
      for (std::size_t ii = 0; ii < neumannEvaluations.size(); ++ii)
        delete neumannEvaluations[ii];
      for (std::size_t ii = 0; ii < forceFunctionals.size(); ++ii)
        delete forceFunctionals[ii];
      for (std::size_t ii = 0; ii < forceEvaluations.size(); ++ii)
        delete forceEvaluations[ii];
      for (std::size_t ii = 0; ii < diffusionOperators.size(); ++ii)
        delete diffusionOperators[ii];
      for (std::size_t ii = 0; ii < diffusionEvaluations.size(); ++ii)
        delete diffusionEvaluations[ii];

      // done
      initialized_ = true;
    } // if !(initialized_)
  } // void init()

  Dune::shared_ptr< const AnsatzSpaceType > ansatzSpace() const
  {
    assert(initialized_);
    return ansatzSpace_;
  }

  Dune::shared_ptr< const TestSpaceType > testSpace() const
  {
    assert(initialized_);
    return testSpace_;
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

  /**
   *  \todo In order to only compute the system matrix and the right hand side once, ceratin members need to be
   *        mutable. Maybe, someone knows a better solution? Since neither the right hand side nor the system matrix
   *        is going to change, the overall scenario of more than one solve might be stupid anyway...
   */
  void solve(Dune::shared_ptr< VectorType > solutionVector,
             const std::string linearSolverType = "eigen.iterative.bicgstab.diagonal",
             const unsigned int linearSolverMaxIter = 5000,
             const double linearSolverPrecision = 1e-12,
             const std::string prefix = "",
             std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
  {
    if (!initialized_)
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " call init() before calling solve()!");
    // check, that we are really in the nonparametric setting!
    assert(matrices_.find("diffusion") != matrices_.end());
    assert(vectors_.find("force") != vectors_.end());
    assert(vectors_.find("neumann") != vectors_.end());
    assert(vectors_.find("dirichlet") != vectors_.end());
    const SeparableMatrixType& diffusionMatrices = *(matrices_.find("diffusion")->second);
    const SeparableVectorType& forceVectors = *(vectors_.find("force")->second);
    const SeparableVectorType& neumannVectors = *(vectors_.find("neumann")->second);
    const SeparableVectorType& dirichletVectors = *(vectors_.find("dirichlet")->second);
    if (model_->parametric())
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " nonparametric solve() called for a parametric model!");
    // since we are in a nonparametric setting, we need to compute the system matrix and the right hand side only once
    Dune::Timer timer;
    if (!systemComputed_) {
      out << prefix << "computing system matrix...   " << std::flush;
      systemMatrix_ = Dune::shared_ptr< MatrixType >(new MatrixType(*(diffusionMatrices.components()[0])));
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;
      out << prefix << "computing right hand side... " << std::flush;
      timer.reset();
      rightHandSide_ = ContainerFactory::createDenseVector(*testSpace_);
      rightHandSide_->backend() = forceVectors.components()[0]->backend()
          + neumannVectors.components()[0]->backend()
          - systemMatrix_->backend() * dirichletVectors.components()[0]->backend();
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;
      out << prefix << "applying constraints...      " << std::flush;
      timer.reset();
      typedef Dune::Detailed::Discretizations::Assembler::System< TestSpaceType, AnsatzSpaceType > SystemAssemblerType;
      const SystemAssemblerType systemAssembler(*testSpace_, *ansatzSpace_);
      systemAssembler.applyConstraints(*systemMatrix_, *rightHandSide_);
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;
      systemComputed_ = true;
    } // since we are in a nonparametric setting, we need to compute the system matrix and the right hand side only once
    out << prefix << "solving linear system (of size " << systemMatrix_->rows()
        << "x" << systemMatrix_->cols() << ")" << std::endl;
    out << prefix << "  using '" << linearSolverType << "'... " << std::flush;
    timer.reset();
    typedef typename Dune::Stuff::LA::Solver::Interface< MatrixType, VectorType > SolverType;
    Dune::shared_ptr< SolverType > solver = Dune::Stuff::LA::Solver::create< MatrixType, VectorType >(linearSolverType);
    const unsigned int failure = solver->apply(*systemMatrix_,
                                               *rightHandSide_,
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
    if (solutionVector->size() != ansatzSpace_->map().size())
      DUNE_THROW(Dune::MathError,
                 "\n"
                 << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " linear solver '" << linearSolverType << "' produced a solution of wrong size (is "
                 << solutionVector->size() << ", should be " << ansatzSpace_->map().size() << ")!");
    solutionVector->backend() += dirichletVectors.components()[0]->backend();
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

//  void visualizeFunction(const Dune::shared_ptr< const DiscreteTestFunctionType > discreteFunction,
//                         const std::string filename = id() + ".discreteTestFunction",
//                         const std::string prefix = "",
//                         std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
//  {
//    // preparations
//    assert(initialized_ && "Please call init() before calling visualize()");
//    Dune::Timer timer;
//    out << prefix << "writing '" << discreteFunction->name() << "' to '" << filename;
//    if (dimDomain == 1)
//      out << ".vtp";
//    else
//      out << ".vtu";
//    out << "'... " << std::flush;
//    typedef Dune::VTKWriter< typename DiscreteAnsatzFunctionType::DiscreteFunctionSpaceType::GridViewType > VTKWriterType;
//    VTKWriterType vtkWriter(discreteFunction->space().gridView());
//    vtkWriter.addVertexData(discreteFunction);
//    vtkWriter.write(filename);
//    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
//  } // void visualizeFunction(...)

private:
  const Dune::shared_ptr< const GridPartType > gridPart_;
  const Dune::shared_ptr< const BoundaryInfoType > boundaryInfo_;
  const Dune::shared_ptr< const ModelType > model_;
  bool initialized_;
  mutable bool systemComputed_;
  Dune::shared_ptr< const LagrangeSpaceType > lagrangeSpace_;
  Dune::shared_ptr< const TestSpaceType > testSpace_;
  Dune::shared_ptr< const AnsatzSpaceType > ansatzSpace_;
  std::map< const std::string, Dune::shared_ptr< const PatternType > > patterns_;
  std::map< const std::string, Dune::shared_ptr< SeparableMatrixType > > matrices_;
  std::map< const std::string, Dune::shared_ptr< SeparableVectorType > > vectors_;
  mutable Dune::shared_ptr< MatrixType > systemMatrix_;
  mutable Dune::shared_ptr< VectorType > rightHandSide_;
}; // class DetailedDiscretizations


#define ENABLED 0
#if ENABLED
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
    const unsigned int failure = solver->apply(*systemMatrix,
                                       *rhsVector,
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
    if (solutionVector->size() != ansatzSpace_->map().size())
      DUNE_THROW(Dune::MathError,
                 "\n"
                 << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " linear solver '" << linearSolverType << "' produced a solution of wrong size (is "
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
#endif // ENABLED

} // namespace CG
} // namespace Elliptic
} // namespace Linear
} // namespace Stationary
} // namespace Solvers
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_CG_DETAILED_DISCRETIZATIONS_HH
