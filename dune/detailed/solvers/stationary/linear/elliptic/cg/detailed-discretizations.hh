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


template< class GridPartImp,
          int polynomialOrder,
          class RangeFieldImp, int rangeDim >
class DetailedDiscretizations;


template< class GridPartImp, int polynomialOrder, class RangeFieldImp >
class DetailedDiscretizations< GridPartImp, polynomialOrder, RangeFieldImp, 1 >
{
public:
  typedef DetailedDiscretizations<  GridPartImp,
                                    polynomialOrder,
                                    RangeFieldImp, 1 >                ThisType;
  typedef Dune::grid::Part::Interface< typename GridPartImp::Traits > GridPartType;
  static const int                                                    polOrder = polynomialOrder;

  typedef typename GridPartType::ctype                                DomainFieldType;
  static const int                                                    dimDomain = GridPartType::dimension;
  typedef RangeFieldImp                                               RangeFieldType;
  static const int                                                    dimRange = 1;

  typedef Dune::Stuff::Grid::BoundaryInfo::Interface< typename GridPartType::GridViewType > BoundaryInfoType;
  typedef typename Model::Interface< DomainFieldType, dimDomain, RangeFieldType, dimRange > ModelType;

  typedef typename ModelType::ParamType ParamType;

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
  typedef typename ContainerFactory::RowMajorSparseMatrixType MatrixType;
  typedef Dune::Stuff::LA::Container::Separable< MatrixType > SeparableMatrixType;
  typedef typename ContainerFactory::DenseVectorType          VectorType;
  typedef Dune::Stuff::LA::Container::Separable< VectorType > SeparableVectorType;
  typedef typename TestSpaceType::PatternType                 PatternType;

  typedef Dune::Detailed::Discretizations::DiscreteFunction::Default< TestSpaceType, VectorType >
                                                                                          DiscreteTestFunctionType;
  typedef Dune::Detailed::Discretizations::DiscreteFunction::DefaultConst< TestSpaceType, VectorType >
                                                                                          DiscreteTestFunctionConstType;
  typedef Dune::Detailed::Discretizations::DiscreteFunction::Default< AnsatzSpaceType, VectorType >
                                                                                          DiscreteAnsatzFunctionType;
  typedef Dune::Detailed::Discretizations::DiscreteFunction::DefaultConst< AnsatzSpaceType, VectorType >
                                                                                          DiscreteAnsatzFunctionConstType;

  typedef Dune::Stuff::Common::ExtendedParameterTree DescriptionType;

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
    if (model_->parametric() && !model_->separable()) {
      msg << "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
          << " only implemented for nonparametric or separable-parametric models!";
      ++throw_up;
    }
    if (throw_up)
      DUNE_THROW(Dune::InvalidStateException, msg.str());
  } // DetailedDiscretizations

  static DescriptionType createSampleDescription(const std::string /*subName*/ = "")
  {
    return DescriptionType();
  } // ... createSampleDescription(...)

  static ThisType* createFromDescription(const Dune::shared_ptr< const GridPartType > /*_gridPart*/,
                                         const Dune::shared_ptr< const ModelType > /*_model*/,
                                         const Dune::shared_ptr< const BoundaryInfoType > /*_boundaryInfo*/,
                                         const DescriptionType& _description,
                                         const std::string _subName = id())
  {
    // get correct description
    DescriptionType description;
    if (_description.hasSub(_subName))
      description = _description.sub(_subName);
    else
      description = _description;
    assert(false);
  } // ... createFromParamTree(...)

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

  bool parametric() const
  {
    return model_->parametric();
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
      if (!model_->dirichlet()->parametric()) {
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
        for (size_t qq = 0; qq < model_->dirichlet()->numComponents(); ++qq) {
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
        out << prefix << "  1 diffusion operator...    " << std::flush;
        diffusionEvaluations.push_back(new EllipticEvaluationType(model_->diffusion(), model_->diffusion()->order()));
        diffusionOperators.push_back(new EllipticOperatorType(*(diffusionEvaluations[0])));
      } else {
        // we are separable (see constructor), loop over all components
        out << prefix << "  " << model_->diffusion()->numComponents() << " diffusion operators...   " << std::flush;
        for (size_t qq = 0; qq < model_->diffusion()->numComponents(); ++qq) {
          diffusionEvaluations.push_back(new EllipticEvaluationType(model_->diffusion()->components()[qq],
                                                                    model_->diffusion()->order()));
          diffusionOperators.push_back(new EllipticOperatorType(*(diffusionEvaluations[qq])));
        } // loop over all components
      } // if (!model_->diffusion()->parametric())
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
      if (!model_->force()->parametric()) {
        out << prefix << "  1 force     functional...  " << std::flush;
        forceEvaluations.push_back(new ProductEvaluationType(model_->force(), model_->force()->order()));
        forceFunctionals.push_back(new L2VolumeFunctionalType(*(forceEvaluations[0])));
      } else {
        // we are separable (see constructor), loop over all components
        out << prefix << "  " << model_->force()->numComponents() << " force     functionals... " << std::flush;
        for (size_t qq = 0; qq < model_->force()->numComponents(); ++qq) {
          forceEvaluations.push_back(new ProductEvaluationType(model_->force()->components()[qq],
                                                               model_->force()->order()));
          forceFunctionals.push_back(new L2VolumeFunctionalType(*(forceEvaluations[qq])));
        } // loop over all components
      } // if (!model_->force()->parametric())
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;
      timer.reset();
      //   * L2 neumann functional
      std::vector< ProductEvaluationType* > neumannEvaluations;
      std::vector< L2VolumeFunctionalType* > neumannFunctionals;
      if (!model_->neumann()->parametric())
      {
        out << prefix << "  1 neumann   functional...  " << std::flush;
        neumannEvaluations.push_back(new ProductEvaluationType(model_->neumann(), model_->neumann()->order()));
        neumannFunctionals.push_back(new L2VolumeFunctionalType(*(neumannEvaluations[0])));
      } else {
        // we are separable (see constructor), loop over all components
        out << prefix << "  " << model_->neumann()->numComponents() << " neumann   functionals... " << std::flush;
        for (size_t qq = 0; qq < model_->neumann()->numComponents(); ++qq) {
          neumannEvaluations.push_back(new ProductEvaluationType(model_->neumann()->components()[qq],
                                                                 model_->neumann()->order()));
          neumannFunctionals.push_back(new L2VolumeFunctionalType(*(neumannEvaluations[qq])));
        } // loop over all components
      }
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      out << prefix << "initializing matrices and vectors:" << std::endl;
      timer.reset();
      // create the left hand side matrices for the diffusion
      // * therefore create the pattern
      Dune::shared_ptr< const PatternType > diffusionPattern = testSpace_->computePattern(*ansatzSpace_);
      patterns_.insert(std::make_pair("diffusion", diffusionPattern));
      // * and the matrices
      Dune::shared_ptr< SeparableMatrixType > diffusionMatrix;
      if (!model_->diffusion()->parametric()) {
        out << prefix << "  1 diffusion matrix   (of size "
            << testSpace_->map().size() << "x" << ansatzSpace_->map().size() << ")... " << std::flush;
        diffusionMatrix
            = Dune::make_shared< SeparableMatrixType >(ContainerFactory::createRowMajorSparseMatrix(*testSpace_,
                                                                                                    *ansatzSpace_,
                                                                                                    *diffusionPattern));
      } else {
        // create one matrix for each component
        out << prefix << "  " << model_->diffusion()->numComponents() <<   " diffusion matrices (of size "
            << testSpace_->map().size() << "x" << ansatzSpace_->map().size() << ")... " << std::flush;
        std::vector< Dune::shared_ptr< MatrixType > > diffusionMatrices;
        for (size_t qq = 0; qq < model_->diffusion()->numComponents(); ++qq)
          diffusionMatrices.push_back(ContainerFactory::createRowMajorSparseMatrix(*testSpace_,
                                                                                   *ansatzSpace_,
                                                                                   *diffusionPattern));
        diffusionMatrix = Dune::make_shared< SeparableMatrixType >(model_->diffusion()->paramSize(),
                                                                   diffusionMatrices,
                                                                   model_->diffusion()->coefficients());
      } // if (!model_->diffusion()->parametric())
      matrices_.insert(std::make_pair("diffusion", diffusionMatrix));
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;
      // create the right hand side vectors
      // * for the force
      Dune::shared_ptr< SeparableVectorType > forceVector;
      if (!model_->force()->parametric()) {
        out << prefix << "  1 force     vector   (of size " << testSpace_->map().size() << ")... "
            << " " << Dune::Stuff::Common::whitespaceify(ansatzSpace_->map().size()) << std::flush;
        forceVector
            = Dune::make_shared< SeparableVectorType >(ContainerFactory::createDenseVector(*testSpace_));
      } else {
        out << prefix << "  " << model_->force()->numComponents()
            << " force     vectors  (of size " << testSpace_->map().size() << ")... "
            << " " << Dune::Stuff::Common::whitespaceify(ansatzSpace_->map().size()) << std::flush;
        std::vector< Dune::shared_ptr< VectorType > > forceVectors;
        for (size_t qq = 0; qq < model_->force()->numComponents(); ++qq)
          forceVectors.push_back(ContainerFactory::createDenseVector(*testSpace_));
        forceVector = Dune::make_shared< SeparableVectorType >(model_->force()->paramSize(),
                                                               forceVectors,
                                                               model_->force()->coefficients());
      }
      vectors_.insert(std::make_pair("force", forceVector));
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;
      // * for the neumann values
      Dune::shared_ptr< SeparableVectorType > neumannVector;
      if (!model_->neumann()->parametric()) {
        out << prefix << "  1 neumann   vector   (of size " << testSpace_->map().size() << ")... "
            << " " << Dune::Stuff::Common::whitespaceify(ansatzSpace_->map().size()) << std::flush;
        neumannVector
            = Dune::make_shared< SeparableVectorType >(ContainerFactory::createDenseVector(*testSpace_));
      } else {
        out << prefix << "  " << model_->neumann()->numComponents()
            << " neumann   vectors  (of size " << testSpace_->map().size() << ")... "
            << " " << Dune::Stuff::Common::whitespaceify(ansatzSpace_->map().size()) << std::flush;
        std::vector< Dune::shared_ptr< VectorType > > neumannVectors;
        for (size_t qq = 0; qq < model_->neumann()->numComponents(); ++qq)
          neumannVectors.push_back(ContainerFactory::createDenseVector(*testSpace_));
        neumannVector = Dune::make_shared< SeparableVectorType >(model_->neumann()->paramSize(),
                                                                 neumannVectors,
                                                                 model_->neumann()->coefficients());
      }
      vectors_.insert(std::make_pair("neumann", neumannVector));
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      out << prefix << "initialiting assemblers:" << std::endl;
      timer.reset();
      // * local matrix assembler
      typedef Dune::Detailed::Discretizations::Assembler::Local::Codim0::Matrix< EllipticOperatorType >
          LocalDiffusionMatrixAssemblerType;
      std::vector< Dune::shared_ptr< const LocalDiffusionMatrixAssemblerType > > localDiffusionMatrixAssembler;
      out << prefix << "  " << diffusionOperators.size() << " diffusion matrix assembler... " << std::flush;
      for (size_t qq = 0; qq < diffusionOperators.size(); ++qq)
        localDiffusionMatrixAssembler.push_back(
              Dune::make_shared< LocalDiffusionMatrixAssemblerType >(*(diffusionOperators[qq])));
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;
      timer.reset();
      // * local vector assemblers
      //   * force vector
      typedef Dune::Detailed::Discretizations::Assembler::Local::Codim0::Vector< L2VolumeFunctionalType >
          LocalFunctionalVectorAssemblerType;
      std::vector< Dune::shared_ptr< const LocalFunctionalVectorAssemblerType > > localForceVectorAssembler;
      out << prefix << "  " << forceFunctionals.size() << " force     vector assembler... " << std::flush;
      for (size_t qq = 0; qq < forceFunctionals.size(); ++qq)
        localForceVectorAssembler.push_back(
              Dune::make_shared< LocalFunctionalVectorAssemblerType >(*(forceFunctionals[qq])));
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;
      timer.reset();
      //   * neumann vector
      std::vector< Dune::shared_ptr< const LocalFunctionalVectorAssemblerType > > localNeumannVectorAssembler;
      out << prefix << "  " << neumannFunctionals.size() << " neumann   vector assembler... " << std::flush;
      for (size_t qq = 0; qq < neumannFunctionals.size(); ++qq)
        localNeumannVectorAssembler.push_back(
              Dune::make_shared< LocalFunctionalVectorAssemblerType >(*(neumannFunctionals[qq])));
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;
      // * system assembler

      out << prefix << "assembling system... " << std::flush;
      typedef Dune::Detailed::Discretizations::Assembler::System< TestSpaceType, AnsatzSpaceType > SystemAssemblerType;
      SystemAssemblerType systemAssembler(*testSpace_, *ansatzSpace_);
      for (size_t qq = 0; qq < localDiffusionMatrixAssembler.size(); ++qq)
        systemAssembler.addLocalMatrixAssembler(localDiffusionMatrixAssembler[qq],
                                                diffusionMatrix->components()[qq]);
      for (size_t qq = 0; qq < localForceVectorAssembler.size(); ++qq)
        systemAssembler.addLocalVectorAssembler(localForceVectorAssembler[qq],
                                                forceVector->components()[qq]);
      for (size_t qq = 0; qq < localNeumannVectorAssembler.size(); ++qq)
        systemAssembler.addLocalVectorAssembler(localNeumannVectorAssembler[qq],
                                                neumannVector->components()[qq]);
      systemAssembler.assemble();
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      // clean up
      for (size_t ii = 0; ii < neumannFunctionals.size(); ++ii)
        delete neumannFunctionals[ii];
      for (size_t ii = 0; ii < neumannEvaluations.size(); ++ii)
        delete neumannEvaluations[ii];
      for (size_t ii = 0; ii < forceFunctionals.size(); ++ii)
        delete forceFunctionals[ii];
      for (size_t ii = 0; ii < forceEvaluations.size(); ++ii)
        delete forceEvaluations[ii];
      for (size_t ii = 0; ii < diffusionOperators.size(); ++ii)
        delete diffusionOperators[ii];
      for (size_t ii = 0; ii < diffusionEvaluations.size(); ++ii)
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

  Dune::shared_ptr< const PatternType > systemPattern(const std::string type = "diffusion") const
  {
    if (!initialized_)
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " please call init() before calling pattern()!");
    assert(patterns_.find(type) != patterns_.end());
    return patterns_.find(type)->second;
  } // Dune::shared_ptr< const PatternType > pattern() const

  Dune::shared_ptr< const SeparableMatrixType > systemMatrix(const std::string type) const
  {
    if (!initialized_)
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " please call init() before calling systemMatrix()!");
    assert(matrices_.find(type) != matrices_.end());
    return matrices_.find(type)->second;
  } // Dune::shared_ptr< const PatternType > pattern() const

  Dune::shared_ptr< SeparableMatrixType > matrix(const std::string type)
  {
    if (!initialized_)
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " please call init() before calling systemMatrix()!");
    assert(matrices_.find(type) != matrices_.end());
    return matrices_.find(type)->second;
  } // Dune::shared_ptr< const PatternType > pattern() const

  Dune::shared_ptr< MatrixType > systemMatrix(const ParamType mu = ParamType()) const
  {
    if (!initialized_)
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " please call init() before calling systemMatrix()!");
    // first of all, get the corect parameters (the model returns empty ones for nonparametric functions)
    const ParamType muDiffusion = model_->mapParam(mu, "diffusion");
    assert(matrices_.find("diffusion") != matrices_.end());
    const SeparableMatrixType& diffusionMatrix = *(matrices_.find("diffusion")->second);
    Dune::shared_ptr< MatrixType > systemMatrix = diffusionMatrix.fix(muDiffusion);
    typedef Dune::Detailed::Discretizations::Assembler::System< TestSpaceType, AnsatzSpaceType > SystemAssemblerType;
    const SystemAssemblerType systemAssembler(*testSpace_, *ansatzSpace_);
    systemAssembler.applyMatrixConstraints(*systemMatrix);
    return systemMatrix;
  }

  Dune::shared_ptr< const SeparableVectorType > systemVector(const std::string type) const
  {
    if (!initialized_)
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " please call init() before calling systemVector()!");
    assert(vectors_.find(type) != vectors_.end());
    return vectors_.find(type)->second;
  } // Dune::shared_ptr< const PatternType > pattern() const

  Dune::shared_ptr< SeparableVectorType > vector(const std::string type)
  {
    if (!initialized_)
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " please call init() before calling systemVector()!");
    assert(vectors_.find(type) != vectors_.end());
    return vectors_.find(type)->second;
  }

  Dune::shared_ptr< VectorType > rightHandSide(const ParamType mu = ParamType()) const
  {
    if (!initialized_)
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " please call init() before calling rightHandSide()!");
    // first of all, get the corect parameters (the model returns empty ones for nonparametric functions)
    const ParamType muDiffusion = model_->mapParam(mu, "diffusion");
    const ParamType muForce = model_->mapParam(mu, "force");
    const ParamType muDirichlet = model_->mapParam(mu, "dirichlet");
    const ParamType muNeumann = model_->mapParam(mu, "neumann");
    assert(matrices_.find("diffusion") != matrices_.end());
    assert(vectors_.find("force") != vectors_.end());
    assert(vectors_.find("neumann") != vectors_.end());
    assert(vectors_.find("dirichlet") != vectors_.end());
    const SeparableMatrixType& diffusionMatrix = *(matrices_.find("diffusion")->second);
    const SeparableVectorType& forceVector = *(vectors_.find("force")->second);
    const SeparableVectorType& neumannVector = *(vectors_.find("neumann")->second);
    const SeparableVectorType& dirichletVector = *(vectors_.find("dirichlet")->second);
    Dune::shared_ptr< MatrixType > systemMatrix = diffusionMatrix.fix(muDiffusion);
    Dune::shared_ptr< VectorType > rightHandSide(new VectorType());
    rightHandSide->backend() = forceVector.fix(muForce)->backend()
        + neumannVector.fix(muNeumann)->backend()
        - systemMatrix->backend() * dirichletVector.fix(muDirichlet)->backend();
    typedef Dune::Detailed::Discretizations::Assembler::System< TestSpaceType, AnsatzSpaceType > SystemAssemblerType;
    const SystemAssemblerType systemAssembler(*testSpace_, *ansatzSpace_);
    systemAssembler.applyVectorConstraints(*rightHandSide);
    return rightHandSide;
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

  Dune::shared_ptr< DiscreteAnsatzFunctionConstType > createConstAnsatzFunction(const Dune::shared_ptr< const VectorType > vector,
                                                                                const std::string name = "ansatzFunction") const
  {
    assert(initialized_ && "Please call init() before calling createAnsatzFunction()!");
    Dune::shared_ptr< DiscreteAnsatzFunctionConstType > ansatzFunction(new DiscreteAnsatzFunctionConstType(
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

  Dune::shared_ptr< DiscreteTestFunctionConstType > createConstTestFunction(const Dune::shared_ptr< const VectorType > vector,
                                                                            const std::string name = "testFunction") const
  {
    assert(initialized_ && "Please call init() before calling createTestFunction()!");
    Dune::shared_ptr< DiscreteTestFunctionConstType > testFunction(new DiscreteTestFunctionConstType(
                                                                     *testSpace_,
                                                                     vector,
                                                                     name));
    return testFunction;
  } // Dune::shared_ptr< DiscreteTestFunctionType > createTestFunction(...) const

private:
  void generic_solve(Dune::shared_ptr< VectorType >& solutionVector,
                     const ParamType& mu,
                     const std::string& linearSolverType,
                     const size_t& linearSolverMaxIter,
                     const double linearSolverPrecision,
                     const std::string prefix,
                     std::ostream& out) const
  {
    // first of all, get the corect parameters (the model returns empty ones for nonparametric functions)
    const ParamType muDiffusion = model_->mapParam(mu, "diffusion");
    const ParamType muForce = model_->mapParam(mu, "force");
    const ParamType muDirichlet = model_->mapParam(mu, "dirichlet");
    const ParamType muNeumann = model_->mapParam(mu, "neumann");
    assert(matrices_.find("diffusion") != matrices_.end());
    assert(vectors_.find("force") != vectors_.end());
    assert(vectors_.find("neumann") != vectors_.end());
    assert(vectors_.find("dirichlet") != vectors_.end());
    const SeparableMatrixType& diffusionMatrix = *(matrices_.find("diffusion")->second);
    const SeparableVectorType& forceVector = *(vectors_.find("force")->second);
    const SeparableVectorType& neumannVector = *(vectors_.find("neumann")->second);
    const SeparableVectorType& dirichletVector = *(vectors_.find("dirichlet")->second);
    Dune::Timer timer;
    out << prefix << "computing system matrix...   " << std::flush;
    Dune::shared_ptr< MatrixType > systemMatrix = diffusionMatrix.fix(muDiffusion);
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
    out << prefix << "computing right hand side... " << std::flush;
    timer.reset();
    VectorType rightHandSide;
    rightHandSide.backend() = forceVector.fix(muForce)->backend()
        + neumannVector.fix(muNeumann)->backend()
        - systemMatrix->backend() * dirichletVector.fix(muDirichlet)->backend();
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
    out << prefix << "applying constraints...      " << std::flush;
    timer.reset();
    typedef Dune::Detailed::Discretizations::Assembler::System< TestSpaceType, AnsatzSpaceType > SystemAssemblerType;
    const SystemAssemblerType systemAssembler(*testSpace_, *ansatzSpace_);
    systemAssembler.applyConstraints(*systemMatrix, rightHandSide);
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
    out << prefix << "solving linear system (of size " << systemMatrix->rows()
        << "x" << systemMatrix->cols() << ")" << std::endl;
    out << prefix << "  using '" << linearSolverType << "'... " << std::flush;
    timer.reset();
    typedef typename Dune::Stuff::LA::Solver::Interface< MatrixType, VectorType > SolverType;
    const Dune::shared_ptr< const SolverType > solver(Dune::Stuff::LA::Solver::create< MatrixType, VectorType >(linearSolverType));
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
    if (solutionVector->size() != ansatzSpace_->map().size())
      DUNE_THROW(Dune::MathError,
                 "\n"
                 << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " linear solver '" << linearSolverType << "' produced a solution of wrong size (is "
                 << solutionVector->size() << ", should be " << ansatzSpace_->map().size() << ")!");
    solutionVector->backend() += dirichletVector.components()[0]->backend();
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void solve(...)

public:
  void solve(Dune::shared_ptr< VectorType > solutionVector,
             const std::string linearSolverType = "bicgstab.diagonal",
             const size_t linearSolverMaxIter = 5000,
             const double linearSolverPrecision = 1e-12,
             const std::string prefix = "",
             std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
  {
    if (!initialized_)
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " call init() before calling solve()!");
    // check, that we are really in the nonparametric setting!
    if (parametric())
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " nonparametric solve() called for a parametric model!");
    generic_solve(solutionVector,
                  ParamType(),
                  linearSolverType, linearSolverMaxIter, linearSolverPrecision,
                  prefix, out);
  }

  void solve(Dune::shared_ptr< VectorType > solutionVector,
             const ParamType& mu,
             const std::string linearSolverType = "bicgstab.diagonal",
             const size_t linearSolverMaxIter = 5000,
             const double linearSolverPrecision = 1e-12,
             const std::string prefix = "",
             std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
  {
    if (!initialized_)
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " call init() before calling solve()!");
    // check, that we are really in the nonparametric setting!
    if (!parametric())
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " parametric solve() called for a nonparametric model!");
    generic_solve(solutionVector,
                  mu,
                  linearSolverType, linearSolverMaxIter, linearSolverPrecision,
                  prefix, out);
  }

  void solveFullNeumann(Dune::shared_ptr< VectorType >& solutionVector,
                        const std::string& linearSolverType,
                        const size_t& linearSolverMaxIter,
                        const double linearSolverPrecision,
                        const std::string prefix,
                        std::ostream& out) const
  {
    if (!initialized_)
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " call init() before calling solveFullNeumann()!");
    // check, that we are really in the nonparametric setting!
    if (parametric())
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " nonparametric solveFullNeumann() called for a parametric model!");
    // first of all, get the corect parameters (the model returns empty ones for nonparametric functions)
    const ParamType mu;
    const ParamType muDiffusion = model_->mapParam(mu, "diffusion");
    const ParamType muForce = model_->mapParam(mu, "force");
    const ParamType muNeumann = model_->mapParam(mu, "neumann");
    assert(matrices_.find("diffusion") != matrices_.end());
    assert(vectors_.find("force") != vectors_.end());
    assert(vectors_.find("neumann") != vectors_.end());
    assert(vectors_.find("dirichlet") != vectors_.end());
    const SeparableMatrixType& diffusionMatrix = *(matrices_.find("diffusion")->second);
    const SeparableVectorType& forceVector = *(vectors_.find("force")->second);
    const SeparableVectorType& neumannVector = *(vectors_.find("neumann")->second);
    Dune::Timer timer;
    out << prefix << "computing system matrix...   " << std::flush;
    Dune::shared_ptr< MatrixType > systemMatrix = diffusionMatrix.fix(muDiffusion);
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
    out << prefix << "computing right hand side... " << std::flush;
    timer.reset();
    VectorType rightHandSide;
    rightHandSide.backend() = forceVector.fix(muForce)->backend()
        + neumannVector.fix(muNeumann)->backend();
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
    out << prefix << "fixing first dof...      " << std::flush;
    timer.reset();
    assert(patterns_.find("diffusion") != patterns_.end());
    for (size_t jj : patterns_.find("diffusion")->second->set(0))
      systemMatrix->set(0, jj, 0.0);
    systemMatrix->set(0, 0, 1.0);
    rightHandSide.set(0, 0.0);
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
    out << prefix << "solving linear system (of size " << systemMatrix->rows()
        << "x" << systemMatrix->cols() << ")" << std::endl;
    out << prefix << "  using '" << linearSolverType << "'... " << std::flush;
    timer.reset();
    typedef typename Dune::Stuff::LA::Solver::Interface< MatrixType, VectorType > SolverType;
    const Dune::shared_ptr< const SolverType > solver(Dune::Stuff::LA::Solver::create< MatrixType, VectorType >(linearSolverType));
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
    if (solutionVector->size() != ansatzSpace_->map().size())
      DUNE_THROW(Dune::MathError,
                 "\n"
                 << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " linear solver '" << linearSolverType << "' produced a solution of wrong size (is "
                 << solutionVector->size() << ", should be " << ansatzSpace_->map().size() << ")!");
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void solve(...)

  void visualizeAnsatzVector(const Dune::shared_ptr< const VectorType > vector,
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
    const Dune::shared_ptr< const DiscreteAnsatzFunctionConstType > discreteFunction
        = createConstAnsatzFunction(vector, name);
    visualizeFunction(discreteFunction, filename, "", Dune::Stuff::Common::Logger().devnull());
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void visualizeAnsatzVector(...)

  void visualizeTestVector(const Dune::shared_ptr< const VectorType > vector,
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
    const Dune::shared_ptr< const DiscreteTestFunctionConstType > discreteFunction
        = createConstTestFunction(vector, name);
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

  void visualizeFunction(const Dune::shared_ptr< const DiscreteAnsatzFunctionConstType > discreteFunction,
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
    typedef Dune::VTKWriter< typename DiscreteAnsatzFunctionConstType::DiscreteFunctionSpaceType::GridViewType > VTKWriterType;
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
  mutable bool systemComputed_;
  Dune::shared_ptr< const LagrangeSpaceType > lagrangeSpace_;
  Dune::shared_ptr< const TestSpaceType > testSpace_;
  Dune::shared_ptr< const AnsatzSpaceType > ansatzSpace_;
  std::map< const std::string, Dune::shared_ptr< const PatternType > > patterns_;
  std::map< const std::string, Dune::shared_ptr< SeparableMatrixType > > matrices_;
  std::map< const std::string, Dune::shared_ptr< SeparableVectorType > > vectors_;
}; // class DetailedDiscretizations


} // namespace CG
} // namespace Elliptic
} // namespace Linear
} // namespace Stationary
} // namespace Solvers
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_CG_DETAILED_DISCRETIZATIONS_HH
