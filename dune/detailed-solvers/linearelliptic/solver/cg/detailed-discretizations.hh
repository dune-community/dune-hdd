#ifndef DUNE_DETAILED_SOLVERS_LINEARELLIPTIC_SOLVER_CG_DETAILED_DISCRETIZATIONS_HH
#define DUNE_DETAILED_SOLVERS_LINEARELLIPTIC_SOLVER_CG_DETAILED_DISCRETIZATIONS_HH

#include <memory>
#include <sstream>

#include <dune/common/exceptions.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/timer.hh>

#include <dune/grid/part/interface.hh>

#include <dune/fem/misc/mpimanager.hh>

#include <dune/stuff/common/logging.hh>
#include <dune/stuff/discretefunction/projection/dirichlet.hh>
#include <dune/stuff/la/solver.hh>
#include <dune/stuff/la/container/affineparametric.hh>
#include <dune/stuff/common/color.hh>
#include <dune/stuff/common/parameter/tree.hh>

#include <dune/detailed/discretizations/space/continuouslagrange/fem.hh>
#include <dune/detailed/discretizations/la/containerfactory/eigen.hh>
#include <dune/detailed/discretizations/localevaluation/elliptic.hh>
#include <dune/detailed/discretizations/localevaluation/product.hh>
#include <dune/detailed/discretizations/localoperator/codim0.hh>
#include <dune/detailed/discretizations/localfunctional/codim0.hh>
#include <dune/detailed/discretizations/localfunctional/codim1.hh>
#include <dune/detailed/discretizations/assembler/local/codim0.hh>
#include <dune/detailed/discretizations/assembler/local/codim1.hh>
#include <dune/detailed/discretizations/space/constraints.hh>
#include <dune/detailed/discretizations/assembler/system.hh>
#include <dune/detailed/discretizations/discretefunction/default.hh>

#include "../../model/interface.hh"
#include "../interface.hh"

namespace Sane = Dune::Detailed::Discretizations;

namespace Dune {
namespace DetailedSolvers {
namespace LinearElliptic {


//// forward of the multiscale solver, to allow for some friendlyness
//template< class GridImp, class RangeFieldImp, int rangeDim, int polynomialOrder >
//class MultiscaleSolverSemiContinuousGalerkinDD;


// forward of the solver, to be used in the traits and allow for specialization
template< class GridPartImp, class RangeFieldImp, int rangeDim, int polynomialOrder, bool scalarDiffusion = true >
class SolverContinuousGalerkinDD
{
public:
  SolverContinuousGalerkinDD() = delete;
};


/**
 *  \brief  Traits for SolverContinuousGalerkinDD
 */
template< class GridPartImp, class RangeFieldImp, int rangeDim, int polynomialOrder, bool scalarDiffusion = true >
class SolverContinuousGalerkinDDTraits
{
public:
  typedef SolverContinuousGalerkinDD< GridPartImp, RangeFieldImp, rangeDim, polynomialOrder, scalarDiffusion > derived_type;
  typedef typename GridPartImp::Traits                  GridPartTraits;
  typedef Dune::grid::Part::Interface< GridPartTraits > GridPartType;
  typedef typename GridPartType::ctype  DomainFieldType;
  static const unsigned int             dimDomain = GridPartType::dimension;
  static const unsigned int             polOrder = polynomialOrder;
  typedef RangeFieldImp                 RangeFieldType;
  static const unsigned int             dimRange = rangeDim;
public:
  typedef Dune::Stuff::GridboundaryInterface< typename GridPartType::GridViewType >                 BoundaryInfoType;
  typedef ModelInterface< DomainFieldType, dimDomain, RangeFieldType, dimRange, scalarDiffusion >   ModelType;
  typedef typename Sane::ContainerFactoryEigen< RangeFieldImp >                                     ContainerFactory;
  typedef typename ContainerFactory::RowMajorSparseMatrixType                                       MatrixType;
  typedef typename ContainerFactory::DenseVectorType                                                VectorType;
}; // class ContinuousGalerkinDDTraits


/**
 *  \brief  Solver of linear elliptic pdes using a continuous galerkin discretization provided by dune-detailed-discretizations
 */
template< class GridPartImp, class RangeFieldImp, int polynomialOrder >
class SolverContinuousGalerkinDD< GridPartImp, RangeFieldImp, 1, polynomialOrder, true >
    : public SolverInterface< SolverContinuousGalerkinDDTraits< GridPartImp, RangeFieldImp, 1, polynomialOrder, true > >
//    , public SolverParametricInterface< SolverContinuousGalerkinDDTraits< GridPartImp, RangeFieldImp, 1, polynomialOrder > >
{
  typedef SolverInterface< SolverContinuousGalerkinDDTraits< GridPartImp, RangeFieldImp, 1, polynomialOrder, true > > BaseType;
//  typedef SolverParametricInterface< Traits >                                                 ParametricBaseType;
public:
  typedef SolverContinuousGalerkinDDTraits< GridPartImp, RangeFieldImp, 1, polynomialOrder, true >  Traits;

  typedef typename Traits::GridPartType GridPartType;
  static const int polOrder = Traits::polOrder;

  typedef typename Traits::DomainFieldType  DomainFieldType;
  static const unsigned int                 dimDomain = GridPartType::dimension;
  typedef typename Traits::RangeFieldType   RangeFieldType;
  static const unsigned int                 dimRange = Traits::dimRange;

  typedef typename Traits::ModelType        ModelType;
  typedef typename Traits::BoundaryInfoType BoundaryInfoType;

  typedef typename ModelType::ParamType ParamType;

private:
  typedef typename Traits::ContainerFactory ContainerFactory;
public:
  typedef typename Traits::MatrixType       MatrixType;
  typedef typename Traits::VectorType       VectorType;

private:
  typedef Dune::Stuff::LA::AffineParametricContainer< MatrixType > AffineParametricMatrixType;
  typedef Dune::Stuff::LA::AffineParametricContainer< VectorType > AffineParametricVectorType;

public:
  typedef Sane::ContinuousLagrangeSpace::FemWrapper< GridPartType, polOrder, RangeFieldType, dimRange > TestSpaceType;
  typedef TestSpaceType                         AnsatzSpaceType;
  typedef typename TestSpaceType::PatternType   PatternType;

private:
  typedef Sane::DiscreteFunctionDefault< AnsatzSpaceType, VectorType >      DiscreteFunctionType;
  typedef Sane::DiscreteFunctionDefaultConst< AnsatzSpaceType, VectorType > ConstDiscreteFunctionType;

public:
  typedef Dune::Stuff::Common::ExtendedParameterTree SettingsType;

  static const std::string id()
  {
    return BaseType::id() + ".cg.dd";
  }

  SolverContinuousGalerkinDD(const std::shared_ptr< const GridPartType > _gridPart,
                             const std::shared_ptr< const BoundaryInfoType > _boundaryInfo,
                             const std::shared_ptr< const ModelType > _model)
    : gridPart_(_gridPart)
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
    if (model_->dirichlet()->parametric()) {
      msg << "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
          << " not implemented for parametric dirichlet values!";
      ++throw_up;
    }
    if (throw_up)
      DUNE_THROW(Dune::InvalidStateException, msg.str());
    // function spaces
    space_ = std::make_shared< const TestSpaceType >(*gridPart_);
  } // SolverContinuousGalerkinDD

  std::shared_ptr< const GridPartType > gridPart() const
  {
    return gridPart_;
  }

  std::shared_ptr< const BoundaryInfoType > boundaryInfo() const
  {
    return boundaryInfo_;
  }

  std::shared_ptr< const ModelType > model() const
  {
    return model_;
  }

  std::shared_ptr< VectorType > createVector() const
  {
    return std::shared_ptr< VectorType >(ContainerFactory::createDenseVector(*space_));
  }

  void visualize(const std::shared_ptr< const VectorType > vector,
                 const std::string filename = id() + ".vector",
                 const std::string name = id() + ".vector",
                 std::ostream& out = Dune::Stuff::Common::Logger().devnull(),
                 const std::string prefix = "") const
  {
    // preparations
    assert(vector->size() == space_->mapper().size() && "Given vector has wrong size!");
    Dune::Timer timer;
    out << prefix << "writing '" << name << "'" << std::endl;
    out << prefix << "     to '" << filename;
    if (dimDomain == 1)
      out << ".vtp";
    else
      out << ".vtu";
    out << "'... " << std::flush;
    const std::shared_ptr< const ConstDiscreteFunctionType > discreteFunction = createConstAnsatzFunction(vector, name);
    visualizeFunction(discreteFunction, filename, Dune::Stuff::Common::Logger().devnull());
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // ... visualize(...)

  /**
   *  \defgroup multiscale ´´These methods are needed by the multiscale solver.``
   *  @{
   */
  std::shared_ptr< const AnsatzSpaceType > ansatzSpace() const
  {
    return space_;
  }

  std::shared_ptr< const TestSpaceType > testSpace() const
  {
    return space_;
  }
  /** @} */

  void init(std::ostream& out = Dune::Stuff::Common::Logger().devnull(),
            const std::string prefix = "")
  {
    if (!initialized_) {
      Dune::Timer timer;

      out << prefix << "projecting dirichlet boundary values... " << std::flush;
//      // if the dirichlet values are not parametric
//      if (!model_->dirichlet()->parametric()) {
        // project them
        dirichletVector_ = std::shared_ptr< VectorType >(createVector());
        DiscreteFunctionType dirichlet(*space_, dirichletVector_);
        assert(!model_->dirichlet()->parametric());
        Dune::Stuff::DiscreteFunction::project(*boundaryInfo_, *(model_->dirichlet()), dirichlet);
//      } else {
//        // we can assume they are separable (see constructor), so we project each component
//        std::vector< std::shared_ptr< VectorType > > dirichletComponents;
//        for (size_t qq = 0; qq < model_->dirichlet()->components().size(); ++qq) {
//          DiscreteFunctionType dirichletComponent(*space_); // <- this is supposed to be here and not before for!
//          Dune::Stuff::DiscreteFunction::project(*boundaryInfo_,
//                                                 *(model_->dirichlet()->components()[qq]),
//                                                 dirichletComponent);
//          dirichletComponents.emplace_back(dirichletComponent.vector());
//        }
//        if (!model_->dirichlet()->hasAffinePart()) {
//          dirichletVector_ = std::make_shared< AffineParametricVectorType >(model_->dirichlet()->paramSize(),
//                                                                            dirichletComponents,
//                                                                            model_->dirichlet()->coefficients());
//        } else {
//          DiscreteFunctionType dirichletAffineShift(*space_);
//          Dune::Stuff::DiscreteFunction::project(*boundaryInfo_,
//                                                 *(model_->dirichlet()->affinePart()),
//                                                 dirichletAffineShift);
//          dirichletVector_ = std::make_shared< AffineParametricVectorType >(model_->dirichlet()->paramSize(),
//                                                                            dirichletComponents,
//                                                                            model_->dirichlet()->coefficients(),
//                                                                            dirichletAffineShift.vector());
//        }
//      } // if the dirichlet values are not parametric
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      out << prefix << "assembing system... " << std::flush;
      timer.reset();
      typedef typename ModelType::DiffusionType DiffusionType;
      typedef typename ModelType::ForceType     ForceType;
      typedef typename ModelType::NeumannType   NeumannType;
      // prepare matrix, pattern, vector and assembler
      std::shared_ptr< AffineParametricMatrixType > diffusionMatrix;
      pattern_ = std::shared_ptr< const PatternType >(space_->computePattern());
      // we need the affine part in any case (for the dirichlet row identification)
      typedef Sane::SystemAssembler< TestSpaceType > SystemAssemblerType;
      SystemAssemblerType systemAssembler(*space_);
      // * elliptic diffusion operator
      typedef Sane::LocalOperator::Codim0Integral< Sane::LocalEvaluation::Elliptic< DiffusionType > > EllipticOperatorType;
      typedef Sane::LocalAssembler::Codim0Matrix< EllipticOperatorType > LocalMatrixAssemblerType;
      std::vector< EllipticOperatorType* > diffusionOperators;
      std::vector< LocalMatrixAssemblerType* > diffusionMatrixAssemblers;
      if (!model_->diffusion()->parametric()) {
        diffusionOperators.emplace_back(new EllipticOperatorType(*(model_->diffusion())));
        diffusionMatrixAssemblers.emplace_back(new LocalMatrixAssemblerType(*(diffusionOperators[0])));
        diffusionMatrix = std::make_shared< AffineParametricMatrixType >(
                            std::make_shared< MatrixType >(space_->mapper().size(),
                                                           space_->mapper().size(),
                                                           *pattern_));
        systemAssembler.addLocalAssembler(*(diffusionMatrixAssemblers[0]), *(diffusionMatrix->affinePart()));
      } else {
        // we are affine parametrix (see constructor)
        std::vector< std::shared_ptr< MatrixType > > diffusionMatrixComponents;
        size_t qq = 0;
        for (; qq < model_->diffusion()->components().size(); ++qq) {
          diffusionOperators.emplace_back(new EllipticOperatorType(*(model_->diffusion()->components()[qq])));
          diffusionMatrixAssemblers.emplace_back(new LocalMatrixAssemblerType(*(diffusionOperators[qq])));
          diffusionMatrixComponents.emplace_back(std::make_shared< MatrixType >(space_->mapper().size(),
                                                                                space_->mapper().size(),
                                                                                *pattern_));
          systemAssembler.addLocalAssembler(*(diffusionMatrixAssemblers[qq]), *(diffusionMatrixComponents[qq]));
        }
        if (model_->diffusion()->hasAffinePart()) {
          ++qq;
          diffusionOperators.emplace_back(new EllipticOperatorType(*(model_->diffusion()->affinePart())));
          diffusionMatrixAssemblers.emplace_back(new LocalMatrixAssemblerType(*(diffusionOperators[qq])));
          diffusionMatrix = std::make_shared< AffineParametricMatrixType >(model_->diffusion()->paramSize(),
                                                                           diffusionMatrixComponents,
                                                                           model_->diffusion()->coefficients(),
                                                                           std::make_shared< MatrixType >(space_->mapper().size(),
                                                                                                          space_->mapper().size(),
                                                                                                          *pattern_));
          systemAssembler.addLocalAssembler(*(diffusionMatrixAssemblers[qq]), *(diffusionMatrix->affinePart()));
        } else {
          diffusionMatrix = std::make_shared< AffineParametricMatrixType >(model_->diffusion()->paramSize(),
                                                                           diffusionMatrixComponents,
                                                                           model_->diffusion()->coefficients());
        }
      } // if (!model_->diffusion()->parametric())
      //   * L2 force functional
      typedef Sane::LocalFunctional::Codim0Integral< Sane::LocalEvaluation::Product< ForceType > > L2VolumeFunctionalType;
      typedef Sane::LocalAssembler::Codim0Vector< L2VolumeFunctionalType > LocalVolumeVectorAssemblerType;
      std::shared_ptr< AffineParametricVectorType > forceVector;
      std::vector< L2VolumeFunctionalType* > forceFunctionals;
      std::vector< LocalVolumeVectorAssemblerType* > forceVectorAssemblers;
      if (!model_->force()->parametric()) {
        forceFunctionals.emplace_back(new L2VolumeFunctionalType(*(model_->force())));
        forceVectorAssemblers.emplace_back(new LocalVolumeVectorAssemblerType(*(forceFunctionals[0])));
        auto forceVectorAffinePart = std::make_shared< VectorType >(space_->mapper().size());
        systemAssembler.addLocalAssembler(*(forceVectorAssemblers[0]), *forceVectorAffinePart);
        forceVector = std::make_shared< AffineParametricVectorType >(forceVectorAffinePart);
      } else {
        // we are affine parametric (see constructor)
        std::vector< std::shared_ptr< VectorType > > forceVectorComponents;
        size_t qq = 0;
        for (; qq < model_->force()->components().size(); ++qq) {
          forceFunctionals.emplace_back(new L2VolumeFunctionalType(*(model_->force()->components()[qq])));
          forceVectorAssemblers.emplace_back(new LocalVolumeVectorAssemblerType(*(forceFunctionals[qq])));
          forceVectorComponents.emplace_back(std::make_shared< VectorType >(space_->mapper().size()));
          systemAssembler.addLocalAssembler(*(forceVectorAssemblers[qq]), *(forceVectorComponents[qq]));
        }
        if (model_->force()->hasAffinePart()) {
          auto forceVectorAffinePart = std::make_shared< VectorType >(space_->mapper().size());
          ++qq;
          forceFunctionals.emplace_back(new L2VolumeFunctionalType(*(model_->force()->affinePart())));
          forceVectorAssemblers.emplace_back(new LocalVolumeVectorAssemblerType(*(forceFunctionals[qq])));
          systemAssembler.addLocalAssembler(*(forceVectorAssemblers[qq]), *forceVectorAffinePart);
          forceVector = std::make_shared< AffineParametricVectorType >(model_->force()->paramSize(),
                                                                       forceVectorComponents,
                                                                       model_->force()->coefficients(),
                                                                       forceVectorAffinePart);
        } else {
          forceVector = std::make_shared< AffineParametricVectorType >(model_->force()->paramSize(),
                                                                       forceVectorComponents,
                                                                       model_->force()->coefficients());
        }
      } // if (!model_->force()->parametric())
      //   * L2 neumann functional
      typedef Sane::LocalFunctional::Codim1Integral< Sane::LocalEvaluation::Product< NeumannType > > L2FaceFunctionalType;
      typedef Sane::LocalAssembler::Codim1Vector< L2FaceFunctionalType > LocalFaceVectorAssemblerType;
      std::shared_ptr< AffineParametricVectorType > neumannVector;
      std::vector< L2FaceFunctionalType* > neumannFunctionals;
      std::vector< LocalFaceVectorAssemblerType* > neumannVectorAssemblers;
      if (!model_->neumann()->parametric()) {
        neumannFunctionals.emplace_back(new L2FaceFunctionalType(*(model_->neumann())));
        neumannVectorAssemblers.emplace_back(new LocalFaceVectorAssemblerType(*(neumannFunctionals[0])));
        auto neumannVectorAffinePart = std::make_shared< VectorType >(space_->mapper().size());
        systemAssembler.addLocalAssembler(*(neumannVectorAssemblers[0]),
                                          typename SystemAssemblerType::AssembleOnNeumann(*boundaryInfo_),
                                          *neumannVectorAffinePart);
        neumannVector = std::make_shared< AffineParametricVectorType >(neumannVectorAffinePart);
      } else {
        // we are affine parametric (see constructor)
        std::vector< std::shared_ptr< VectorType > > neumannVectorComponents;
        size_t qq = 0;
        for (; qq < model_->neumann()->components().size(); ++qq) {
          neumannFunctionals.emplace_back(new L2FaceFunctionalType(*(model_->neumann()->components()[qq])));
          neumannVectorAssemblers.emplace_back(new LocalFaceVectorAssemblerType(*(neumannFunctionals[qq])));
          neumannVectorComponents.emplace_back(std::make_shared< VectorType >(space_->mapper().size()));
          systemAssembler.addLocalAssembler(*(neumannVectorAssemblers[qq]),
                                            typename SystemAssemblerType::AssembleOnNeumann(*boundaryInfo_),
                                            *(neumannVectorComponents[qq]));
        }
        if (model_->neumann()->hasAffinePart()) {
          auto neumannVectorAffinePart = std::make_shared< VectorType >(space_->mapper().size());
          ++qq;
          neumannFunctionals.emplace_back(new L2FaceFunctionalType(*(model_->neumann()->affinePart())));
          neumannVectorAssemblers.emplace_back(new LocalFaceVectorAssemblerType(*(neumannFunctionals[qq])));
          systemAssembler.addLocalAssembler(*(neumannVectorAssemblers[qq]),
                                            typename SystemAssemblerType::AssembleOnNeumann(*boundaryInfo_),
                                            *neumannVectorAffinePart);
          neumannVector = std::make_shared< AffineParametricVectorType >(model_->neumann()->paramSize(),
                                                                         neumannVectorComponents,
                                                                         model_->neumann()->coefficients(),
                                                                         neumannVectorAffinePart);
        } else {
          neumannVector = std::make_shared< AffineParametricVectorType >(model_->neumann()->paramSize(),
                                                                         neumannVectorComponents,
                                                                         model_->neumann()->coefficients());
        }
      } // if (!model_->force()->parametric())
      systemAssembler.assemble();

      // build the system matrix and prepare constraints
      Sane::Constraints::Dirichlet< typename GridPartType::GridViewType,
                                    RangeFieldType, true > clearAndSetRows(*boundaryInfo_,
                                                                           space_->mapper().maxNumDofs(),
                                                                           space_->mapper().maxNumDofs());
      Sane::Constraints::Dirichlet< typename GridPartType::GridViewType,
                                    RangeFieldType, false > clearRows(*boundaryInfo_,
                                                                      space_->mapper().maxNumDofs(),
                                                                      space_->mapper().maxNumDofs());
      if (!diffusionMatrix->parametric()) {
        systemMatrix_ = diffusionMatrix;
        systemAssembler.addLocalConstraints(clearAndSetRows, *(systemMatrix_->affinePart()));
      } else {
        auto systemMatrixComponents = diffusionMatrix->components();
        if (diffusionMatrix->hasAffinePart()) {
          systemMatrix_ = std::make_shared< AffineParametricMatrixType >(diffusionMatrix->paramSize(),
                                                                         systemMatrixComponents,
                                                                         diffusionMatrix->coefficients(),
                                                                         diffusionMatrix->affinePart());
        } else {
          auto systemMatrixAffinePart = std::shared_ptr< MatrixType >(
                                          ContainerFactory::createRowMajorSparseMatrix(*space_, *space_));
          systemMatrix_ = std::make_shared< AffineParametricMatrixType >(diffusionMatrix->paramSize(),
                                                                         systemMatrixComponents,
                                                                         diffusionMatrix->coefficients(),
                                                                         systemMatrixAffinePart);
        }
        for (size_t qq = 0; qq < systemMatrix_->components().size(); ++qq)
          systemAssembler.addLocalConstraints(clearRows, *(systemMatrix_->components()[qq]));
        systemAssembler.addLocalConstraints(clearAndSetRows, *(systemMatrix_->affinePart()));
      }

      // build the right hand side vectors
//      std::vector< std::shared_ptr< VectorType > > rhsComponents;
      typedef typename ForceType::CoefficientType CoefficientType;
//      std::vector< std::shared_ptr< const CoefficientType > > rhsCoefficients;
//      std::shared_ptr< VectorType > rhsAffinePart;
//      bool rhsHasAffinePart = forceVector->hasAffinePart()
//                              || neumannVector->hasAffinePart()
//                              || (diffusionMatrix->hasAffinePart() /*&& dirichletVector_->hasAffinePart()*/);
//      if (rhsHasAffinePart)
//        rhsAffinePart = std::make_shared< VectorType >(space_->mapper().size());
      // * force
      forceVector_ = forceVector;
//      for (size_t qq = 0; qq < forceVector->components().size(); ++qq) {
//        rhsComponents.push_back(forceVector->components()[qq]);
//        rhsCoefficients.push_back(forceVector->coefficients()[qq]);
//      }
//      if (forceVector->hasAffinePart())
//        rhsAffinePart->backend() += forceVector->affinePart()->backend();
      // * neumann
      neumannVector_ = neumannVector;
//      for (size_t qq = 0; qq < neumannVector->components().size(); ++qq) {
//        rhsComponents.push_back(neumannVector->components()[qq]);
//        rhsCoefficients.push_back(neumannVector->coefficients()[qq]);
//      }
//      if (neumannVector->hasAffinePart())
//        rhsAffinePart->backend() += neumannVector->affinePart()->backend();
//      // * diffusionAffinePart * dirichlet
//      if (diffusionMatrix->hasAffinePart()) {
//        for (size_t qq = 0; qq < dirichletVector_->components().size(); ++qq) {
//          auto component = std::make_shared< VectorType >(space_->mapper().size());
//          component->backend() = diffusionMatrix->affinePart()->backend() * dirichletVector_->components()[qq]->backend();
//          rhsComponents.push_back(component);
//          rhsCoefficients.push_back(std::make_shared< const CoefficientType >(
//                                      "-1.0*(" + dirichletVector_->coefficients()[qq]->expression() + ")"));
//        }
//        if (dirichletVector_->hasAffinePart())
//          rhsAffinePart->backend() += diffusionMatrix->affinePart()->backend() * dirichletVector_->affinePart()->backend();
//      }
//      // * diffusion * dirichletAffinePart
//      if (dirichletVector_->hasAffinePart()) {
//        for (size_t qq = 0; qq < diffusionMatrix->components().size(); ++qq) {
//          auto component = std::make_shared< VectorType >(space_->mapper().size());
//          component->backend() = diffusionMatrix->components()[qq]->backend() * dirichletVector_->affinePart()->backend();
//          rhsComponents.push_back(component);
//          rhsCoefficients.push_back(std::make_shared< const CoefficientType >(
//                                      "-1.0*(" + diffusionMatrix->coefficients()[qq]->expression() + ")"));
//        }
//      }
      // * diffusion * dirichlet
      if (!diffusionMatrix->parametric()) {
        auto tmp = std::make_shared< VectorType >(space_->mapper().size());
        tmp->backend() = RangeFieldType(-1) * (diffusionMatrix->affinePart()->backend() * dirichletVector_->backend());
        diffusionDirichletVector_ = std::make_shared< AffineParametricVectorType >(tmp);
      } else {
        std::vector< std::shared_ptr< VectorType > > components;
        for (size_t pp = 0; pp < diffusionMatrix->components().size(); ++pp) {
//        for (size_t qq = 0; qq < dirichletVector_->components().size(); ++qq) {
          auto component = std::make_shared< VectorType >(space_->mapper().size());
          component->backend() = RangeFieldType(-1) * (diffusionMatrix->components()[pp]->backend() * dirichletVector_/*->components()[qq]*/->backend());
          components.push_back(component);
//        }
        }
        if (!diffusionMatrix->hasAffinePart()) {
          diffusionDirichletVector_ = std::make_shared< AffineParametricVectorType >(diffusionMatrix->paramSize(),
                                                                                     components,
                                                                                     diffusionMatrix->coefficients());
        } else {
          auto affinePart = std::make_shared< VectorType >(space_->mapper().size());
          affinePart->backend() += diffusionMatrix->affinePart()->backend() * dirichletVector_->backend();
          diffusionDirichletVector_ = std::make_shared< AffineParametricVectorType >(diffusionMatrix->paramSize(),
                                                                                     components,
                                                                                     diffusionMatrix->coefficients(),
                                                                                     affinePart);
        }
      }
//      // actually create the rhs vector
//      if (rhsComponents.size() == 0) {
//        assert(rhsCoefficients.size() == 0);
//        assert(rhsHasAffinePart);
//        rhsVector_ = std::make_shared< AffineParametricVectorType >(rhsAffinePart);
//      } else {
//        assert(rhsComponents.size() == rhsCoefficients.size());
//        if (!rhsHasAffinePart)
//          rhsVector_ = std::make_shared< AffineParametricVectorType >(forceVector->paramSize()
//                                                                      + neumannVector->paramSize()
//                                                                      + diffusionMatrix->paramSize(),
//                                                                      rhsComponents,
//                                                                      rhsCoefficients);
//        else
//          rhsVector_ = std::make_shared< AffineParametricVectorType >(forceVector->paramSize()
//                                                                      + neumannVector->paramSize()
//                                                                      + diffusionMatrix->paramSize(),
//                                                                      rhsComponents,
//                                                                      rhsCoefficients,
//                                                                      rhsAffinePart);
//      } // actually create the rhs vector
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      out << prefix << "applying constraints... " << std::flush;
      // adding the rhs constraints here, the matrix constraints were added above
      if (forceVector_->hasAffinePart())
        systemAssembler.addLocalConstraints(clearRows, *(forceVector_->affinePart()));
      if (forceVector_->parametric()) {
        for (size_t qq = 0; qq < forceVector_->components().size(); ++qq)
          systemAssembler.addLocalConstraints(clearRows, *(forceVector_->components()[qq]));
      }
      if (neumannVector_->hasAffinePart())
        systemAssembler.addLocalConstraints(clearRows, *(neumannVector_->affinePart()));
      if (neumannVector_->parametric()) {
        for (size_t qq = 0; qq < neumannVector_->components().size(); ++qq)
          systemAssembler.addLocalConstraints(clearRows, *(neumannVector_->components()[qq]));
      }
      if (diffusionDirichletVector_->hasAffinePart())
        systemAssembler.addLocalConstraints(clearRows, *(diffusionDirichletVector_->affinePart()));
      if (diffusionDirichletVector_->parametric()) {
        for (size_t qq = 0; qq < diffusionDirichletVector_->components().size(); ++qq)
          systemAssembler.addLocalConstraints(clearRows, *(diffusionDirichletVector_->components()[qq]));
      }
      systemAssembler.applyConstraints();
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      // clean up
      for (auto& element : neumannVectorAssemblers)
        delete element;
      for (auto& element : neumannFunctionals)
        delete element;
      for (auto& element : forceVectorAssemblers)
        delete element;
      for (auto& element : forceFunctionals)
        delete element;
      for (auto& element : diffusionMatrixAssemblers)
        delete element;
      for (auto& element : diffusionOperators)
        delete element;

      // done
      initialized_ = true;
    } // if !(initialized_)
  } // void init(...)

  bool initialized() const
  {
    return initialized_;
  }

private:
  void generic_solve(std::shared_ptr< VectorType >& solutionVector,
                     const ParamType& mu,
                     const SettingsType& linearSolverSettings,
                     const std::string prefix,
                     std::ostream& out) const
  {
    // first of all, get the corect parameters (the model returns empty ones for nonparametric functions)
    const ParamType muDiffusion = model_->mapParam(mu, "diffusion");
    const ParamType muForce = model_->mapParam(mu, "force");
//    const ParamType muDirichlet = model_->mapParam(mu, "dirichlet");
    const ParamType muNeumann = model_->mapParam(mu, "neumann");
    Dune::Timer timer;
    out << prefix << "computing system matrix...   " << std::flush;
    std::shared_ptr< const MatrixType > systemMatrix;
    if (systemMatrix_->parametric())
      systemMatrix = systemMatrix_->fix(muDiffusion);
    else
      systemMatrix = systemMatrix_->affinePart();
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
    out << prefix << "computing right hand side... " << std::flush;
    timer.reset();
    std::shared_ptr< VectorType > rhsVector(createVector());
    if (forceVector_->parametric()) {
      const auto forceVector = forceVector_->fix(muForce);
      rhsVector->backend() = forceVector->backend();
    } else {
      rhsVector->backend() = forceVector_->affinePart()->backend();
    }
    if (neumannVector_->parametric()) {
      const auto neumannVector = neumannVector_->fix(muNeumann);
      rhsVector->backend() += neumannVector->backend();
    } else {
      rhsVector->backend() += neumannVector_->affinePart()->backend();
    }
    if (diffusionDirichletVector_->parametric()) {
      const auto diffusionDirichletVector = diffusionDirichletVector_->fix(muDiffusion);
      rhsVector->backend() += diffusionDirichletVector->backend();
    } else {
      rhsVector->backend() += diffusionDirichletVector_->affinePart()->backend();
    }
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;

    out << prefix << "solving linear system (of size " << systemMatrix->rows()
        << "x" << systemMatrix->cols() << ")" << std::endl;
    const std::string linearSolverType = linearSolverSettings.get< std::string >("type");
    out << prefix << "  using '" << linearSolverType << "'... " << std::flush;
    timer.reset();
    typedef typename Dune::Stuff::LA::SolverInterface< MatrixType, VectorType > SolverType;
    const std::shared_ptr< const SolverType > solver(Dune::Stuff::LA::createSolver< MatrixType, VectorType >(linearSolverType));
    const size_t failure = solver->apply(*systemMatrix,
                                         *rhsVector,
                                         *solutionVector,
                                         linearSolverSettings);
    if (failure)
      DUNE_THROW(Dune::MathError,
                 "\n"
                 << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " linear solver '" << linearSolverType << "' reported error code " << failure << "!\n"
                 << "  1: did not converge\n"
                 << "  2: had numerical issues\n"
                 << "  3: dude, I have no idea");
    if (solutionVector->size() != space_->mapper().size())
      DUNE_THROW(Dune::MathError,
                 "\n"
                 << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " linear solver '" << linearSolverType << "' produced a solution of wrong size (is "
                 << solutionVector->size() << ", should be " << space_->mapper().size() << ")!");
//    if (dirichletVector_->parametric()) {
//      const auto dirichletVector = dirichletVector_->fix(muDirichlet);
//      solutionVector->backend() += dirichletVector->backend();
//    } else
      solutionVector->backend() += dirichletVector_/*->affinePart()*/->backend();
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void generic_solve(...)

public:
  void solve(std::shared_ptr< VectorType > solutionVector,
             const SettingsType& linearSolverSettings,
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
                  linearSolverSettings,
                  prefix, out);
  } // ... solve(...)

  void solve(std::shared_ptr< VectorType > solutionVector,
             const ParamType& mu,
             const SettingsType& linearSolverSettings,
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
                  linearSolverSettings,
                  prefix, out);
  } // ... solve(..., mu, ...)

//  std::shared_ptr< const PatternType > pattern(const std::string type = "diffusion") const
//  {
//    if (!initialized_)
//      DUNE_THROW(Dune::InvalidStateException,
//                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
//                 << " please call init() before calling pattern()!");
//    const auto result = patterns_.find(type);
//    if (result == patterns_.end())
//      DUNE_THROW(Dune::InvalidStateException,
//                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
//                 << " wrong type given (is '" << type << "', has to be 'diffusion')!");
//    return result->second;
//  } // ... pattern(...) const

//  std::shared_ptr< const AffineParametricMatrixType > matrix(const std::string type = "diffusion") const
//  {
//    if (!initialized_)
//      DUNE_THROW(Dune::InvalidStateException,
//                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
//                 << " please call init() before calling matrix()!");
//    const auto result = matrices_.find(type);
//    if (result == matrices_.end())
//      DUNE_THROW(Dune::InvalidStateException,
//                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
//                 << " wrong type given (is '" << type << "', has to be 'diffusion')!");
//    return result->second;
//  } // ... matrix(...) const

//  std::shared_ptr< const AffineParametricVectorType > vector(const std::string type) const
//  {
//    if (!initialized_)
//      DUNE_THROW(Dune::InvalidStateException,
//                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
//                 << " please call init() before calling matrix()!");
//    const auto result = vectors_.find(type);
//    if (result == vectors_.end())
//      DUNE_THROW(Dune::InvalidStateException,
//                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
//                 << " wrong type given (is '" << type << "', has to be one of 'force', 'dirichlet', 'neumann')!");
//    return result->second;
//  } // ... vector(...) const

private:
  std::shared_ptr< ConstDiscreteFunctionType > createConstAnsatzFunction(const std::shared_ptr< const VectorType > vector,
                                                                         const std::string name) const
  {
    return std::make_shared< ConstDiscreteFunctionType >(*space_, vector, name);
  }

  void visualizeFunction(const std::shared_ptr< const ConstDiscreteFunctionType > discreteFunction,
                         const std::string filename = id() + ".function",
                         std::ostream& out = Dune::Stuff::Common::Logger().devnull(),
                         const std::string prefix = "") const
  {
    // preparations
    Dune::Timer timer;
    out << prefix << "writing '" << discreteFunction->name() << "' to '" << filename;
    if (dimDomain == 1)
      out << ".vtp";
    else
      out << ".vtu";
    out << "'... " << std::flush;
    typedef Dune::VTKWriter< typename GridPartType::GridViewType > VTKWriterType;
    VTKWriterType vtkWriter(gridPart_->gridView());
    vtkWriter.addVertexData(discreteFunction);
    vtkWriter.write(filename);
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void visualizeFunction(...)

//  template< class G, class R, int r, int p >
//  friend class MultiscaleSolverSemiContinuousGalerkinDD;

  const std::shared_ptr< const GridPartType > gridPart_;
  const std::shared_ptr< const BoundaryInfoType > boundaryInfo_;
  const std::shared_ptr< const ModelType > model_;
  bool initialized_;
  std::shared_ptr< const TestSpaceType > space_;
  std::shared_ptr< const PatternType > pattern_;
  std::shared_ptr< AffineParametricMatrixType > systemMatrix_;
  std::shared_ptr< AffineParametricVectorType > forceVector_;
  std::shared_ptr< AffineParametricVectorType > neumannVector_;
  std::shared_ptr< /*AffineParametric*/VectorType > dirichletVector_;
  std::shared_ptr< AffineParametricVectorType > diffusionDirichletVector_;
}; // class SolverContinuousGalerkinDD


} // namespace LinearElliptic
} // namespace DetailedSolver
} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_LINEARELLIPTIC_SOLVER_CG_DETAILED_DISCRETIZATIONS_HH
