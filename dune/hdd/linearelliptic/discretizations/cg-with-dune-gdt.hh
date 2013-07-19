// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_DISCRETIZATION_CG_GDT_HH
#define DUNE_HDD_LINEARELLIPTIC_DISCRETIZATION_CG_GDT_HH

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

#include <dune/gdt/space/continuouslagrange/fem.hh>
#include <dune/gdt/la/containerfactory/eigen.hh>
#include <dune/gdt/localevaluation/elliptic.hh>
#include <dune/gdt/localevaluation/product.hh>
#include <dune/gdt/localoperator/codim0.hh>
#include <dune/gdt/localfunctional/codim0.hh>
#include <dune/gdt/localfunctional/codim1.hh>
#include <dune/gdt/assembler/local/codim0.hh>
#include <dune/gdt/assembler/local/codim1.hh>
#include <dune/gdt/space/constraints.hh>
#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/discretefunction/default.hh>

#include <dune/pymor/la/container/eigen.hh>
#include <dune/pymor/operators/eigen.hh>
#include <dune/pymor/operators/affine.hh>

#include "../model/interface.hh"
#include "interface.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {


// forward, needed in the Traits
template< class GridPartImp, class RangeFieldImp, int rangeDim, int polynomialOrder, bool scalarDiffusion = true >
class ContinuousGalerkinDiscretizationGDT;


/**
 *  \brief  Traits for SolverContinuousGalerkinDD
 */
template< class GridPartImp, class RangeFieldImp, int rangeDim, int polynomialOrder, bool scalarDiffusion = true >
class ContinuousGalerkinDiscretizationGDTTraits
{
public:
  typedef ContinuousGalerkinDiscretizationGDT<  GridPartImp,
                                                RangeFieldImp,
                                                rangeDim,
                                                polynomialOrder,
                                                scalarDiffusion > derived_type;
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
  typedef typename Dune::GDT::ContainerFactoryEigen< RangeFieldImp >                                     ContainerFactory;
  typedef typename ContainerFactory::RowMajorSparseMatrixType                                       MatrixType;
  typedef typename ContainerFactory::DenseVectorType                                                VectorType;
}; // class ContinuousGalerkinDiscretizationGDTTraits


template< class GridPartImp, class RangeFieldImp, int rangeDim, int polynomialOrder, bool scalarDiffusion >
class ContinuousGalerkinDiscretizationGDT
    : public DiscretizationInterface< ContinuousGalerkinDiscretizationGDTTraits<  GridPartImp,
                                                                                  RangeFieldImp,
                                                                                  rangeDim,
                                                                                  polynomialOrder,
                                                                                  scalarDiffusion > >
//    , public SolverParametricInterface< SolverContinuousGalerkinDDTraits< GridPartImp, RangeFieldImp, 1, polynomialOrder > >
{
  typedef ContinuousGalerkinDiscretizationGDTTraits<  GridPartImp,
                                                      RangeFieldImp,
                                                      rangeDim,
                                                      polynomialOrder,
                                                      scalarDiffusion > Traits;
//  typedef SolverInterface< SolverContinuousGalerkinDDTraits< GridPartImp, RangeFieldImp, 1, polynomialOrder, true > > BaseType;
//  typedef SolverParametricInterface< Traits >                                                 ParametricBaseType;
public:
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
//  typedef typename Traits::MatrixType       MatrixType;
//  typedef typename Traits::VectorType       VectorType;
  typedef Dune::Pymor::Operators::EigenRowMajorSparseMatrix MatrixType;
  typedef Dune::Pymor::LA::EigenDenseVector                 VectorType;

private:
//  typedef Dune::Stuff::LA::AffineParametricContainer< MatrixType > AffineParametricMatrixType;
  typedef Dune::Stuff::LA::AffineParametricContainer< VectorType > AffineParametricVectorType;
  typedef Dune::Pymor::Operators::AffinelyDecomposedEigenRowMajorSparseMatrix AffineParametricMatrixType;

public:
  typedef Dune::GDT::ContinuousLagrangeSpace::FemWrapper< GridPartType, polOrder, RangeFieldType, dimRange > TestSpaceType;
  typedef TestSpaceType                         AnsatzSpaceType;
  typedef typename TestSpaceType::PatternType   PatternType;

private:
  typedef Dune::GDT::DiscreteFunctionDefault< AnsatzSpaceType, VectorType >      DiscreteFunctionType;
  typedef Dune::GDT::DiscreteFunctionDefaultConst< AnsatzSpaceType, VectorType > ConstDiscreteFunctionType;

public:
  typedef Dune::Stuff::Common::ExtendedParameterTree SettingsType;

  static const std::string id()
  {
    return typename DiscretizationInterface< ContinuousGalerkinDiscretizationGDTTraits< GridPartType,
                                                                                        RangeFieldType,
                                                                                        dimRange,
                                                                                        polOrder,
                                                                                        scalarDiffusion >
                                            >::id() + ".cg.gdt";
  }

  ContinuousGalerkinDiscretizationGDT(const std::shared_ptr< const GridPartType > gP,
                                      const std::shared_ptr< const BoundaryInfoType > bI,
                                      const std::shared_ptr< const ModelType > mod)
    : gridPart_(gP)
    , boundaryInfo_(bI)
    , model_(mod)
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
  } // ContinuousGalerkinDiscretizationGDT

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
    return std::make_shared< VectorType >(space_->mapper().size());
  }

//  void visualize(const std::shared_ptr< const VectorType > vector,
//                 const std::string filename = id() + ".vector",
//                 const std::string name = id() + ".vector",
//                 std::ostream& out = Dune::Stuff::Common::Logger().devnull(),
//                 const std::string prefix = "") const
//  {
//    // preparations
//    assert(vector->size() == space_->mapper().size() && "Given vector has wrong size!");
//    Dune::Timer timer;
//    out << prefix << "writing '" << name << "'" << std::endl;
//    out << prefix << "     to '" << filename;
//    if (dimDomain == 1)
//      out << ".vtp";
//    else
//      out << ".vtu";
//    out << "'... " << std::flush;
//    const std::shared_ptr< const ConstDiscreteFunctionType > discreteFunction = createConstAnsatzFunction(vector, name);
//    visualizeFunction(discreteFunction, filename, Dune::Stuff::Common::Logger().devnull());
//    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
//  } // ... visualize(...)

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
      typedef Dune::GDT::SystemAssembler< TestSpaceType > SystemAssemblerType;
      SystemAssemblerType systemAssembler(*space_);
      // * elliptic diffusion operator
      typedef Dune::GDT::LocalOperator::Codim0Integral< Dune::GDT::LocalEvaluation::Elliptic< DiffusionType > > EllipticOperatorType;
      typedef Dune::GDT::LocalAssembler::Codim0Matrix< EllipticOperatorType > LocalMatrixAssemblerType;
      std::vector< EllipticOperatorType* > diffusionOperators;
      std::vector< LocalMatrixAssemblerType* > diffusionMatrixAssemblers;
      if (!model_->diffusion()->parametric()) {
        diffusionOperators.emplace_back(new EllipticOperatorType(*(model_->diffusion())));
        diffusionMatrixAssemblers.emplace_back(new LocalMatrixAssemblerType(*(diffusionOperators[0])));

        diffusionMatrix = std::make_shared< AffineParametricMatrixType >(new MatrixType(space_->mapper().size(),
                                                                                        space_->mapper().size(),
                                                                                        *pattern_));
        systemAssembler.addLocalAssembler(*(diffusionMatrixAssemblers[0]), *(diffusionMatrix->affinePart()));
      } else {
//        // we are affine parametric (see constructor)
//        diffusionMatrix = std::make_shared< AffineParametricMatrixType >(
//                            Dune::Pymor::ParameterType("diffusion", model_->diffusion()->paramSize()));
////        std::vector< std::shared_ptr< MatrixType > > diffusionMatrixComponents;
//        size_t qq = 0;
//        for (; qq < model_->diffusion()->components().size(); ++qq) {
//          diffusionOperators.emplace_back(new EllipticOperatorType(*(model_->diffusion()->components()[qq])));
//          diffusionMatrixAssemblers.emplace_back(new LocalMatrixAssemblerType(*(diffusionOperators[qq])));
//          diffusionMatrixComponents.emplace_back(std::make_shared< MatrixType >(space_->mapper().size(),
//                                                                                space_->mapper().size(),
//                                                                                *pattern_));
//          systemAssembler.addLocalAssembler(*(diffusionMatrixAssemblers[qq]), *(diffusionMatrixComponents[qq]));
//        }
//        if (model_->diffusion()->hasAffinePart()) {
//          ++qq;
//          diffusionOperators.emplace_back(new EllipticOperatorType(*(model_->diffusion()->affinePart())));
//          diffusionMatrixAssemblers.emplace_back(new LocalMatrixAssemblerType(*(diffusionOperators[qq])));
//          diffusionMatrix = std::make_shared< AffineParametricMatrixType >(model_->diffusion()->paramSize(),
//                                                                           diffusionMatrixComponents,
//                                                                           model_->diffusion()->coefficients(),
//                                                                           std::make_shared< MatrixType >(space_->mapper().size(),
//                                                                                                          space_->mapper().size(),
//                                                                                                          *pattern_));
//          systemAssembler.addLocalAssembler(*(diffusionMatrixAssemblers[qq]), *(diffusionMatrix->affinePart()));
//        } else {
//          diffusionMatrix = std::make_shared< AffineParametricMatrixType >(model_->diffusion()->paramSize(),
//                                                                           diffusionMatrixComponents,
//                                                                           model_->diffusion()->coefficients());
//        }
      } // if (!model_->diffusion()->parametric())
//      //   * L2 force functional
//      typedef Dune::GDT::LocalFunctional::Codim0Integral< Dune::GDT::LocalEvaluation::Product< ForceType > > L2VolumeFunctionalType;
//      typedef Dune::GDT::LocalAssembler::Codim0Vector< L2VolumeFunctionalType > LocalVolumeVectorAssemblerType;
//      std::shared_ptr< AffineParametricVectorType > forceVector;
//      std::vector< L2VolumeFunctionalType* > forceFunctionals;
//      std::vector< LocalVolumeVectorAssemblerType* > forceVectorAssemblers;
//      if (!model_->force()->parametric()) {
//        forceFunctionals.emplace_back(new L2VolumeFunctionalType(*(model_->force())));
//        forceVectorAssemblers.emplace_back(new LocalVolumeVectorAssemblerType(*(forceFunctionals[0])));
//        auto forceVectorAffinePart = std::make_shared< VectorType >(space_->mapper().size());
//        systemAssembler.addLocalAssembler(*(forceVectorAssemblers[0]), *forceVectorAffinePart);
//        forceVector = std::make_shared< AffineParametricVectorType >(forceVectorAffinePart);
//      } else {
//        // we are affine parametric (see constructor)
//        std::vector< std::shared_ptr< VectorType > > forceVectorComponents;
//        size_t qq = 0;
//        for (; qq < model_->force()->components().size(); ++qq) {
//          forceFunctionals.emplace_back(new L2VolumeFunctionalType(*(model_->force()->components()[qq])));
//          forceVectorAssemblers.emplace_back(new LocalVolumeVectorAssemblerType(*(forceFunctionals[qq])));
//          forceVectorComponents.emplace_back(std::make_shared< VectorType >(space_->mapper().size()));
//          systemAssembler.addLocalAssembler(*(forceVectorAssemblers[qq]), *(forceVectorComponents[qq]));
//        }
//        if (model_->force()->hasAffinePart()) {
//          auto forceVectorAffinePart = std::make_shared< VectorType >(space_->mapper().size());
//          ++qq;
//          forceFunctionals.emplace_back(new L2VolumeFunctionalType(*(model_->force()->affinePart())));
//          forceVectorAssemblers.emplace_back(new LocalVolumeVectorAssemblerType(*(forceFunctionals[qq])));
//          systemAssembler.addLocalAssembler(*(forceVectorAssemblers[qq]), *forceVectorAffinePart);
//          forceVector = std::make_shared< AffineParametricVectorType >(model_->force()->paramSize(),
//                                                                       forceVectorComponents,
//                                                                       model_->force()->coefficients(),
//                                                                       forceVectorAffinePart);
//        } else {
//          forceVector = std::make_shared< AffineParametricVectorType >(model_->force()->paramSize(),
//                                                                       forceVectorComponents,
//                                                                       model_->force()->coefficients());
//        }
//      } // if (!model_->force()->parametric())
//      //   * L2 neumann functional
//      typedef Dune::GDT::LocalFunctional::Codim1Integral< Dune::GDT::LocalEvaluation::Product< NeumannType > > L2FaceFunctionalType;
//      typedef Dune::GDT::LocalAssembler::Codim1Vector< L2FaceFunctionalType > LocalFaceVectorAssemblerType;
//      std::shared_ptr< AffineParametricVectorType > neumannVector;
//      std::vector< L2FaceFunctionalType* > neumannFunctionals;
//      std::vector< LocalFaceVectorAssemblerType* > neumannVectorAssemblers;
//      if (!model_->neumann()->parametric()) {
//        neumannFunctionals.emplace_back(new L2FaceFunctionalType(*(model_->neumann())));
//        neumannVectorAssemblers.emplace_back(new LocalFaceVectorAssemblerType(*(neumannFunctionals[0])));
//        auto neumannVectorAffinePart = std::make_shared< VectorType >(space_->mapper().size());
//        systemAssembler.addLocalAssembler(*(neumannVectorAssemblers[0]),
//                                          typename SystemAssemblerType::AssembleOnNeumann(*boundaryInfo_),
//                                          *neumannVectorAffinePart);
//        neumannVector = std::make_shared< AffineParametricVectorType >(neumannVectorAffinePart);
//      } else {
//        // we are affine parametric (see constructor)
//        std::vector< std::shared_ptr< VectorType > > neumannVectorComponents;
//        size_t qq = 0;
//        for (; qq < model_->neumann()->components().size(); ++qq) {
//          neumannFunctionals.emplace_back(new L2FaceFunctionalType(*(model_->neumann()->components()[qq])));
//          neumannVectorAssemblers.emplace_back(new LocalFaceVectorAssemblerType(*(neumannFunctionals[qq])));
//          neumannVectorComponents.emplace_back(std::make_shared< VectorType >(space_->mapper().size()));
//          systemAssembler.addLocalAssembler(*(neumannVectorAssemblers[qq]),
//                                            typename SystemAssemblerType::AssembleOnNeumann(*boundaryInfo_),
//                                            *(neumannVectorComponents[qq]));
//        }
//        if (model_->neumann()->hasAffinePart()) {
//          auto neumannVectorAffinePart = std::make_shared< VectorType >(space_->mapper().size());
//          ++qq;
//          neumannFunctionals.emplace_back(new L2FaceFunctionalType(*(model_->neumann()->affinePart())));
//          neumannVectorAssemblers.emplace_back(new LocalFaceVectorAssemblerType(*(neumannFunctionals[qq])));
//          systemAssembler.addLocalAssembler(*(neumannVectorAssemblers[qq]),
//                                            typename SystemAssemblerType::AssembleOnNeumann(*boundaryInfo_),
//                                            *neumannVectorAffinePart);
//          neumannVector = std::make_shared< AffineParametricVectorType >(model_->neumann()->paramSize(),
//                                                                         neumannVectorComponents,
//                                                                         model_->neumann()->coefficients(),
//                                                                         neumannVectorAffinePart);
//        } else {
//          neumannVector = std::make_shared< AffineParametricVectorType >(model_->neumann()->paramSize(),
//                                                                         neumannVectorComponents,
//                                                                         model_->neumann()->coefficients());
//        }
//      } // if (!model_->force()->parametric())
//      systemAssembler.assemble();

//      // build the system matrix and prepare constraints
//      Dune::GDT::Constraints::Dirichlet< typename GridPartType::GridViewType,
//                                    RangeFieldType, true > clearAndSetRows(*boundaryInfo_,
//                                                                           space_->mapper().maxNumDofs(),
//                                                                           space_->mapper().maxNumDofs());
//      Dune::GDT::Constraints::Dirichlet< typename GridPartType::GridViewType,
//                                    RangeFieldType, false > clearRows(*boundaryInfo_,
//                                                                      space_->mapper().maxNumDofs(),
//                                                                      space_->mapper().maxNumDofs());
//      if (!diffusionMatrix->parametric()) {
//        systemMatrix_ = diffusionMatrix;
//        systemAssembler.addLocalConstraints(clearAndSetRows, *(systemMatrix_->affinePart()));
//      } else {
//        auto systemMatrixComponents = diffusionMatrix->components();
//        if (diffusionMatrix->hasAffinePart()) {
//          systemMatrix_ = std::make_shared< AffineParametricMatrixType >(diffusionMatrix->paramSize(),
//                                                                         systemMatrixComponents,
//                                                                         diffusionMatrix->coefficients(),
//                                                                         diffusionMatrix->affinePart());
//        } else {
//          auto systemMatrixAffinePart = std::shared_ptr< MatrixType >(
//                                          ContainerFactory::createRowMajorSparseMatrix(*space_, *space_));
//          systemMatrix_ = std::make_shared< AffineParametricMatrixType >(diffusionMatrix->paramSize(),
//                                                                         systemMatrixComponents,
//                                                                         diffusionMatrix->coefficients(),
//                                                                         systemMatrixAffinePart);
//        }
//        for (size_t qq = 0; qq < systemMatrix_->components().size(); ++qq)
//          systemAssembler.addLocalConstraints(clearRows, *(systemMatrix_->components()[qq]));
//        systemAssembler.addLocalConstraints(clearAndSetRows, *(systemMatrix_->affinePart()));
//      }

//      // build the right hand side vectors
////      std::vector< std::shared_ptr< VectorType > > rhsComponents;
//      typedef typename ForceType::CoefficientType CoefficientType;
////      std::vector< std::shared_ptr< const CoefficientType > > rhsCoefficients;
////      std::shared_ptr< VectorType > rhsAffinePart;
////      bool rhsHasAffinePart = forceVector->hasAffinePart()
////                              || neumannVector->hasAffinePart()
////                              || (diffusionMatrix->hasAffinePart() /*&& dirichletVector_->hasAffinePart()*/);
////      if (rhsHasAffinePart)
////        rhsAffinePart = std::make_shared< VectorType >(space_->mapper().size());
//      // * force
//      forceVector_ = forceVector;
////      for (size_t qq = 0; qq < forceVector->components().size(); ++qq) {
////        rhsComponents.push_back(forceVector->components()[qq]);
////        rhsCoefficients.push_back(forceVector->coefficients()[qq]);
////      }
////      if (forceVector->hasAffinePart())
////        rhsAffinePart->backend() += forceVector->affinePart()->backend();
//      // * neumann
//      neumannVector_ = neumannVector;
////      for (size_t qq = 0; qq < neumannVector->components().size(); ++qq) {
////        rhsComponents.push_back(neumannVector->components()[qq]);
////        rhsCoefficients.push_back(neumannVector->coefficients()[qq]);
////      }
////      if (neumannVector->hasAffinePart())
////        rhsAffinePart->backend() += neumannVector->affinePart()->backend();
////      // * diffusionAffinePart * dirichlet
////      if (diffusionMatrix->hasAffinePart()) {
////        for (size_t qq = 0; qq < dirichletVector_->components().size(); ++qq) {
////          auto component = std::make_shared< VectorType >(space_->mapper().size());
////          component->backend() = diffusionMatrix->affinePart()->backend() * dirichletVector_->components()[qq]->backend();
////          rhsComponents.push_back(component);
////          rhsCoefficients.push_back(std::make_shared< const CoefficientType >(
////                                      "-1.0*(" + dirichletVector_->coefficients()[qq]->expression() + ")"));
////        }
////        if (dirichletVector_->hasAffinePart())
////          rhsAffinePart->backend() += diffusionMatrix->affinePart()->backend() * dirichletVector_->affinePart()->backend();
////      }
////      // * diffusion * dirichletAffinePart
////      if (dirichletVector_->hasAffinePart()) {
////        for (size_t qq = 0; qq < diffusionMatrix->components().size(); ++qq) {
////          auto component = std::make_shared< VectorType >(space_->mapper().size());
////          component->backend() = diffusionMatrix->components()[qq]->backend() * dirichletVector_->affinePart()->backend();
////          rhsComponents.push_back(component);
////          rhsCoefficients.push_back(std::make_shared< const CoefficientType >(
////                                      "-1.0*(" + diffusionMatrix->coefficients()[qq]->expression() + ")"));
////        }
////      }
//      // * diffusion * dirichlet
//      if (!diffusionMatrix->parametric()) {
//        auto tmp = std::make_shared< VectorType >(space_->mapper().size());
//        tmp->backend() = RangeFieldType(-1) * (diffusionMatrix->affinePart()->backend() * dirichletVector_->backend());
//        diffusionDirichletVector_ = std::make_shared< AffineParametricVectorType >(tmp);
//      } else {
//        std::vector< std::shared_ptr< VectorType > > components;
//        for (size_t pp = 0; pp < diffusionMatrix->components().size(); ++pp) {
////        for (size_t qq = 0; qq < dirichletVector_->components().size(); ++qq) {
//          auto component = std::make_shared< VectorType >(space_->mapper().size());
//          component->backend() = RangeFieldType(-1) * (diffusionMatrix->components()[pp]->backend() * dirichletVector_/*->components()[qq]*/->backend());
//          components.push_back(component);
////        }
//        }
//        if (!diffusionMatrix->hasAffinePart()) {
//          diffusionDirichletVector_ = std::make_shared< AffineParametricVectorType >(diffusionMatrix->paramSize(),
//                                                                                     components,
//                                                                                     diffusionMatrix->coefficients());
//        } else {
//          auto affinePart = std::make_shared< VectorType >(space_->mapper().size());
//          affinePart->backend() += diffusionMatrix->affinePart()->backend() * dirichletVector_->backend();
//          diffusionDirichletVector_ = std::make_shared< AffineParametricVectorType >(diffusionMatrix->paramSize(),
//                                                                                     components,
//                                                                                     diffusionMatrix->coefficients(),
//                                                                                     affinePart);
//        }
//      }
////      // actually create the rhs vector
////      if (rhsComponents.size() == 0) {
////        assert(rhsCoefficients.size() == 0);
////        assert(rhsHasAffinePart);
////        rhsVector_ = std::make_shared< AffineParametricVectorType >(rhsAffinePart);
////      } else {
////        assert(rhsComponents.size() == rhsCoefficients.size());
////        if (!rhsHasAffinePart)
////          rhsVector_ = std::make_shared< AffineParametricVectorType >(forceVector->paramSize()
////                                                                      + neumannVector->paramSize()
////                                                                      + diffusionMatrix->paramSize(),
////                                                                      rhsComponents,
////                                                                      rhsCoefficients);
////        else
////          rhsVector_ = std::make_shared< AffineParametricVectorType >(forceVector->paramSize()
////                                                                      + neumannVector->paramSize()
////                                                                      + diffusionMatrix->paramSize(),
////                                                                      rhsComponents,
////                                                                      rhsCoefficients,
////                                                                      rhsAffinePart);
////      } // actually create the rhs vector
//      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

//      out << prefix << "applying constraints... " << std::flush;
//      // adding the rhs constraints here, the matrix constraints were added above
//      if (forceVector_->hasAffinePart())
//        systemAssembler.addLocalConstraints(clearRows, *(forceVector_->affinePart()));
//      if (forceVector_->parametric()) {
//        for (size_t qq = 0; qq < forceVector_->components().size(); ++qq)
//          systemAssembler.addLocalConstraints(clearRows, *(forceVector_->components()[qq]));
//      }
//      if (neumannVector_->hasAffinePart())
//        systemAssembler.addLocalConstraints(clearRows, *(neumannVector_->affinePart()));
//      if (neumannVector_->parametric()) {
//        for (size_t qq = 0; qq < neumannVector_->components().size(); ++qq)
//          systemAssembler.addLocalConstraints(clearRows, *(neumannVector_->components()[qq]));
//      }
//      if (diffusionDirichletVector_->hasAffinePart())
//        systemAssembler.addLocalConstraints(clearRows, *(diffusionDirichletVector_->affinePart()));
//      if (diffusionDirichletVector_->parametric()) {
//        for (size_t qq = 0; qq < diffusionDirichletVector_->components().size(); ++qq)
//          systemAssembler.addLocalConstraints(clearRows, *(diffusionDirichletVector_->components()[qq]));
//      }
//      systemAssembler.applyConstraints();
//      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

//      // clean up
//      for (auto& element : neumannVectorAssemblers)
//        delete element;
//      for (auto& element : neumannFunctionals)
//        delete element;
//      for (auto& element : forceVectorAssemblers)
//        delete element;
//      for (auto& element : forceFunctionals)
//        delete element;
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
//    // first of all, get the corect parameters (the model returns empty ones for nonparametric functions)
//    const ParamType muDiffusion = model_->mapParam(mu, "diffusion");
//    const ParamType muForce = model_->mapParam(mu, "force");
////    const ParamType muDirichlet = model_->mapParam(mu, "dirichlet");
//    const ParamType muNeumann = model_->mapParam(mu, "neumann");
//    Dune::Timer timer;
//    out << prefix << "computing system matrix...   " << std::flush;
//    std::shared_ptr< const MatrixType > systemMatrix;
//    if (systemMatrix_->parametric())
//      systemMatrix = systemMatrix_->fix(muDiffusion);
//    else
//      systemMatrix = systemMatrix_->affinePart();
//    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
//    out << prefix << "computing right hand side... " << std::flush;
//    timer.reset();
//    std::shared_ptr< VectorType > rhsVector(createVector());
//    if (forceVector_->parametric()) {
//      const auto forceVector = forceVector_->fix(muForce);
//      rhsVector->backend() = forceVector->backend();
//    } else {
//      rhsVector->backend() = forceVector_->affinePart()->backend();
//    }
//    if (neumannVector_->parametric()) {
//      const auto neumannVector = neumannVector_->fix(muNeumann);
//      rhsVector->backend() += neumannVector->backend();
//    } else {
//      rhsVector->backend() += neumannVector_->affinePart()->backend();
//    }
//    if (diffusionDirichletVector_->parametric()) {
//      const auto diffusionDirichletVector = diffusionDirichletVector_->fix(muDiffusion);
//      rhsVector->backend() += diffusionDirichletVector->backend();
//    } else {
//      rhsVector->backend() += diffusionDirichletVector_->affinePart()->backend();
//    }
//    out << "done (took " << timer.elapsed() << " sec)" << std::endl;

//    out << prefix << "solving linear system (of size " << systemMatrix->rows()
//        << "x" << systemMatrix->cols() << ")" << std::endl;
//    const std::string linearSolverType = linearSolverSettings.get< std::string >("type");
//    out << prefix << "  using '" << linearSolverType << "'... " << std::flush;
//    timer.reset();
//    typedef typename Dune::Stuff::LA::SolverInterface< MatrixType, VectorType > SolverType;
//    const std::shared_ptr< const SolverType > solver(Dune::Stuff::LA::createSolver< MatrixType, VectorType >(linearSolverType));
//    const size_t failure = solver->apply(*systemMatrix,
//                                         *rhsVector,
//                                         *solutionVector,
//                                         linearSolverSettings);
//    if (failure)
//      DUNE_THROW(Dune::MathError,
//                 "\n"
//                 << Dune::Stuff::Common::colorStringRed("ERROR:")
//                 << " linear solver '" << linearSolverType << "' reported error code " << failure << "!\n"
//                 << "  1: did not converge\n"
//                 << "  2: had numerical issues\n"
//                 << "  3: dude, I have no idea");
//    if (solutionVector->size() != space_->mapper().size())
//      DUNE_THROW(Dune::MathError,
//                 "\n"
//                 << Dune::Stuff::Common::colorStringRed("ERROR:")
//                 << " linear solver '" << linearSolverType << "' produced a solution of wrong size (is "
//                 << solutionVector->size() << ", should be " << space_->mapper().size() << ")!");
////    if (dirichletVector_->parametric()) {
////      const auto dirichletVector = dirichletVector_->fix(muDirichlet);
////      solutionVector->backend() += dirichletVector->backend();
////    } else
//      solutionVector->backend() += dirichletVector_/*->affinePart()*/->backend();
//    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
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

private:
//  std::shared_ptr< ConstDiscreteFunctionType > createConstAnsatzFunction(const std::shared_ptr< const VectorType > vector,
//                                                                         const std::string name) const
//  {
//    return std::make_shared< ConstDiscreteFunctionType >(*space_, vector, name);
//  }

//  void visualizeFunction(const std::shared_ptr< const ConstDiscreteFunctionType > discreteFunction,
//                         const std::string filename = id() + ".function",
//                         std::ostream& out = Dune::Stuff::Common::Logger().devnull(),
//                         const std::string prefix = "") const
//  {
//    // preparations
//    Dune::Timer timer;
//    out << prefix << "writing '" << discreteFunction->name() << "' to '" << filename;
//    if (dimDomain == 1)
//      out << ".vtp";
//    else
//      out << ".vtu";
//    out << "'... " << std::flush;
//    typedef Dune::VTKWriter< typename GridPartType::GridViewType > VTKWriterType;
//    VTKWriterType vtkWriter(gridPart_->gridView());
//    vtkWriter.addVertexData(discreteFunction);
//    vtkWriter.write(filename);
//    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
//  } // void visualizeFunction(...)

  const std::shared_ptr< const GridPartType > gridPart_;
  const std::shared_ptr< const BoundaryInfoType > boundaryInfo_;
  const std::shared_ptr< const ModelType > model_;
  bool initialized_;
  std::shared_ptr< const TestSpaceType > space_;
  std::shared_ptr< const PatternType > pattern_;
  AffineParametricMatrixType* systemMatrix_;
  std::shared_ptr< AffineParametricVectorType > forceVector_;
  std::shared_ptr< AffineParametricVectorType > neumannVector_;
  std::shared_ptr< /*AffineParametric*/VectorType > dirichletVector_;
  std::shared_ptr< AffineParametricVectorType > diffusionDirichletVector_;
}; // class ContinuousGalerkinDiscretizationGDT


} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_DISCRETIZATION_CG_GDT_HH
