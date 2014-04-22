// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_CG_WITH_DUNE_GDT_HH
#define DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_CG_WITH_DUNE_GDT_HH

#include <memory>
#include <sstream>

#include <dune/common/exceptions.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/timer.hh>

#include <dune/grid/part/interface.hh>

#include <dune/fem/misc/mpimanager.hh>

#include <dune/stuff/common/logging.hh>
#include <dune/stuff/discretefunction/projection/dirichlet.hh>
#include <dune/stuff/common/parameter/tree.hh>

#include <dune/gdt/space/continuouslagrange/fem.hh>
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

#include <dune/pymor/common/exceptions.hh>
#include <dune/pymor/la/container/eigen.hh>
#include <dune/pymor/la/container/affine.hh>
#include <dune/pymor/functionals/default.hh>
#include <dune/pymor/functionals/affine.hh>
#include <dune/pymor/operators/eigen.hh>
#include <dune/pymor/operators/affine.hh>

#include "../problems/interfaces.hh"
#include "interfaces.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Discretization {


// forward, needed in the Traits
template< class GridPartImp, class RangeFieldImp, int rangeDim, int polynomialOrder, bool scalarDiffusion = true >
class ContinuousGalerkinWithDuneGDT;


/**
 *  \brief  Traits for SolverContinuousGalerkinDD
 */
template< class GridPartImp, class RangeFieldImp, int rangeDim, int polynomialOrder, bool scalarDiffusion = true >
class ContinuousGalerkinWithDuneGDTTraits
{
public:
  typedef ContinuousGalerkinWithDuneGDT<  GridPartImp, RangeFieldImp, rangeDim,
                                          polynomialOrder, scalarDiffusion >  derived_type;
  typedef typename GridPartImp::Traits                                        GridPartTraits;
  typedef Dune::grid::Part::Interface< GridPartTraits >                       GridPartType;
  typedef typename GridPartType::ctype  DomainFieldType;
  static const unsigned int             dimDomain = GridPartType::dimension;
  static const unsigned int             polOrder = polynomialOrder;
  typedef RangeFieldImp                 RangeFieldType;
  static const unsigned int             dimRange = rangeDim;
  typedef Dune::Stuff::GridboundaryInterface< typename GridPartType::GridViewType >                 BoundaryInfoType;
  typedef ProblemInterface< DomainFieldType, dimDomain, RangeFieldType, dimRange, scalarDiffusion > ProblemType;
  typedef Dune::Pymor::LA::EigenRowMajorSparseMatrix< double >  MatrixType;
  typedef Dune::Pymor::LA::EigenDenseVector< double >           VectorType;
private:
  typedef Dune::Pymor::Operators::EigenRowMajorSparse< double > OperatorComponentType;
  typedef Dune::Pymor::Functionals::VectorBased< VectorType >   FunctionalComponentType;
public:
  typedef Dune::Pymor::Functionals::LinearAffinelyDecomposedVectorBased< VectorType >             FunctionalType;
  typedef Dune::Pymor::Operators::LinearAffinelyDecomposedContainerBased< OperatorComponentType > OperatorType;
}; // class ContinuousGalerkinWithDuneGDTTraits


template< class GridPartImp, class RangeFieldImp, int rangeDim, int polynomialOrder, bool scalarDiffusion >
class ContinuousGalerkinWithDuneGDT
    : public DiscretizationInterface< ContinuousGalerkinWithDuneGDTTraits<  GridPartImp,
                                                                            RangeFieldImp,
                                                                            rangeDim,
                                                                            polynomialOrder,
                                                                            scalarDiffusion > >
{
public:
  typedef ContinuousGalerkinWithDuneGDTTraits<  GridPartImp, RangeFieldImp, rangeDim,
                                                polynomialOrder, scalarDiffusion > Traits;
private:
  typedef DiscretizationInterface< Traits > BaseType;
public:
  typedef typename Traits::GridPartType GridPartType;
  static const int polOrder = Traits::polOrder;

  typedef typename Traits::DomainFieldType  DomainFieldType;
  static const unsigned int                 dimDomain = GridPartType::dimension;
  typedef typename Traits::RangeFieldType   RangeFieldType;
  static const unsigned int                 dimRange = Traits::dimRange;

  typedef typename Traits::ProblemType      ProblemType;
  typedef typename Traits::BoundaryInfoType BoundaryInfoType;

  typedef GDT::ContinuousLagrangeSpace::FemWrapper< GridPartType, polOrder, RangeFieldType, dimRange > TestSpaceType;
  typedef TestSpaceType                         AnsatzSpaceType;
  typedef typename TestSpaceType::PatternType   PatternType;

  typedef typename Traits::OperatorType   OperatorType;
  typedef typename Traits::FunctionalType FunctionalType;

private:
  typedef typename Traits::MatrixType       MatrixType;
  typedef typename Traits::VectorType       VectorType;
  typedef typename MatrixType::BackendType  MatrixBackendType;
  typedef typename VectorType::BackendType  VectorBackendType;
  typedef Pymor::LA::AffinelyDecomposedContainer< MatrixType > AffinelyDecomposedMatrixType;
  typedef Pymor::LA::AffinelyDecomposedContainer< VectorType > AffinelyDecomposedVectorType;
  typedef GDT::DiscreteFunctionDefaultConst< AnsatzSpaceType, VectorBackendType > ConstDiscreteFunctionType;

public:
//  typedef Stuff::Common::ExtendedParameterTree SettingsType;

  static std::string static_id()
  {
    return typename DiscretizationInterface< Traits >::static_id() + ".cg-with-dune-gdt";
  }

  ContinuousGalerkinWithDuneGDT(const std::shared_ptr< const GridPartType > gP,
                                const std::shared_ptr< const BoundaryInfoType > bI,
                                const std::shared_ptr< const ProblemType > prob)
    : BaseType(*prob)
    , gridPart_(gP)
    , boundaryInfo_(bI)
    , problem_(prob)
    , initialized_(false)
  {
    // sanity checks
    std::stringstream msg;
    size_t failure = 0;
    // * integration orders
    if (problem_->diffusion()->order() < 0) {
      msg << "negative integration order given for the diffusion!\n";
      ++failure;
    }
    if (problem_->force()->order() < 0) {
      msg << "negative integration order given for the force!\n";
      ++failure;
    }
    if (problem_->neumann()->order() < 0) {
      msg << "negative integration order given for the neumann values!\n";
      ++failure;
    }
    // * parametrization
    if (problem_->parametric()) {
      if (!problem_->diffusion()->affinely_decomposable()
          || !problem_->force()->affinely_decomposable()
          || !problem_->dirichlet()->affinely_decomposable()
          || !problem_->neumann()->affinely_decomposable()) {
        msg << "only implemented for nonparametric or affinely decomposable functions!\n";
        ++failure;
      }
    }
    if (failure)
      DUNE_PYMOR_THROW(Pymor::Exception::wrong_input, msg.str());
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

  const ProblemType& problem() const
  {
    return problem_;
  }

  VectorType create_vector() const
  {
    return VectorType(space_->mapper().size());
  }

  void visualize(const VectorType& vector, const std::string filename, const std::string name) const
  {
    // preparations
    if (vector.dim() != space_->mapper().size())
      DUNE_PYMOR_THROW(Pymor::Exception::sizes_do_not_match,
                       "the dim of vector (" << vector.dim() << ") does not match the size of the ansatz space ("
                       << space_->mapper().size() << ")!");
    const std::shared_ptr< const ConstDiscreteFunctionType > discreteFunction = createConstAnsatzFunction(vector, name);
    visualizeFunction(discreteFunction, filename, Dune::Stuff::Common::Logger().devnull());
  }

  void initialize(std::ostream& out = Dune::Stuff::Common::Logger().devnull(),
                  const std::string prefix = "")
  {
    if (!initialized_) {
      Dune::Timer timer;
      systemMatrix_ = std::make_shared< AffinelyDecomposedMatrixType >();
      rhsVector_ = std::make_shared< AffinelyDecomposedVectorType >();
      dirichletVector_ = std::make_shared< AffinelyDecomposedVectorType >();

      out << prefix << "projecting dirichlet boundary values... " << std::flush;
      typedef GDT::DiscreteFunctionDefault< AnsatzSpaceType, VectorBackendType > DiscreteFunctionType;
      DiscreteFunctionType discreteFunction(*space_, std::make_shared< VectorBackendType >(space_->mapper().size()));
      for (size_t qq = 0; qq < problem_->dirichlet()->num_components(); ++qq) {
        Dune::Stuff::DiscreteFunction::project(*boundaryInfo_,
                                               *(problem_->dirichlet()->component(qq)),
                                               discreteFunction);
        dirichletVector_->register_component(new VectorType(new VectorBackendType(*(discreteFunction.vector()))),
                                             problem_->dirichlet()->coefficient(qq));
      }
      if (!problem_->dirichlet()->parametric() || problem_->dirichlet()->has_affine_part()) {
        Dune::Stuff::DiscreteFunction::project(*boundaryInfo_,
                                               *(problem_->dirichlet()->affine_part()),
                                               discreteFunction);
        dirichletVector_->register_affine_part(new VectorType(new VectorBackendType(*(discreteFunction.vector()))));
      }
      // clearing tmp storage, do not use discreteFunction after this point!
      discreteFunction.vector().reset();
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      out << prefix << "assembing system... " << std::flush;
      timer.reset();
      typedef typename ProblemType::DiffusionType::NonparametricType  DiffusionType;
      typedef typename ProblemType::ForceType::NonparametricType      ForceType;
      typedef typename ProblemType::NeumannType::NonparametricType    NeumannType;
      // prepare matrix, pattern, vector and assembler
      pattern_ = std::shared_ptr< const PatternType >(space_->computePattern());
      typedef GDT::SystemAssembler< TestSpaceType > SystemAssemblerType;
      SystemAssemblerType systemAssembler(*space_);
      // * elliptic diffusion operator
      typedef GDT::LocalOperator::Codim0Integral< GDT::LocalEvaluation::Elliptic< DiffusionType > >
          EllipticOperatorType;
      typedef GDT::LocalAssembler::Codim0Matrix< EllipticOperatorType > LocalMatrixAssemblerType;
      std::vector< EllipticOperatorType* > diffusionOperators;
      std::vector< LocalMatrixAssemblerType* > systemMatrix_Assemblers;
      for (size_t qq = 0; qq < problem_->diffusion()->num_components(); ++qq) {
        diffusionOperators.push_back(new EllipticOperatorType(*(problem_->diffusion()->component(qq))));
        systemMatrix_Assemblers.push_back(new LocalMatrixAssemblerType(*(diffusionOperators[qq])));
        systemMatrix_->register_component(new MatrixType(new MatrixBackendType(space_->mapper().size(),
                                                                              space_->mapper().size(),
                                                                              *pattern_)),
                                         problem_->diffusion()->coefficient(qq));
        systemAssembler.addLocalAssembler(*(systemMatrix_Assemblers[qq]), systemMatrix_->component(qq)->backend());
      }
      if (!problem_->diffusion()->parametric() || problem_->diffusion()->has_affine_part()) {
        diffusionOperators.push_back(new EllipticOperatorType(*(problem_->diffusion()->affine_part())));
        systemMatrix_Assemblers.push_back(new LocalMatrixAssemblerType(*(diffusionOperators[diffusionOperators.size() - 1])));
        systemMatrix_->register_affine_part(new MatrixType(new MatrixBackendType(space_->mapper().size(),
                                                                                space_->mapper().size(),
                                                                                *pattern_)));
        systemAssembler.addLocalAssembler(*(systemMatrix_Assemblers[systemMatrix_Assemblers.size() - 1]),
                                          systemMatrix_->affine_part()->backend());
      }
      //   * L2 force functional
      typedef GDT::LocalFunctional::Codim0Integral< GDT::LocalEvaluation::Product< ForceType > > L2VolumeFunctionalType;
      typedef GDT::LocalAssembler::Codim0Vector< L2VolumeFunctionalType > LocalVolumeVectorAssemblerType;
      std::vector< L2VolumeFunctionalType* > forceFunctionals;
      std::vector< LocalVolumeVectorAssemblerType* > forceVectorAssemblers;
      for (size_t qq = 0; qq < problem_->force()->num_components(); ++qq) {
        forceFunctionals.push_back(new L2VolumeFunctionalType(*(problem_->force()->component(qq))));
        forceVectorAssemblers.push_back(new LocalVolumeVectorAssemblerType(*(forceFunctionals[qq])));
        const size_t ind = rhsVector_->register_component(new VectorType(space_->mapper().size()),
                                                         problem_->force()->coefficient(qq));
        systemAssembler.addLocalAssembler(*(forceVectorAssemblers[qq]), rhsVector_->component(ind)->backend());
      }
      if (!problem_->force()->parametric() || problem_->force()->has_affine_part()) {
        forceFunctionals.push_back(new L2VolumeFunctionalType(*(problem_->force()->affine_part())));
        forceVectorAssemblers.push_back(new LocalVolumeVectorAssemblerType(*(forceFunctionals[forceFunctionals.size() - 1])));
        if (!rhsVector_->has_affine_part())
          rhsVector_->register_affine_part(new VectorType(space_->mapper().size()));
        systemAssembler.addLocalAssembler(*(forceVectorAssemblers[forceVectorAssemblers.size() - 1]),
                                          rhsVector_->affine_part()->backend());
      }
      //   * L2 neumann functional
      typedef GDT::LocalFunctional::Codim1Integral< GDT::LocalEvaluation::Product< NeumannType > > L2FaceFunctionalType;
      typedef GDT::LocalAssembler::Codim1Vector< L2FaceFunctionalType > LocalFaceVectorAssemblerType;
      std::vector< L2FaceFunctionalType* > neumannFunctionals;
      std::vector< LocalFaceVectorAssemblerType* > neumannVectorAssemblers;
      for (size_t qq = 0; qq < problem_->neumann()->num_components(); ++qq) {
        neumannFunctionals.push_back(new L2FaceFunctionalType(*(problem_->neumann()->component(qq))));
        neumannVectorAssemblers.push_back(new LocalFaceVectorAssemblerType(*(neumannFunctionals[qq])));
        const size_t ind = rhsVector_->register_component(new VectorType(space_->mapper().size()),
                                                         problem_->neumann()->coefficient(qq));
        systemAssembler.addLocalAssembler(*(neumannVectorAssemblers[qq]),
                                          typename SystemAssemblerType::AssembleOnNeumann(*boundaryInfo_),
                                          rhsVector_->component(ind)->backend());
      }
      if (!problem_->neumann()->parametric() || problem_->neumann()->has_affine_part()) {
        neumannFunctionals.push_back(new L2FaceFunctionalType(*(problem_->neumann()->affine_part())));
        neumannVectorAssemblers.push_back(new LocalFaceVectorAssemblerType(*(neumannFunctionals[neumannFunctionals.size() - 1])));
        if (!rhsVector_->has_affine_part())
          rhsVector_->register_affine_part(new VectorType(space_->mapper().size()));
        systemAssembler.addLocalAssembler(*(neumannVectorAssemblers[neumannVectorAssemblers.size() - 1]),
                                          typename SystemAssemblerType::AssembleOnNeumann(*boundaryInfo_),
                                          rhsVector_->affine_part()->backend());
      }
      // do the actual work
      systemAssembler.assemble();

      // compute the dirichlet shift
      if (systemMatrix_->has_affine_part() && dirichletVector_->has_affine_part()) {
        if (!rhsVector_->has_affine_part())
          rhsVector_->register_affine_part(new VectorType(space_->mapper().size()));
        rhsVector_->affine_part()->backend().backend() -= systemMatrix_->affine_part()->backend().backend()
                                                         * dirichletVector_->affine_part()->backend().backend();
      }
      if (systemMatrix_->has_affine_part()) {
        for (size_t qq = 0; qq < dirichletVector_->num_components(); ++qq) {
          const size_t ind = rhsVector_->register_component(new VectorType(space_->mapper().size()),
                                                           dirichletVector_->coefficient(qq));
          rhsVector_->component(ind)->backend().backend() = -1.0
                                                           * systemMatrix_->affine_part()->backend().backend()
                                                           * dirichletVector_->component(qq)->backend().backend();
        }
      }
      if (dirichletVector_->has_affine_part()) {
        for (size_t qq = 0; qq < systemMatrix_->num_components(); ++qq) {
          const size_t ind = rhsVector_->register_component(new VectorType(space_->mapper().size()),
                                                           systemMatrix_->coefficient(qq));
          rhsVector_->component(ind)->backend().backend() = -1.0
                                                           * systemMatrix_->component(qq)->backend().backend()
                                                           * dirichletVector_->affine_part()->backend().backend();
        }
      }
      Pymor::ParameterType diffusionDirichletMu;
      for (auto key : systemMatrix_->parameter_type().keys())
        diffusionDirichletMu.set(key, systemMatrix_->parameter_type().get(key));
      for (auto key : dirichletVector_->parameter_type().keys())
        diffusionDirichletMu.set(key, dirichletVector_->parameter_type().get(key));
      for (size_t pp = 0; pp < systemMatrix_->num_components(); ++ pp) {
        for (size_t qq = 0; qq < dirichletVector_->num_components(); ++qq) {
          const std::string expression = "(" + systemMatrix_->coefficient(pp)->expression()
                                         + ")*(" + dirichletVector_->coefficient(qq)->expression() + ")";
          const size_t ind = rhsVector_->register_component(new VectorType(space_->mapper().size()),
                                                           new Pymor::ParameterFunctional(diffusionDirichletMu,
                                                                                          expression));
          rhsVector_->component(ind)->backend().backend() = -1.0
                                                           * systemMatrix_->component(pp)->backend().backend()
                                                           * dirichletVector_->component(qq)->backend().backend();
        }
      }
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      out << prefix << "applying constraints... " << std::flush;
      GDT::Constraints::Dirichlet<  typename GridPartType::GridViewType,
                                    RangeFieldType, true > clearAndSetRows(*boundaryInfo_,
                                                                           space_->mapper().maxNumDofs(),
                                                                           space_->mapper().maxNumDofs());
      GDT::Constraints::Dirichlet<  typename GridPartType::GridViewType,
                                    RangeFieldType, false > clearRows(*boundaryInfo_,
                                                                      space_->mapper().maxNumDofs(),
                                                                      space_->mapper().maxNumDofs());
      // we always need an affine shift in the system matrix for the dirichlet rows
      if (!systemMatrix_->has_affine_part())
        systemMatrix_->register_affine_part(new MatrixType(new MatrixBackendType(space_->mapper().size(),
                                                                                space_->mapper().size(),
                                                                                *pattern_)));
      systemAssembler.addLocalConstraints(clearAndSetRows, systemMatrix_->affine_part()->backend());
      for (size_t qq = 0; qq < systemMatrix_->num_components(); ++qq)
        systemAssembler.addLocalConstraints(clearRows, systemMatrix_->component(qq)->backend());
      if (rhsVector_->has_affine_part())
        systemAssembler.addLocalConstraints(clearRows, rhsVector_->affine_part()->backend());
      for (size_t qq = 0; qq < systemMatrix_->num_components(); ++qq)
        systemAssembler.addLocalConstraints(clearRows, rhsVector_->component(qq)->backend());
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
      for (auto& element : systemMatrix_Assemblers)
        delete element;
      for (auto& element : diffusionOperators)
        delete element;

      // inherit parameter type
      inherit_parameter_type(systemMatrix_->parameter_type(), "lhs");
      inherit_parameter_type(rhsVector_->parameter_type(), "rhs");
      inherit_parameter_type(dirichletVector_->parameter_type(), "dirichlet");

      // done
      initialized_ = true;
    } // if !(initialized_)
  } // void initialize(...)

  bool initialized() const
  {
    return initialized_;
  }

  std::vector< std::string > available_operators() const
  {
    return { "lhs" };
  }

  OperatorType get_operator(const std::string id) const
  {
    if (id != "lhs") DUNE_PYMOR_THROW(Pymor::Exception::key_is_not_valid, "id has to be 'lhs' (is '" << id << "')!");
    return OperatorType(*systemMatrix_);
  }

  std::vector< std::string > available_functionals() const
  {
    return { "rhs" };
  }

  FunctionalType get_functional(const std::string id) const
  {
    if (id != "rhs") DUNE_PYMOR_THROW(Pymor::Exception::key_is_not_valid, "id has to be 'rhs' (is '" << id << "')!");
    return FunctionalType(*rhsVector_);
  }

  std::vector< std::string > solver_options() const
  {
    return { "problem" };
  }

  std::string solver_options(const std::string context) const
  {
    if (context != "problem")
      DUNE_PYMOR_THROW(Pymor::Exception::key_is_not_valid, "context has to be 'problem' (is '" << context << "')!");
    return OperatorType::invert_options()[0];
  }

  void solve(VectorType& vector,
             const Pymor::Parameter mu = Pymor::Parameter(),
             std::ostream& out = Dune::Stuff::Common::Logger().devnull(),
             const std::string prefix = "") const
  {
    if (mu.type() != Pymor::Parametric::parameter_type())
      DUNE_PYMOR_THROW(Pymor::Exception::wrong_parameter_type,
                       "the type of mu (" << mu.type() << ") does not match the parameter_type of this ("
                       << Pymor::Parametric::parameter_type() << ")!");
    Dune::Timer timer;
    // compute right hand side vector
    std::shared_ptr< const VectorType > rhs;
    if (!rhsVector_->parametric())
      rhs = rhsVector_->affine_part();
    else {
      out << prefix << "computing rhs... " << std::flush;
      timer.reset();
      Pymor::Parameter muRhs = Pymor::Parametric::map_parameter(mu, "rhs");
      rhs = std::make_shared< const VectorType >(rhsVector_->freeze_parameter(muRhs));
      out << "done (took " << timer.elapsed() << "s)" << std::endl;
    }
    const OperatorType lhsOperator(*systemMatrix_);
    if (lhsOperator.parametric()) {
      out << prefix << "computing lhs... " << std::flush;
      timer.reset();
      Pymor::Parameter muLhs = Pymor::Parametric::map_parameter(mu, "lhs");
      const auto frozenOperator = lhsOperator.freeze_parameter(muLhs);
      out << "done (took " << timer.elapsed() << "s)" << std::endl;
      const std::string option = frozenOperator.invert_options()[0];
      out << prefix << "solving with '" << option << "' option... " << std::flush;
      timer.reset();
      frozenOperator.apply_inverse(*rhs, vector, option);
      out << "done (took " << timer.elapsed() << "s)" << std::endl;
    } else {
      const auto nonparametricOperator = lhsOperator.affine_part();
      const std::string option = nonparametricOperator.invert_options()[0];
      out << prefix << "solving with '" << option << "' option... " << std::flush;
      timer.reset();
      nonparametricOperator.apply_inverse(*rhs, vector, option);
      out << "done (took " << timer.elapsed() << "s)" << std::endl;
    }
    out << prefix << "applying dirichlet correction... " << std::flush;
    timer.reset();
    std::shared_ptr< const VectorType > dirichletCorrection;
    if (dirichletVector_->parametric()) {
      Pymor::Parameter muDirichlet = Pymor::Parametric::map_parameter(mu, "dirichlet");
      dirichletCorrection = std::make_shared< const VectorType >(dirichletVector_->freeze_parameter(muDirichlet));
    } else
      dirichletCorrection = dirichletVector_->affine_part();
    vector.backend().backend() += dirichletCorrection->backend().backend();
    out << "done (took " << timer.elapsed() << "s)" << std::endl;
  }

private:
  std::shared_ptr< ConstDiscreteFunctionType > createConstAnsatzFunction(const VectorType& vector,
                                                                         const std::string name) const
  {
    return std::make_shared< ConstDiscreteFunctionType >(*space_,
                                                         std::make_shared< VectorBackendType >(vector.backend()),
                                                         name);
  }

  void visualizeFunction(const std::shared_ptr< const ConstDiscreteFunctionType > discreteFunction,
                         const std::string filename,
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

  const std::shared_ptr< const GridPartType > gridPart_;
  const std::shared_ptr< const BoundaryInfoType > boundaryInfo_;
  const std::shared_ptr< const ProblemType > problem_;
  bool initialized_;
  std::shared_ptr< const TestSpaceType > space_;
  std::shared_ptr< const PatternType > pattern_;
  std::shared_ptr< AffinelyDecomposedVectorType > dirichletVector_;
  std::shared_ptr< AffinelyDecomposedMatrixType > systemMatrix_;
  std::shared_ptr< AffinelyDecomposedVectorType > rhsVector_;
}; // class ContinuousGalerkinWithDuneGDT


} // namespace Discretization
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_CG_WITH_DUNE_GDT_HH
