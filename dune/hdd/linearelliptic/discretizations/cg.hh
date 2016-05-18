// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_CG_HH
#define DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_CG_HH

#include <memory>
#include <vector>

#include <dune/common/timer.hh>

#if HAVE_DUNE_GRID_MULTISCALE
# include <dune/grid/multiscale/provider/interface.hh>
#endif

#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/timedlogging.hh>
#include <dune/stuff/grid/layers.hh>
#include <dune/stuff/grid/provider.hh>
#include <dune/stuff/la/container.hh>
#include <dune/stuff/la/solver.hh>

#include <dune/gdt/spaces/cg.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/playground/operators/elliptic-cg.hh>
#include <dune/gdt/functionals/l2.hh>

#include "base.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Discretizations {


// forward, for friendlyness
template< class GridImp, class RangeFieldImp, int rangeDim, int polynomialOrder, Stuff::LA::ChooseBackend la_backend >
class BlockSWIPDG;


// forward, needed in the Traits
template< class GridImp, Stuff::Grid::ChooseLayer layer, class RangeFieldImp, int rangeDim, int polynomialOrder = 1,
#if HAVE_DUNE_FEM
          GDT::ChooseSpaceBackend space_backend = GDT::ChooseSpaceBackend::fem,
#elif HAVE_DUNE_PDELAB
          GDT::ChooseSpaceBackend space_backend = GDT::ChooseSpaceBackend::pdelab,
#else
# error No suitable space backend available!
#endif
          Stuff::LA::ChooseBackend la_backend = Stuff::LA::default_sparse_backend >
class CG;


namespace internal {


template< class GridImp, Stuff::Grid::ChooseLayer layer, class RangeFieldImp, int rangeDim, int polynomialOrder,
          GDT::ChooseSpaceBackend space_backend,
          Stuff::LA::ChooseBackend la_backend >
class CGTraits
  : public internal::ContainerBasedDefaultTraits< typename Stuff::LA::Container< RangeFieldImp, la_backend >::MatrixType,
                                                  typename Stuff::LA::Container< RangeFieldImp, la_backend >::VectorType >
{
  typedef internal::ContainerBasedDefaultTraits< typename Stuff::LA::Container< RangeFieldImp, la_backend >::MatrixType,
  typename Stuff::LA::Container< RangeFieldImp, la_backend >::VectorType > BaseType;
public:
  typedef CG< GridImp, layer, RangeFieldImp, rangeDim, polynomialOrder, space_backend, la_backend > derived_type;
  typedef GridImp GridType;
  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRange = rangeDim;
  static const unsigned int polOrder = polynomialOrder;

private:
  typedef GDT::Spaces::CGProvider< GridType, layer, space_backend, polOrder, RangeFieldType, dimRange > SpaceProvider;

  friend class CG< GridImp, layer, RangeFieldImp, rangeDim, polynomialOrder, space_backend, la_backend >;

public:
  typedef typename SpaceProvider::Type TestSpaceType;
  typedef TestSpaceType AnsatzSpaceType;
  typedef typename TestSpaceType::GridViewType GridViewType;
}; // class CGTraits


} // namespace internal


template< class GridImp, Stuff::Grid::ChooseLayer layer, class RangeFieldImp, int rangeDim, int polynomialOrder,
          GDT::ChooseSpaceBackend space_backend,
          Stuff::LA::ChooseBackend la_backend >
class CG
  : public ContainerBasedDefault< internal::CGTraits< GridImp, layer, RangeFieldImp, rangeDim, polynomialOrder, space_backend, la_backend > >
{
  typedef ContainerBasedDefault< internal::CGTraits< GridImp, layer, RangeFieldImp, rangeDim, polynomialOrder, space_backend, la_backend > > BaseType;
public:
  typedef internal::CGTraits< GridImp, layer, RangeFieldImp, rangeDim, polynomialOrder, space_backend, la_backend > Traits;
  typedef GridImp GridType;
  using typename BaseType::ProblemType;
  using typename BaseType::BoundaryInfoType;
  using typename BaseType::TestSpaceType;
  using typename BaseType::AnsatzSpaceType;
  using typename BaseType::MatrixType;
  using typename BaseType::VectorType;
  using typename BaseType::GridViewType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::PatternType;

protected:
  typedef typename Traits::SpaceProvider SpaceProvider;

  typedef Stuff::Grid::ProviderInterface< GridType >      GridProviderType;
#if HAVE_DUNE_GRID_MULTISCALE
  typedef grid::Multiscale::ProviderInterface< GridType > MsGridProviderType;
#endif
  using typename BaseType::AffinelyDecomposedMatrixType;
  using typename BaseType::AffinelyDecomposedVectorType;

  typedef typename ProblemType::DiffusionFactorType::NonparametricType DiffusionFactorType;
  typedef typename ProblemType::DiffusionTensorType::NonparametricType DiffusionTensorType;
  typedef GDT::Operators::EllipticCG< DiffusionFactorType, MatrixType, TestSpaceType,
                                      AnsatzSpaceType, GridViewType, DiffusionTensorType > EllipticOperatorType;

public:
  static std::string static_id()
  {
    return DiscretizationInterface< Traits >::static_id() + ".cg";
  }

  CG(GridProviderType& grid_provider,
     Stuff::Common::Configuration bound_inf_cfg,
     const ProblemType& prob,
     const int level = 0,
     const std::vector< std::string >& only_these_products = {})
    : BaseType(SpaceProvider::create(grid_provider, level),
               SpaceProvider::create(grid_provider, level),
               bound_inf_cfg,
               prob)
    , grid_provider_(grid_provider.copy())
    , only_these_products_(only_these_products)
  {
    // in that case we would have to build the elliptic operators like the dirichlet shift
    if (this->problem_.diffusion_factor()->parametric() && this->problem_.diffusion_tensor()->parametric())
      DUNE_THROW(NotImplemented, "Both parametric diffusion factor and tensor not supported!");
  } // CG(...)

#if HAVE_DUNE_GRID_MULTISCALE

  CG(MsGridProviderType& grid_provider,
     const Stuff::Common::Configuration& bound_inf_cfg,
     const ProblemType& prob,
     const int level_or_subdomain = 0,
     const std::vector< std::string >& only_these_products = {})
    : BaseType(SpaceProvider::create(grid_provider, level_or_subdomain),
               SpaceProvider::create(grid_provider, level_or_subdomain),
               bound_inf_cfg,
               prob)
    , grid_provider_(grid_provider.copy())
    , only_these_products_(only_these_products)
  {
    // in that case we would have to build the elliptic operators like the dirichlet shift
    if (this->problem_.diffusion_factor()->parametric() && this->problem_.diffusion_tensor()->parametric())
      DUNE_THROW(NotImplemented, "Both parametric diffusion factor and tensor not supported!");
  } // CG(...)

#endif // HAVE_DUNE_GRID_MULTISCALE

  void init(const bool prune = false)
  {
    if (this->container_based_initialized_)
      return;

    using namespace Dune::GDT;
    Dune::Timer timer;
    auto logger = Stuff::Common::TimedLogger().get("hdd.linearelliptic.discretizations.block-swipdg.init");
    auto dirichlet_vector = std::make_shared< AffinelyDecomposedVectorType >();

    auto& matrix = *(this->matrix_);
    auto& rhs = *(this->rhs_);
    const auto& space = this->test_space_;
    const auto& grid_view = space.grid_view();
    const auto& boundary_info = *(this->boundary_info_);

    logger.info() << "assembling... " << std::flush;
    timer.reset();
    GDT::SystemAssembler< TestSpaceType > system_assembler(space);

    pattern_ = std::make_shared< PatternType >(EllipticOperatorType::pattern(this->test_space(), this->test_space()));

    // project dirichlet boundary values
    typedef typename ProblemType::FunctionType::NonparametricType DirichletType;
    const auto& dirichlet = *(this->problem_.dirichlet());
    typedef GDT::DiscreteFunction< AnsatzSpaceType, VectorType > DiscreteFunctionType;
    typedef GDT::Operators::DirichletProjectionLocalizable< GridViewType, DirichletType, DiscreteFunctionType >
        DirichletProjectionOperator;
    std::vector< std::unique_ptr< DiscreteFunctionType > > dirichlet_projections;
    std::vector< std::unique_ptr< DirichletProjectionOperator > > dirichlet_projection_operators;
    for (DUNE_STUFF_SSIZE_T qq = 0; qq < dirichlet.num_components(); ++qq) {
      const size_t id = dirichlet_vector->register_component(new VectorType(space.mapper().size()),
                                                             dirichlet.coefficient(qq));
      dirichlet_projections.emplace_back(new DiscreteFunctionType(space, *(dirichlet_vector->component(id))));
      dirichlet_projection_operators.emplace_back(
            new DirichletProjectionOperator(grid_view,
                                            boundary_info,
                                            *(dirichlet.component(qq)),
                                            *(dirichlet_projections[dirichlet_projections.size() - 1])));
    }
    if (dirichlet.has_affine_part()) {
      dirichlet_vector->register_affine_part(new VectorType(space.mapper().size()));
      dirichlet_projections.emplace_back(new DiscreteFunctionType(space, *(dirichlet_vector->affine_part())));
      dirichlet_projection_operators.emplace_back(
            new DirichletProjectionOperator(grid_view,
                                            boundary_info,
                                            *(dirichlet.affine_part()),
                                            *(dirichlet_projections[dirichlet_projections.size() - 1])));
    }
    for (auto& projection_operator : dirichlet_projection_operators)
      system_assembler.add(*projection_operator, new Stuff::Grid::ApplyOn::BoundaryEntities< GridViewType >());

    // lhs operator
    const auto& diffusion_factor = *(this->problem_.diffusion_factor());
    const auto& diffusion_tensor = *(this->problem_.diffusion_tensor());
    std::vector< std::unique_ptr< EllipticOperatorType > > elliptic_operators;
    if (diffusion_factor.has_affine_part() && diffusion_tensor.has_affine_part()) {
      matrix.register_affine_part(new MatrixType(space.mapper().size(), space.mapper().size(), *pattern_));
      elliptic_operators.emplace_back(new EllipticOperatorType(*(diffusion_factor.affine_part()),
                                                               *(diffusion_tensor.affine_part()),
                                                               *(matrix.affine_part()),
                                                               space));
    }
    if (diffusion_factor.parametric() && !diffusion_tensor.parametric()) {
      for (DUNE_STUFF_SSIZE_T qq = 0; qq < diffusion_factor.num_components(); ++qq) {
        const size_t id = matrix.register_component(new MatrixType(space.mapper().size(),
                                                                   space.mapper().size(),
                                                                   *pattern_),
                                                    diffusion_factor.coefficient(qq));
        elliptic_operators.emplace_back(new EllipticOperatorType(*(diffusion_factor.component(qq)),
                                                                 *(diffusion_tensor.affine_part()),
                                                                 *(matrix.component(id)),
                                                                 space));
      }
    } else if (!diffusion_factor.parametric() && diffusion_tensor.parametric()) {
      for (DUNE_STUFF_SSIZE_T qq = 0; qq < diffusion_tensor.num_components(); ++qq) {
        const size_t id = matrix.register_component(new MatrixType(space.mapper().size(),
                                                                   space.mapper().size(),
                                                                   *pattern_),
                                                    diffusion_tensor.coefficient(qq));
        elliptic_operators.emplace_back(new EllipticOperatorType(*(diffusion_factor.affine_part()),
                                                                 *(diffusion_tensor.component(qq)),
                                                                 *(matrix.component(id)),
                                                                 space));
      }
    } else if (diffusion_factor.parametric() && diffusion_tensor.parametric()) {
      DUNE_THROW(Stuff::Exceptions::internal_error, "This should not happen!");
    }
    for (auto& elliptic_operator : elliptic_operators)
      system_assembler.add(*elliptic_operator);

    // rhs functional
    // * force
    typedef typename ProblemType::FunctionType::NonparametricType ForceType;
    const auto& force = *(this->problem_.force());
    typedef GDT::Functionals::L2Volume< ForceType, VectorType, TestSpaceType > L2VolumeFunctionalType;
    std::vector< std::unique_ptr< L2VolumeFunctionalType > > force_functionals;
    for (DUNE_STUFF_SSIZE_T qq = 0; qq < force.num_components(); ++qq) {
      const size_t id = rhs.register_component(new VectorType(space.mapper().size()), force.coefficient(qq));
      force_functionals.emplace_back(new L2VolumeFunctionalType(*(force.component(qq)),
                                                                *(rhs.component(id)),
                                                                space));
    }
    if (force.has_affine_part()) {
      rhs.register_affine_part(new VectorType(space.mapper().size()));
      force_functionals.emplace_back(new L2VolumeFunctionalType(*(force.affine_part()),
                                                                *(rhs.affine_part()),
                                                                space));
    }
    for (auto& force_functional : force_functionals)
      system_assembler.add(*force_functional);
    // * neumann
    typedef typename ProblemType::FunctionType::NonparametricType NeumannType;
    const auto& neumann = *(this->problem_.neumann());
    typedef GDT::Functionals::L2Face< NeumannType, VectorType, TestSpaceType > L2FaceFunctionalType;
    std::vector< std::unique_ptr< L2FaceFunctionalType > > neumann_functionals;
    for (DUNE_STUFF_SSIZE_T qq = 0; qq < neumann.num_components(); ++qq) {
      const size_t id = rhs.register_component(new VectorType(space.mapper().size()), neumann.coefficient(qq));
      neumann_functionals.emplace_back(new L2FaceFunctionalType(*(neumann.component(qq)),
                                                                *(rhs.component(id)),
                                                                space));
    }
    if (neumann.has_affine_part()) {
      if (!rhs.has_affine_part())
        rhs.register_affine_part(new VectorType(space.mapper().size()));
      neumann_functionals.emplace_back(new L2FaceFunctionalType(*(neumann.affine_part()),
                                                                *(rhs.affine_part()),
                                                                space));
    }
    for (auto& neumann_functional : neumann_functionals)
      system_assembler.add(*neumann_functional,
                           new Stuff::Grid::ApplyOn::NeumannIntersections< GridViewType >(boundary_info));

    // constraints
    Spaces::DirichletConstraints< typename GridViewType::Intersection >
        clear_and_set_dirichlet_rows(boundary_info, space.mapper().size());
    Spaces::DirichletConstraints< typename GridViewType::Intersection >
        clear_dirichlet_rows(boundary_info, space.mapper().size(), false);
    system_assembler.add(clear_and_set_dirichlet_rows, new Stuff::Grid::ApplyOn::BoundaryEntities< GridViewType >());
    system_assembler.add(clear_dirichlet_rows, new Stuff::Grid::ApplyOn::BoundaryEntities< GridViewType >());
    // do the actual assembling
    system_assembler.walk();

    logger.info() << "done (took " << timer.elapsed() << "s)" << std::endl;

    logger.info() << "computing dirichlet shift... " << std::flush;
    VectorType tmp(space.mapper().size());
    if (matrix.has_affine_part() && dirichlet_vector->has_affine_part()) {
      if (!rhs.has_affine_part())
        rhs.register_affine_part(new VectorType(space.mapper().size()));
      matrix.affine_part()->mv(*(dirichlet_vector->affine_part()), tmp);
      *(rhs.affine_part()) -= tmp;
    }
    if (matrix.has_affine_part()) {
      for (DUNE_STUFF_SSIZE_T qq = 0; qq < dirichlet_vector->num_components(); ++qq) {
        const size_t ind = rhs.register_component(new VectorType(space.mapper().size()),
                                                  dirichlet_vector->coefficient(qq));
        matrix.affine_part()->mv(*(dirichlet_vector->component(qq)), tmp);
        *(rhs.component(ind)) -= tmp;
      }
    }
    if (dirichlet_vector->has_affine_part()) {
      for (DUNE_STUFF_SSIZE_T qq = 0; qq < matrix.num_components(); ++qq) {
        const size_t ind = rhs.register_component(new VectorType(space.mapper().size()), matrix.coefficient(qq));
        matrix.component(qq)->mv(*(dirichlet_vector->affine_part()), tmp);
        *(rhs.component(ind)) -= tmp;
      }
    }
    Pymor::ParameterType param;
    for (auto key : matrix.parameter_type().keys())
      param.set(key, matrix.parameter_type().get(key));
    for (auto key : dirichlet_vector->parameter_type().keys())
      param.set(key, dirichlet_vector->parameter_type().get(key));
    for (DUNE_STUFF_SSIZE_T pp = 0; pp < matrix.num_components(); ++ pp) {
      for (DUNE_STUFF_SSIZE_T qq = 0; qq < dirichlet_vector->num_components(); ++qq) {
        const std::string expression = "(" + matrix.coefficient(pp)->expression()
                                       + ")*(" + dirichlet_vector->coefficient(qq)->expression() + ")";
        const size_t ind = rhs.register_component(new VectorType(space.mapper().size()),
                                                  new Pymor::ParameterFunctional(param, expression));
        const auto& matrix_component = *matrix.component(pp);
        const auto& dirichlet_component = *dirichlet_vector->component(qq);
        auto& rhs_component = *rhs.component(ind);
        rhs_component -= matrix_component * dirichlet_component;
      }
    }
    logger.info() << "done (took " << timer.elapsed() << " sec)" << std::endl;

    logger.info() << "assembling products... " << std::flush;
    timer.reset();
    this->assemble_products(only_these_products_,
                            clear_and_set_dirichlet_rows,
                            clear_dirichlet_rows,
                            2);
    logger.info() << "done (took " << timer.elapsed() << " sec)" << std::endl;

    logger.info() << "applying constraints... " << std::flush;
    timer.reset();
    // we always need an affine part in the system matrix for the dirichlet rows
    if (!matrix.has_affine_part())
      matrix.register_affine_part(new MatrixType(space.mapper().size(), space.mapper().size(), *pattern_));
    clear_and_set_dirichlet_rows.apply(*(matrix.affine_part()));
    for (DUNE_STUFF_SSIZE_T qq = 0; qq < matrix.num_components(); ++qq)
      clear_dirichlet_rows.apply(*(matrix.component(qq)));
    if (rhs.has_affine_part())
      clear_dirichlet_rows.apply(*(rhs.affine_part()));
    for (DUNE_STUFF_SSIZE_T qq = 0; qq < rhs.num_components(); ++qq)
      clear_dirichlet_rows.apply(*(rhs.component(qq)));
    logger.info() << "done (took " << timer.elapsed() << " sec)" << std::endl;

    // build parameter type
    this->inherit_parameter_type(matrix.parameter_type(), "lhs");
    this->inherit_parameter_type(rhs.parameter_type(), "rhs");
    this->inherit_parameter_type(dirichlet_vector->parameter_type(), "dirichlet");

    // finalize
    this->vectors_.insert(std::make_pair("dirichlet", dirichlet_vector));
    this->finalize_init(prune);
  } // ... init(...)

private:
  friend class BlockSWIPDG< GridImp, RangeFieldImp, rangeDim, polynomialOrder, la_backend >;

  template< class SetConstraints, class ClearConstraints >
  std::shared_ptr< AffinelyDecomposedMatrixType >
  make_zero_dirichlet_product(const std::shared_ptr< AffinelyDecomposedMatrixType >& prod,
                              const SetConstraints& set,
                              const ClearConstraints& clear)
  {
    auto ret = std::make_shared<AffinelyDecomposedMatrixType>(prod->copy());
    for (ssize_t qq = 0; qq < ret->num_components(); ++qq) {
      // clear rows
      clear.apply(*ret->component(qq));
      // clear columns
      for (const auto& column : clear.dirichlet_DoFs())
        ret->component(qq)->unit_col(column);
    }
    if (!ret->has_affine_part()) {
      // we can assume that there is at least one component
      assert(ret->num_components() > 0);
      ret->register_affine_part(new MatrixType(ret->component(0)->copy()));
      *ret->affine_part() *= 0;
    }
    // set rows
    set.apply(*ret->affine_part());
    // set columns
    for (const auto& column : set.dirichlet_DoFs())
      ret->affine_part()->unit_col(column);
    return ret;
  } // ... make_zero_dirichlet_product(...)

  using BaseType::pattern_;
  std::shared_ptr< GridProviderType > grid_provider_;
  const std::vector< std::string > only_these_products_;
}; // class CG


} // namespace Discretizations
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_CG_HH
