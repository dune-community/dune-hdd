// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_CG_HXX
#define DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_CG_HXX

#include "cg.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Discretizations {


template< class G, Stuff::Grid::ChooseLayer gl, class R, int r, int p, GDT::ChooseSpaceBackend s, Stuff::LA::ChooseBackend la >
    std::string CG< G, gl, R, r, p, s, la >::
static_id()
{
  return DiscretizationInterface< Traits >::static_id() + ".cg";
}

template< class G, Stuff::Grid::ChooseLayer gl, class R, int r, int p, GDT::ChooseSpaceBackend s, Stuff::LA::ChooseBackend la >
    CG< G, gl, R, r, p, s, la >::
CG(typename CG< G, gl, R, r, p, s, la >::GridProviderType& grid_provider,
   Stuff::Common::Configuration bound_inf_cfg,
   const typename CG< G, gl, R, r, p, s, la >::ProblemType& prob,
   const int level,
   const std::vector< std::string >& only_these_products)
  : BaseType(SpaceProvider::create(grid_provider, level),
             SpaceProvider::create(grid_provider, level),
             bound_inf_cfg,
             prob)
  , only_these_products_(only_these_products)
{
  // in that case we would have to build the elliptic operators like the dirichlet shift
  if (this->problem_.diffusion_factor()->parametric() && this->problem_.diffusion_tensor()->parametric())
    DUNE_THROW(NotImplemented, "Both parametric diffusion factor and tensor not supported!");
} // CG(...)

#if HAVE_DUNE_GRID_MULTISCALE

template< class G, Stuff::Grid::ChooseLayer gl, class R, int r, int p, GDT::ChooseSpaceBackend s, Stuff::LA::ChooseBackend la >
    CG< G, gl, R, r, p, s, la >::
CG(typename CG< G, gl, R, r, p, s, la >::MsGridProviderType& grid_provider,
   const Stuff::Common::Configuration& bound_inf_cfg,
   const typename CG< G, gl, R, r, p, s, la >::ProblemType& prob,
   const int level_or_subdomain,
   const std::vector< std::string >& only_these_products)
  : BaseType(SpaceProvider::create(grid_provider, level_or_subdomain),
             SpaceProvider::create(grid_provider, level_or_subdomain),
             bound_inf_cfg,
             prob)
  , only_these_products_(only_these_products)
{
  // in that case we would have to build the elliptic operators like the dirichlet shift
  if (this->problem_.diffusion_factor()->parametric() && this->problem_.diffusion_tensor()->parametric())
    DUNE_THROW(NotImplemented, "Both parametric diffusion factor and tensor not supported!");
} // CG(...)

#endif // HAVE_DUNE_GRID_MULTISCALE

template< class G, Stuff::Grid::ChooseLayer gl, class R, int r, int p, GDT::ChooseSpaceBackend s, Stuff::LA::ChooseBackend la >
    void CG< G, gl, R, r, p, s, la >::
init(const bool prune)
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

  // build parameter type (matrix and rhs are done in base)
  this->inherit_parameter_type(dirichlet_vector->parameter_type(), "dirichlet");

  // finalize
  this->vectors_.insert(std::make_pair("dirichlet", dirichlet_vector));
  this->finalize_init(prune);
} // ... init(...)


} // namespace Discretizations
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_CG_HXX
