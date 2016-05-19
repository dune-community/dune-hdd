// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_SWIPDG_HXX
#define DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_SWIPDG_HXX

#include "swipdg.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Discretizations {


template< class G, Stuff::Grid::ChooseLayer gl, class R, int r, int p, GDT::ChooseSpaceBackend s,
          Stuff::LA::ChooseBackend la >
    std::string SWIPDG< G, gl, R, r, p, s, la >::
static_id()
{
  return DiscretizationInterface< Traits >::static_id() + ".swipdg";
}

template< class G, Stuff::Grid::ChooseLayer gl, class R, int r, int p, GDT::ChooseSpaceBackend s,
          Stuff::LA::ChooseBackend la >
    SWIPDG< G, gl, R, r, p, s, la >::
SWIPDG(typename SWIPDG< G, gl, R, r, p, s, la >::GridProviderType& grid_provider,
       const Stuff::Common::Configuration& bound_inf_cfg,
       const typename SWIPDG< G, gl, R, r, p, s, la >::ProblemType& prob,
       const int level_or_subdomain,
       const std::vector< std::string >& only_these_products)
  : BaseType(SpaceProvider::create(grid_provider, level_or_subdomain),
             SpaceProvider::create(grid_provider, level_or_subdomain),
             bound_inf_cfg,
             prob)
  , beta_(GDT::LocalEvaluation::SIPDG::internal::default_beta(dimDomain))
  , only_these_products_(only_these_products)
{
  // in case of parametric diffusion tensor this discretization is not affinely decomposable any more
  if (this->problem_.diffusion_tensor()->parametric())
    DUNE_THROW(NotImplemented, "The diffusion tensor must not be parametric!");
  if (!this->problem_.diffusion_tensor()->has_affine_part())
    DUNE_THROW(Stuff::Exceptions::wrong_input_given, "The diffusion tensor must not be empty!");
} // SWIPDG(...)

#if HAVE_DUNE_GRID_MULTISCALE

template< class G, Stuff::Grid::ChooseLayer gl, class R, int r, int p, GDT::ChooseSpaceBackend s,
          Stuff::LA::ChooseBackend la >
    SWIPDG< G, gl, R, r, p, s, la >::
SWIPDG(const typename SWIPDG< G, gl, R, r, p, s, la >::MsGridProviderType& grid_provider,
       const Stuff::Common::Configuration& bound_inf_cfg,
       const typename SWIPDG< G, gl, R, r, p, s, la >::ProblemType& prob,
       const int level_or_subdomain,
       const std::vector< std::string >& only_these_products)
  : BaseType(SpaceProvider::create(grid_provider, level_or_subdomain),
             SpaceProvider::create(grid_provider, level_or_subdomain),
             bound_inf_cfg,
             prob)
  , beta_(GDT::LocalEvaluation::SIPDG::internal::default_beta(dimDomain))
  , only_these_products_(only_these_products)
{
  // in case of parametric diffusion tensor this discretization is not affinely decomposable any more
  if (this->problem_.diffusion_tensor()->parametric())
    DUNE_THROW(NotImplemented, "The diffusion tensor must not be parametric!");
  if (!this->problem_.diffusion_tensor()->has_affine_part())
    DUNE_THROW(Stuff::Exceptions::wrong_input_given, "The diffusion tensor must not be empty!");
} // SWIPDG(...)

#endif // HAVE_DUNE_GRID_MULTISCALE

template< class G, Stuff::Grid::ChooseLayer gl, class R, int r, int p, GDT::ChooseSpaceBackend s,
          Stuff::LA::ChooseBackend la >
  void SWIPDG< G, gl, R, r, p, s, la >::
init(const bool prune)
{
  if (!this->container_based_initialized_) {
    using namespace GDT;

    auto logger = Stuff::Common::TimedLogger().get("hdd.linearelliptic.discretizations.swipdg.init");

    auto& matrix = *(this->matrix_);
    auto& rhs = *(this->rhs_);
    const auto& space = this->test_space_;
    const auto& boundary_info = *(this->boundary_info_);

    logger.info() << "assembling... " << std::flush;
    Dune::Timer timer;
    SystemAssembler< TestSpaceType > system_assembler(space);
    Stuff::Grid::Functor::DirichletDetector< GridViewType > dirichlet_detector(this->boundary_info());
    system_assembler.add(dirichlet_detector);

    pattern_ = std::make_shared< PatternType >(EllipticOperatorType::pattern(BaseType::test_space(),
                                                                             BaseType::ansatz_space()));

    // lhs operator
    const auto& diffusion_factor = *(this->problem_.diffusion_factor());
    const auto& diffusion_tensor = *(this->problem_.diffusion_tensor());
    assert(!diffusion_tensor.parametric());
    assert(diffusion_tensor.has_affine_part());
    std::vector< std::unique_ptr< EllipticOperatorType > > elliptic_operators;
    for (size_t qq = 0; qq < boost::numeric_cast< size_t >(diffusion_factor.num_components()); ++qq) {
      const size_t id = matrix.register_component(diffusion_factor.coefficient(qq),
                                                  space.mapper().size(), space.mapper().size(), *pattern_);
      elliptic_operators.emplace_back(new EllipticOperatorType(
          *(diffusion_factor.component(qq)),
          *(diffusion_tensor.affine_part()),
          boundary_info,
          *(matrix.component(id)),
          space));
    }
    if (diffusion_factor.has_affine_part()) {
      if (!matrix.has_affine_part())
        matrix.register_affine_part(space.mapper().size(), space.mapper().size(), *pattern_);
      elliptic_operators.emplace_back(new EllipticOperatorType(
          *(diffusion_factor.affine_part()),
          *(diffusion_tensor.affine_part()),
          boundary_info,
          *(matrix.affine_part()),
          space));
    }
    for (auto& elliptic_operator : elliptic_operators)
      system_assembler.add(*elliptic_operator);

    // rhs functional
    // * volume
    typedef typename ProblemType::FunctionType::NonparametricType FunctionType;
    const auto& force = *(this->problem_.force());
    typedef Functionals::L2Volume< FunctionType, VectorType, TestSpaceType > L2VolumeFunctionalType;
    std::vector< std::unique_ptr< L2VolumeFunctionalType > > force_functionals;
    for (size_t qq = 0; qq < boost::numeric_cast< size_t >(force.num_components()); ++qq) {
      const size_t id = rhs.register_component(force.coefficient(qq), space.mapper().size());
      force_functionals.emplace_back(new L2VolumeFunctionalType(*(force.component(qq)),
                                                                *(rhs.component(id)),
                                                                space));
    }
    if (force.has_affine_part()) {
      if (!rhs.has_affine_part())
        rhs.register_affine_part(space.mapper().size());
      force_functionals.emplace_back(new L2VolumeFunctionalType(*(force.affine_part()),
                                                                *(rhs.affine_part()),
                                                                space));
    }
    for (auto& force_functional : force_functionals)
      system_assembler.add(*force_functional);
    // * dirichlet boundary
    const auto& dirichlet = *(this->problem_.dirichlet());
    typedef Functionals::DirichletBoundarySWIPDG< DiffusionFactorType, FunctionType, VectorType, TestSpaceType,
                                                  GridViewType, DiffusionTensorType > DirichletBoundaryFunctionalType;
    std::vector< std::unique_ptr< DirichletBoundaryFunctionalType > > dirichlet_boundary_functionals;
    if (diffusion_factor.has_affine_part() && dirichlet.has_affine_part()) {
      if (!rhs.has_affine_part())
        rhs.register_affine_part(space.mapper().size());
      dirichlet_boundary_functionals.emplace_back(new DirichletBoundaryFunctionalType(
          *(diffusion_factor.affine_part()),
          *(diffusion_tensor.affine_part()),
          *(dirichlet.affine_part()),
          boundary_info,
          *(rhs.affine_part()),
          space));
    }
    if (diffusion_factor.has_affine_part()) {
      for (size_t qq = 0; qq < boost::numeric_cast< size_t >(dirichlet.num_components()); ++qq) {
        const size_t id = rhs.register_component(dirichlet.coefficient(qq), space.mapper().size());
        dirichlet_boundary_functionals.emplace_back(new DirichletBoundaryFunctionalType(
            *(diffusion_factor.affine_part()),
            *(diffusion_tensor.affine_part()),
            *(dirichlet.component(qq)),
            boundary_info,
            *(rhs.component(id)),
            space));
      }
    }
    if (dirichlet.has_affine_part()) {
      for (size_t qq = 0; qq < boost::numeric_cast< size_t >(diffusion_factor.num_components()); ++qq) {
        const size_t id = rhs.register_component(diffusion_factor.coefficient(qq), space.mapper().size());
        dirichlet_boundary_functionals.emplace_back(new DirichletBoundaryFunctionalType(
            *(diffusion_factor.component(qq)),
            *(diffusion_tensor.affine_part()),
            *(dirichlet.affine_part()),
            boundary_info,
            *(rhs.component(id)),
            space));
      }
    }
    Pymor::ParameterType param;
    for (const auto& key : diffusion_factor.parameter_type().keys())
      param.set(key, diffusion_factor.parameter_type().get(key));
    for (const auto& key : dirichlet.parameter_type().keys())
      param.set(key, dirichlet.parameter_type().get(key));
    for (size_t pp = 0; pp < boost::numeric_cast< size_t >(diffusion_factor.num_components()); ++ pp) {
      for (size_t qq = 0; qq < boost::numeric_cast< size_t >(dirichlet.num_components()); ++qq) {
        const std::string expression = "(" + diffusion_factor.coefficient(pp)->expression()
                                           + ")*(" + dirichlet.coefficient(qq)->expression() + ")";
        const size_t id = rhs.register_component(param, expression, space.mapper().size());
        dirichlet_boundary_functionals.emplace_back(new DirichletBoundaryFunctionalType(
            *(diffusion_factor.component(pp)),
            *(diffusion_tensor.affine_part()),
            *(dirichlet.component(qq)),
            boundary_info,
            *(rhs.component(id)),
            space));
      }
    }
    for (auto& dirichlet_boundary_functional : dirichlet_boundary_functionals)
      system_assembler.add(*dirichlet_boundary_functional);

    // * neumann boundary
    const auto& neumann = *(this->problem_.neumann());
    typedef Functionals::L2Face< FunctionType, VectorType, TestSpaceType > L2FaceFunctionalType;
    std::vector< std::unique_ptr< L2FaceFunctionalType > > neumann_boundary_functionals;
    for (size_t qq = 0; qq < boost::numeric_cast< size_t >(neumann.num_components()); ++qq) {
      const size_t id = rhs.register_component(neumann.coefficient(qq), space.mapper().size());
      neumann_boundary_functionals.emplace_back(new L2FaceFunctionalType(
          *(neumann.component(qq)),
          *(rhs.component(id)),
          space,
          new Stuff::Grid::ApplyOn::NeumannIntersections< GridViewType >(boundary_info)));
    }
    if (neumann.has_affine_part()) {
      if (!rhs.has_affine_part())
        rhs.register_affine_part(space.mapper().size());
      neumann_boundary_functionals.emplace_back(new L2FaceFunctionalType(
          *(neumann.affine_part()),
          *(rhs.affine_part()),
          space,
          new Stuff::Grid::ApplyOn::NeumannIntersections< GridViewType >(boundary_info)));
    }
    for (auto& neumann_boundary_functional : neumann_boundary_functionals)
      system_assembler.add(*neumann_boundary_functional);

    // swipdg penalty product
    const size_t over_integrate = 2;
    typedef Products::SwipdgPenaltyAssemblable
        < MatrixType, DiffusionFactorType, DiffusionTensorType, TestSpaceType > PenaltyProductType;
    std::vector< std::unique_ptr< PenaltyProductType > > penalty_products;
    auto penalty_product_matrix = std::make_shared< AffinelyDecomposedMatrixType >();
    if (std::find(only_these_products_.begin(), only_these_products_.end(), "penalty") != only_these_products_.end()) {
      for (DUNE_STUFF_SSIZE_T qq = 0; qq < diffusion_factor.num_components(); ++qq) {
        const auto id = penalty_product_matrix->register_component(diffusion_factor.coefficient(qq),
                                                                   this->test_space().mapper().size(),
                                                                   this->ansatz_space().mapper().size(),
                                                                   *pattern_);
        penalty_products.emplace_back(new PenaltyProductType(*penalty_product_matrix->component(id),
                                                             this->test_space(),
                                                             this->grid_view(),
                                                             this->ansatz_space(),
                                                             *diffusion_factor.component(qq),
                                                             *diffusion_tensor.affine_part(),
                                                             over_integrate));
      }
      if (diffusion_factor.has_affine_part()) {
        penalty_product_matrix->register_affine_part(this->test_space().mapper().size(),
                                                     this->ansatz_space().mapper().size(),
                                                     *pattern_);
        penalty_products.emplace_back(new PenaltyProductType(*penalty_product_matrix->affine_part(),
                                                             this->test_space(),
                                                             this->grid_view(),
                                                             this->ansatz_space(),
                                                             *diffusion_factor.affine_part(),
                                                             *diffusion_tensor.affine_part(),
                                                             over_integrate));
      }
      for (auto& product : penalty_products)
        system_assembler.add(*product);
    }
    // do the actual assembling
    system_assembler.walk();
    logger.info() << "done (took " << timer.elapsed() << "s)" << std::endl;

    logger.info() << "assembling products... " << std::flush;
    timer.reset();
    this->assemble_products(only_these_products_, 2);
    logger.info() << "done (took " << timer.elapsed() << " sec)" << std::endl;

    if (!dirichlet_detector.found())
      this->purely_neumann_ = true;

    if (std::find(only_these_products_.begin(), only_these_products_.end(), "penalty") != only_these_products_.end())
      this->products_.insert(std::make_pair("penalty", penalty_product_matrix));

    this->finalize_init(prune);
  } // if (!this->container_based_initialized_)
} // ... init(...)


} // namespace Discretizations
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_SWIPDG_HXX
