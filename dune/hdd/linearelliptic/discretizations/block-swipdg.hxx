// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_BLOCK_SWIPDG_HXX
#define DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_BLOCK_SWIPDG_HXX

#include "block-swipdg.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Discretizations {


template< class G, class R, int r, int p, Stuff::LA::ChooseBackend la >
    std::string BlockSWIPDG< G, R, r, p, la >::
static_id()
{
  return DiscretizationInterface< Traits >::static_id() + ".block-swipdg";
}

template< class G, class R, int r, int p, Stuff::LA::ChooseBackend la >
    BlockSWIPDG< G, R, r, p, la >::
BlockSWIPDG(const typename BlockSWIPDG< G, R, r, p, la >::GridProviderType& grid_provider,
            const Stuff::Common::Configuration& bound_inf_cfg,
            const typename BlockSWIPDG< G, R, r, p, la >::ProblemType& prob,
            const std::vector< std::string >& only_these_products)
  : LocalDiscretizationsBaseType(grid_provider, prob, only_these_products)
  , BaseType(TestSpaceType(grid_provider.ms_grid(), this->local_test_spaces_),
             AnsatzSpaceType(grid_provider.ms_grid(), this->local_ansatz_spaces_),
             bound_inf_cfg,
             this->zero_boundary_problem_)
  , grid_provider_(grid_provider)
  , ms_grid_(grid_provider.ms_grid())
  , only_these_products_(only_these_products)
  , local_matrices_(ms_grid_->size())
  , local_vectors_(ms_grid_->size())
  , inside_outside_patterns_(ms_grid_->size())
  , outside_inside_patterns_(ms_grid_->size())
  , inside_outside_matrices_(ms_grid_->size())
  , outside_inside_matrices_(ms_grid_->size())
{
  // in case of parametric diffusion tensor everything is too complicated
  if (this->problem_.diffusion_tensor()->parametric())
    DUNE_THROW(NotImplemented, "The diffusion tensor must not be parametric!");
  if (!this->problem_.diffusion_tensor()->has_affine_part())
    DUNE_THROW(Stuff::Exceptions::wrong_input_given, "The diffusion tensor must not be empty!");
} // BlockSWIPDG(...)

template< class G, class R, int r, int p, Stuff::LA::ChooseBackend la >
    const std::vector< std::shared_ptr< typename BlockSWIPDG< G, R, r, p, la >::LocalDiscretizationType > >& BlockSWIPDG< G, R, r, p, la >::
local_discretizations() const
{
  return this->local_discretizations_;
}

template< class G, class R, int r, int p, Stuff::LA::ChooseBackend la >
    void BlockSWIPDG< G, R, r, p, la >::
init(const bool prune)
{
  if (this->container_based_initialized_)
    return;

  auto logger = Stuff::Common::TimedLogger().get("hdd.linearelliptic.discretizations.block-swipdg.init");

  pattern_ = std::make_shared< PatternType >(BaseType::test_space().mapper().size());

  const size_t subdomains = ms_grid_->size();
  logger.info() << "discretizing on " << subdomains << " subdomains..." << std::endl;
  // walk the subdomains for the first time
  //   * to initialize the coupling pattern,
  //   * to finalize the global sparsity pattern
  logger.info() << "  computing patterns and local contributions... " << std::endl;
  for (size_t ss = 0; ss < subdomains; ++ss) {
    // init the local discretizations (assembles matrices and patterns)
    this->local_discretizations_[ss]->init(false);
    // and create the local containers
    // * the matrices
    //   * just copy those from the local discretizations
    const auto local_operator = this->local_discretizations_[ss]->get_operator();
    local_matrices_[ss] = std::make_shared< AffinelyDecomposedMatrixType >();
    //   * we take the affine part only if the diffusion has one, otherwise it contains only the dirichlet rows,
    //     thus it is empty, since the local problems are purely neumann
    if (this->problem().diffusion_factor()->has_affine_part()) {
      if (!local_operator.has_affine_part())
        DUNE_THROW(Stuff::Exceptions::internal_error, "The local operator is missing the affine part!");
      local_matrices_[ss]->register_affine_part(new MatrixType(*(local_operator.affine_part().container())));
    }
    if (local_operator.num_components() < this->problem().diffusion_factor()->num_components())
      DUNE_THROW(Stuff::Exceptions::requirements_not_met,
                 "The local operator should have " << this->problem().diffusion_factor()->num_components()
                 << " components (but has only " << local_operator.num_components() << ")!");
    for (ssize_t qq = 0; qq < this->problem().diffusion_factor()->num_components(); ++qq) {
      local_matrices_[ss]->register_component(new MatrixType(
          local_operator.component(qq).container()->backend()),
          this->problem().diffusion_factor()->coefficient(qq));
    }
    // * and the vectors
    const auto local_functional = this->local_discretizations_[ss]->get_rhs();
    local_vectors_[ss] = std::make_shared< AffinelyDecomposedVectorType >();
    for (size_t qq = 0; qq < boost::numeric_cast< size_t >(local_functional.num_components()); ++qq)
      local_vectors_[ss]->register_component(new VectorType(*(local_functional.component(qq).container())),
                                             new Pymor::ParameterFunctional(local_functional.coefficient(qq)));
    if (local_functional.has_affine_part())
      local_vectors_[ss]->register_affine_part(new VectorType(*(local_functional.affine_part().container())));

    // create and copy the local patterns
    add_local_to_global_pattern(this->local_discretizations_[ss]->pattern(), ss, ss, *pattern_);
    const auto& inner_test_space = this->local_discretizations_[ss]->test_space();
    const auto& inner_ansatz_space = this->local_discretizations_[ss]->ansatz_space();
    // walk the neighbors
    for (const size_t& nn : ms_grid_->neighborsOf(ss)) {
      // visit each coupling only once (assemble primally)
      if (ss < nn) {
        const auto& outer_test_space = this->local_discretizations_[nn]->test_space();
        const auto& outer_ansatz_space = this->local_discretizations_[nn]->ansatz_space();
        const auto inside_outside_grid_part = ms_grid_->couplingGridPart(ss, nn);
        const auto outside_inside_grid_part = ms_grid_->couplingGridPart(nn, ss);
        // create the coupling patterns
        auto inside_outside_pattern = std::make_shared< PatternType >(
              inner_test_space.compute_face_pattern(inside_outside_grid_part, outer_ansatz_space));
        inside_outside_patterns_[ss].insert(std::make_pair(nn, inside_outside_pattern));
        auto outside_inside_pattern = std::make_shared< PatternType >(
              outer_test_space.compute_face_pattern(outside_inside_grid_part, inner_ansatz_space));
        outside_inside_patterns_[nn].insert(std::make_pair(ss, outside_inside_pattern));
        // and copy them
        add_local_to_global_pattern(*inside_outside_pattern,  ss, nn, *pattern_);
        add_local_to_global_pattern(*outside_inside_pattern,  nn, ss, *pattern_);
      } // visit each coupling only once (assemble primaly)
    } // walk the neighbors
  } // walk the subdomains for the first time

  // walk the subdomains for the second time
  //   * to assemble the boundary matrices and vectors and
  //   * to assemble the coupling matrices
  logger.info() << "computing coupling and boundary contributions... " << std::endl;
  for (size_t ss = 0; ss < ms_grid_->size(); ++ss) {
    const auto& inner_test_mapper = this->local_discretizations_[ss]->test_space().mapper();
    const auto& inner_ansatz_mapper = this->local_discretizations_[ss]->ansatz_space().mapper();
    if (ms_grid_->boundary(ss))
      assemble_boundary_contributions(ss);
    // walk the neighbors
    for (const size_t& nn : ms_grid_->neighborsOf(ss)) {
      // visit each coupling only once (assemble primaly)
      if (ss < nn) {
        const auto& outer_test_mapper = this->local_discretizations_[nn]->test_space().mapper();
        const auto& outer_ansatz_mapper = this->local_discretizations_[nn]->ansatz_space().mapper();
        // get the patterns
        const auto in_out_result = inside_outside_patterns_[ss].find(nn);
        if (in_out_result == inside_outside_patterns_[ss].end())
          DUNE_THROW(Stuff::Exceptions::internal_error, "subdomain " << ss << ", neighbour " << nn);
        const auto& inside_outside_pattern = *(in_out_result->second);
        const auto out_in_result = outside_inside_patterns_[nn].find(ss);
        if (out_in_result == outside_inside_patterns_[nn].end())
          DUNE_THROW(Stuff::Exceptions::internal_error, "subdomain " << ss << ", neighbour " << nn);
        const auto& outside_inside_pattern = *(out_in_result->second);
        // create the coupling matrices
        auto inside_outside_matrix = std::make_shared< AffinelyDecomposedMatrixType >();
        auto outside_inside_matrix = std::make_shared< AffinelyDecomposedMatrixType >();
        if (this->problem().diffusion_factor()->has_affine_part()) {
          inside_outside_matrix->register_affine_part(new MatrixType(inner_test_mapper.size(),
                                                                     outer_ansatz_mapper.size(),
                                                                     inside_outside_pattern));
          outside_inside_matrix->register_affine_part(new MatrixType(outer_test_mapper.size(),
                                                                     inner_ansatz_mapper.size(),
                                                                     outside_inside_pattern));
        }
        for (ssize_t qq = 0; qq < this->problem().diffusion_factor()->num_components(); ++qq) {
          inside_outside_matrix->register_component(new MatrixType(inner_test_mapper.size(),
                                                                   outer_ansatz_mapper.size(),
                                                                   inside_outside_pattern),
                                                    this->problem().diffusion_factor()->coefficient(qq));
          outside_inside_matrix->register_component(new MatrixType(outer_test_mapper.size(),
                                                                   inner_ansatz_mapper.size(),
                                                                   outside_inside_pattern),
                                                    this->problem().diffusion_factor()->coefficient(qq));
        }
        // and assemble them
        assemble_coupling_contributions(ss, nn,
                                        *(local_matrices_[ss]),
                                        *(inside_outside_matrix),
                                        *(outside_inside_matrix),
                                        *(local_matrices_[nn]));
        inside_outside_matrices_[ss].insert(std::make_pair(nn, inside_outside_matrix));
        outside_inside_matrices_[nn].insert(std::make_pair(ss, outside_inside_matrix));
      } // visit each coupling only once
    } // walk the neighbors
  } // walk the subdomains for the second time

  // build global containers
  pattern_->sort();
  build_global_containers();

  logger.info() << "assembling products... " << std::endl;
  this->assemble_products(only_these_products_, 2);

  // finalize
  this->finalize_init(prune);

  logger.info() << "finished!" << std::endl;
} // ... init(...)

template< class G, class R, int r, int p, Stuff::LA::ChooseBackend la >
    ssize_t BlockSWIPDG< G, R, r, p, la >::
num_subdomains() const
{
  return ms_grid_->size();
}

template< class G, class R, int r, int p, Stuff::LA::ChooseBackend la >
    std::vector< ssize_t > BlockSWIPDG< G, R, r, p, la >::
neighbouring_subdomains(const ssize_t ss) const
{
  if (ss < 0 || ss >= num_subdomains())
    DUNE_THROW(Stuff::Exceptions::index_out_of_range,
               "0 <= ss < num_subdomains() = " << num_subdomains() << " is not true for ss = " << ss << "!");
  const auto set_of_neighbours = ms_grid_->neighborsOf(ss);
  return std::vector< ssize_t >(set_of_neighbours.begin(), set_of_neighbours.end());
}

template< class G, class R, int r, int p, Stuff::LA::ChooseBackend la >
    typename BlockSWIPDG< G, R, r, p, la >::VectorType BlockSWIPDG< G, R, r, p, la >::
localize_vector(const typename BlockSWIPDG< G, R, r, p, la >::VectorType& global_vector, const size_t ss) const
{
  if ((std::make_signed< size_t >::type)(ss) >= num_subdomains())
    DUNE_THROW(Stuff::Exceptions::index_out_of_range,
               "0 <= ss < num_subdomains() = " << num_subdomains() << " is not true for ss = " << ss << "!");
  if (global_vector.size() != this->ansatz_space().mapper().size())
    DUNE_THROW(Stuff::Exceptions::index_out_of_range,
               "The size() of global_vector (" << global_vector.dim()
               << ") does not match the size() of the ansatz space (" << this->ansatz_space().mapper().size() << ")!");
  assert(ss < this->local_discretizations_.size());
  VectorType local_vector = this->local_discretizations_[ss]->create_vector();
  for (size_t ii = 0; ii < local_vector.size(); ++ii)
    local_vector.set_entry(ii, global_vector.get_entry(this->ansatz_space().mapper().mapToGlobal(ss, ii)));
  return local_vector;
} // ... localize_vetor(...)

template< class G, class R, int r, int p, Stuff::LA::ChooseBackend la >
    typename BlockSWIPDG< G, R, r, p, la >::VectorType BlockSWIPDG< G, R, r, p, la >::
globalize_vectors(const std::vector< typename BlockSWIPDG< G, R, r, p, la >::VectorType >& local_vectors) const
{
  if (local_vectors.size() != boost::numeric_cast< size_t >(num_subdomains()))
    DUNE_THROW(Stuff::Exceptions::wrong_input_given,
               "Given local_vectors has wrong size (is " << local_vectors.size() << ", should be "
               << num_subdomains() << ")!");
  VectorType ret(this->ansatz_space().mapper().size());
  for (size_t ss = 0; ss < boost::numeric_cast< size_t >(num_subdomains()); ++ss) {
    const auto& local_vector = local_vectors[ss];
    if (local_vector.size() != this->local_discretizations_[ss]->ansatz_space().mapper().size())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given,
                 "Given local_vectors[" << ss << "] has wrong size (is "
                 << local_vector.size() << ", should be "
                 << this->local_discretizations_[ss]->ansatz_space().mapper().size() << ")!");
    copy_local_to_global_vector(local_vector, ss, ret);
  }
  return ret;
}

template< class G, class R, int r, int p, Stuff::LA::ChooseBackend la >
    typename BlockSWIPDG< G, R, r, p, la >::VectorType* BlockSWIPDG< G, R, r, p, la >::
globalize_vectors_and_return_ptr(const std::vector< typename BlockSWIPDG< G, R, r, p, la >::VectorType >& local_vectors) const
{
  return new VectorType(globalize_vectors(local_vectors));
}

template< class G, class R, int r, int p, Stuff::LA::ChooseBackend la >
    typename BlockSWIPDG< G, R, r, p, la >::VectorType* BlockSWIPDG< G, R, r, p, la >::
localize_vector_and_return_ptr(const typename BlockSWIPDG< G, R, r, p, la >::VectorType& global_vector,
                               const ssize_t ss) const
{
  return new VectorType(localize_vector(global_vector, boost::numeric_cast< size_t >(ss)));
}

template< class G, class R, int r, int p, Stuff::LA::ChooseBackend la >
    typename BlockSWIPDG< G, R, r, p, la >::ProductType BlockSWIPDG< G, R, r, p, la >::
get_local_product(const size_t ss, const std::string id) const
{
  if (boost::numeric_cast< ssize_t >(ss) >= num_subdomains())
    DUNE_THROW(Stuff::Exceptions::index_out_of_range,
               "0 <= ss < num_subdomains() = " << num_subdomains() << " is not true for ss = " << ss << "!");
  return this->local_discretizations_[ss]->get_product(id);
}

template< class G, class R, int r, int p, Stuff::LA::ChooseBackend la >
    typename BlockSWIPDG< G, R, r, p, la >::ProductType* BlockSWIPDG< G, R, r, p, la >::
get_local_product_and_return_ptr(const ssize_t ss, const std::string id) const
{
  return new ProductType(get_local_product(boost::numeric_cast< size_t >(ss), id));
}

template< class G, class R, int r, int p, Stuff::LA::ChooseBackend la >
    typename BlockSWIPDG< G, R, r, p, la >::OperatorType BlockSWIPDG< G, R, r, p, la >::
get_local_operator(const size_t ss) const
{
  if (boost::numeric_cast< ssize_t >(ss) >= num_subdomains())
    DUNE_THROW(Stuff::Exceptions::index_out_of_range,
               "0 <= ss < num_subdomains() = " << num_subdomains() << " is not true for ss = " << ss << "!");
  assert(ss < local_matrices_.size());
  return OperatorType(*(local_matrices_[ss]));
}

template< class G, class R, int r, int p, Stuff::LA::ChooseBackend la >
    typename BlockSWIPDG< G, R, r, p, la >::OperatorType* BlockSWIPDG< G, R, r, p, la >::
get_local_operator_and_return_ptr(const ssize_t ss) const
{
  return new OperatorType(get_local_operator(boost::numeric_cast< size_t >(ss)));
}

template< class G, class R, int r, int p, Stuff::LA::ChooseBackend la >
    typename BlockSWIPDG< G, R, r, p, la >::OperatorType BlockSWIPDG< G, R, r, p, la >::
get_coupling_operator(const size_t ss, const size_t nn) const
{
  if (ss >= boost::numeric_cast< size_t >(num_subdomains()))
    DUNE_THROW(Stuff::Exceptions::index_out_of_range,
               "0 <= ss < num_subdomains() = " << num_subdomains() << " is not true for ss = " << ss << "!");
  const auto neighbours = ms_grid_->neighborsOf(ss);
  if (neighbours.count(nn) == 0)
    DUNE_THROW(Stuff::Exceptions::index_out_of_range,
               "Subdomain " << nn << " is not a neighbour of subdomain " << ss
               << " (call neighbouring_subdomains(" << ss << ") to find out)!");
  if (ss < nn) {
    // we need to look for this coupling operator in the inside/outside context
    const auto result_inside_outside_matrix = inside_outside_matrices_[ss].find(nn);
    if (result_inside_outside_matrix == inside_outside_matrices_[ss].end())
      DUNE_THROW(Stuff::Exceptions::internal_error,
                 "The coupling matrix for subdomain " << ss << " and neighbour " << nn << " is missing!");
    const auto inside_outside_matrix = result_inside_outside_matrix->second;
    return OperatorType(*inside_outside_matrix);
  } else if (nn < ss) {
    // we need to look for this coupling operator in the outside/inside context
    const auto result_outside_inside_matrix = outside_inside_matrices_[ss].find(nn);
    if (result_outside_inside_matrix == outside_inside_matrices_[ss].end())
      DUNE_THROW(Stuff::Exceptions::internal_error,
                 "The coupling matrix for neighbour " << nn << " and subdomain " << ss << " is missing!");
    const auto outside_inside_matrix = result_outside_inside_matrix->second;
    return OperatorType(*outside_inside_matrix);
  } else {
    // the above exception should have cought this
    DUNE_THROW(Stuff::Exceptions::internal_error,
               "The multiscale grid is corrupted! Subdomain " << ss << " must not be its own neighbour!");
  }
} // ... get_coupling_operator(...)

template< class G, class R, int r, int p, Stuff::LA::ChooseBackend la >
    typename BlockSWIPDG< G, R, r, p, la >::OperatorType* BlockSWIPDG< G, R, r, p, la >::
get_coupling_operator_and_return_ptr(const ssize_t ss, const ssize_t nn) const
{
  return new OperatorType(get_coupling_operator(boost::numeric_cast< size_t >(ss),
                                                boost::numeric_cast< size_t >(nn)));
}

template< class G, class R, int r, int p, Stuff::LA::ChooseBackend la >
    typename BlockSWIPDG< G, R, r, p, la >::FunctionalType BlockSWIPDG< G, R, r, p, la >::
get_local_functional(const size_t ss) const
{
  if (ss >= boost::numeric_cast< size_t >(num_subdomains()))
    DUNE_THROW(Stuff::Exceptions::index_out_of_range,
               "0 <= ss < num_subdomains() = " << num_subdomains() << " is not true for ss = " << ss << "!");
  assert(ss < local_vectors_.size());
  return FunctionalType(*(local_vectors_[ss]));
}

template< class G, class R, int r, int p, Stuff::LA::ChooseBackend la >
    typename BlockSWIPDG< G, R, r, p, la >::FunctionalType* BlockSWIPDG< G, R, r, p, la >::
get_local_functional_and_return_ptr(const ssize_t ss) const
{
  return new FunctionalType(get_local_functional(boost::numeric_cast< size_t >(ss)));
}

template< class G, class R, int r, int p, Stuff::LA::ChooseBackend la >
    typename BlockSWIPDG< G, R, r, p, la >::VectorType BlockSWIPDG< G, R, r, p, la >::
solve_for_local_correction(const std::vector< typename BlockSWIPDG< G, R, r, p, la >::VectorType >& local_vectors,
                           const size_t subdomain,
                           const Pymor::Parameter mu) const
{
  DUNE_THROW(Stuff::Exceptions::internal_error, "Do not call this method, I do not trust it!");
  using namespace GDT;

  if (mu.type() != this->parameter_type())
    DUNE_THROW(Pymor::Exceptions::wrong_parameter_type,
               "mu is " << mu.type() << ", should be " << this->parameter_type() << "!");
  if (local_vectors.size() != ms_grid_->size())
    DUNE_THROW(Stuff::Exceptions::wrong_input_given,
               "local_vectors is of size " << local_vectors.size() << " and should be of size " << ms_grid_->size());
  VectorType vector(this->ansatz_space().mapper().size());
  for (size_t ss = 0; ss < ms_grid_->size(); ++ss) {
    if (local_vectors[ss].size() != this->local_discretizations_[ss]->ansatz_space().mapper().size())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given,
                 "local_vectors[" << ss << "] is of size " << local_vectors[ss].size() << " and should be of size "
                 << this->local_discretizations_[ss]->ansatz_space().mapper().size());
    if (!local_vectors[ss].valid())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given, "local_vectors[" << ss << "] contains NaN or INF!");
    copy_local_to_global_vector(local_vectors[ss], ss, vector);
  }

//    const std::string prefix = "subdomain_" + DSC::toString(subdomain) + "_";

  const ConstDiscreteFunction< AnsatzSpaceType, VectorType > current_global_solution(this->ansatz_space(), vector);
//    current_global_solution.visualize(prefix + "current_solution_global");

  typedef SWIPDG< typename MsGridType::GridType, Stuff::Grid::ChooseLayer::local_oversampled, RangeFieldType, dimRange
                , 1, GDT::ChooseSpaceBackend::fem, la > OversampledDiscretizationType;
  OversampledDiscretizationType oversampled_discretization(grid_provider_,
                                                           this->multiscale_boundary_config_,
                                                           this->problem(),
                                                           boost::numeric_cast< int >(subdomain));
  oversampled_discretization.init();

  DiscreteFunction< typename OversampledDiscretizationType::AnsatzSpaceType, VectorType >
      current_oversampled_solution(oversampled_discretization.ansatz_space());
  const Operators::Projection< typename OversampledDiscretizationType::GridViewType >
      oversampled_projection_operator(oversampled_discretization.grid_view());
  oversampled_projection_operator.apply(current_global_solution, current_oversampled_solution);
//    current_oversampled_solution.visualize(prefix + "current_solution_oversampled");

  if (!oversampled_discretization.rhs()->has_affine_part())
    oversampled_discretization.rhs()->register_affine_part(oversampled_discretization.test_space().mapper().size());
  if (oversampled_discretization.system_matrix()->parametric()) {
    const auto oversampled_system_matrix = oversampled_discretization.system_matrix()->freeze_parameter(mu);
    *(oversampled_discretization.rhs()->affine_part())
        -= oversampled_system_matrix * current_oversampled_solution.vector();
  } else {
    const auto& oversampled_system_matrix = *(oversampled_discretization.system_matrix()->affine_part());
    *(oversampled_discretization.rhs()->affine_part())
        -= oversampled_system_matrix * current_oversampled_solution.vector();
  }

  oversampled_discretization.solve(current_oversampled_solution.vector(), mu);
//    current_oversampled_solution.visualize(prefix + "correction_oversampled");

  DiscreteFunction< typename LocalDiscretizationType::AnsatzSpaceType, VectorType >
      local_solution(this->local_discretizations_[subdomain]->ansatz_space());
  const Operators::Projection< typename LocalDiscretizationType::GridViewType >
      local_projection_operator(this->local_discretizations_[subdomain]->grid_view());
  local_projection_operator.apply(current_oversampled_solution, local_solution);
//    local_solution.visualize(prefix + "correction_local");

  return local_solution.vector();
} // ... solve_for_local_correction(...)

template< class G, class R, int r, int p, Stuff::LA::ChooseBackend la >
    typename BlockSWIPDG< G, R, r, p, la >::LocalDiscretizationType BlockSWIPDG< G, R, r, p, la >::
get_local_discretization(const size_t subdomain) const
{
  if (subdomain >= this->grid_provider_.num_subdomains())
    DUNE_THROW(Stuff::Exceptions::index_out_of_range,
               "Given subdomain " << subdomain << " too large (has to be smaller than "
               << this->grid_provider_.num_subdomains() << "!");
  return *(this->local_discretizations_[subdomain]);
} // ... get_local_discretization(...)

template< class G, class R, int r, int p, Stuff::LA::ChooseBackend la >
    typename BlockSWIPDG< G, R, r, p, la >::LocalDiscretizationType* BlockSWIPDG< G, R, r, p, la >::
pb_get_local_discretization(const ssize_t subdomain) const
{
  size_t ss = std::numeric_limits< size_t >::max();
  try {
    ss = boost::numeric_cast< size_t >(subdomain);
  } catch (boost::bad_numeric_cast& ee) {
    DUNE_THROW(Stuff::Exceptions::index_out_of_range,
               "There was an error in boost converting " << subdomain << " to "
               << Stuff::Common::Typename< size_t >::value() << ":\n\n" << ee.what());
  }
  return new LocalDiscretizationType(get_local_discretization(ss));
} // ... pb_get_local_discretization(...)

template< class G, class R, int r, int p, Stuff::LA::ChooseBackend la >
    typename BlockSWIPDG< G, R, r, p, la >::OversampledDiscretizationType BlockSWIPDG< G, R, r, p, la >::
get_oversampled_discretization(const size_t subdomain, const std::string boundary_value_type) const
{
  if (subdomain >= this->grid_provider_.num_subdomains())
    DUNE_THROW(Stuff::Exceptions::index_out_of_range,
               "Given subdomain " << subdomain << " too large (has to be smaller than "
               << this->grid_provider_.num_subdomains() << "!");
  if (boundary_value_type == "dirichlet") {
    if (!this->oversampled_discretizations_dirichlet_[subdomain]) {
      this->oversampled_discretizations_dirichlet_[subdomain] = std::make_shared< OversampledDiscretizationType >(
                                                                  grid_provider_,
                                                                  this->all_dirichlet_boundary_config_,
                                                                  this->zero_boundary_problem_,
                                                                  subdomain,
                                                                  only_these_products_);
      this->oversampled_discretizations_dirichlet_[subdomain]->init();
    }
    return *(this->oversampled_discretizations_dirichlet_[subdomain]);
  } else if (boundary_value_type == "neumann") {
    if (!this->oversampled_discretizations_neumann_[subdomain]) {
      this->oversampled_discretizations_neumann_[subdomain] = std::make_shared< OversampledDiscretizationType >(
                                                                grid_provider_,
                                                                this->all_neumann_boundary_config_,
                                                                this->zero_boundary_problem_,
                                                                subdomain,
                                                                only_these_products_);
      this->oversampled_discretizations_neumann_[subdomain]->init();
    }
    return *(this->oversampled_discretizations_neumann_[subdomain]);
  } else {
    DUNE_THROW(Stuff::Exceptions::wrong_input_given,
               "Unknown boundary_value_type given (has to be dirichlet or neumann): " << boundary_value_type);
    return *(this->oversampled_discretizations_dirichlet_[subdomain]);
  }
} // ... get_oversampled_discretization(...)

template< class G, class R, int r, int p, Stuff::LA::ChooseBackend la >
    typename BlockSWIPDG< G, R, r, p, la >::OversampledDiscretizationType* BlockSWIPDG< G, R, r, p, la >::
pb_get_oversampled_discretization(const ssize_t subdomain, const std::string boundary_value_type) const
{
  size_t ss = std::numeric_limits< size_t >::max();
  try {
    ss = boost::numeric_cast< size_t >(subdomain);
  } catch (boost::bad_numeric_cast& ee) {
    DUNE_THROW(Stuff::Exceptions::index_out_of_range,
               "There was an error in boost converting " << subdomain << " to "
               << Stuff::Common::Typename< size_t >::value() << ":\n\n" << ee.what());
  }
  return new OversampledDiscretizationType(get_oversampled_discretization(ss, boundary_value_type));
} // ... pb_get_local_discretization(...)

//template< class G, class R, int r, int p, Stuff::LA::ChooseBackend la >
//    typename BlockSWIPDG< G, R, r, p, la >::OversampledDiscretizationType* BlockSWIPDG< G, R, r, p, la >::
//pb_get_oversampled_discretization(const ssize_t subdomain,
//                                  const std::string boundary_value_type,
//                                  const typename BlockSWIPDG< G, R, r, p, la >::VectorType& boundary_values) const
//{
//  size_t ss = std::numeric_limits< size_t >::max();
//  try {
//    ss = boost::numeric_cast< size_t >(subdomain);
//  } catch (boost::bad_numeric_cast& ee) {
//    DUNE_THROW(Stuff::Exceptions::index_out_of_range,
//               "There was an error in boost converting " << subdomain << " to "
//               << Stuff::Common::Typename< size_t >::value() << ":\n\n" << ee.what());
//  }
//  return new OversampledDiscretizationType(get_oversampled_discretization(ss, boundary_value_type, boundary_values));
//} // ... pb_get_local_discretization(...)

template< class G, class R, int r, int p, Stuff::LA::ChooseBackend la >
    void BlockSWIPDG< G, R, r, p, la >::
add_local_to_global_pattern(const BlockSWIPDG< G, R, r, p, la >::PatternType& local,
                            const size_t test_subdomain,
                            const size_t ansatz_subdomain,
                            typename BlockSWIPDG< G, R, r, p, la >::PatternType& global) const
{
  for (size_t local_ii = 0; local_ii < local.size(); ++local_ii) {
    const size_t global_ii = this->test_space().mapper().mapToGlobal(test_subdomain, local_ii);
    const auto& local_rows = local.inner(local_ii);
    for (const auto& local_jj : local_rows) {
      const size_t global_jj = this->ansatz_space().mapper().mapToGlobal(ansatz_subdomain, local_jj);
      global.insert(global_ii, global_jj);
    }
  }
} // ... add_local_to_global_pattern(...)

template< class G, class R, int r, int p, Stuff::LA::ChooseBackend la >
    void BlockSWIPDG< G, R, r, p, la >::
copy_local_to_global_matrix(const typename BlockSWIPDG< G, R, r, p, la >::AffinelyDecomposedConstMatrixType& local_matrix,
                            const typename BlockSWIPDG< G, R, r, p, la >::PatternType& local_pattern,
                            const size_t subdomain,
                            const size_t neighbor,
                            BlockSWIPDG< G, R, r, p, la >::AffinelyDecomposedMatrixType& global_matrix) const
{
  for (size_t qq = 0; qq < boost::numeric_cast< size_t >(local_matrix.num_components()); ++qq) {
    auto comp = find_add_component(global_matrix,
                                   *local_matrix.coefficient(qq),
                                   this->test_space(),
                                   this->ansatz_space(),
                                   *pattern_);
    copy_local_to_global_matrix(*(local_matrix.component(qq)),
                                local_pattern,
                                subdomain,
                                neighbor,
                                *(global_matrix.component(comp)));
  }
  if (local_matrix.has_affine_part()) {
    if (!global_matrix.has_affine_part())
      global_matrix.register_affine_part(this->test_space_.mapper().size(),
                                         this->ansatz_space_.mapper().size(),
                                         *pattern_);
    copy_local_to_global_matrix(*(local_matrix.affine_part()),
                                local_pattern,
                                subdomain,
                                neighbor,
                                *(global_matrix.affine_part()));
  }
} // copy_local_to_global_matrix(...)

template< class G, class R, int r, int p, Stuff::LA::ChooseBackend la >
    void BlockSWIPDG< G, R, r, p, la >::
copy_local_to_global_vector(const typename BlockSWIPDG< G, R, r, p, la >::AffinelyDecomposedConstVectorType& local_vector,
                            const size_t subdomain,
                            typename BlockSWIPDG< G, R, r, p, la >::AffinelyDecomposedVectorType& global_vector) const
{
  for (size_t qq = 0; qq < boost::numeric_cast< size_t >(local_vector.num_components()); ++qq) {
    auto comp = find_add_component(global_vector, *local_vector.coefficient(qq), this->test_space());
    copy_local_to_global_vector(*(local_vector.component(qq)),
                                subdomain,
                                *(global_vector.component(comp)));
  }
  if (local_vector.has_affine_part()) {
    if (!global_vector.has_affine_part())
      global_vector.register_affine_part(this->test_space().mapper().size());
    copy_local_to_global_vector(*(local_vector.affine_part()),
                                subdomain,
                                *(global_vector.affine_part()));
  }
} // copy_local_to_global_vector(...)

template< class G, class R, int r, int p, Stuff::LA::ChooseBackend la >
    void BlockSWIPDG< G, R, r, p, la >::
assemble_boundary_contributions(const size_t subdomain) const
{
  typedef typename MsGridType::BoundaryGridPartType BoundaryGridPartType;
  typedef typename LocalDiscretizationType::TestSpaceType   LocalTestSpaceType;
  typedef typename LocalDiscretizationType::AnsatzSpaceType LocalAnsatzSpaceType;
  const LocalTestSpaceType&   local_test_space   = this->local_discretizations_[subdomain]->test_space();
  const LocalAnsatzSpaceType& local_ansatz_space = this->local_discretizations_[subdomain]->ansatz_space();
  const auto& local_pattern = this->local_discretizations_[subdomain]->pattern();
  typedef GDT::SystemAssembler< LocalTestSpaceType, BoundaryGridPartType, LocalAnsatzSpaceType > BoundaryAssemblerType;
  BoundaryAssemblerType boundary_assembler(local_test_space,
                                           local_ansatz_space,
                                           ms_grid_->boundaryGridPart(subdomain));

  auto& local_matrix = *(local_matrices_[subdomain]);
  auto& local_vector = *(local_vectors_[subdomain]);

  // lhs
  // * dirichlet boundary terms
  typedef typename ProblemType::DiffusionFactorType::NonparametricType DiffusionFactorType;
  typedef typename ProblemType::DiffusionTensorType::NonparametricType DiffusionTensorType;
  const auto& diffusion_tensor = *(this->problem().diffusion_tensor());
  assert(!diffusion_tensor.parametric());
  assert(diffusion_tensor.has_affine_part());
  typedef GDT::LocalOperator::Codim1BoundaryIntegral< GDT::LocalEvaluation::SWIPDG::BoundaryLHS< DiffusionFactorType, DiffusionTensorType > >
      DirichletOperatorType;
  typedef GDT::LocalAssembler::Codim1BoundaryMatrix< DirichletOperatorType > DirichletMatrixAssemblerType;
  std::vector< std::unique_ptr< DirichletOperatorType > > dirichlet_operators;
  std::vector< std::unique_ptr< DirichletMatrixAssemblerType > > dirichlet_matrix_assemblers;
  for (ssize_t qq = 0; qq < this->problem().diffusion_factor()->num_components(); ++qq) {
    auto mat_comp = find_add_component(local_matrix,
                                       *this->problem().diffusion_factor()->coefficient(qq),
                                       local_test_space,
                                       local_ansatz_space,
                                       local_pattern);
    dirichlet_operators.emplace_back(new DirichletOperatorType(*(this->problem().diffusion_factor()->component(qq)),
                                                               *(diffusion_tensor.affine_part())));
    dirichlet_matrix_assemblers.emplace_back(new DirichletMatrixAssemblerType(*dirichlet_operators.back()));
    boundary_assembler.add(*dirichlet_matrix_assemblers.back(),
                           *(local_matrix.component(mat_comp)),
                           new Stuff::Grid::ApplyOn::DirichletIntersections< BoundaryGridPartType >(this->boundary_info()));
  }
  if (this->problem().diffusion_factor()->has_affine_part()) {
    if (!local_matrix.has_affine_part())
      local_matrix.register_affine_part(local_test_space.mapper().size(),
                                        local_ansatz_space.mapper().size(),
                                        local_pattern);
    dirichlet_operators.emplace_back(new DirichletOperatorType(*(this->problem().diffusion_factor()->affine_part()),
                                                               *(diffusion_tensor.affine_part())));
    dirichlet_matrix_assemblers.emplace_back(new DirichletMatrixAssemblerType(*dirichlet_operators.back()));
    boundary_assembler.add(*dirichlet_matrix_assemblers.back(),
                           *(local_matrix.affine_part()),
                           new Stuff::Grid::ApplyOn::DirichletIntersections< BoundaryGridPartType >(this->boundary_info()));
  }

  // rhs
  // * neumann boundary terms
  typedef typename ProblemType::FunctionType::NonparametricType NeumannType;
  typedef GDT::LocalFunctional::Codim1Integral< GDT::LocalEvaluation::Product< NeumannType > > NeumannFunctionalType;
  typedef GDT::LocalAssembler::Codim1Vector< NeumannFunctionalType > NeumannVectorAssemblerType;
  std::vector< std::unique_ptr< NeumannFunctionalType > > neumann_functionals;
  std::vector< std::unique_ptr< NeumannVectorAssemblerType > > neumann_vector_assemblers;
  for (ssize_t qq = 0; qq < this->problem().neumann()->num_components(); ++qq) {
    auto vec_comp = find_add_component(local_vector,
                                       *this->problem().neumann()->coefficient(qq),
                                       local_test_space);
    neumann_functionals.emplace_back(new NeumannFunctionalType(*(this->problem().neumann()->component(qq))));
    neumann_vector_assemblers.emplace_back(new NeumannVectorAssemblerType(*(neumann_functionals.back())));
    boundary_assembler.add(*(neumann_vector_assemblers.back()),
                           *(local_vector.component(vec_comp)),
                           new Stuff::Grid::ApplyOn::NeumannIntersections< BoundaryGridPartType >(this->boundary_info()));
  }
  if (this->problem().neumann()->has_affine_part()) {
    if (!local_vector.has_affine_part())
      local_vector.register_affine_part(local_test_space.mapper().size());
    neumann_functionals.emplace_back(new NeumannFunctionalType(*(this->problem().neumann()->affine_part())));
    neumann_vector_assemblers.emplace_back(new NeumannVectorAssemblerType(*(neumann_functionals.back())));
    boundary_assembler.add(*(neumann_vector_assemblers.back()),
                           *(local_vector.affine_part()),
                           new Stuff::Grid::ApplyOn::NeumannIntersections< BoundaryGridPartType >(this->boundary_info()));
  }

  // * dirichlet boundary terms
  typedef typename ProblemType::FunctionType::NonparametricType DirichletType;
  typedef GDT::LocalFunctional::Codim1Integral< GDT::LocalEvaluation::SWIPDG::BoundaryRHS< DiffusionFactorType,
                                                                                           DirichletType,
                                                                                           DiffusionTensorType > >
      DirichletFunctionalType;
  typedef GDT::LocalAssembler::Codim1Vector< DirichletFunctionalType > DirichletVectorAssemblerType;
  std::vector< std::unique_ptr< DirichletFunctionalType > > dirichlet_functionals;
  std::vector< std::unique_ptr< DirichletVectorAssemblerType > > dirichlet_vector_assemblers;
  if (this->problem().diffusion_factor()->has_affine_part()) {
    for (ssize_t qq = 0; qq < this->problem().dirichlet()->num_components(); ++qq) {
      auto vec_comp = find_add_component(local_vector,
                                         *this->problem().dirichlet()->coefficient(qq),
                                         local_test_space);
      dirichlet_functionals.emplace_back(new DirichletFunctionalType(*(this->problem().diffusion_factor()->affine_part()),
                                                                     *(diffusion_tensor.affine_part()),
                                                                     *(this->problem().dirichlet()->component(qq))));
      dirichlet_vector_assemblers.emplace_back(new DirichletVectorAssemblerType(*dirichlet_functionals.back()));
      boundary_assembler.add(*dirichlet_vector_assemblers.back(),
                             *local_vector.component(vec_comp),
                             new Stuff::Grid::ApplyOn::DirichletIntersections< BoundaryGridPartType >(this->boundary_info()));
    }
  }
  if (this->problem().dirichlet()->has_affine_part()) {
    for (ssize_t qq = 0; qq < this->problem().diffusion_factor()->num_components(); ++qq) {
      auto vec_comp = find_add_component(local_vector,
                                         *this->problem().diffusion_factor()->coefficient(qq),
                                         local_test_space);
      dirichlet_functionals.emplace_back(new DirichletFunctionalType(*(this->problem().diffusion_factor()->component(qq)),
                                                                     *(diffusion_tensor.affine_part()),
                                                                     *(this->problem().dirichlet()->affine_part())));
      dirichlet_vector_assemblers.emplace_back(new DirichletVectorAssemblerType(*dirichlet_functionals.back()));
      boundary_assembler.add(*dirichlet_vector_assemblers.back(),
                             *local_vector.component(vec_comp),
                             new Stuff::Grid::ApplyOn::DirichletIntersections< BoundaryGridPartType >(this->boundary_info()));
    }
  }
  Pymor::ParameterType combined_param_type;
  for (auto key : this->problem().diffusion_factor()->parameter_type().keys())
    combined_param_type.set(key, this->problem().diffusion_factor()->parameter_type().get(key));
  for (auto key : this->problem().dirichlet()->parameter_type().keys())
    combined_param_type.set(key, this->problem().dirichlet()->parameter_type().get(key));
  for (ssize_t pp = 0; pp < this->problem().diffusion_factor()->num_components(); ++ pp) {
    for (ssize_t qq = 0; qq < this->problem().dirichlet()->num_components(); ++qq) {
      const Pymor::ParameterFunctional coefficient(combined_param_type,
                                                   "("
                                                   + this->problem().diffusion_factor()->coefficient(pp)->expression()
                                                   + ")*(" + this->problem().dirichlet()->coefficient(qq)->expression()
                                                   + ")");
      auto vec_comp = find_add_component(local_vector, coefficient, local_test_space);
      dirichlet_functionals.emplace_back(new DirichletFunctionalType(*(this->problem().diffusion_factor()->component(pp)),
                                                                     *(diffusion_tensor.affine_part()),
                                                                     *(this->problem().dirichlet()->component(qq))));
      dirichlet_vector_assemblers.emplace_back(new DirichletVectorAssemblerType(*(dirichlet_functionals.back())));
      boundary_assembler.add(*(dirichlet_vector_assemblers.back()),
                             *(local_vector.component(vec_comp)),
                             new Stuff::Grid::ApplyOn::DirichletIntersections< BoundaryGridPartType >(this->boundary_info()));
    }
  }
  if (this->problem().dirichlet()->has_affine_part() && this->problem().diffusion_factor()->has_affine_part()) {
    if (!local_vector.has_affine_part())
      local_vector.register_affine_part(new VectorType(local_test_space.mapper().size()));
    dirichlet_functionals.emplace_back(new DirichletFunctionalType(*(this->problem().diffusion_factor()->affine_part()),
                                                                   *(diffusion_tensor.affine_part()),
                                                                   *(this->problem().dirichlet()->affine_part())));
    dirichlet_vector_assemblers.emplace_back(new DirichletVectorAssemblerType(*(dirichlet_functionals.back())));
    boundary_assembler.add(*(dirichlet_vector_assemblers.back()),
                           *(local_vector.affine_part()),
                           new Stuff::Grid::ApplyOn::DirichletIntersections< BoundaryGridPartType >(this->boundary_info()));
  } // dirichlet boundary terms

  // do the actual work
  boundary_assembler.assemble();
} // ... assemble_boundary_contributions(...)

template< class G, class R, int r, int p, Stuff::LA::ChooseBackend la >
    void BlockSWIPDG< G, R, r, p, la >::
assemble_coupling_contributions(const size_t subdomain,
                                const size_t neighbour,
                                typename BlockSWIPDG< G, R, r, p, la >::AffinelyDecomposedMatrixType& inside_inside_matrix,
                                typename BlockSWIPDG< G, R, r, p, la >::AffinelyDecomposedMatrixType& inside_outside_matrix,
                                typename BlockSWIPDG< G, R, r, p, la >::AffinelyDecomposedMatrixType& outside_inside_matrix,
                                typename BlockSWIPDG< G, R, r, p, la >::AffinelyDecomposedMatrixType& outside_outside_matrix) const
{
  typedef typename LocalDiscretizationType::TestSpaceType   LocalTestSpaceType;
  typedef typename LocalDiscretizationType::AnsatzSpaceType LocalAnsatzSpaceType;
  const LocalTestSpaceType&   inner_test_space   = this->local_discretizations_[subdomain]->test_space();
  const LocalAnsatzSpaceType& inner_ansatz_space = this->local_discretizations_[subdomain]->ansatz_space();
  const LocalTestSpaceType&   outer_test_space   = this->local_discretizations_[neighbour]->test_space();
  const LocalAnsatzSpaceType& outer_ansatz_space = this->local_discretizations_[neighbour]->ansatz_space();
  CouplingAssembler coupling_assembler(inner_test_space, inner_ansatz_space,
                                       outer_test_space, outer_ansatz_space,
                                       ms_grid_->couplingGridPart(subdomain, neighbour));

  typedef typename ProblemType::DiffusionFactorType::NonparametricType DiffusionFactorType;
  typedef typename ProblemType::DiffusionTensorType::NonparametricType DiffusionTensorType;
  const auto& diffusion_tensor = *(this->problem().diffusion_tensor());
  assert(!diffusion_tensor.parametric());
  assert(diffusion_tensor.has_affine_part());
  typedef GDT::LocalOperator::Codim1CouplingIntegral< GDT::LocalEvaluation::SWIPDG::Inner< DiffusionFactorType,
                                                                                           DiffusionTensorType > >
      CouplingOperatorType;
  typedef GDT::LocalAssembler::Codim1CouplingMatrix< CouplingOperatorType > CouplingMatrixAssemblerType;
  std::vector< std::unique_ptr< CouplingOperatorType > > coupling_operators;
  std::vector< std::unique_ptr< CouplingMatrixAssemblerType > > coupling_matrix_assemblers;
  for (ssize_t qq = 0; qq < this->problem().diffusion_factor()->num_components(); ++qq) {
    coupling_operators.emplace_back(new CouplingOperatorType(*(this->problem().diffusion_factor()->component(qq)),
                                                             *(diffusion_tensor.affine_part())));
    coupling_matrix_assemblers.emplace_back(new CouplingMatrixAssemblerType(*(coupling_operators[qq])));
    coupling_assembler.addLocalAssembler(*(coupling_matrix_assemblers[qq]),
                                         *(inside_inside_matrix.component(qq)),
                                         *(inside_outside_matrix.component(qq)),
                                         *(outside_inside_matrix.component(qq)),
                                         *(outside_outside_matrix.component(qq)));
  }
  if (this->problem().diffusion_factor()->has_affine_part()) {
    coupling_operators.emplace_back(new CouplingOperatorType(*(this->problem().diffusion_factor()->affine_part()),
                                                             *(diffusion_tensor.affine_part())));
    coupling_matrix_assemblers.emplace_back(new CouplingMatrixAssemblerType(*(
        coupling_operators[coupling_operators.size() - 1])));
    coupling_assembler.addLocalAssembler(*(coupling_matrix_assemblers[coupling_matrix_assemblers.size() - 1]),
                                         *(inside_inside_matrix.affine_part()),
                                         *(inside_outside_matrix.affine_part()),
                                         *(outside_inside_matrix.affine_part()),
                                         *(outside_outside_matrix.affine_part()));
  }

  // do the actual work
  coupling_assembler.assemble();
} // ... assemble_coupling_contributions(...)

template< class G, class R, int r, int p, Stuff::LA::ChooseBackend la >
    void BlockSWIPDG< G, R, r, p, la >::
build_global_containers()
{
  // walk the subdomains
  for (size_t ss = 0; ss < ms_grid_->size(); ++ss) {

    copy_local_to_global_matrix(*(local_matrices_[ss]),
                                this->local_discretizations_[ss]->pattern(),
                                ss,
                                ss,
                                *(this->matrix_));
    copy_local_to_global_vector(*(local_vectors_[ss]),
                                ss,
                                *(this->rhs_));

    // walk the neighbours
    for (const size_t& nn : ms_grid_->neighborsOf(ss)) {
      if (ss < nn) {
        // get the coupling patterns
        const auto result_inside_outside_pattern = inside_outside_patterns_[ss].find(nn);
        if (result_inside_outside_pattern == inside_outside_patterns_[ss].end())
          DUNE_THROW(Stuff::Exceptions::internal_error,
                     "The coupling pattern for subdomain " << ss << " and neighbour " << nn << "is missing!");
        const auto& inside_outside_pattern = *(result_inside_outside_pattern->second);
        const auto result_outside_inside_pattern = outside_inside_patterns_[nn].find(ss);
        if (result_outside_inside_pattern == outside_inside_patterns_[nn].end())
          DUNE_THROW(Stuff::Exceptions::internal_error,
                     "The coupling pattern for neighbour " << nn << " and subdomain " << ss << "is missing!");
        const auto& outside_inside_pattern = *(result_outside_inside_pattern->second);
        // and the coupling matrices
        auto result_inside_outside_matrix = inside_outside_matrices_[ss].find(nn);
        if (result_inside_outside_matrix == inside_outside_matrices_[ss].end())
          DUNE_THROW(Stuff::Exceptions::internal_error,
                     "The coupling matrix for subdomain " << ss << " and neighbour " << nn << "is missing!");
        auto& inside_outside_matrix = *(result_inside_outside_matrix->second);
        auto result_outside_inside_matrix = outside_inside_matrices_[nn].find(ss);
        if (result_outside_inside_matrix == outside_inside_matrices_[nn].end())
          DUNE_THROW(Stuff::Exceptions::internal_error,
                     "The coupling matrix for neighbour " << nn << " and subdomain " << ss << "is missing!");
        auto& outside_inside_matrix = *(result_outside_inside_matrix->second);
        // and copy them into the global matrix
        copy_local_to_global_matrix(inside_outside_matrix,
                                    inside_outside_pattern,
                                    ss, nn,
                                    *(this->matrix_));
        copy_local_to_global_matrix(outside_inside_matrix,
                                    outside_inside_pattern,
                                    nn, ss,
                                    *(this->matrix_));
      }
    } // walk the neighbours
  } // walk the subdomains
} // ... build_global_containers(...)


} // namespace Discretizations
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_BLOCK_SWIPDG_HXX
