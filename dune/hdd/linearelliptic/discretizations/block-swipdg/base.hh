// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_BLOCK_SWIPDG_BASE_HH
#define DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_BLOCK_SWIPDG_BASE_HH

#include <memory>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <limits>
#include <type_traits>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/stuff/common/disable_warnings.hh>
# if HAVE_EIGEN
#   include <Eigen/Eigenvalues>
# endif
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/common/timer.hh>
#include <dune/common/dynmatrix.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/timedlogging.hh>
#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/float_cmp.hh>
#include <dune/stuff/common/fixed_map.hh>
#include <dune/stuff/grid/layers.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/la/container.hh>
#include <dune/stuff/la/solver.hh>
#include <dune/stuff/grid/walker.hh>

#include <dune/pymor/common/exceptions.hh>

#include <dune/gdt/spaces/dg.hh>
#include <dune/gdt/playground/localevaluation/swipdg.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/oswaldinterpolation.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/spaces/fv/default.hh>
#include <dune/gdt/spaces/rt/pdelab.hh>
#include <dune/gdt/playground/operators/fluxreconstruction.hh>
#include <dune/gdt/playground/products/swipdgpenalty.hh>
#include <dune/gdt/products/boundaryl2.hh>
#include <dune/gdt/products/l2.hh>
#include <dune/gdt/products/h1.hh>
#include <dune/gdt/products/elliptic.hh>
#include <dune/gdt/assembler/system.hh>

#include <dune/hdd/linearelliptic/problems/default.hh>
#include <dune/hdd/linearelliptic/problems/zero-boundary.hh>

#include "../base.hh"
#include "../swipdg.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Discretizations {
namespace BlockSwipdg {
namespace internal {


template <class LocalGridsProviderType, Stuff::Grid::ChooseLayer local_layer, Stuff::Grid::ChooseLayer overlap_layer,
          class R, int r, int p, Stuff::LA::ChooseBackend la>
class LocalDiscretizationsContainer
{
public:
  typedef typename LocalGridsProviderType::LocalGridType LocalGridType;
  typedef SWIPDG<LocalGridType, local_layer, R, r, p, GDT::ChooseSpaceBackend::fem, la>    DiscretizationType;
  typedef SWIPDG<LocalGridType, overlap_layer, R, r, p, GDT::ChooseSpaceBackend::fem, la > OversampledDiscretizationType;
  typedef typename DiscretizationType::ProblemType     LocalProblemType;
  typedef typename DiscretizationType::TestSpaceType   TestSpaceType;
  typedef typename DiscretizationType::AnsatzSpaceType AnsatzSpaceType;

private:
  typedef Problems::ZeroBoundary< LocalProblemType > FakeProblemType;
  typedef typename DiscretizationType::GridViewType::Intersection IntersectionType;

public:
  LocalDiscretizationsContainer(LocalGridsProviderType& local_grids_provider,
                                const LocalProblemType& prob,
                                const std::vector< std::string >& only_these_products)
    : zero_boundary_problem_(prob)
    , all_dirichlet_boundary_config_(Stuff::Grid::BoundaryInfos::AllDirichlet< IntersectionType >::default_config())
    , all_neumann_boundary_config_(Stuff::Grid::BoundaryInfos::AllNeumann< IntersectionType >::default_config())
    , multiscale_boundary_config_(Stuff::Grid::BoundaryInfoConfigs::IdBased::default_config())
    , local_discretizations_(local_grids_provider.num_subdomains(), nullptr)
    , oversampled_discretizations_dirichlet_(local_grids_provider.num_subdomains(), nullptr)
    , oversampled_discretizations_neumann_(local_grids_provider.num_subdomains(), nullptr)
    , local_test_spaces_(local_grids_provider.num_subdomains(), nullptr)
    , local_ansatz_spaces_(local_grids_provider.num_subdomains(), nullptr)
  {
    multiscale_boundary_config_["neumann"] = "7";
    for (size_t ss = 0; ss < local_grids_provider.num_subdomains(); ++ss) {
      local_discretizations_[ss] = std::make_shared< DiscretizationType >(local_grids_provider.local_grid(ss),
                                                                          all_neumann_boundary_config_,
                                                                          zero_boundary_problem_,
                                                                          local_grids_provider.local_level(ss),
                                                                          only_these_products);
      local_test_spaces_[ss]   = std::make_shared< TestSpaceType >(  local_discretizations_[ss]->test_space());
      local_ansatz_spaces_[ss] = std::make_shared< AnsatzSpaceType >(local_discretizations_[ss]->ansatz_space());
    }
  }

protected:
  const FakeProblemType zero_boundary_problem_;
  const Stuff::Common::Configuration all_dirichlet_boundary_config_;
  const Stuff::Common::Configuration all_neumann_boundary_config_;
  Stuff::Common::Configuration multiscale_boundary_config_;
  std::vector< std::shared_ptr< DiscretizationType > > local_discretizations_;
  mutable std::vector< std::shared_ptr< OversampledDiscretizationType > > oversampled_discretizations_dirichlet_;
  mutable std::vector< std::shared_ptr< OversampledDiscretizationType > > oversampled_discretizations_neumann_;
  std::vector< std::shared_ptr< const TestSpaceType > > local_test_spaces_;
  std::vector< std::shared_ptr< const AnsatzSpaceType > > local_ansatz_spaces_;
}; // class LocalDiscretizationsContainer


template< class ImpTraits, class R, int r, int p, Stuff::LA::ChooseBackend la >
class BaseTraits
  : public Discretizations::internal::ContainerBasedDefaultTraits<typename Stuff::LA::Container<R, la>::MatrixType,
                                                                  typename Stuff::LA::Container<R, la>::VectorType>
{
public:
  static const constexpr Stuff::Grid::ChooseLayer local_layer = ImpTraits::local_layer;
  static const constexpr Stuff::Grid::ChooseLayer overlap_layer = ImpTraits::overlap_layer;
  typedef typename ImpTraits::derived_type derived_type;
  typedef R RangeFieldType;
  static const unsigned int dimRange = r;
//  static const unsigned int polOrder = p;
  typedef typename ImpTraits::LocalGridsProviderType LocalGridsProviderType;
  typedef LocalDiscretizationsContainer<LocalGridsProviderType, local_layer, overlap_layer, R, r, p, la>
      LocalDiscretizationsContainerType;
  typedef typename LocalDiscretizationsContainerType::TestSpaceType   LocalTestSpaceType;
  typedef typename LocalDiscretizationsContainerType::AnsatzSpaceType LocalAnsatzSpaceType;

  typedef typename ImpTraits::template BlockSpace<LocalTestSpaceType>::type   TestSpaceType;
  typedef typename ImpTraits::template BlockSpace<LocalAnsatzSpaceType>::type AnsatzSpaceType;
  typedef typename TestSpaceType::GridViewType GridViewType;
}; // class BaseTraits


/**
 * \attention The given problem is replaced by a Problems::ZeroBoundary.
 * \attention The given boundary info config is replaced by a Stuff::Grid::BoundaryInfos::AllDirichlet.
 * \attention The boundary info for the local oversampled discretizations is hardwired to dirichlet zero atm!
 */
template< class ImpTraits, class R, int r, int p, Stuff::LA::ChooseBackend la >
class Base
  : protected internal::LocalDiscretizationsContainer
        <typename ImpTraits::LocalGridsProviderType, ImpTraits::local_layer, ImpTraits::overlap_layer, R, r, p, la>
  , public ContainerBasedDefault<internal::BaseTraits<ImpTraits, R, r, p, la>>

{
  typedef internal::LocalDiscretizationsContainer
      <typename ImpTraits::LocalGridsProviderType, ImpTraits::local_layer, ImpTraits::overlap_layer, R, r, p, la>
                                                                              LocalDiscretizationsBaseType;
  typedef ContainerBasedDefault<internal::BaseTraits<ImpTraits, R, r, p, la>> BaseType;
//  typedef Base< ImpTraits, GridImp, R, r, p, la >        ThisType;
public:
  typedef internal::BaseTraits<ImpTraits, R, r, p, la> Traits;
  typedef typename Traits::LocalGridsProviderType LocalGridsProviderType;

  using typename BaseType::ProblemType;
//  using typename BaseType::GridViewType    GridViewType;
  using typename BaseType::TestSpaceType;
  using typename BaseType::AnsatzSpaceType;
  using typename BaseType::PatternType;
//  using typename BaseType::EntityType      EntityType;
//  using typename BaseType::DomainFieldType DomainFieldType;
//  using typename BaseType::RangeFieldType  RangeFieldType;
  using typename BaseType::MatrixType;
  using typename BaseType::VectorType;
//  using typename BaseType::OperatorType    OperatorType;
//  using typename BaseType::ProductType     ProductType;
//  using typename BaseType::FunctionalType  FunctionalType;

//  static const unsigned int dimDomain = BaseType::dimDomain;
//  static const unsigned int dimRange  = BaseType::dimRange;

//  typedef grid::Multiscale::ProviderInterface< GridImp > GridProviderType;
//  typedef typename GridProviderType::GridType   GridType;
//  typedef typename GridProviderType::MsGridType MsGridType;

//  typedef typename Traits::LocalDiscretizationsContainerType::DiscretizationType            LocalDiscretizationType;
//  typedef typename Traits::LocalDiscretizationsContainerType::OversampledDiscretizationType OversampledDiscretizationType;
//  typedef typename TestSpaceType::PatternType PatternType;

protected:
  using typename BaseType::AffinelyDecomposedMatrixType;
  using typename BaseType::AffinelyDecomposedVectorType;

  typedef typename Traits::LocalTestSpaceType   LocalTestSpaceType;
  typedef typename Traits::LocalAnsatzSpaceType LocalAnsatzSpaceType;
//  typedef Pymor::LA::AffinelyDecomposedConstContainer< MatrixType > AffinelyDecomposedConstMatrixType;
//  typedef Pymor::LA::AffinelyDecomposedConstContainer< VectorType > AffinelyDecomposedConstVectorType;

public:
//  typedef typename LocalDiscretizationType::ProblemType LocalProblemType;
//  typedef typename LocalDiscretizationType::ProductType LocalProductType;

  static std::string static_id()
  {
    return DiscretizationInterface< Traits >::static_id() + ".block-swipdg";
  }

  Base(LocalGridsProviderType& local_grids_provider,
       const Stuff::Common::Configuration& bound_inf_cfg,
       const ProblemType& prob,
       const std::vector< std::string >& only_these_products = {})
    : LocalDiscretizationsBaseType(local_grids_provider, prob, only_these_products)
    , BaseType(TestSpaceType(local_grids_provider.ms_grid(), this->local_test_spaces_),
               AnsatzSpaceType(local_grids_provider.ms_grid(), this->local_ansatz_spaces_),
               bound_inf_cfg,
               this->zero_boundary_problem_)
    , local_grids_provider_(local_grids_provider)
    , only_these_products_(only_these_products)
    , local_matrices_(local_grids_provider_.num_subdomains())
    , local_vectors_(local_grids_provider_.num_subdomains())
    , inside_outside_patterns_(local_grids_provider_.num_subdomains())
    , outside_inside_patterns_(local_grids_provider_.num_subdomains())
    , inside_outside_matrices_(local_grids_provider_.num_subdomains())
    , outside_inside_matrices_(local_grids_provider_.num_subdomains())
  {
    // in case of parametric diffusion tensor everything is too complicated
    if (this->problem_.diffusion_tensor()->parametric())
      DUNE_THROW(NotImplemented, "The diffusion tensor must not be parametric!");
    if (!this->problem_.diffusion_tensor()->has_affine_part())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given, "The diffusion tensor must not be empty!");

    this->pattern_ = std::make_shared<PatternType>(this->test_space().mapper().size());
  } // Base(...)

//  const std::vector< std::shared_ptr< LocalDiscretizationType > >& local_discretizations() const
//  {
//    return this->local_discretizations_;
//  }

//  void init(const bool prune = false)
//  {
//    if (this->container_based_initialized_)
//      return;

//    auto logger = Stuff::Common::TimedLogger().get("hdd.linearelliptic.discretizations.block-swipdg.init");

//    pattern_ = std::make_shared< PatternType >(BaseType::test_space().mapper().size());

//    const size_t subdomains = local_grids_provider_.num_subdomains();
//    logger.info() << "discretizing on " << subdomains << " subdomains..." << std::endl;
//    // walk the subdomains for the first time
//    //   * to initialize the coupling pattern,
//    //   * to finalize the global sparsity pattern
//    logger.info() << "  computing patterns and local contributions... " << std::endl;
//    for (size_t ss = 0; ss < subdomains; ++ss) {
//      // init the local discretizations (assembles matrices and patterns)
//      this->local_discretizations_[ss]->init(false);
//      // and create the local containers
//      // * the matrices
//      //   * just copy those from the local discretizations
//      const auto local_operator = this->local_discretizations_[ss]->get_operator();
//      local_matrices_[ss] = std::make_shared< AffinelyDecomposedMatrixType >();
//      //   * we take the affine part only if the diffusion has one, otherwise it contains only the dirichlet rows,
//      //     thus it is empty, since the local problems are purely neumann
//      if (this->problem().diffusion_factor()->has_affine_part()) {
//        if (!local_operator.has_affine_part())
//          DUNE_THROW(Stuff::Exceptions::internal_error, "The local operator is missing the affine part!");
//        local_matrices_[ss]->register_affine_part(new MatrixType(*(local_operator.affine_part().container())));
//      }
//      if (local_operator.num_components() < this->problem().diffusion_factor()->num_components())
//        DUNE_THROW(Stuff::Exceptions::requirements_not_met,
//                   "The local operator should have " << this->problem().diffusion_factor()->num_components()
//                   << " components (but has only " << local_operator.num_components() << ")!");
//      for (ssize_t qq = 0; qq < this->problem().diffusion_factor()->num_components(); ++qq) {
//        local_matrices_[ss]->register_component(new MatrixType(
//            local_operator.component(qq).container()->backend()),
//            this->problem().diffusion_factor()->coefficient(qq));
//      }
//      // * and the vectors
//      const auto local_functional = this->local_discretizations_[ss]->get_rhs();
//      local_vectors_[ss] = std::make_shared< AffinelyDecomposedVectorType >();
//      for (size_t qq = 0; qq < boost::numeric_cast< size_t >(local_functional.num_components()); ++qq)
//        local_vectors_[ss]->register_component(new VectorType(*(local_functional.component(qq).container())),
//                                               new Pymor::ParameterFunctional(local_functional.coefficient(qq)));
//      if (local_functional.has_affine_part())
//        local_vectors_[ss]->register_affine_part(new VectorType(*(local_functional.affine_part().container())));

//      // create and copy the local patterns
//      add_local_to_global_pattern(this->local_discretizations_[ss]->pattern(), ss, ss, *pattern_);
//      const auto& inner_test_space = this->local_discretizations_[ss]->test_space();
//      const auto& inner_ansatz_space = this->local_discretizations_[ss]->ansatz_space();
//      // walk the neighbors
//      for (const size_t& nn : ms_grid_->neighborsOf(ss)) {
//        // visit each coupling only once (assemble primally)
//        if (ss < nn) {
//          const auto& outer_test_space = this->local_discretizations_[nn]->test_space();
//          const auto& outer_ansatz_space = this->local_discretizations_[nn]->ansatz_space();
//          const auto inside_outside_grid_part = ms_grid_->couplingGridPart(ss, nn);
//          const auto outside_inside_grid_part = ms_grid_->couplingGridPart(nn, ss);
//          // create the coupling patterns
//          auto inside_outside_pattern = std::make_shared< PatternType >(
//                inner_test_space.compute_face_pattern(inside_outside_grid_part, outer_ansatz_space));
//          inside_outside_patterns_[ss].insert(std::make_pair(nn, inside_outside_pattern));
//          auto outside_inside_pattern = std::make_shared< PatternType >(
//                outer_test_space.compute_face_pattern(outside_inside_grid_part, inner_ansatz_space));
//          outside_inside_patterns_[nn].insert(std::make_pair(ss, outside_inside_pattern));
//          // and copy them
//          add_local_to_global_pattern(*inside_outside_pattern,  ss, nn, *pattern_);
//          add_local_to_global_pattern(*outside_inside_pattern,  nn, ss, *pattern_);
//        } // visit each coupling only once (assemble primaly)
//      } // walk the neighbors
//    } // walk the subdomains for the first time

//    // walk the subdomains for the second time
//    //   * to assemble the boundary matrices and vectors and
//    //   * to assemble the coupling matrices
//    logger.info() << "computing coupling and boundary contributions... " << std::endl;
//    for (size_t ss = 0; ss < local_grids_provider_.num_subdomains(); ++ss) {
//      const auto& inner_test_mapper = this->local_discretizations_[ss]->test_space().mapper();
//      const auto& inner_ansatz_mapper = this->local_discretizations_[ss]->ansatz_space().mapper();
//      if (ms_grid_->boundary(ss))
//        assemble_boundary_contributions(ss);
//      // walk the neighbors
//      for (const size_t& nn : ms_grid_->neighborsOf(ss)) {
//        // visit each coupling only once (assemble primaly)
//        if (ss < nn) {
//          const auto& outer_test_mapper = this->local_discretizations_[nn]->test_space().mapper();
//          const auto& outer_ansatz_mapper = this->local_discretizations_[nn]->ansatz_space().mapper();
//          // get the patterns
//          const auto in_out_result = inside_outside_patterns_[ss].find(nn);
//          if (in_out_result == inside_outside_patterns_[ss].end())
//            DUNE_THROW(Stuff::Exceptions::internal_error, "subdomain " << ss << ", neighbour " << nn);
//          const auto& inside_outside_pattern = *(in_out_result->second);
//          const auto out_in_result = outside_inside_patterns_[nn].find(ss);
//          if (out_in_result == outside_inside_patterns_[nn].end())
//            DUNE_THROW(Stuff::Exceptions::internal_error, "subdomain " << ss << ", neighbour " << nn);
//          const auto& outside_inside_pattern = *(out_in_result->second);
//          // create the coupling matrices
//          auto inside_outside_matrix = std::make_shared< AffinelyDecomposedMatrixType >();
//          auto outside_inside_matrix = std::make_shared< AffinelyDecomposedMatrixType >();
//          if (this->problem().diffusion_factor()->has_affine_part()) {
//            inside_outside_matrix->register_affine_part(new MatrixType(inner_test_mapper.size(),
//                                                                       outer_ansatz_mapper.size(),
//                                                                       inside_outside_pattern));
//            outside_inside_matrix->register_affine_part(new MatrixType(outer_test_mapper.size(),
//                                                                       inner_ansatz_mapper.size(),
//                                                                       outside_inside_pattern));
//          }
//          for (ssize_t qq = 0; qq < this->problem().diffusion_factor()->num_components(); ++qq) {
//            inside_outside_matrix->register_component(new MatrixType(inner_test_mapper.size(),
//                                                                     outer_ansatz_mapper.size(),
//                                                                     inside_outside_pattern),
//                                                      this->problem().diffusion_factor()->coefficient(qq));
//            outside_inside_matrix->register_component(new MatrixType(outer_test_mapper.size(),
//                                                                     inner_ansatz_mapper.size(),
//                                                                     outside_inside_pattern),
//                                                      this->problem().diffusion_factor()->coefficient(qq));
//          }
//          // and assemble them
//          assemble_coupling_contributions(ss, nn,
//                                          *(local_matrices_[ss]),
//                                          *(inside_outside_matrix),
//                                          *(outside_inside_matrix),
//                                          *(local_matrices_[nn]));
//          inside_outside_matrices_[ss].insert(std::make_pair(nn, inside_outside_matrix));
//          outside_inside_matrices_[nn].insert(std::make_pair(ss, outside_inside_matrix));
//        } // visit each coupling only once
//      } // walk the neighbors
//    } // walk the subdomains for the second time

//    // build global containers
//    pattern_->sort();
//    build_global_containers();

//    logger.info() << "assembling products... " << std::endl;
//    this->assemble_products(only_these_products_, 2);

//    // finalize
//    this->finalize_init(prune);

//    logger.info() << "finished!" << std::endl;
//  } // ... init(...)

//  ssize_t num_subdomains() const
//  {
//    return local_grids_provider_.num_subdomains();
//  }

//  std::vector< ssize_t > neighbouring_subdomains(const ssize_t ss) const
//  {
//    if (ss < 0 || ss >= num_subdomains())
//      DUNE_THROW(Stuff::Exceptions::index_out_of_range,
//                 "0 <= ss < num_subdomains() = " << num_subdomains() << " is not true for ss = " << ss << "!");
//    const auto set_of_neighbours = ms_grid_->neighborsOf(ss);
//    return std::vector< ssize_t >(set_of_neighbours.begin(), set_of_neighbours.end());
//  }

//  VectorType localize_vector(const VectorType& global_vector, const size_t ss) const
//  {
//    if ((std::make_signed< size_t >::type)(ss) >= num_subdomains())
//      DUNE_THROW(Stuff::Exceptions::index_out_of_range,
//                 "0 <= ss < num_subdomains() = " << num_subdomains() << " is not true for ss = " << ss << "!");
//    if (global_vector.size() != this->ansatz_space().mapper().size())
//      DUNE_THROW(Stuff::Exceptions::index_out_of_range,
//                 "The size() of global_vector (" << global_vector.dim()
//                 << ") does not match the size() of the ansatz space (" << this->ansatz_space().mapper().size() << ")!");
//    assert(ss < this->local_discretizations_.size());
//    VectorType local_vector = this->local_discretizations_[ss]->create_vector();
//    for (size_t ii = 0; ii < local_vector.size(); ++ii)
//      local_vector.set_entry(ii, global_vector.get_entry(this->ansatz_space().mapper().mapToGlobal(ss, ii)));
//    return local_vector;
//  } // ... localize_vetor(...)

//  VectorType globalize_vectors(const std::vector< VectorType >& local_vectors) const
//  {
//    if (local_vectors.size() != boost::numeric_cast< size_t >(num_subdomains()))
//      DUNE_THROW(Stuff::Exceptions::wrong_input_given,
//                 "Given local_vectors has wrong size (is " << local_vectors.size() << ", should be "
//                 << num_subdomains() << ")!");
//    VectorType ret(this->ansatz_space().mapper().size());
//    for (size_t ss = 0; ss < boost::numeric_cast< size_t >(num_subdomains()); ++ss) {
//      const auto& local_vector = local_vectors[ss];
//      if (local_vector.size() != this->local_discretizations_[ss]->ansatz_space().mapper().size())
//        DUNE_THROW(Stuff::Exceptions::wrong_input_given,
//                   "Given local_vectors[" << ss << "] has wrong size (is "
//                   << local_vector.size() << ", should be "
//                   << this->local_discretizations_[ss]->ansatz_space().mapper().size() << ")!");
//      copy_local_to_global_vector(local_vector, ss, ret);
//    }
//    return ret;
//  }

//  VectorType* globalize_vectors_and_return_ptr(const std::vector< VectorType >& local_vectors) const
//  {
//    return new VectorType(globalize_vectors(local_vectors));
//  }

//  VectorType* localize_vector_and_return_ptr(const VectorType& global_vector, const ssize_t ss) const
//  {
//    return new VectorType(localize_vector(global_vector, boost::numeric_cast< size_t >(ss)));
//  }

//  ProductType get_local_product(const size_t ss, const std::string id) const
//  {
//    if (boost::numeric_cast< ssize_t >(ss) >= num_subdomains())
//      DUNE_THROW(Stuff::Exceptions::index_out_of_range,
//                 "0 <= ss < num_subdomains() = " << num_subdomains() << " is not true for ss = " << ss << "!");
//    return this->local_discretizations_[ss]->get_product(id);
//  }

//  ProductType* get_local_product_and_return_ptr(const ssize_t ss, const std::string id) const
//  {
//    return new ProductType(get_local_product(boost::numeric_cast< size_t >(ss), id));
//  }

//  OperatorType get_local_operator(const size_t ss) const
//  {
//    if (boost::numeric_cast< ssize_t >(ss) >= num_subdomains())
//      DUNE_THROW(Stuff::Exceptions::index_out_of_range,
//                 "0 <= ss < num_subdomains() = " << num_subdomains() << " is not true for ss = " << ss << "!");
//    assert(ss < local_matrices_.size());
//    return OperatorType(*(local_matrices_[ss]));
//  }

//  OperatorType* get_local_operator_and_return_ptr(const ssize_t ss) const
//  {
//    return new OperatorType(get_local_operator(boost::numeric_cast< size_t >(ss)));
//  }

//  OperatorType get_coupling_operator(const size_t ss, const size_t nn) const
//  {
//    if (ss >= boost::numeric_cast< size_t >(num_subdomains()))
//      DUNE_THROW(Stuff::Exceptions::index_out_of_range,
//                 "0 <= ss < num_subdomains() = " << num_subdomains() << " is not true for ss = " << ss << "!");
//    const auto neighbours = ms_grid_->neighborsOf(ss);
//    if (neighbours.count(nn) == 0)
//      DUNE_THROW(Stuff::Exceptions::index_out_of_range,
//                 "Subdomain " << nn << " is not a neighbour of subdomain " << ss
//                 << " (call neighbouring_subdomains(" << ss << ") to find out)!");
//    if (ss < nn) {
//      // we need to look for this coupling operator in the inside/outside context
//      const auto result_inside_outside_matrix = inside_outside_matrices_[ss].find(nn);
//      if (result_inside_outside_matrix == inside_outside_matrices_[ss].end())
//        DUNE_THROW(Stuff::Exceptions::internal_error,
//                   "The coupling matrix for subdomain " << ss << " and neighbour " << nn << " is missing!");
//      const auto inside_outside_matrix = result_inside_outside_matrix->second;
//      return OperatorType(*inside_outside_matrix);
//    } else if (nn < ss) {
//      // we need to look for this coupling operator in the outside/inside context
//      const auto result_outside_inside_matrix = outside_inside_matrices_[ss].find(nn);
//      if (result_outside_inside_matrix == outside_inside_matrices_[ss].end())
//        DUNE_THROW(Stuff::Exceptions::internal_error,
//                   "The coupling matrix for neighbour " << nn << " and subdomain " << ss << " is missing!");
//      const auto outside_inside_matrix = result_outside_inside_matrix->second;
//      return OperatorType(*outside_inside_matrix);
//    } else {
//      // the above exception should have cought this
//      DUNE_THROW(Stuff::Exceptions::internal_error,
//                 "The multiscale grid is corrupted! Subdomain " << ss << " must not be its own neighbour!");
//    }
//  } // ... get_coupling_operator(...)

//  OperatorType* get_coupling_operator_and_return_ptr(const ssize_t ss, const ssize_t nn) const
//  {
//    return new OperatorType(get_coupling_operator(boost::numeric_cast< size_t >(ss),
//                                                  boost::numeric_cast< size_t >(nn)));
//  }

//  FunctionalType get_local_functional(const size_t ss) const
//  {
//    if (ss >= boost::numeric_cast< size_t >(num_subdomains()))
//      DUNE_THROW(Stuff::Exceptions::index_out_of_range,
//                 "0 <= ss < num_subdomains() = " << num_subdomains() << " is not true for ss = " << ss << "!");
//    assert(ss < local_vectors_.size());
//    return FunctionalType(*(local_vectors_[ss]));
//  }

//  FunctionalType* get_local_functional_and_return_ptr(const ssize_t ss) const
//  {
//    return new FunctionalType(get_local_functional(boost::numeric_cast< size_t >(ss)));
//  }

//  VectorType solve_for_local_correction(const std::vector< VectorType >& local_vectors,
//                                        const size_t subdomain,
//                                        const Pymor::Parameter mu = Pymor::Parameter()) const
//  {
//    DUNE_THROW(Stuff::Exceptions::internal_error, "Do not call this method, I do not trust it!");
//    using namespace GDT;

//    if (mu.type() != this->parameter_type())
//      DUNE_THROW(Pymor::Exceptions::wrong_parameter_type,
//                 "mu is " << mu.type() << ", should be " << this->parameter_type() << "!");
//    if (local_vectors.size() != local_grids_provider_.num_subdomains())
//      DUNE_THROW(Stuff::Exceptions::wrong_input_given,
//                 "local_vectors is of size " << local_vectors.size() << " and should be of size " << local_grids_provider_.num_subdomains());
//    VectorType vector(this->ansatz_space().mapper().size());
//    for (size_t ss = 0; ss < local_grids_provider_.num_subdomains(); ++ss) {
//      if (local_vectors[ss].size() != this->local_discretizations_[ss]->ansatz_space().mapper().size())
//        DUNE_THROW(Stuff::Exceptions::wrong_input_given,
//                   "local_vectors[" << ss << "] is of size " << local_vectors[ss].size() << " and should be of size "
//                   << this->local_discretizations_[ss]->ansatz_space().mapper().size());
//      if (!local_vectors[ss].valid())
//        DUNE_THROW(Stuff::Exceptions::wrong_input_given, "local_vectors[" << ss << "] contains NaN or INF!");
//      copy_local_to_global_vector(local_vectors[ss], ss, vector);
//    }

//  //    const std::string prefix = "subdomain_" + DSC::toString(subdomain) + "_";

//    const ConstDiscreteFunction< AnsatzSpaceType, VectorType > current_global_solution(this->ansatz_space(), vector);
//  //    current_global_solution.visualize(prefix + "current_solution_global");

//    typedef SWIPDG< typename MsGridType::GridType, Stuff::Grid::ChooseLayer::local_oversampled, RangeFieldType, dimRange
//                  , 1, GDT::ChooseSpaceBackend::fem, la > OversampledDiscretizationType;
//    OversampledDiscretizationType oversampled_discretization(grid_provider_,
//                                                             this->multiscale_boundary_config_,
//                                                             this->problem(),
//                                                             boost::numeric_cast< int >(subdomain));
//    oversampled_discretization.init();

//    DiscreteFunction< typename OversampledDiscretizationType::AnsatzSpaceType, VectorType >
//        current_oversampled_solution(oversampled_discretization.ansatz_space());
//    const Operators::Projection< typename OversampledDiscretizationType::GridViewType >
//        oversampled_projection_operator(oversampled_discretization.grid_view());
//    oversampled_projection_operator.apply(current_global_solution, current_oversampled_solution);
//  //    current_oversampled_solution.visualize(prefix + "current_solution_oversampled");

//    if (!oversampled_discretization.rhs()->has_affine_part())
//      oversampled_discretization.rhs()->register_affine_part(oversampled_discretization.test_space().mapper().size());
//    if (oversampled_discretization.system_matrix()->parametric()) {
//      const auto oversampled_system_matrix = oversampled_discretization.system_matrix()->freeze_parameter(mu);
//      *(oversampled_discretization.rhs()->affine_part())
//          -= oversampled_system_matrix * current_oversampled_solution.vector();
//    } else {
//      const auto& oversampled_system_matrix = *(oversampled_discretization.system_matrix()->affine_part());
//      *(oversampled_discretization.rhs()->affine_part())
//          -= oversampled_system_matrix * current_oversampled_solution.vector();
//    }

//    oversampled_discretization.solve(current_oversampled_solution.vector(), mu);
//  //    current_oversampled_solution.visualize(prefix + "correction_oversampled");

//    DiscreteFunction< typename LocalDiscretizationType::AnsatzSpaceType, VectorType >
//        local_solution(this->local_discretizations_[subdomain]->ansatz_space());
//    const Operators::Projection< typename LocalDiscretizationType::GridViewType >
//        local_projection_operator(this->local_discretizations_[subdomain]->grid_view());
//    local_projection_operator.apply(current_oversampled_solution, local_solution);
//  //    local_solution.visualize(prefix + "correction_local");

//    return local_solution.vector();
//  } // ... solve_for_local_correction(...)

//  LocalDiscretizationType get_local_discretization(const size_t subdomain) const
//  {
//    if (subdomain >= this->grid_provider_.num_subdomains())
//      DUNE_THROW(Stuff::Exceptions::index_out_of_range,
//                 "Given subdomain " << subdomain << " too large (has to be smaller than "
//                 << this->grid_provider_.num_subdomains() << "!");
//    return *(this->local_discretizations_[subdomain]);
//  } // ... get_local_discretization(...)

//  LocalDiscretizationType* pb_get_local_discretization(const ssize_t subdomain) const
//  {
//    size_t ss = std::numeric_limits< size_t >::max();
//    try {
//      ss = boost::numeric_cast< size_t >(subdomain);
//    } catch (boost::bad_numeric_cast& ee) {
//      DUNE_THROW(Stuff::Exceptions::index_out_of_range,
//                 "There was an error in boost converting " << subdomain << " to "
//                 << Stuff::Common::Typename< size_t >::value() << ":\n\n" << ee.what());
//    }
//    return new LocalDiscretizationType(get_local_discretization(ss));
//  } // ... pb_get_local_discretization(...)

//  OversampledDiscretizationType get_oversampled_discretization(const size_t subdomain,
//                                                               const std::string boundary_value_type) const
//  {
//    if (subdomain >= this->grid_provider_.num_subdomains())
//      DUNE_THROW(Stuff::Exceptions::index_out_of_range,
//                 "Given subdomain " << subdomain << " too large (has to be smaller than "
//                 << this->grid_provider_.num_subdomains() << "!");
//    if (boundary_value_type == "dirichlet") {
//      if (!this->oversampled_discretizations_dirichlet_[subdomain]) {
//        this->oversampled_discretizations_dirichlet_[subdomain] = std::make_shared< OversampledDiscretizationType >(
//                                                                    grid_provider_,
//                                                                    this->all_dirichlet_boundary_config_,
//                                                                    this->zero_boundary_problem_,
//                                                                    subdomain,
//                                                                    only_these_products_);
//        this->oversampled_discretizations_dirichlet_[subdomain]->init();
//      }
//      return *(this->oversampled_discretizations_dirichlet_[subdomain]);
//    } else if (boundary_value_type == "neumann") {
//      if (!this->oversampled_discretizations_neumann_[subdomain]) {
//        this->oversampled_discretizations_neumann_[subdomain] = std::make_shared< OversampledDiscretizationType >(
//                                                                  grid_provider_,
//                                                                  this->all_neumann_boundary_config_,
//                                                                  this->zero_boundary_problem_,
//                                                                  subdomain,
//                                                                  only_these_products_);
//        this->oversampled_discretizations_neumann_[subdomain]->init();
//      }
//      return *(this->oversampled_discretizations_neumann_[subdomain]);
//    } else {
//      DUNE_THROW(Stuff::Exceptions::wrong_input_given,
//                 "Unknown boundary_value_type given (has to be dirichlet or neumann): " << boundary_value_type);
//      return *(this->oversampled_discretizations_dirichlet_[subdomain]);
//    }
//  } // ... get_oversampled_discretization(...)

//  OversampledDiscretizationType* pb_get_oversampled_discretization(const ssize_t subdomain,
//                                                                   const std::string boundary_value_type) const
//  {
//    size_t ss = std::numeric_limits< size_t >::max();
//    try {
//      ss = boost::numeric_cast< size_t >(subdomain);
//    } catch (boost::bad_numeric_cast& ee) {
//      DUNE_THROW(Stuff::Exceptions::index_out_of_range,
//                 "There was an error in boost converting " << subdomain << " to "
//                 << Stuff::Common::Typename< size_t >::value() << ":\n\n" << ee.what());
//    }
//    return new OversampledDiscretizationType(get_oversampled_discretization(ss, boundary_value_type));
//  } // ... pb_get_local_discretization(...)

protected:
  void assemble_local_contributions(const size_t subdomain)
  {
    // init the local discretizations (assembles matrices and patterns)
    this->local_discretizations_[subdomain]->init(false);
    // and create the local containers
    // * the matrices
    local_matrices_[subdomain] = std::make_shared<AffinelyDecomposedMatrixType>();
    //   * just copy those from the local discretizations
    const auto local_operator = this->local_discretizations_[subdomain]->get_operator();
    //   * we take the affine part only if the diffusion has one, otherwise it contains only the dirichlet rows,
    //     thus it is empty, since the local problems are purely neumann
    if (this->problem().diffusion_factor()->has_affine_part()) {
      if (!local_operator.has_affine_part())
        DUNE_THROW(Stuff::Exceptions::internal_error, "The local operator is missing the affine part!");
      local_matrices_[subdomain]->register_affine_part(new MatrixType(*(local_operator.affine_part().container())));
    }
    if (local_operator.num_components() < this->problem().diffusion_factor()->num_components())
      DUNE_THROW(Stuff::Exceptions::requirements_not_met,
                 "The local operator should have " << this->problem().diffusion_factor()->num_components()
                 << " components (but has only " << local_operator.num_components() << ")!");
    for (ssize_t qq = 0; qq < this->problem().diffusion_factor()->num_components(); ++qq) {
      local_matrices_[subdomain]->register_component(new MatrixType(
          local_operator.component(qq).container()->backend()),
          this->problem().diffusion_factor()->coefficient(qq));
    }
    // * and the vectors
    const auto local_functional = this->local_discretizations_[subdomain]->get_rhs();
    local_vectors_[subdomain] = std::make_shared< AffinelyDecomposedVectorType >();
    for (size_t qq = 0; qq < boost::numeric_cast< size_t >(local_functional.num_components()); ++qq)
      local_vectors_[subdomain]->register_component(new VectorType(*(local_functional.component(qq).container())),
                                             new Pymor::ParameterFunctional(local_functional.coefficient(qq)));
    if (local_functional.has_affine_part())
      local_vectors_[subdomain]->register_affine_part(new VectorType(*(local_functional.affine_part().container())));
  } // ... assemble_local_contributions(...)

  template <class BoundaryGridViewType, class BoundaryAssemblerType>
  void assemble_boundary_contributions(const size_t subdomain, BoundaryAssemblerType& boundary_assembler) const
  {
    const auto& local_test_space   = this->local_discretizations_[subdomain]->test_space();
    const auto& local_ansatz_space = this->local_discretizations_[subdomain]->ansatz_space();
    const auto& local_pattern = this->local_discretizations_[subdomain]->pattern();
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
                                         local_test_space.mapper().size(),
                                         local_ansatz_space.mapper().size(),
                                         local_pattern);
      dirichlet_operators.emplace_back(new DirichletOperatorType(*(this->problem().diffusion_factor()->component(qq)),
                                                                 *(diffusion_tensor.affine_part())));
      dirichlet_matrix_assemblers.emplace_back(new DirichletMatrixAssemblerType(*dirichlet_operators.back()));
      boundary_assembler.add(*dirichlet_matrix_assemblers.back(),
                             *(local_matrix.component(mat_comp)),
                             new Stuff::Grid::ApplyOn::DirichletIntersections< BoundaryGridViewType >(this->boundary_info()));
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
                             new Stuff::Grid::ApplyOn::DirichletIntersections< BoundaryGridViewType >(this->boundary_info()));
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
                                         local_test_space.mapper().size());
      neumann_functionals.emplace_back(new NeumannFunctionalType(*(this->problem().neumann()->component(qq))));
      neumann_vector_assemblers.emplace_back(new NeumannVectorAssemblerType(*(neumann_functionals.back())));
      boundary_assembler.add(*(neumann_vector_assemblers.back()),
                             *(local_vector.component(vec_comp)),
                             new Stuff::Grid::ApplyOn::NeumannIntersections< BoundaryGridViewType >(this->boundary_info()));
    }
    if (this->problem().neumann()->has_affine_part()) {
      if (!local_vector.has_affine_part())
        local_vector.register_affine_part(local_test_space.mapper().size());
      neumann_functionals.emplace_back(new NeumannFunctionalType(*(this->problem().neumann()->affine_part())));
      neumann_vector_assemblers.emplace_back(new NeumannVectorAssemblerType(*(neumann_functionals.back())));
      boundary_assembler.add(*(neumann_vector_assemblers.back()),
                             *(local_vector.affine_part()),
                             new Stuff::Grid::ApplyOn::NeumannIntersections< BoundaryGridViewType >(this->boundary_info()));
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
                                           local_test_space.mapper().size());
        dirichlet_functionals.emplace_back(new DirichletFunctionalType(*(this->problem().diffusion_factor()->affine_part()),
                                                                       *(diffusion_tensor.affine_part()),
                                                                       *(this->problem().dirichlet()->component(qq))));
        dirichlet_vector_assemblers.emplace_back(new DirichletVectorAssemblerType(*dirichlet_functionals.back()));
        boundary_assembler.add(*dirichlet_vector_assemblers.back(),
                               *local_vector.component(vec_comp),
                               new Stuff::Grid::ApplyOn::DirichletIntersections< BoundaryGridViewType >(this->boundary_info()));
      }
    }
    if (this->problem().dirichlet()->has_affine_part()) {
      for (ssize_t qq = 0; qq < this->problem().diffusion_factor()->num_components(); ++qq) {
        auto vec_comp = find_add_component(local_vector,
                                           *this->problem().diffusion_factor()->coefficient(qq),
                                           local_test_space.mapper().size());
        dirichlet_functionals.emplace_back(new DirichletFunctionalType(*(this->problem().diffusion_factor()->component(qq)),
                                                                       *(diffusion_tensor.affine_part()),
                                                                       *(this->problem().dirichlet()->affine_part())));
        dirichlet_vector_assemblers.emplace_back(new DirichletVectorAssemblerType(*dirichlet_functionals.back()));
        boundary_assembler.add(*dirichlet_vector_assemblers.back(),
                               *local_vector.component(vec_comp),
                               new Stuff::Grid::ApplyOn::DirichletIntersections< BoundaryGridViewType >(this->boundary_info()));
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
        auto vec_comp = find_add_component(local_vector, coefficient, local_test_space.mapper().size());
        dirichlet_functionals.emplace_back(new DirichletFunctionalType(*(this->problem().diffusion_factor()->component(pp)),
                                                                       *(diffusion_tensor.affine_part()),
                                                                       *(this->problem().dirichlet()->component(qq))));
        dirichlet_vector_assemblers.emplace_back(new DirichletVectorAssemblerType(*(dirichlet_functionals.back())));
        boundary_assembler.add(*(dirichlet_vector_assemblers.back()),
                               *(local_vector.component(vec_comp)),
                               new Stuff::Grid::ApplyOn::DirichletIntersections< BoundaryGridViewType >(this->boundary_info()));
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
                             new Stuff::Grid::ApplyOn::DirichletIntersections< BoundaryGridViewType >(this->boundary_info()));
    } // dirichlet boundary terms

    // do the actual work
    boundary_assembler.assemble();
  } // ... assemble_boundary_contributions(...)

  template<class IntersectionType>
  class CouplingAssemblerBase
  {
  protected:
    typedef Dune::DynamicMatrix<R> LocalMatrixType;
    typedef Dune::DynamicVector<R> LocalVectorType;
    typedef std::vector<std::vector<LocalMatrixType>> LocalMatricesContainerType;
    typedef std::vector<std::vector<LocalVectorType>> LocalVectorsContainerType;
    typedef std::vector<Dune::DynamicVector<size_t>> IndicesContainer;

    typedef typename Traits::LocalTestSpaceType LocalTestSpaceType;
    typedef typename Traits::LocalAnsatzSpaceType LocalAnsatzSpaceType;

    class LocalCodim1MatrixAssemblerApplication
    {
    public:
      virtual ~LocalCodim1MatrixAssemblerApplication(){}

      virtual void apply(const LocalTestSpaceType& /*inner_test_space*/,
                         const LocalAnsatzSpaceType& /*inner_ansatz_space*/,
                         const LocalTestSpaceType& /*outer_test_space*/,
                         const LocalAnsatzSpaceType& /*outer_ansatz_space*/,
                         const IntersectionType& /*_intersection*/,
                         LocalMatricesContainerType& /*_localMatricesContainer*/,
                         IndicesContainer& /*indicesContainer*/) const = 0;

      virtual std::vector< size_t > numTmpObjectsRequired() const = 0;
    };

    template< class LocalAssemblerType, class M >
    class LocalCodim1MatrixAssemblerWrapper
      : public LocalCodim1MatrixAssemblerApplication
    {
    public:
      LocalCodim1MatrixAssemblerWrapper(const LocalAssemblerType& localAssembler,
                                        Dune::Stuff::LA::MatrixInterface< M >& in_in_matrix,
                                        Dune::Stuff::LA::MatrixInterface< M >& in_out_matrix,
                                        Dune::Stuff::LA::MatrixInterface< M >& out_in_matrix,
                                        Dune::Stuff::LA::MatrixInterface< M >& out_out_matrix)
        : localMatrixAssembler_(localAssembler)
        , in_in_matrix_(in_in_matrix)
        , out_out_matrix_(out_out_matrix)
        , in_out_matrix_(in_out_matrix)
        , out_in_matrix_(out_in_matrix)
      {}

      virtual void apply(const LocalTestSpaceType& inner_test_space,
                         const LocalAnsatzSpaceType& inner_ansatz_space,
                         const LocalTestSpaceType& outer_test_space,
                         const LocalAnsatzSpaceType& outer_ansatz_space,
                         const IntersectionType& intersection,
                         LocalMatricesContainerType& localMatricesContainer,
                         IndicesContainer& indicesContainer) const
      {
        localMatrixAssembler_.assembleLocal(inner_test_space, inner_ansatz_space,
                                            outer_test_space, outer_ansatz_space,
                                            intersection,
                                            in_in_matrix_, out_out_matrix_, in_out_matrix_, out_in_matrix_,
                                            localMatricesContainer, indicesContainer);
      }

      virtual std::vector< size_t > numTmpObjectsRequired() const
      {
        return localMatrixAssembler_.numTmpObjectsRequired();
      }

    private:
      const LocalAssemblerType& localMatrixAssembler_;
      Dune::Stuff::LA::MatrixInterface< M >& in_in_matrix_;
      Dune::Stuff::LA::MatrixInterface< M >& out_out_matrix_;
      Dune::Stuff::LA::MatrixInterface< M >& in_out_matrix_;
      Dune::Stuff::LA::MatrixInterface< M >& out_in_matrix_;
    }; // class LocalCodim1MatrixAssemblerWrapper

  public:
    CouplingAssemblerBase(const LocalTestSpaceType& inner_test_space,
                          const LocalAnsatzSpaceType& inner_ansatz_space,
                          const LocalTestSpaceType& outer_test_space,
                          const LocalAnsatzSpaceType& outer_ansatz_space)
      : innerTestSpace_(inner_test_space)
      , innerAnsatzSpace_(inner_ansatz_space)
      , outerTestSpace_(outer_test_space)
      , outerAnsatzSpace_(outer_ansatz_space)
    {}

    ~CouplingAssemblerBase()
    {
      clearLocalAssemblers();
    }

    void clearLocalAssemblers()
    {
      for (auto& element: localCodim1MatrixAssemblers_)
        delete element;
    }

    template< class L, class M >
    void addLocalAssembler(const GDT::LocalAssembler::Codim1CouplingMatrix< L >& localAssembler,
                           Dune::Stuff::LA::MatrixInterface< M >& in_in_matrix,
                           Dune::Stuff::LA::MatrixInterface< M >& in_out_matrix,
                           Dune::Stuff::LA::MatrixInterface< M >& out_in_matrix,
                           Dune::Stuff::LA::MatrixInterface< M >& out_out_matrix)
    {
      assert(in_in_matrix.rows() == innerTestSpace_.mapper().size());
      assert(in_in_matrix.cols() == innerAnsatzSpace_.mapper().size());
      assert(in_out_matrix.rows() == innerTestSpace_.mapper().size());
      assert(in_out_matrix.cols() == outerAnsatzSpace_.mapper().size());
      assert(out_in_matrix.rows() == outerTestSpace_.mapper().size());
      assert(out_in_matrix.cols() == innerAnsatzSpace_.mapper().size());
      assert(out_out_matrix.rows() == outerTestSpace_.mapper().size());
      assert(out_out_matrix.cols() == outerAnsatzSpace_.mapper().size());
      localCodim1MatrixAssemblers_.push_back(
            new LocalCodim1MatrixAssemblerWrapper< GDT::LocalAssembler::Codim1CouplingMatrix< L >, M >(
              localAssembler, in_in_matrix, in_out_matrix, out_in_matrix, out_out_matrix));
    }

  protected:
    const LocalTestSpaceType& innerTestSpace_;
    const LocalAnsatzSpaceType& innerAnsatzSpace_;
    const LocalTestSpaceType& outerTestSpace_;
    const LocalAnsatzSpaceType& outerAnsatzSpace_;
    std::vector< LocalCodim1MatrixAssemblerApplication* > localCodim1MatrixAssemblers_;
  }; // class CouplingAssemblerBase

  template <class CouplingAssemblerType>
  void assemble_coupling_contributions(const size_t subdomain,
                                       const size_t neighboring_subomain,
                                       CouplingAssemblerType& coupling_assembler)
  {
    const auto& inner_test_space   = this->local_discretizations_[subdomain]->test_space();
    const auto& inner_ansatz_space = this->local_discretizations_[subdomain]->ansatz_space();
    const auto& outer_test_space   = this->local_discretizations_[neighboring_subomain]->test_space();
    const auto& outer_ansatz_space = this->local_discretizations_[neighboring_subomain]->ansatz_space();
    // get matrices and patterns
    inside_outside_matrices_[subdomain].insert(std::make_pair(neighboring_subomain,
                                                              std::make_shared<AffinelyDecomposedMatrixType>()));
    outside_inside_matrices_[neighboring_subomain].insert(std::make_pair(subdomain,
                                                                         std::make_shared<AffinelyDecomposedMatrixType>()));
    auto& inside_inside_matrix = *local_matrices_[subdomain];
    assert(local_matrices_[neighboring_subomain]);
    auto& outside_outside_matrix = *local_matrices_[neighboring_subomain];
    auto& inside_outside_matrix = *inside_outside_matrices_[subdomain][neighboring_subomain];
    auto& outside_inside_matrix = *outside_inside_matrices_[neighboring_subomain][subdomain];
    const auto& inside_inside_pattern = this->local_discretizations_[subdomain]->pattern();
    const auto& outside_outside_pattern = this->local_discretizations_[neighboring_subomain]->pattern();
    const auto& inside_outside_pattern = *inside_outside_patterns_[subdomain][neighboring_subomain];
    const auto& outside_inside_pattern = *outside_inside_patterns_[neighboring_subomain][subdomain];
    // assemble
    typedef typename ProblemType::DiffusionFactorType::NonparametricType DiffusionFactorType;
    typedef typename ProblemType::DiffusionTensorType::NonparametricType DiffusionTensorType;
    const auto& diffusion_factor = *(this->problem().diffusion_factor());
    const auto& diffusion_tensor = *(this->problem().diffusion_tensor());
    typedef GDT::LocalOperator::Codim1CouplingIntegral
        <GDT::LocalEvaluation::SWIPDG::Inner<DiffusionFactorType, DiffusionTensorType>> CouplingOperatorType;
    typedef GDT::LocalAssembler::Codim1CouplingMatrix<CouplingOperatorType>             CouplingMatrixAssemblerType;
    std::vector<std::unique_ptr<CouplingOperatorType>> coupling_operators;
    std::vector<std::unique_ptr<CouplingMatrixAssemblerType>> coupling_matrix_assemblers;
    for (ssize_t qq = 0; qq < diffusion_factor.num_components(); ++qq) {
      auto id_in_in = find_add_component(inside_inside_matrix, *diffusion_factor.coefficient(qq),
                                         inner_test_space.mapper().size(), inner_ansatz_space.mapper().size(),
                                         inside_inside_pattern);
      auto id_out_out = find_add_component(outside_outside_matrix, *diffusion_factor.coefficient(qq),
                                           outer_test_space.mapper().size(), outer_ansatz_space.mapper().size(),
                                           outside_outside_pattern);
      auto id_in_out = find_add_component(inside_outside_matrix, *diffusion_factor.coefficient(qq),
                                          inner_test_space.mapper().size(), outer_ansatz_space.mapper().size(),
                                          inside_outside_pattern);
      auto id_out_in = find_add_component(outside_inside_matrix, *diffusion_factor.coefficient(qq),
                                          outer_test_space.mapper().size(), inner_ansatz_space.mapper().size(),
                                          outside_inside_pattern);
      coupling_operators.emplace_back(new CouplingOperatorType(*diffusion_factor.component(qq),
                                                               *diffusion_tensor.affine_part()));
      coupling_matrix_assemblers.emplace_back(new CouplingMatrixAssemblerType(*coupling_operators.back()));
      coupling_assembler.addLocalAssembler(*coupling_matrix_assemblers.back(),
                                           *inside_inside_matrix.component(id_in_in),
                                           *inside_outside_matrix.component(id_in_out),
                                           *outside_inside_matrix.component(id_out_in),
                                           *outside_outside_matrix.component(id_out_out));
    }
    if (this->problem().diffusion_factor()->has_affine_part()) {
      ensure_affine_part(inside_inside_matrix,
                         inner_test_space.mapper().size(), inner_ansatz_space.mapper().size(), inside_inside_pattern);
      ensure_affine_part(inside_outside_matrix,
                         inner_test_space.mapper().size(), outer_ansatz_space.mapper().size(), inside_outside_pattern);
      ensure_affine_part(outside_inside_matrix,
                         outer_test_space.mapper().size(), inner_ansatz_space.mapper().size(), outside_inside_pattern);
      ensure_affine_part(outside_outside_matrix,
                         outer_test_space.mapper().size(), outer_ansatz_space.mapper().size(), inside_inside_pattern);
      coupling_operators.emplace_back(new CouplingOperatorType(*diffusion_factor.affine_part(),
                                                               *diffusion_tensor.affine_part()));
      coupling_matrix_assemblers.emplace_back(new CouplingMatrixAssemblerType(*coupling_operators.back()));
      coupling_assembler.addLocalAssembler(*coupling_matrix_assemblers.back(),
                                           *inside_inside_matrix.affine_part(),
                                           *inside_outside_matrix.affine_part(),
                                           *outside_inside_matrix.affine_part(),
                                           *outside_outside_matrix.affine_part());
    }

    // do the actual work
    coupling_assembler.assemble();
  } // ... assemble_coupling_contributions(...)

  void build_global_containers()
  {
    // walk the subdomains to build the global pattern
    for (size_t subdomain = 0; subdomain < local_grids_provider_.num_subdomains(); ++subdomain) {
      // copy the local patterns
      add_local_to_global_pattern(this->local_discretizations_[subdomain]->pattern(), subdomain, subdomain);
      // walk the neighbors
      for (const auto& element : inside_outside_patterns_[subdomain]) {
        const size_t neighbor = element.first;
        // copy the coupling patterns
        const auto& inside_outside_pattern = *element.second;
        add_local_to_global_pattern(inside_outside_pattern, subdomain, neighbor);
        const auto& outside_inside_pattern = *outside_inside_patterns_[neighbor].at(subdomain);
        add_local_to_global_pattern(outside_inside_pattern, neighbor, subdomain);
      } // walk the neighbors
    } // walk the subdomains to build the global pattern
    this->pattern_->sort();

    // walk the subdomains to build the global container
    for (size_t subdomain = 0; subdomain < local_grids_provider_.num_subdomains(); ++subdomain) {
      // local container
      copy_local_to_global_matrix(*local_matrices_[subdomain],
                                  this->local_discretizations_[subdomain]->pattern(),
                                  subdomain,
                                  subdomain);
      copy_local_to_global_vector(*local_vectors_[subdomain], subdomain);
      // walk the neighbors
      for (const auto& element : inside_outside_matrices_[subdomain]) {
        const size_t neighbor = element.first;
        // copy the coupling matrices
        const auto& inside_outside_matrix = *element.second;
        copy_local_to_global_matrix(inside_outside_matrix,
                                    *inside_outside_patterns_[subdomain].at(neighbor),
                                    subdomain,
                                    neighbor);
        const auto& outside_inside_matrix = *outside_inside_matrices_[neighbor].at(subdomain);
        copy_local_to_global_matrix(outside_inside_matrix,
                                    *outside_inside_patterns_[neighbor].at(subdomain),
                                    neighbor,
                                    subdomain);
      } // walk the neighbors
    } // walk the subdomains to build the global container
  } // ... build_global_containers(...)

  void copy_local_to_global_matrix(const AffinelyDecomposedMatrixType& local_matrix,
                                   const PatternType& local_pattern,
                                   const size_t subdomain,
                                   const size_t neighbor)
  {
    for (ssize_t qq = 0; qq < local_matrix.num_components(); ++qq) {
      auto comp = find_add_component(*this->matrix_,
                                     *local_matrix.coefficient(qq),
                                     this->test_space().mapper().size(),
                                     this->ansatz_space().mapper().size(),
                                     *this->pattern_);
      copy_local_to_global_matrix(*local_matrix.component(qq),
                                  local_pattern,
                                  subdomain,
                                  neighbor,
                                  *this->matrix_->component(comp));
    }
    if (local_matrix.has_affine_part()) {
      ensure_affine_part(*this->matrix_,
                         this->test_space_.mapper().size(),
                         this->ansatz_space_.mapper().size(),
                         *this->pattern_);
      copy_local_to_global_matrix(*local_matrix.affine_part(),
                                  local_pattern,
                                  subdomain,
                                  neighbor,
                                  *this->matrix_->affine_part());
    }
  } // copy_local_to_global_matrix(...)

  void copy_local_to_global_matrix(const MatrixType& local_matrix,
                                   const PatternType& local_pattern,
                                   const size_t test_subdomain,
                                   const size_t ansatz_subdomain,
                                   MatrixType& global_matrix) const
  {
    const auto& test_mapper = this->test_space_.mapper();
    const auto& ansatz_mapper = this->ansatz_space_.mapper();
    for (size_t local_ii = 0; local_ii < local_pattern.size(); ++local_ii) {
      const size_t global_ii = test_mapper.mapToGlobal(test_subdomain, local_ii);
      for (const size_t& local_jj : local_pattern.inner(local_ii)) {
        const size_t global_jj = ansatz_mapper.mapToGlobal(ansatz_subdomain, local_jj);
        global_matrix.add_to_entry(global_ii, global_jj, local_matrix.get_entry(local_ii, local_jj));
      }
    }
  } // ... copy_local_to_global_matrix(...)

  void copy_local_to_global_vector(const AffinelyDecomposedVectorType& local_vector, const size_t subdomain)
  {
    for (size_t qq = 0; qq < boost::numeric_cast< size_t >(local_vector.num_components()); ++qq) {
      auto comp = find_add_component(*this->rhs_, *local_vector.coefficient(qq), this->test_space().mapper().size());
      copy_local_to_global_vector(*local_vector.component(qq),
                                  subdomain,
                                  *this->rhs_->component(comp));
    }
    if (local_vector.has_affine_part()) {
      ensure_affine_part(*this->rhs_, this->test_space_.mapper().size());
      copy_local_to_global_vector(*local_vector.affine_part(),
                                  subdomain,
                                  *this->rhs()->affine_part());
    }
  } // copy_local_to_global_vector(...)

  void copy_local_to_global_vector(const VectorType& local_vector,
                                   const size_t subdomain,
                                   VectorType& global_vector) const
  {
    const auto& test_mapper = this->test_space_.mapper();
    for (size_t local_ii = 0; local_ii < local_vector.size(); ++local_ii) {
      const size_t global_ii = test_mapper.mapToGlobal(subdomain, local_ii);
      global_vector.add_to_entry(global_ii, local_vector.get_entry(local_ii));
    }
  } // ... copy_local_to_global_vector(...)

  ssize_t find_add_component(AffinelyDecomposedMatrixType& matrix,
                             const Pymor::ParameterFunctional& coefficient,
                             const size_t rows,
                             const size_t cols,
                             const PatternType& ptrn) const
  {
    for (size_t qq = 0; qq < boost::numeric_cast<size_t>(matrix.num_components()); ++qq)
      if (*(matrix.coefficient(qq)) == coefficient)
        return qq;
    return matrix.register_component(coefficient, rows, cols, ptrn);
  } // ... find_add_component(...)

  ssize_t find_add_component(AffinelyDecomposedVectorType& vector,
                             const Pymor::ParameterFunctional& coefficient,
                             const size_t size) const
  {
    for (size_t qq = 0; qq < boost::numeric_cast<size_t>(vector.num_components()); ++qq)
      if (*(vector.coefficient(qq)) == coefficient)
        return qq;
    return vector.register_component(coefficient, size);
  } // ... find_add_component(...)

  void ensure_affine_part(AffinelyDecomposedMatrixType& matrix,
                          const size_t rows,
                          const size_t cols,
                          const PatternType& ptrn) const
  {
    if (!matrix.has_affine_part())
      matrix.register_affine_part(rows, cols, ptrn);
  }

  void ensure_affine_part(AffinelyDecomposedVectorType& vec, const size_t size) const
  {
    if (!vec.has_affine_part())
      vec.register_affine_part(size);
  }

  void add_local_to_global_pattern(const PatternType& local,
                                   const size_t test_subdomain,
                                   const size_t ansatz_subdomain) const
  {
    const auto& test_mapper = this->test_space_.mapper();
    const auto& ansatz_mapper = this->ansatz_space_.mapper();
    for (size_t local_ii = 0; local_ii < local.size(); ++local_ii) {
      const size_t global_ii = test_mapper.mapToGlobal(test_subdomain, local_ii);
      const auto& local_rows = local.inner(local_ii);
      for (const auto& local_jj : local_rows) {
        const size_t global_jj = ansatz_mapper.mapToGlobal(ansatz_subdomain, local_jj);
        this->pattern_->insert(global_ii, global_jj);
      }
    }
  } // ... add_local_to_global_pattern(...)

  LocalGridsProviderType& local_grids_provider_;
  const std::vector<std::string> only_these_products_;
  std::vector<std::shared_ptr<AffinelyDecomposedMatrixType>> local_matrices_;
  std::vector<std::shared_ptr<AffinelyDecomposedVectorType>> local_vectors_;
  std::vector<std::map<size_t,std::shared_ptr<PatternType>>> inside_outside_patterns_;
  std::vector<std::map<size_t,std::shared_ptr<PatternType>>> outside_inside_patterns_;
  std::vector<std::map<size_t,std::shared_ptr<AffinelyDecomposedMatrixType>>> inside_outside_matrices_;
  std::vector<std::map<size_t,std::shared_ptr<AffinelyDecomposedMatrixType>>> outside_inside_matrices_;
}; // Base


} // namespace internal
} // namespace BlockSwipdg
} // namespace Discretizations
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_BLOCK_SWIPDG_BASE_HH
