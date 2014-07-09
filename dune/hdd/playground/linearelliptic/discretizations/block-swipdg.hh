// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_BLOCK_SWIPDG_HH
#define DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_BLOCK_SWIPDG_HH

#include <memory>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <limits>

#if HAVE_EIGEN
# include <Eigen/Eigenvalues>
#endif

#include <dune/common/timer.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/multiscale/provider.hh>

#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/configtree.hh>
#include <dune/stuff/common/float_cmp.hh>
#include <dune/stuff/grid/layers.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/la/container.hh>
#include <dune/stuff/la/solver.hh>

#include <dune/gdt/spaces/discontinuouslagrange.hh>
#include <dune/gdt/playground/spaces/block.hh>
#include <dune/gdt/playground/localevaluation/swipdg.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/oswaldinterpolation.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/playground/spaces/finitevolume/default.hh>
#include <dune/gdt/playground/spaces/raviartthomas/pdelab.hh>
#include <dune/gdt/playground/operators/fluxreconstruction.hh>
#include <dune/gdt/products/l2.hh>
#include <dune/gdt/products/h1.hh>
#include <dune/gdt/assembler/system.hh>

#include <dune/hdd/linearelliptic/problems/default.hh>
#include <dune/hdd/playground/linearelliptic/problems/zero-boundary.hh>

#include "../../../linearelliptic/discretizations/base.hh"
#include "../../../linearelliptic/discretizations/swipdg.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Discretizations {


// forward, needed in the Traits
template< class GridImp, class RangeFieldImp, int rangeDim, int polynomialOrder = 1
        , Stuff::LA::ChooseBackend la_backend = Stuff::LA::ChooseBackend::istl_sparse >
class BlockSWIPDG;


namespace internal {


template< class GridType, class RangeFieldType, int dimRange, int polOrder, Stuff::LA::ChooseBackend la_backend >
class LocalDiscretizationsContainer
{
  typedef grid::Multiscale::ProviderInterface< GridType > GridProviderType;
public:
  typedef SWIPDG< GridType, Stuff::Grid::ChooseLayer::local, RangeFieldType, dimRange
                , polOrder, GDT::ChooseSpaceBackend::fem, la_backend > DiscretizationType;
  typedef typename DiscretizationType::ProblemType     ProblemType;
  typedef typename DiscretizationType::TestSpaceType   TestSpaceType;
  typedef typename DiscretizationType::AnsatzSpaceType AnsatzSpaceType;

private:
  typedef Problems::ZeroBoundary< ProblemType > FakeProblemType;
  typedef typename DiscretizationType::GridViewType::Intersection IntersectionType;

public:
  LocalDiscretizationsContainer(const GridProviderType& grid_provider,
                                const ProblemType& prob)
    : zero_boundary_problem_(prob)
    , all_dirichlet_boundary_config_(Stuff::Grid::BoundaryInfos::AllDirichlet< IntersectionType >::default_config())
    , all_neumann_boundary_config_(Stuff::Grid::BoundaryInfos::AllNeumann< IntersectionType >::default_config())
    , multiscale_boundary_config_(Stuff::Grid::BoundaryInfoConfigs::IdBased::default_config())
    , local_discretizations_(grid_provider.num_subdomains(), nullptr)
    , local_test_spaces_(grid_provider.num_subdomains(), nullptr)
    , local_ansatz_spaces_(grid_provider.num_subdomains(), nullptr)
  {
    multiscale_boundary_config_["neumann"] = "7";
    for (size_t ss = 0; ss < grid_provider.num_subdomains(); ++ss) {
      local_discretizations_[ss] = std::make_shared< DiscretizationType >(grid_provider,
                                                                          all_neumann_boundary_config_,
                                                                          zero_boundary_problem_,
                                                                          ss);
      local_test_spaces_[ss] = local_discretizations_[ss]->test_space();
      local_ansatz_spaces_[ss] = local_discretizations_[ss]->ansatz_space();
    }
  }

protected:
  const FakeProblemType zero_boundary_problem_;
  const Stuff::Common::ConfigTree all_dirichlet_boundary_config_;
  const Stuff::Common::ConfigTree all_neumann_boundary_config_;
  Stuff::Common::ConfigTree multiscale_boundary_config_;
  std::vector< std::shared_ptr< DiscretizationType > > local_discretizations_;
  std::vector< std::shared_ptr< const TestSpaceType > > local_test_spaces_;
  std::vector< std::shared_ptr< const AnsatzSpaceType > > local_ansatz_spaces_;
}; // class LocalDiscretizationsContainer


template< class GridImp, class RangeFieldImp, int rangeDim, int polynomialOrder, Stuff::LA::ChooseBackend la_backend >
class BlockSWIPDGTraits
  : public ContainerBasedDefaultTraits< typename Stuff::LA::Container< RangeFieldImp, la_backend >::MatrixType,
                                        typename Stuff::LA::Container< RangeFieldImp, la_backend >::VectorType >
{
public:
  typedef BlockSWIPDG< GridImp, RangeFieldImp, rangeDim, polynomialOrder, la_backend > derived_type;
  typedef GridImp GridType;
  typedef RangeFieldImp     RangeFieldType;
  static const unsigned int dimRange = rangeDim;
  static const unsigned int polOrder = polynomialOrder;
private:
  friend class BlockSWIPDG< GridImp, RangeFieldImp, rangeDim, polynomialOrder, la_backend >;
  typedef grid::Multiscale::ProviderInterface< GridType > GridProviderType;
  typedef LocalDiscretizationsContainer< GridType, RangeFieldType, dimRange, polOrder, la_backend >
      LocalDiscretizationsContainerType;
  typedef typename LocalDiscretizationsContainerType::TestSpaceType   LocalTestSpaceType;
  typedef typename LocalDiscretizationsContainerType::AnsatzSpaceType LocalAnsatzSpaceType;
public:
  typedef GDT::Spaces::Block< LocalTestSpaceType >   TestSpaceType;
  typedef GDT::Spaces::Block< LocalAnsatzSpaceType > AnsatzSpaceType;
  typedef typename TestSpaceType::GridViewType GridViewType;
}; // class BlockSWIPDGTraits


template< class FunctionType, class EntityType, int dimRange, int dimRangeCols >
class Minimum
{
  static_assert(AlwaysFalse< FunctionType >::value, "Not implemented for these dimensions!");
};

template< class FunctionType, class EntityType >
class Minimum< FunctionType, EntityType, 1, 1 >
{
  typedef typename FunctionType::RangeFieldType RangeFieldType;

public:
  /**
   * We try to find the minimum of a polynomial of given order by evaluating it at the points of a quadrature that
   * would integrate this polynomial exactly.
   * \todo These are just some heuristics and should be replaced by something proper.
   */
  static RangeFieldType compute(const FunctionType& function, const EntityType& entity)
  {
    typename FunctionType::RangeType tmp_value(0);
    RangeFieldType ret = std::numeric_limits< RangeFieldType >::max();
    const auto local_function = function.local_function(entity);
    const size_t ord = local_function->order();
    assert(ord < std::numeric_limits< int >::max());
    const auto& quadrature = QuadratureRules< typename FunctionType::DomainFieldType
                                            , FunctionType::dimDomain >::rule(entity.type(), int(ord));
    const auto quad_point_it_end = quadrature.end();
    for (auto quad_point_it = quadrature.begin(); quad_point_it != quad_point_it_end; ++quad_point_it) {
      local_function->evaluate(quad_point_it->position(), tmp_value);
      ret = std::min(ret, tmp_value[0]);
    }
    return ret;
  } // ... compute(...)
}; // class Minimum< ..., 1, 1 >


template< class FunctionType, class EntityType, int dimDomain >
class Minimum< FunctionType, EntityType, dimDomain, dimDomain >
{
  typedef typename FunctionType::RangeFieldType RangeFieldType;

public:
  static RangeFieldType compute(const FunctionType& function, const EntityType& entity)
  {
#if !HAVE_EIGEN
    static_assert(AlwaysFalse< DT >::value, "You are missing eigen!");
#endif
    const auto local_function = function.local_function(entity);
    assert(local_function->order() == 0);
    const auto& reference_element = ReferenceElements< typename FunctionType::DomainFieldType
                                                     , FunctionType::dimDomain >::general(entity.type());
    const Stuff::LA::EigenDenseMatrix< RangeFieldType >
        tensor = local_function->evaluate(reference_element.position(0, 0));
    ::Eigen::EigenSolver< typename Stuff::LA::EigenDenseMatrix< RangeFieldType >::BackendType >
        eigen_solver(tensor.backend());
    assert(eigen_solver.info() == ::Eigen::Success);
    const auto eigenvalues = eigen_solver.eigenvalues(); // <- this should be an Eigen vector of std::complex
    RangeFieldType min_ev = std::numeric_limits< RangeFieldType >::max();
    for (size_t ii = 0; ii < eigenvalues.size(); ++ii) {
      // assert this is real
      assert(std::abs(eigenvalues[ii].imag()) < 1e-15);
      // assert that this eigenvalue is positive
      const RangeFieldType eigenvalue = eigenvalues[ii].real();
      assert(eigenvalue > 1e-15);
      min_ev = std::min(min_ev, eigenvalue);
    }
    return min_ev;
  } // ... compute(...)
}; // class Minimum< ..., dimDomain, dimDomain >


template< class FunctionType, class EntityType >
static typename FunctionType::RangeFieldType compute_minimum(const FunctionType& function, const EntityType& entity)
{
  static const int dimRange = FunctionType::dimRange;
  static const int dimRangeCols = FunctionType::dimRangeCols;
  return Minimum< FunctionType, EntityType, dimRange, dimRangeCols >::compute(function, entity);
} // ... compute_minimum(...)


} // namespace internal


/**
 * \attention The given problem is replaced by a Problems::ZeroBoundary.
 * \attention The given boundary info config is replaced by a Stuff::Grid::BoundaryInfos::AllDirichlet.
 * \attention The boundary info for the local oversampled discretizations is hardwired to dirichlet zero atm!
 */
template< class GridImp, class RangeFieldImp, int rangeDim, int polynomialOrder, Stuff::LA::ChooseBackend la_backend >
class BlockSWIPDG
  : internal::LocalDiscretizationsContainer< GridImp, RangeFieldImp, rangeDim, polynomialOrder, la_backend >
  , public ContainerBasedDefault< internal::BlockSWIPDGTraits< GridImp, RangeFieldImp, rangeDim
                                                             , polynomialOrder, la_backend > >

{
  typedef internal::LocalDiscretizationsContainer< GridImp, RangeFieldImp, rangeDim, polynomialOrder, la_backend >
    LocalDiscretizationsBaseType;
  typedef ContainerBasedDefault< internal::BlockSWIPDGTraits< GridImp, RangeFieldImp, rangeDim
                                                            , polynomialOrder, la_backend > > BaseType;
public:
  typedef internal::BlockSWIPDGTraits< GridImp, RangeFieldImp, rangeDim, polynomialOrder, la_backend > Traits;
  using typename BaseType::ProblemType;
  using typename BaseType::TestSpaceType;
  using typename BaseType::AnsatzSpaceType;
  using typename BaseType::EntityType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::MatrixType;
  using typename BaseType::VectorType;
  using typename BaseType::OperatorType;
  using typename BaseType::ProductType;
  using typename BaseType::FunctionalType;

  static const unsigned int dimDomain = BaseType::dimDomain;
  static const unsigned int dimRange = BaseType::dimRange;

  typedef grid::Multiscale::ProviderInterface< GridImp > GridProviderType;
  typedef typename GridProviderType::MsGridType MsGridType;

  typedef typename Traits::LocalDiscretizationsContainerType::DiscretizationType LocalDiscretizationType;
  typedef typename TestSpaceType::PatternType PatternType;

private:
  using typename BaseType::AffinelyDecomposedMatrixType;
  using typename BaseType::AffinelyDecomposedVectorType;
  typedef Pymor::LA::AffinelyDecomposedConstContainer< MatrixType > AffinelyDecomposedConstMatrixType;
  typedef Pymor::LA::AffinelyDecomposedConstContainer< VectorType > AffinelyDecomposedConstVectorType;

public:
  typedef typename LocalDiscretizationType::ProblemType LocalProblemType;
  typedef typename LocalDiscretizationType::ProductType LocalProductType;

  static std::string static_id()
  {
    return DiscretizationInterface< Traits >::static_id() + ".block-swipdg";
  }

  BlockSWIPDG(const GridProviderType& grid_provider,
              const Stuff::Common::ConfigTree& /*bound_inf_cfg*/,
              const ProblemType& prob)
    : LocalDiscretizationsBaseType(grid_provider, prob)
    , BaseType(std::make_shared< TestSpaceType >(grid_provider.ms_grid(), this->local_test_spaces_),
               std::make_shared< AnsatzSpaceType >(grid_provider.ms_grid(), this->local_ansatz_spaces_),
               this->all_dirichlet_boundary_config_,
               this->zero_boundary_problem_)
    , grid_provider_(grid_provider)
    , ms_grid_(grid_provider.ms_grid())
    , pattern_(BaseType::test_space()->mapper().size())
    , local_matrices_(ms_grid_->size())
    , local_vectors_(ms_grid_->size())
    , inside_outside_patterns_(ms_grid_->size())
    , outside_inside_patterns_(ms_grid_->size())
    , inside_outside_matrices_(ms_grid_->size())
    , outside_inside_matrices_(ms_grid_->size())
//    , local_boundary_patterns_(ms_grid_->size())
//    , local_coupling_patterns_(ms_grid_->size())
//    , inside_outside_coupling_patterns_(ms_grid_->size())
//    , outside_inside_coupling_patterns_(ms_grid_->size())
//    , boundary_matrices_(ms_grid_->size(), nullptr)
//    , boundary_vectors_(ms_grid_->size(), nullptr)
//    , local_coupling_matrices_(ms_grid_->size(), nullptr)
//    , inside_outside_coupling_matrices_(ms_grid_->size())
//    , outside_inside_coupling_matrices_(ms_grid_->size())
  {
    // in case of parametric diffusion tensor everything is too complicated
    if (this->problem_.diffusion_tensor()->parametric())
      DUNE_THROW_COLORFULLY(NotImplemented, "The diffusion tensor must not be parametric!");
    if (!this->problem_.diffusion_tensor()->has_affine_part())
      DUNE_THROW_COLORFULLY(Stuff::Exceptions::wrong_input_given, "The diffusion tensor must not be empty!");
  } // BlockSWIPDG(...)

  const std::vector< std::shared_ptr< LocalDiscretizationType > >& local_discretizations() const
  {
    return this->local_discretizations_;
  }

  void init(std::ostream& out = Stuff::Common::Logger().devnull(), const std::string prefix = "")
  {
    if (!this->container_based_initialized_) {
      const size_t subdomains = ms_grid_->size();
      // walk the subdomains for the first time
      //   * to initialize the coupling pattern,
      //   * to finalize the global sparsity pattern
      out << prefix << "walking subdomains for the first time... " << std::flush;
      Dune::Timer timer;
      for (size_t ss = 0; ss < subdomains; ++ss) {
        // init the local discretizations (assembles matrices and patterns)
        this->local_discretizations_[ss]->init();
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
          DUNE_PYMOR_THROW(Pymor::Exception::requirements_not_met,
                           "The local operator should have " << this->problem().diffusion_factor()->num_components()
                           << " components (but has only " << local_operator.num_components() << ")!");
        for (DUNE_STUFF_SSIZE_T qq = 0; qq < this->problem().diffusion_factor()->num_components(); ++qq)
          local_matrices_[ss]->register_component(new MatrixType(
              local_operator.component(qq).container()->backend()),
              this->problem().diffusion_factor()->coefficient(qq));
        // * and the vectors
        const auto local_functional = this->local_discretizations_[ss]->get_rhs();
        local_vectors_[ss] = std::make_shared< AffinelyDecomposedVectorType >();
        for (size_t qq = 0; qq < local_functional.num_components(); ++qq)
          local_vectors_[ss]->register_component(new VectorType(*(local_functional.component(qq).container())),
                                                 new Pymor::ParameterFunctional(local_functional.coefficient(qq)));
        if (local_functional.has_affine_part())
          local_vectors_[ss]->register_affine_part(new VectorType(*(local_functional.affine_part().container())));

        // create and copy the local patterns
        add_local_to_global_pattern(this->local_discretizations_[ss]->pattern(), ss, ss, pattern_);
        const auto& inner_test_space = *(this->local_discretizations_[ss]->test_space());
        const auto& inner_ansatz_space = *(this->local_discretizations_[ss]->ansatz_space());
        // walk the neighbors
        for (const size_t& nn : ms_grid_->neighborsOf(ss)) {
          // visit each coupling only once (assemble primally)
          if (ss < nn) {
            const auto& outer_test_space = *(this->local_discretizations_[nn]->test_space());
            const auto& outer_ansatz_space = *(this->local_discretizations_[nn]->ansatz_space());
            const auto& inside_outside_grid_part = *(ms_grid_->couplingGridPart(ss, nn));
            const auto& outside_inside_grid_part = *(ms_grid_->couplingGridPart(nn, ss));
            // create the coupling patterns
            auto inside_outside_pattern = std::make_shared< PatternType >(
                  inner_test_space.compute_face_pattern(inside_outside_grid_part, outer_ansatz_space));
            inside_outside_patterns_[ss].insert(std::make_pair(nn, inside_outside_pattern));
            auto outside_inside_pattern = std::make_shared< PatternType >(
                  outer_test_space.compute_face_pattern(outside_inside_grid_part, inner_ansatz_space));
            outside_inside_patterns_[nn].insert(std::make_pair(ss, outside_inside_pattern));
            // and copy them
            add_local_to_global_pattern(*inside_outside_pattern,  ss, nn, pattern_);
            add_local_to_global_pattern(*outside_inside_pattern,  nn, ss, pattern_);
          } // visit each coupling only once (assemble primaly)
        } // walk the neighbors
      } // walk the subdomains for the first time
      out<< "done (took " << timer.elapsed() << " sek)" << std::endl;

      // walk the subdomains for the second time
      //   * to assemble the boundary matrices and vectors and
      //   * to assemble the coupling matrices
      out << prefix << "walking subdomains for the second time... " << std::flush;
      for (size_t ss = 0; ss < ms_grid_->size(); ++ss) {
        const auto& inner_test_mapper = this->local_discretizations_[ss]->test_space()->mapper();
        const auto& inner_ansatz_mapper = this->local_discretizations_[ss]->ansatz_space()->mapper();
        if (ms_grid_->boundary(ss))
          assemble_boundary_contributions(ss);
        // walk the neighbors
        for (const size_t& nn : ms_grid_->neighborsOf(ss)) {
          // visit each coupling only once (assemble primaly)
          if (ss < nn) {
            const auto& outer_test_mapper = this->local_discretizations_[nn]->test_space()->mapper();
            const auto& outer_ansatz_mapper = this->local_discretizations_[nn]->ansatz_space()->mapper();
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
            for (DUNE_STUFF_SSIZE_T qq = 0; qq < this->problem().diffusion_factor()->num_components(); ++qq) {
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
      out<< "done (took " << timer.elapsed() << " sek)" << std::endl;

      // build global containers
      build_global_containers();

      // products
      typedef typename MsGridType::GlobalGridPartType::GridViewType GlobalGridViewType;
      const auto& global_grid_view = ms_grid_->globalGridPart()->gridView();
      GDT::SystemAssembler< TestSpaceType, GlobalGridViewType, AnsatzSpaceType > system_assembler(*this->test_space(),
                                                                                                  *this->ansatz_space(),
                                                                                                  global_grid_view);

      // * L2
      typedef GDT::Products::L2Assemblable< MatrixType, TestSpaceType, GlobalGridViewType, AnsatzSpaceType >
          L2ProductType;
      auto l2_product_matrix = std::make_shared< AffinelyDecomposedMatrixType >();
      l2_product_matrix->register_affine_part(this->test_space()->mapper().size(),
                                              this->ansatz_space()->mapper().size(),
                                              L2ProductType::pattern(*this->test_space(),
                                                                     *this->ansatz_space()));
      L2ProductType l2_product(*(l2_product_matrix->affine_part()),
                               *this->test_space(),
                               global_grid_view,
                               *this->ansatz_space());
      system_assembler.add(l2_product);
      // * H1 semi
      typedef GDT::Products::H1SemiAssemblable< MatrixType, TestSpaceType, GlobalGridViewType, AnsatzSpaceType >
          H1ProductType;
      auto h1_product_matrix = std::make_shared< AffinelyDecomposedMatrixType >();
      h1_product_matrix->register_affine_part(this->test_space()->mapper().size(),
                                              this->ansatz_space()->mapper().size(),
                                              H1ProductType::pattern(*this->test_space(),
                                                                     *this->ansatz_space()));
      H1ProductType h1_product(*(h1_product_matrix->affine_part()),
                               *this->test_space(),
                               global_grid_view,
                               *this->ansatz_space());
      system_assembler.add(h1_product);
      system_assembler.assemble();

      // finalize
      this->inherit_parameter_type(*(this->matrix_), "lhs");
      this->inherit_parameter_type(*(this->rhs_), "rhs");
      this->products_.insert(std::make_pair("l2", l2_product_matrix));
      this->products_.insert(std::make_pair("h1_semi", h1_product_matrix));

      this->container_based_initialized_ = true;
    } // if (!this->container_based_initialized_)
  } // ... init(...)

  DUNE_STUFF_SSIZE_T num_subdomains() const
  {
    return ms_grid_->size();
  }

  std::vector< DUNE_STUFF_SSIZE_T > neighbouring_subdomains(const DUNE_STUFF_SSIZE_T ss) const
  {
    if (ss < 0 || ss >= num_subdomains())
      DUNE_THROW(Stuff::Exceptions::index_out_of_range,
                 "0 <= ss < num_subdomains() = " << num_subdomains() << " is not true for ss = " << ss << "!");
    const auto set_of_neighbours = ms_grid_->neighborsOf(ss);
    return std::vector< DUNE_STUFF_SSIZE_T >(set_of_neighbours.begin(), set_of_neighbours.end());
  }

  VectorType localize_vector(const VectorType& global_vector, const size_t ss) const
  {
    if (ss >= num_subdomains())
      DUNE_THROW(Stuff::Exceptions::index_out_of_range,
                 "0 <= ss < num_subdomains() = " << num_subdomains() << " is not true for ss = " << ss << "!");
    if (global_vector.size() != this->ansatz_space()->mapper().size())
      DUNE_THROW(Stuff::Exceptions::index_out_of_range,
                 "The size() of global_vector (" << global_vector.dim()
                 << ") does not match the size() of the ansatz space (" << this->ansatz_space()->mapper().size() << ")!");
    assert(ss < this->local_discretizations_.size());
    VectorType local_vector = this->local_discretizations_[ss]->create_vector();
    for (size_t ii = 0; ii < local_vector.size(); ++ii)
      local_vector.set_entry(ii, global_vector.get_entry(this->ansatz_space()->mapper().mapToGlobal(ss, ii)));
    return local_vector;
  } // ... localize_vetor(...)

  VectorType* localize_vector_and_return_ptr(const VectorType& global_vector, const DUNE_STUFF_SSIZE_T ss) const
  {
    assert(ss >= 0);
    assert(ss < std::numeric_limits< size_t >::max());
    return new VectorType(localize_vector(global_vector, size_t(ss)));
  }

  ProductType get_local_product(const size_t ss, const std::string id) const
  {
    if (ss >= num_subdomains())
      DUNE_THROW(Stuff::Exceptions::index_out_of_range,
                 "0 <= ss < num_subdomains() = " << num_subdomains() << " is not true for ss = " << ss << "!");
    return this->local_discretizations_[ss]->get_product(id);
  }

  ProductType* get_local_product_and_return_ptr(const DUNE_STUFF_SSIZE_T ss, const std::string id) const
  {
    assert(ss >= 0);
    assert(ss < std::numeric_limits< size_t >::max());
    return new ProductType(get_local_product(size_t(ss), id));
  }

  OperatorType get_local_operator(const size_t ss) const
  {
    if (ss >= num_subdomains())
      DUNE_THROW(Stuff::Exceptions::index_out_of_range,
                 "0 <= ss < num_subdomains() = " << num_subdomains() << " is not true for ss = " << ss << "!");
    assert(ss < int(local_matrices_.size()));
    return OperatorType(*(local_matrices_[ss]));
  }

  OperatorType* get_local_operator_and_return_ptr(const DUNE_STUFF_SSIZE_T ss) const
  {
    assert(ss >= 0);
    assert(ss < std::numeric_limits< size_t >::max());
    return new OperatorType(get_local_operator(size_t(ss)));
  }

  OperatorType get_coupling_operator(const size_t ss, const size_t nn) const
  {
    if (ss >= num_subdomains())
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

  OperatorType* get_coupling_operator_and_return_ptr(const DUNE_STUFF_SSIZE_T ss, const DUNE_STUFF_SSIZE_T nn) const
  {
    assert(ss >= 0);
    assert(ss < std::numeric_limits< size_t >::max());
    assert(nn >= 0);
    assert(nn < std::numeric_limits< size_t >::max());
    return new OperatorType(get_coupling_operator(size_t(ss), size_t(nn)));
  }

  FunctionalType get_local_functional(const size_t ss) const
  {
    if (ss >= num_subdomains())
      DUNE_PYMOR_THROW(Pymor::Exception::index_out_of_range,
                       "0 <= ss < num_subdomains() = " << num_subdomains() << " is not true for ss = " << ss << "!");
    assert(ss < local_vectors_.size());
    return FunctionalType(*(local_vectors_[ss]));
  }

  FunctionalType* get_local_functional_and_return_ptr(const DUNE_STUFF_SSIZE_T ss) const
  {
    assert(ss >= 0);
    assert(ss < std::numeric_limits< size_t >::max());
    return new FunctionalType(get_local_functional(size_t(ss)));
  }

  static std::vector< std::string > available_estimators()
  {
    return { "eta" } ;
  }

  RangeFieldType estimate(const VectorType& vector, const std::string type) const
  {
    if (type == "eta") {
      if (this->parametric())
        DUNE_THROW(NotImplemented, "Yet");
      return estimate(vector);
    } else
      DUNE_THROW(Stuff::Exceptions::you_are_using_this_wrong,
                 "type '" << type << "' is not one of available_estimators()!");
  } // ... estimate(...)

  RangeFieldType estimate(const VectorType& vector,
                          const Pymor::Parameter mu = Pymor::Parameter(),
                          const Pymor::Parameter mu_fixed = Pymor::Parameter(),
                          const Pymor::Parameter mu_minimizing = Pymor::Parameter()) const
  {
    if (mu.type() != this->parameter_type())
      DUNE_THROW(Pymor::Exceptions::wrong_parameter_type,
                 "mu is " << mu.type() << ", should be " << this->parameter_type() << "!");
    if (mu_fixed.type() != this->parameter_type())
      DUNE_THROW(Pymor::Exceptions::wrong_parameter_type,
                 "mu_fixed is " << mu_fixed.type() << ", should be " << this->parameter_type() << "!");
    if (mu_minimizing.type() != this->parameter_type())
      DUNE_THROW(Pymor::Exceptions::wrong_parameter_type,
                 "mu_minimizing is " << mu_minimizing.type() << ", should be " << this->parameter_type() << "!");

    using namespace GDT;

    typedef typename MsGridType::GlobalGridViewType GlobalGridViewType;
    const auto grid_view = ms_grid_->globalGridView();
    typedef ConstDiscreteFunction< AnsatzSpaceType, VectorType > ConstDiscreteFunctionType;
    const ConstDiscreteFunctionType discrete_solution(*this->ansatz_space(), vector);

    VectorType oswald_interpolation_vector(this->ansatz_space()->mapper().size());
    typedef DiscreteFunction< AnsatzSpaceType, VectorType > DiscreteFunctionType;
    DiscreteFunctionType oswald_interpolation(*this->ansatz_space(), oswald_interpolation_vector);
    const Operators::OswaldInterpolation< GlobalGridViewType > oswald_interpolation_operator(*grid_view);
    oswald_interpolation_operator.apply(discrete_solution, oswald_interpolation);

    typedef Spaces::FiniteVolume::Default< GlobalGridViewType, RangeFieldType, 1, 1 > P0SpaceType;
    const P0SpaceType p0_space(grid_view);
    VectorType p0_force_vector(p0_space.mapper().size());
    DiscreteFunction< P0SpaceType, VectorType > p0_force(p0_space, p0_force_vector);
    Operators::Projection< GlobalGridViewType > projection_operator(*grid_view);
    const auto force = this->problem().force()->with_mu(this->map_parameter(mu, "force"));
    projection_operator.apply(*force, p0_force);

    typedef Spaces::RaviartThomas::PdelabBased< GlobalGridViewType, 0, RangeFieldType, dimDomain > RTN0SpaceType;
    const RTN0SpaceType rtn0_space(grid_view);
    VectorType diffusive_flux_vector(rtn0_space.mapper().size());
    typedef DiscreteFunction< RTN0SpaceType, VectorType > RTN0DiscreteFunctionType;
    RTN0DiscreteFunctionType diffusive_flux(rtn0_space, diffusive_flux_vector);

    typedef typename ProblemType::DiffusionFactorType::NonparametricType DiffusionFactorType;
    const auto diffusion_factor = this->problem().diffusion_factor()->with_mu(this->map_parameter(mu,
                                                                                                  "diffusion_factor"));
    const auto diffusion_factor_minimizing = this->problem().diffusion_factor()->with_mu(this->map_parameter(mu_minimizing,
                                                                                                             "diffusion_factor"));
    typedef typename ProblemType::DiffusionTensorType::NonparametricType DiffusionTensorType;
    const auto diffusion_tensor = this->problem().diffusion_tensor()->with_mu(this->map_parameter(mu,
                                                                                                  "diffusion_tensor"));
    const auto diffusion_tensor_minimizing = this->problem().diffusion_tensor()->with_mu(this->map_parameter(mu_minimizing,
                                                                                                             "diffusion_tensor"));
    const Operators::DiffusiveFluxReconstruction< GlobalGridViewType, DiffusionFactorType, DiffusionTensorType >
      diffusive_flux_reconstruction(*grid_view, *diffusion_factor, *diffusion_tensor);
    diffusive_flux_reconstruction.apply(discrete_solution, diffusive_flux);

    const LocalOperator::Codim0Integral< LocalEvaluation::Elliptic< DiffusionFactorType, DiffusionTensorType > >
      local_eta_nc_product(1, *diffusion_factor, *diffusion_tensor);
    const auto eta_nc_difference = discrete_solution - oswald_interpolation;

    typedef Stuff::Functions::Constant< EntityType, DomainFieldType, dimDomain, RangeFieldType, 1 >
        ConstantFunctionType;
    typedef typename ConstantFunctionType::DomainType DomainType;
    const ConstantFunctionType constant_one(1);
    const  LocalOperator::Codim0Integral< LocalEvaluation::Product< ConstantFunctionType > >
      local_eta_r_product(1, constant_one);
    const auto eta_r_difference = *force - p0_force;

    const auto diffusion_factor_fixed = this->problem().diffusion_factor()->with_mu(this->map_parameter(mu_fixed,
                                                                                                        "diffusion_factor"));
    const auto diffusion_tensor_fixed = this->problem().diffusion_tensor()->with_mu(this->map_parameter(mu_fixed,
                                                                                                        "diffusion_tensor"));
    const LocalOperator::Codim0Integral< LocalEvaluation::ESV2007::DiffusiveFluxEstimate< DiffusionFactorType
                                                                                        , RTN0DiscreteFunctionType
                                                                                        , DiffusionTensorType > >
      local_eta_df_product(1, *diffusion_factor_fixed, *diffusion_tensor_fixed, diffusive_flux);

    std::vector< DynamicMatrix< RangeFieldType > >
        tmp_matrices(std::max(std::max(local_eta_nc_product.numTmpObjectsRequired(),
                                       local_eta_r_product.numTmpObjectsRequired()),
                              local_eta_df_product.numTmpObjectsRequired()),
                     DynamicMatrix< RangeFieldType >(1, 1, 0.0));
    DynamicMatrix< RangeFieldType > local_result_matrix(1, 1, 0.0);

    // walk the subdomains
    double eta_nc_squared = 0.0;
    double eta_r_squared = 0.0;
    double eta_df_squared = 0.0;
    for (size_t subdomain = 0; subdomain < ms_grid_->size(); ++subdomain) {
      const auto local_grid_part = ms_grid_->localGridPart(subdomain);

      // some storage
      std::vector< DomainType > vertices;
      RangeFieldType min_diffusion_value = std::numeric_limits< RangeFieldType >::max();

      // walk the local grid
      double eta_nc_T_squared = 0.0;
      double eta_r_T_squared = 0.0;
      double eta_df_T_squared = 0.0;
      const auto entity_it_end = local_grid_part->template end< 0 >();
      for (auto entity_it = local_grid_part->template begin< 0 >(); entity_it != entity_it_end; ++entity_it) {
        const auto& entity = *entity_it;

        // collect all vertices
        for (int cc = 0; cc < entity.template count< dimDomain >(); ++cc)
          vertices.push_back(entity.template subEntity< dimDomain >(cc)->geometry().center());

        // compute minimum diffusion
        const RangeFieldType min_diffusion_value_entity
            =   internal::compute_minimum(*diffusion_factor_minimizing, entity)
              * internal::compute_minimum(*diffusion_tensor_minimizing, entity);
        min_diffusion_value = std::min(min_diffusion_value, min_diffusion_value_entity);

        local_result_matrix *= 0.0;
        const auto local_eta_nc_difference = eta_nc_difference.local_function(entity);
        local_eta_nc_product.apply(*local_eta_nc_difference, *local_eta_nc_difference, local_result_matrix, tmp_matrices);
        assert(local_result_matrix.rows() >= 1);
        assert(local_result_matrix.cols() >= 1);
        eta_nc_T_squared += local_result_matrix[0][0];

        local_result_matrix *= 0.0;
        const auto local_eta_r_difference = eta_r_difference.local_function(entity);
        local_eta_r_product.apply(*local_eta_r_difference, *local_eta_r_difference, local_result_matrix, tmp_matrices);
        assert(local_result_matrix.rows() >= 1);
        assert(local_result_matrix.cols() >= 1);
        eta_r_T_squared += local_result_matrix[0][0];

        local_result_matrix *= 0.0;
        const auto local_discrete_solution = discrete_solution.local_function(entity);
        local_eta_df_product.apply(*local_discrete_solution, *local_discrete_solution, local_result_matrix, tmp_matrices);
        assert(local_result_matrix.rows() >= 1);
        assert(local_result_matrix.cols() >= 1);
        eta_df_T_squared += local_result_matrix[0][0];
      } // walk the local grid

      // compute diameter
      DomainFieldType diameter(0);
      DomainType tmp_diff;
      for (size_t cc = 0; cc < vertices.size(); ++cc) {
        const auto& vertex = vertices[cc];
        for (size_t dd = cc + 1; dd < vertices.size(); ++dd) {
          const auto& other_vertex = vertices[dd];
          tmp_diff = vertex - other_vertex;
          diameter = std::max(diameter, tmp_diff.two_norm());
        }
      }

      const RangeFieldType poincare_constant = 1.0 / (M_PIl * M_PIl);
      assert(min_diffusion_value > 0.0);

      eta_nc_squared += eta_nc_T_squared;
      eta_r_squared += ((poincare_constant * diameter * diameter) / min_diffusion_value) * eta_r_T_squared;
      eta_df_squared += eta_df_T_squared;
    } // walk the subdomains

    assert(!this->problem().diffusion_tensor()->parametric());
    assert(!this->problem().force()->parametric());
    assert(!this->problem().dirichlet()->parametric());
    assert(!this->problem().neumann()->parametric());
    const double alpha = this->problem().diffusion_factor()->alpha(mu, mu_fixed);
    const double gamma = this->problem().diffusion_factor()->gamma(mu, mu_fixed);
    assert(alpha > 0.0);
    assert(gamma > 0.0);
    const double parameter_factor = std::max(std::sqrt(gamma), 1.0/std::sqrt(alpha));

    return std::sqrt(eta_nc_squared) + std::sqrt(eta_r_squared) + parameter_factor * std::sqrt(eta_df_squared);
  } // ... estimate(...)

  /**
   * \return eta_T^2 := eta_nc^T^2 + eta_r^T^2 + eta_df^T^2
   */
  RangeFieldType estimate_local(const VectorType& vector, const size_t subdomain) const
  {
    if (this->parametric()) DUNE_THROW(NotImplemented, "No time for that atm!");

    using namespace GDT;

    typedef typename MsGridType::GlobalGridViewType GlobalGridViewType;
    const auto grid_view = ms_grid_->globalGridView();
    typedef ConstDiscreteFunction< AnsatzSpaceType, VectorType > ConstDiscreteFunctionType;
    const ConstDiscreteFunctionType discrete_solution(*this->ansatz_space(), vector);

    VectorType oswald_interpolation_vector(this->ansatz_space()->mapper().size());
    typedef DiscreteFunction< AnsatzSpaceType, VectorType > DiscreteFunctionType;
    DiscreteFunctionType oswald_interpolation(*this->ansatz_space(), oswald_interpolation_vector);
    const Operators::OswaldInterpolation< GlobalGridViewType > oswald_interpolation_operator(*grid_view);
    oswald_interpolation_operator.apply(discrete_solution, oswald_interpolation);

    typedef Spaces::FiniteVolume::Default< GlobalGridViewType, RangeFieldType, 1, 1 > P0SpaceType;
    const P0SpaceType p0_space(grid_view);
    VectorType p0_force_vector(p0_space.mapper().size());
    DiscreteFunction< P0SpaceType, VectorType > p0_force(p0_space, p0_force_vector);
    Operators::Projection< GlobalGridViewType > projection_operator(*grid_view);
    const auto& force = *(this->problem().force());
    assert(!force.parametric());
    assert(force.has_affine_part());
    projection_operator.apply(*force.affine_part(), p0_force);

    typedef Spaces::RaviartThomas::PdelabBased< GlobalGridViewType, 0, RangeFieldType, dimDomain > RTN0SpaceType;
    const RTN0SpaceType rtn0_space(grid_view);
    VectorType diffusive_flux_vector(rtn0_space.mapper().size());
    typedef DiscreteFunction< RTN0SpaceType, VectorType > RTN0DiscreteFunctionType;
    RTN0DiscreteFunctionType diffusive_flux(rtn0_space, diffusive_flux_vector);

    typedef typename ProblemType::DiffusionFactorType::NonparametricType DiffusionFactorType;
    const auto& diffusion_factor = *(this->problem().diffusion_factor());
    assert(!diffusion_factor.parametric());
    assert(diffusion_factor.has_affine_part());
    typedef typename ProblemType::DiffusionTensorType::NonparametricType DiffusionTensorType;
    const auto& diffusion_tensor = *(this->problem().diffusion_tensor());
    assert(!diffusion_tensor.parametric());
    assert(diffusion_tensor.has_affine_part());
    const Operators::DiffusiveFluxReconstruction< GlobalGridViewType, DiffusionFactorType, DiffusionTensorType >
      diffusive_flux_reconstruction(*grid_view, *diffusion_factor.affine_part(), *diffusion_tensor.affine_part());
    diffusive_flux_reconstruction.apply(discrete_solution, diffusive_flux);

    const LocalOperator::Codim0Integral< LocalEvaluation::Elliptic< DiffusionFactorType, DiffusionTensorType > >
      local_eta_nc_product(1, *diffusion_factor.affine_part(), *diffusion_tensor.affine_part());
    const auto eta_nc_difference = discrete_solution - oswald_interpolation;

    typedef Stuff::Functions::Constant< EntityType, DomainFieldType, dimDomain, RangeFieldType, 1 >
        ConstantFunctionType;
    typedef typename ConstantFunctionType::DomainType DomainType;
    const ConstantFunctionType constant_one(1);
    const  LocalOperator::Codim0Integral< LocalEvaluation::Product< ConstantFunctionType > >
      local_eta_r_product(1, constant_one);
    const auto eta_r_difference = *force.affine_part() - p0_force;

    const LocalOperator::Codim0Integral< LocalEvaluation::ESV2007::DiffusiveFluxEstimate< DiffusionFactorType
                                                                                        , RTN0DiscreteFunctionType
                                                                                        , DiffusionTensorType > >
      local_eta_df_product(1, *diffusion_factor.affine_part(), *diffusion_tensor.affine_part(), diffusive_flux);

    std::vector< DynamicMatrix< RangeFieldType > >
        tmp_matrices(std::max(std::max(local_eta_nc_product.numTmpObjectsRequired(),
                                       local_eta_r_product.numTmpObjectsRequired()),
                              local_eta_df_product.numTmpObjectsRequired()),
                     DynamicMatrix< RangeFieldType >(1, 1, 0.0));
    DynamicMatrix< RangeFieldType > local_result_matrix(1, 1, 0.0);

    // estimate
    assert(subdomain < ms_grid_->size());
    double eta_nc_squared = 0.0;
    double eta_r_squared = 0.0;
    double eta_df_squared = 0.0;
//    for (size_t subdomain = 0; subdomain < ms_grid_->size(); ++subdomain) {
      const auto local_grid_part = ms_grid_->localGridPart(subdomain);

      // some storage
      std::vector< DomainType > vertices;
      RangeFieldType min_diffusion_value = std::numeric_limits< RangeFieldType >::max();

      // walk the local grid
      double eta_nc_T_squared = 0.0;
      double eta_r_T_squared = 0.0;
      double eta_df_T_squared = 0.0;
      const auto entity_it_end = local_grid_part->template end< 0 >();
      for (auto entity_it = local_grid_part->template begin< 0 >(); entity_it != entity_it_end; ++entity_it) {
        const auto& entity = *entity_it;

        // collect all vertices
        for (int cc = 0; cc < entity.template count< dimDomain >(); ++cc)
          vertices.push_back(entity.template subEntity< dimDomain >(cc)->geometry().center());

        // compute minimum diffusion
        const RangeFieldType min_diffusion_value_entity
            =   internal::compute_minimum(*diffusion_factor.affine_part(), entity)
              * internal::compute_minimum(*diffusion_tensor.affine_part(), entity);
        min_diffusion_value = std::min(min_diffusion_value, min_diffusion_value_entity);

        local_result_matrix *= 0.0;
        const auto local_eta_nc_difference = eta_nc_difference.local_function(entity);
        local_eta_nc_product.apply(*local_eta_nc_difference, *local_eta_nc_difference, local_result_matrix, tmp_matrices);
        assert(local_result_matrix.rows() >= 1);
        assert(local_result_matrix.cols() >= 1);
        eta_nc_T_squared += local_result_matrix[0][0];

        local_result_matrix *= 0.0;
        const auto local_eta_r_difference = eta_r_difference.local_function(entity);
        local_eta_r_product.apply(*local_eta_r_difference, *local_eta_r_difference, local_result_matrix, tmp_matrices);
        assert(local_result_matrix.rows() >= 1);
        assert(local_result_matrix.cols() >= 1);
        eta_r_T_squared += local_result_matrix[0][0];

        local_result_matrix *= 0.0;
        const auto local_discrete_solution = discrete_solution.local_function(entity);
        local_eta_df_product.apply(*local_discrete_solution, *local_discrete_solution, local_result_matrix, tmp_matrices);
        assert(local_result_matrix.rows() >= 1);
        assert(local_result_matrix.cols() >= 1);
        eta_df_T_squared += local_result_matrix[0][0];
      } // walk the local grid

      // compute diameter
      DomainFieldType diameter(0);
      DomainType tmp_diff;
      for (size_t cc = 0; cc < vertices.size(); ++cc) {
        const auto& vertex = vertices[cc];
        for (size_t dd = cc + 1; dd < vertices.size(); ++dd) {
          const auto& other_vertex = vertices[dd];
          tmp_diff = vertex - other_vertex;
          diameter = std::max(diameter, tmp_diff.two_norm());
        }
      }

      const RangeFieldType poincare_constant = 1.0 / (M_PIl * M_PIl);
      assert(min_diffusion_value > 0.0);

      eta_nc_squared += eta_nc_T_squared;
      eta_r_squared += ((poincare_constant * diameter * diameter) / min_diffusion_value) * eta_r_T_squared;
      eta_df_squared += eta_df_T_squared;
//    } // walk the subdomains
    return eta_nc_squared + eta_r_squared + eta_df_squared;
  } // ... estimate_local(...)

  VectorType* create_ones() const
  {
    return new VectorType(this->ansatz_space()->mapper().size(), 1.0);
  }

  void visualize_information(const std::vector< VectorType >& local_vectors,
                             const std::string& filename,
                             const Pymor::Parameter mu = Pymor::Parameter(),
                             const Pymor::Parameter mu_fixed = Pymor::Parameter()) const
  {
    using namespace GDT;

    if (mu.type() != this->parameter_type())
      DUNE_THROW(Pymor::Exceptions::wrong_parameter_type,
                 "mu is " << mu.type() << ", should be " << this->parameter_type() << "!");
    if (mu_fixed.type() != this->parameter_type())
      DUNE_THROW(Pymor::Exceptions::wrong_parameter_type,
                 "mu_fixed is " << mu_fixed.type() << ", should be " << this->parameter_type() << "!");
    if (local_vectors.size() != ms_grid_->size())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given,
                 "local_vectors is of size " << local_vectors.size() << " and should be of size " << ms_grid_->size());
    VectorType vector(this->ansatz_space()->mapper().size());
    for (size_t ss = 0; ss < ms_grid_->size(); ++ss) {
      if (local_vectors[ss].size() != this->local_discretizations_[ss]->ansatz_space()->mapper().size())
        DUNE_THROW(Stuff::Exceptions::wrong_input_given,
                   "local_vectors[" << ss << "] is of size " << local_vectors[ss].size() << " and should be of size "
                   << this->local_discretizations_[ss]->ansatz_space()->mapper().size());
      if (!local_vectors[ss].valid())
        DUNE_THROW(Stuff::Exceptions::wrong_input_given, "local_vectors[" << ss << "] contains NaN or INF!");
      copy_local_to_global_vector(local_vectors[ss], ss, vector);
    }

    VectorType solution = this->create_vector();
    this->solve(solution, mu);

    // difference
    VectorType difference_vector = solution - vector;
    typedef DiscreteFunction< AnsatzSpaceType, VectorType > DiscreteFunctionType;
    DiscreteFunctionType difference(*this->ansatz_space(), difference_vector);

    // estimator
    auto local_estimates = estimate_local(local_vectors, mu, mu_fixed);

    // energy error
    typedef typename ProblemType::DiffusionTensorType::NonparametricType DiffusionTensorType;
    const auto& diffusion_tensor = this->problem().diffusion_tensor();
    assert(!diffusion_tensor.parametric());
    typedef typename ProblemType::DiffusionFactorType::NonparametricType DiffusionFactorType;
    const auto& diffusion_factor = this->problem().diffusion_factor();
    const auto nonparametric_diffusion_factor = diffusion_factor->with_mu(mu);
    const LocalOperator::Codim0Integral< LocalEvaluation::Elliptic< DiffusionFactorType, DiffusionTensorType > >
        local_energy_product(1, *nonparametric_diffusion_factor, *diffusion_tensor->affine_part());
    std::vector< DynamicMatrix< RangeFieldType > > tmp_matrices(local_energy_product.numTmpObjectsRequired(),
                                                                DynamicMatrix< RangeFieldType >(1, 1, 0.0));
    DynamicMatrix< RangeFieldType > local_result_matrix(1, 1, 0.0);
    const auto grid_view = ms_grid_->globalGridView();
    std::vector< double > energy_error_visualization(grid_view->indexSet().size(0), 0.0);
    std::vector< double > estimator_visualization(grid_view->indexSet().size(0), 0.0);
    const auto entity_it_end = grid_view->template end< 0 >();
    for (auto entity_it = grid_view->template begin< 0 >(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;

      const auto local_difference = difference.local_function(entity);
      local_energy_product.apply(*local_difference,
                                 *local_difference,
                                 local_result_matrix,
                                 tmp_matrices);
      assert(local_result_matrix.rows() >= 1);
      assert(local_result_matrix.cols() >= 1);
      const size_t index = grid_view->indexSet().index(entity);
      energy_error_visualization[index] = local_result_matrix[0][0];
      estimator_visualization[index] = local_estimates[ms_grid_->subdomainOf(entity)];
    }
    VTKWriter< typename MsGridType::GlobalGridViewType > vtk_writer(*grid_view);
    vtk_writer.addCellData(energy_error_visualization, "energy error");
    vtk_writer.addCellData(estimator_visualization, "estimator");
    vtk_writer.write(filename + ".energy_error", VTK::appendedraw);

    // visualize difference
    for (auto& element : difference.vector())
      element = std::abs(element);
    difference.visualize(filename + ".difference");
  } // ... visualize_information(...)

  std::vector< double > estimate_local(const std::vector< VectorType >& local_vectors,
                                       const Pymor::Parameter mu = Pymor::Parameter(),
                                       const Pymor::Parameter mu_fixed = Pymor::Parameter()) const
  {
    using namespace GDT;
    std::vector< double > ret(ms_grid_->size());

    if (mu.type() != this->parameter_type())
      DUNE_THROW(Pymor::Exceptions::wrong_parameter_type,
                 "mu is " << mu.type() << ", should be " << this->parameter_type() << "!");
    if (mu_fixed.type() != this->parameter_type())
      DUNE_THROW(Pymor::Exceptions::wrong_parameter_type,
                 "mu_fixed is " << mu_fixed.type() << ", should be " << this->parameter_type() << "!");

    if (local_vectors.size() != ms_grid_->size())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given,
                 "local_vectors is of size " << local_vectors.size() << " and should be of size " << ms_grid_->size());
    VectorType vector(this->ansatz_space()->mapper().size());
    for (size_t ss = 0; ss < ms_grid_->size(); ++ss) {
      if (local_vectors[ss].size() != this->local_discretizations_[ss]->ansatz_space()->mapper().size())
        DUNE_THROW(Stuff::Exceptions::wrong_input_given,
                   "local_vectors[" << ss << "] is of size " << local_vectors[ss].size() << " and should be of size "
                   << this->local_discretizations_[ss]->ansatz_space()->mapper().size());
      if (!local_vectors[ss].valid())
        DUNE_THROW(Stuff::Exceptions::wrong_input_given, "local_vectors[" << ss << "] contains NaN or INF!");
      copy_local_to_global_vector(local_vectors[ss], ss, vector);
    }

    typedef typename MsGridType::GlobalGridViewType GlobalGridViewType;
    const auto grid_view = ms_grid_->globalGridView();
    typedef ConstDiscreteFunction< AnsatzSpaceType, VectorType > ConstDiscreteFunctionType;
    const ConstDiscreteFunctionType discrete_solution(*this->ansatz_space(), vector);

    VectorType oswald_interpolation_vector(this->ansatz_space()->mapper().size());
    typedef DiscreteFunction< AnsatzSpaceType, VectorType > DiscreteFunctionType;
    DiscreteFunctionType oswald_interpolation(*this->ansatz_space(), oswald_interpolation_vector);
    const Operators::OswaldInterpolation< GlobalGridViewType > oswald_interpolation_operator(*grid_view);
    oswald_interpolation_operator.apply(discrete_solution, oswald_interpolation);

    typedef Spaces::FiniteVolume::Default< GlobalGridViewType, RangeFieldType, 1, 1 > P0SpaceType;
    const P0SpaceType p0_space(grid_view);
    VectorType p0_force_vector(p0_space.mapper().size());
    DiscreteFunction< P0SpaceType, VectorType > p0_force(p0_space, p0_force_vector);
    Operators::Projection< GlobalGridViewType > projection_operator(*grid_view);
    typedef typename ProblemType::FunctionType::NonparametricType FunctionType;
    const auto& force = *(this->problem().force());
    std::shared_ptr< const FunctionType > nonparametric_force(nullptr);
    if (force.parametric())
      nonparametric_force = force.with_mu(mu);
    else
      nonparametric_force = force.affine_part();
    projection_operator.apply(*nonparametric_force, p0_force);

    typedef Spaces::RaviartThomas::PdelabBased< GlobalGridViewType, 0, RangeFieldType, dimDomain > RTN0SpaceType;
    const RTN0SpaceType rtn0_space(grid_view);
    VectorType diffusive_flux_vector(rtn0_space.mapper().size());
    typedef DiscreteFunction< RTN0SpaceType, VectorType > RTN0DiscreteFunctionType;
    RTN0DiscreteFunctionType diffusive_flux(rtn0_space, diffusive_flux_vector);

    typedef typename ProblemType::DiffusionTensorType::NonparametricType DiffusionTensorType;
    const auto& diffusion_tensor = *(this->problem().diffusion_tensor());
    assert(!diffusion_tensor.parametric());
    typedef typename ProblemType::DiffusionFactorType::NonparametricType DiffusionFactorType;
    const auto& diffusion_factor = *(this->problem().diffusion_factor());
    std::shared_ptr< const DiffusionFactorType > nonparametric_diffusion_factor(nullptr);
    std::shared_ptr< const DiffusionFactorType > fixed_nonparametric_diffusion_factor(nullptr);
    if (diffusion_factor.parametric()) {
      nonparametric_diffusion_factor = diffusion_factor.with_mu(mu);
      fixed_nonparametric_diffusion_factor = diffusion_factor.with_mu(mu_fixed);
    } else {
      nonparametric_diffusion_factor = diffusion_factor.affine_part();
      fixed_nonparametric_diffusion_factor = diffusion_factor.affine_part();
    }

    const double gamma = diffusion_factor.gamma(mu, mu_fixed);
    const double alpha = diffusion_factor.alpha(mu, mu_fixed);
    assert(gamma > 0);
    assert(alpha > 0);
    const double parameter_factor = std::pow(std::max(std::sqrt(gamma), 1.0/std::sqrt(alpha)), 2);

    const Operators::DiffusiveFluxReconstruction< GlobalGridViewType, DiffusionFactorType, DiffusionTensorType >
      diffusive_flux_reconstruction(*grid_view, *nonparametric_diffusion_factor, *diffusion_tensor.affine_part());
    diffusive_flux_reconstruction.apply(discrete_solution, diffusive_flux);

    const LocalOperator::Codim0Integral< LocalEvaluation::Elliptic< DiffusionFactorType, DiffusionTensorType > >
        local_eta_nc_product(1, *nonparametric_diffusion_factor, *diffusion_tensor.affine_part());
    const auto eta_nc_difference = discrete_solution - oswald_interpolation;

    typedef Stuff::Functions::Constant< EntityType, DomainFieldType, dimDomain, RangeFieldType, 1 >
        ConstantFunctionType;
    typedef typename ConstantFunctionType::DomainType DomainType;
    const ConstantFunctionType constant_one(1);
    const  LocalOperator::Codim0Integral< LocalEvaluation::Product< ConstantFunctionType > >
      local_eta_r_product(1, constant_one);
    const auto eta_r_difference = *nonparametric_force - p0_force;

    const LocalOperator::Codim0Integral< LocalEvaluation::ESV2007::DiffusiveFluxEstimate< DiffusionFactorType
                                                                                        , RTN0DiscreteFunctionType
                                                                                        , DiffusionTensorType > >
        local_eta_df_product(1,
                             *fixed_nonparametric_diffusion_factor,
                             *diffusion_tensor.affine_part(),
                             diffusive_flux);

    std::vector< DynamicMatrix< RangeFieldType > >
        tmp_matrices(std::max(std::max(local_eta_nc_product.numTmpObjectsRequired(),
                                       local_eta_r_product.numTmpObjectsRequired()),
                              local_eta_df_product.numTmpObjectsRequired()),
                     DynamicMatrix< RangeFieldType >(1, 1, 0.0));
    DynamicMatrix< RangeFieldType > local_result_matrix(1, 1, 0.0);

    // estimate
    double eta_nc_squared = 0.0;
    double eta_r_squared = 0.0;
    double eta_df_squared = 0.0;
    for (size_t subdomain = 0; subdomain < ms_grid_->size(); ++subdomain) {
      const auto local_grid_part = ms_grid_->localGridPart(subdomain);

      // some storage
      std::vector< DomainType > vertices;
      RangeFieldType min_diffusion_value = std::numeric_limits< RangeFieldType >::max();

      // walk the local grid
      double eta_nc_T_squared = 0.0;
      double eta_r_T_squared = 0.0;
      double eta_df_T_squared = 0.0;
      const auto entity_it_end = local_grid_part->template end< 0 >();
      for (auto entity_it = local_grid_part->template begin< 0 >(); entity_it != entity_it_end; ++entity_it) {
        const auto& entity = *entity_it;

        // collect all vertices
        for (int cc = 0; cc < entity.template count< dimDomain >(); ++cc)
          vertices.push_back(entity.template subEntity< dimDomain >(cc)->geometry().center());

        // compute minimum diffusion
        const RangeFieldType min_diffusion_value_entity
            = internal::compute_minimum(*nonparametric_diffusion_factor, entity)
              * internal::compute_minimum(*diffusion_tensor.affine_part(), entity);
        min_diffusion_value = std::min(min_diffusion_value, min_diffusion_value_entity);

        local_result_matrix *= 0.0;
        const auto local_eta_nc_difference = eta_nc_difference.local_function(entity);
        local_eta_nc_product.apply(*local_eta_nc_difference,
                                   *local_eta_nc_difference,
                                   local_result_matrix,
                                   tmp_matrices);
        assert(local_result_matrix.rows() >= 1);
        assert(local_result_matrix.cols() >= 1);
        eta_nc_T_squared += local_result_matrix[0][0];

        local_result_matrix *= 0.0;
        const auto local_eta_r_difference = eta_r_difference.local_function(entity);
        local_eta_r_product.apply(*local_eta_r_difference,
                                  *local_eta_r_difference,
                                  local_result_matrix,
                                  tmp_matrices);
        assert(local_result_matrix.rows() >= 1);
        assert(local_result_matrix.cols() >= 1);
        eta_r_T_squared += local_result_matrix[0][0];

        local_result_matrix *= 0.0;
        const auto local_discrete_solution = discrete_solution.local_function(entity);
        local_eta_df_product.apply(*local_discrete_solution,
                                   *local_discrete_solution,
                                   local_result_matrix,
                                   tmp_matrices);
        assert(local_result_matrix.rows() >= 1);
        assert(local_result_matrix.cols() >= 1);
        eta_df_T_squared += local_result_matrix[0][0];
      } // walk the local grid

      // compute diameter
      DomainFieldType diameter(0);
      DomainType tmp_diff;
      for (size_t cc = 0; cc < vertices.size(); ++cc) {
        const auto& vertex = vertices[cc];
        for (size_t dd = cc + 1; dd < vertices.size(); ++dd) {
          const auto& other_vertex = vertices[dd];
          tmp_diff = vertex - other_vertex;
          diameter = std::max(diameter, tmp_diff.two_norm());
        }
      }

      const RangeFieldType poincare_constant = 1.0 / (M_PIl * M_PIl);
      assert(min_diffusion_value > 0.0);

      eta_nc_squared += eta_nc_T_squared;
      eta_r_squared += ((poincare_constant * diameter * diameter) / min_diffusion_value) * eta_r_T_squared;
      eta_df_squared += eta_df_T_squared;

      ret[subdomain] = eta_nc_squared + eta_r_squared + parameter_factor * eta_df_squared;
    } // walk the subdomains

    return ret;
  } // ... estimate_local(...)

  VectorType solve_for_local_correction(const std::vector< VectorType >& local_vectors,
                                        const size_t subdomain,
                                        const Pymor::Parameter mu = Pymor::Parameter()) const
  {
    using namespace GDT;

    if (mu.type() != this->parameter_type())
      DUNE_THROW(Pymor::Exceptions::wrong_parameter_type,
                 "mu is " << mu.type() << ", should be " << this->parameter_type() << "!");
    if (local_vectors.size() != ms_grid_->size())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given,
                 "local_vectors is of size " << local_vectors.size() << " and should be of size " << ms_grid_->size());
    VectorType vector(this->ansatz_space()->mapper().size());
    for (size_t ss = 0; ss < ms_grid_->size(); ++ss) {
      if (local_vectors[ss].size() != this->local_discretizations_[ss]->ansatz_space()->mapper().size())
        DUNE_THROW(Stuff::Exceptions::wrong_input_given,
                   "local_vectors[" << ss << "] is of size " << local_vectors[ss].size() << " and should be of size "
                   << this->local_discretizations_[ss]->ansatz_space()->mapper().size());
      if (!local_vectors[ss].valid())
        DUNE_THROW(Stuff::Exceptions::wrong_input_given, "local_vectors[" << ss << "] contains NaN or INF!");
      copy_local_to_global_vector(local_vectors[ss], ss, vector);
    }

    const std::string prefix = "subdomain_" + DSC::toString(subdomain) + "_";

    const ConstDiscreteFunction< AnsatzSpaceType, VectorType > current_global_solution(*this->ansatz_space(), vector);
    current_global_solution.visualize(prefix + "current_solution_global");

    typedef SWIPDG< typename MsGridType::GridType, Stuff::Grid::ChooseLayer::local_oversampled, RangeFieldType, dimRange
                  , 1, GDT::ChooseSpaceBackend::fem, la_backend > OversampledDiscretizationType;
    OversampledDiscretizationType oversampled_discretization(grid_provider_,
                                                             this->multiscale_boundary_config_,
                                                             this->problem(),
                                                             subdomain);
    oversampled_discretization.init();

    DiscreteFunction< typename OversampledDiscretizationType::AnsatzSpaceType, VectorType >
        current_oversampled_solution(*oversampled_discretization.ansatz_space());
    const Operators::Projection< typename OversampledDiscretizationType::GridViewType >
        oversampled_projection_operator(*oversampled_discretization.grid_view());
    oversampled_projection_operator.apply(current_global_solution, current_oversampled_solution);
    current_oversampled_solution.visualize(prefix + "current_solution_oversampled");

    if (!oversampled_discretization.rhs()->has_affine_part())
      oversampled_discretization.rhs()->register_affine_part(oversampled_discretization.test_space()->mapper().size());
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
    current_oversampled_solution.visualize(prefix + "correction_oversampled");

    DiscreteFunction< typename LocalDiscretizationType::AnsatzSpaceType, VectorType >
        local_solution(*this->local_discretizations_[subdomain]->ansatz_space());
    const Operators::Projection< typename LocalDiscretizationType::GridViewType >
        local_projection_operator(*this->local_discretizations_[subdomain]->grid_view());
    local_projection_operator.apply(current_oversampled_solution, local_solution);
    local_solution.visualize(prefix + "correction_local");

    return local_solution.vector();
  } // ... solve_for_local_correction(...)

private:
  class CouplingAssembler
  {
    typedef Dune::DynamicMatrix< RangeFieldType > LocalMatrixType;
    typedef Dune::DynamicVector< RangeFieldType > LocalVectorType;
    typedef std::vector< std::vector< LocalMatrixType > > LocalMatricesContainerType;
    typedef std::vector< std::vector< LocalVectorType > > LocalVectorsContainerType;
    typedef std::vector< Dune::DynamicVector< size_t > > IndicesContainer;

    typedef typename LocalDiscretizationType::TestSpaceType LocalTestSpaceType;
    typedef typename LocalDiscretizationType::AnsatzSpaceType LocalAnsatzSpaceType;
    typedef typename MsGridType::CouplingGridPartType CouplingGridPartType;

    class LocalCodim1MatrixAssemblerApplication
    {
    public:
      virtual ~LocalCodim1MatrixAssemblerApplication(){}

      virtual void apply(const LocalTestSpaceType& /*inner_test_space*/,
                         const LocalAnsatzSpaceType& /*inner_ansatz_space*/,
                         const LocalTestSpaceType& /*outer_test_space*/,
                         const LocalAnsatzSpaceType& /*outer_ansatz_space*/,
                         const typename CouplingGridPartType::IntersectionType& /*_intersection*/,
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
                         const typename CouplingGridPartType::IntersectionType& intersection,
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
    CouplingAssembler(const LocalTestSpaceType& inner_test_space,
                      const LocalAnsatzSpaceType& inner_ansatz_space,
                      const LocalTestSpaceType& outer_test_space,
                      const LocalAnsatzSpaceType& outer_ansatz_space,
                      const CouplingGridPartType& grid_part)
      : innerTestSpace_(inner_test_space)
      , innerAnsatzSpace_(inner_ansatz_space)
      , outerTestSpace_(outer_test_space)
      , outerAnsatzSpace_(outer_ansatz_space)
      , grid_part_(grid_part)
    {}

    ~CouplingAssembler()
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

    void assemble() const
    {
      // only do something, if there are local assemblers
      if (localCodim1MatrixAssemblers_.size() > 0) {
        // common tmp storage for all entities
        // * for the matrix assemblers
        std::vector< size_t > numberOfTmpMatricesNeeded(2, 0);
        for (auto& localCodim1MatrixAssembler : localCodim1MatrixAssemblers_) {
          const auto tmp = localCodim1MatrixAssembler->numTmpObjectsRequired();
          assert(tmp.size() == 2);
          numberOfTmpMatricesNeeded[0] = std::max(numberOfTmpMatricesNeeded[0], tmp[0]);
          numberOfTmpMatricesNeeded[1] = std::max(numberOfTmpMatricesNeeded[1], tmp[1]);
        }
        const size_t maxLocalSize = std::max(innerTestSpace_.mapper().maxNumDofs(),
                                             std::max(innerAnsatzSpace_.mapper().maxNumDofs(),
                                                      std::max(outerTestSpace_.mapper().maxNumDofs(),
                                                               outerAnsatzSpace_.mapper().maxNumDofs())));
        std::vector< LocalMatrixType > tmpLocalAssemblerMatrices( numberOfTmpMatricesNeeded[0],
                                                                  LocalMatrixType(maxLocalSize,
                                                                                  maxLocalSize,
                                                                                  RangeFieldType(0)));
        std::vector< LocalMatrixType > tmpLocalOperatorMatrices(numberOfTmpMatricesNeeded[1],
                                                                LocalMatrixType(maxLocalSize,
                                                                                maxLocalSize,
                                                                                RangeFieldType(0)));
        std::vector< std::vector< LocalMatrixType > > tmpLocalMatricesContainer;
        tmpLocalMatricesContainer.push_back(tmpLocalAssemblerMatrices);
        tmpLocalMatricesContainer.push_back(tmpLocalOperatorMatrices);
        // * for the global indices
        std::vector< Dune::DynamicVector< size_t > > tmpIndices = {
            Dune::DynamicVector< size_t >(maxLocalSize)
          , Dune::DynamicVector< size_t >(maxLocalSize)
          , Dune::DynamicVector< size_t >(maxLocalSize)
          , Dune::DynamicVector< size_t >(maxLocalSize)
        };

        // walk the grid
        const auto entityEndIt = grid_part_.template end< 0 >();
        for(auto entityIt = grid_part_.template begin< 0 >(); entityIt != entityEndIt; ++entityIt ) {
          const auto& entity = *entityIt;
          // walk the intersections
          const auto intersectionEndIt = grid_part_.iend(entity);
          for (auto intersectionIt = grid_part_.ibegin(entity);
               intersectionIt != intersectionEndIt;
               ++intersectionIt) {
            const auto& intersection = *intersectionIt;
            // for a coupling grid part, we can be sure to only get the inner coupling intersetcions
            // so no further check neccesarry than
            assert(intersection.neighbor() && !intersection.boundary());
            // call local matrix assemblers
            for (auto& localCodim1MatrixAssembler : localCodim1MatrixAssemblers_) {
              localCodim1MatrixAssembler->apply(innerTestSpace_, innerAnsatzSpace_,
                                                outerTestSpace_, outerAnsatzSpace_,
                                                intersection,
                                                tmpLocalMatricesContainer, tmpIndices);
            }
          } // walk the intersections
        } // walk the grid
      } // only do something, if there are local assemblers
    } // void assemble() const

  private:
    const LocalTestSpaceType& innerTestSpace_;
    const LocalAnsatzSpaceType& innerAnsatzSpace_;
    const LocalTestSpaceType& outerTestSpace_;
    const LocalAnsatzSpaceType& outerAnsatzSpace_;
    const CouplingGridPartType& grid_part_;
    std::vector< LocalCodim1MatrixAssemblerApplication* > localCodim1MatrixAssemblers_;
  }; // class CouplingAssembler

  void add_local_to_global_pattern(const PatternType& local,
                                   const size_t test_subdomain,
                                   const size_t ansatz_subdomain,
                                   PatternType& global) const
  {
    for (size_t local_ii = 0; local_ii < local.size(); ++local_ii) {
      const size_t global_ii = this->test_space()->mapper().mapToGlobal(test_subdomain, local_ii);
      const auto& local_rows = local.inner(local_ii);
      auto& global_rows = global.inner(global_ii);
      for (const auto& local_jj : local_rows) {
        const size_t global_jj = this->ansatz_space()->mapper().mapToGlobal(ansatz_subdomain, local_jj);
        global_rows.insert(global_jj);
      }
    }
  } // ... add_local_to_global_pattern(...)

  void copy_local_to_global_matrix(const AffinelyDecomposedConstMatrixType& local_matrix,
                                   const PatternType& local_pattern,
                                   const size_t subdomain,
                                   const size_t neighbor,
                                   AffinelyDecomposedMatrixType& global_matrix) const
  {
    for (size_t qq = 0; qq < local_matrix.num_components(); ++qq) {
      const auto coefficient = local_matrix.coefficient(qq);
      ssize_t comp = find_component(global_matrix, *coefficient);
      if (comp < 0)
        comp = global_matrix.register_component(coefficient,
                                                this->test_space()->mapper().size(),
                                                this->ansatz_space()->mapper().size(),
                                                pattern_);
      assert(comp >= 0);
      copy_local_to_global_matrix(*(local_matrix.component(qq)),
                                  local_pattern,
                                  subdomain,
                                  neighbor,
                                  *(global_matrix.component(comp)));
    }
    if (local_matrix.has_affine_part()) {
      if (!global_matrix.has_affine_part())
        global_matrix.register_affine_part(this->test_space_->mapper().size(),
                                    this->ansatz_space_->mapper().size(),
                                    pattern_);
      copy_local_to_global_matrix(*(local_matrix.affine_part()),
                                  local_pattern,
                                  subdomain,
                                  neighbor,
                                  *(global_matrix.affine_part()));
    }
  } // copy_local_to_global_matrix(...)

  template< class ML, class MG >
  void copy_local_to_global_matrix(const Stuff::LA::MatrixInterface< ML >& local_matrix,
                                   const PatternType& local_pattern,
                                   const size_t test_subdomain,
                                   const size_t ansatz_subdomain,
                                   Stuff::LA::MatrixInterface< MG >& global_matrix) const
  {
    for (size_t local_ii = 0; local_ii < local_pattern.size(); ++local_ii) {
      const size_t global_ii = this->test_space()->mapper().mapToGlobal(test_subdomain, local_ii);
      for (const size_t& local_jj : local_pattern.inner(local_ii)) {
        const size_t global_jj = this->ansatz_space()->mapper().mapToGlobal(ansatz_subdomain, local_jj);
        global_matrix.add_to_entry(global_ii, global_jj, local_matrix.get_entry(local_ii, local_jj));
      }
    }
  } // ... copy_local_to_global_matrix(...)

  void copy_local_to_global_vector(const AffinelyDecomposedConstVectorType& local_vector,
                                   const size_t subdomain,
                                   AffinelyDecomposedVectorType& global_vector) const
  {
    for (size_t qq = 0; qq < local_vector.num_components(); ++qq) {
      const auto coefficient = local_vector.coefficient(qq);
      ssize_t comp = find_component(global_vector, *coefficient);
      if (comp < 0)
        comp = global_vector.register_component(coefficient,
                                                this->test_space()->mapper().size());
      assert(comp >= 0);
      copy_local_to_global_vector(*(local_vector.component(qq)),
                                  subdomain,
                                  *(global_vector.component(comp)));
    }
    if (local_vector.has_affine_part()) {
      if (!global_vector.has_affine_part())
        global_vector.register_affine_part(this->test_space()->mapper().size());
      copy_local_to_global_vector(*(local_vector.affine_part()),
                                  subdomain,
                                  *(global_vector.affine_part()));
    }
  } // copy_local_to_global_vector(...)

  template< class VL, class VG >
  void copy_local_to_global_vector(const Stuff::LA::VectorInterface< VL >& local_vector,
                                   const size_t subdomain,
                                   Stuff::LA::VectorInterface< VG >& global_vector) const
  {
    for (size_t local_ii = 0; local_ii < local_vector.size(); ++local_ii) {
      const size_t global_ii = this->test_space()->mapper().mapToGlobal(subdomain, local_ii);
      global_vector.add_to_entry(global_ii, local_vector.get_entry(local_ii));
    }
  } // ... copy_local_to_global_vector(...)

  void assemble_boundary_contributions(const size_t subdomain) const
  {
    typedef typename MsGridType::BoundaryGridPartType BoundaryGridPartType;
    typedef typename LocalDiscretizationType::TestSpaceType   LocalTestSpaceType;
    typedef typename LocalDiscretizationType::AnsatzSpaceType LocalAnsatzSpaceType;
    const LocalTestSpaceType&   local_test_space   = *(this->local_discretizations_[subdomain]->test_space());
    const LocalAnsatzSpaceType& local_ansatz_space = *(this->local_discretizations_[subdomain]->ansatz_space());
    typedef GDT::SystemAssembler< LocalTestSpaceType, BoundaryGridPartType, LocalAnsatzSpaceType > BoundaryAssemblerType;
    BoundaryAssemblerType boundary_assembler(local_test_space,
                                             local_ansatz_space,
                                             *(ms_grid_->boundaryGridPart(subdomain)));

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
    std::vector< DirichletOperatorType* > dirichlet_operators;
    std::vector< DirichletMatrixAssemblerType* > dirichlet_matrix_assemblers;
    for (DUNE_STUFF_SSIZE_T qq = 0; qq < this->problem().diffusion_factor()->num_components(); ++qq) {
      dirichlet_operators.push_back(new DirichletOperatorType(*(this->problem().diffusion_factor()->component(qq)),
                                                              *(diffusion_tensor.affine_part())));
      dirichlet_matrix_assemblers.push_back(new DirichletMatrixAssemblerType(*(dirichlet_operators[qq])));
      boundary_assembler.add(*(dirichlet_matrix_assemblers[qq]),
                             *(local_matrix.component(qq)),
                             new GDT::ApplyOn::DirichletIntersections< BoundaryGridPartType >(this->boundary_info()));
    }
    if (this->problem().diffusion_factor()->has_affine_part()) {
      dirichlet_operators.push_back(new DirichletOperatorType(*(this->problem().diffusion_factor()->affine_part()),
                                                              *(diffusion_tensor.affine_part())));
      dirichlet_matrix_assemblers.push_back(new DirichletMatrixAssemblerType(*(
          dirichlet_operators[dirichlet_operators.size() - 1])));
      boundary_assembler.add(*(dirichlet_matrix_assemblers[dirichlet_matrix_assemblers.size() - 1]),
                             *(local_matrix.affine_part()),
                             new GDT::ApplyOn::DirichletIntersections< BoundaryGridPartType >(this->boundary_info()));
    }

    // rhs
    // * neumann boundary terms
    typedef typename ProblemType::FunctionType::NonparametricType NeumannType;
    typedef GDT::LocalFunctional::Codim1Integral< GDT::LocalEvaluation::Product< NeumannType > > NeumannFunctionalType;
    typedef GDT::LocalAssembler::Codim1Vector< NeumannFunctionalType > NeumannVectorAssemblerType;
    std::vector< NeumannFunctionalType* > neumann_functionals;
    std::vector< NeumannVectorAssemblerType* > neumann_vector_assemblers;
    for (DUNE_STUFF_SSIZE_T qq = 0; qq < this->problem().neumann()->num_components(); ++qq) {
      neumann_functionals.push_back(new NeumannFunctionalType(*(this->problem().neumann()->component(qq))));
      neumann_vector_assemblers.push_back(new NeumannVectorAssemblerType(*(neumann_functionals[qq])));
      boundary_assembler.add(*(neumann_vector_assemblers[qq]),
                             *(local_vector.component(this->problem().force()->num_components() + qq)),
                             new GDT::ApplyOn::NeumannIntersections< BoundaryGridPartType >(this->boundary_info()));
    }
    if (this->problem().neumann()->has_affine_part()) {
      neumann_functionals.push_back(new NeumannFunctionalType(*(this->problem().neumann()->affine_part())));
      neumann_vector_assemblers.push_back(new NeumannVectorAssemblerType(*(
          neumann_functionals[neumann_functionals.size() - 1])));
      boundary_assembler.add(*(neumann_vector_assemblers[neumann_vector_assemblers.size() - 1]),
                             *(local_vector.affine_part()),
                             new GDT::ApplyOn::NeumannIntersections< BoundaryGridPartType >(this->boundary_info()));
    }

    // * dirichlet boundary terms
    typedef typename ProblemType::FunctionType::NonparametricType  DirichletType;
    typedef GDT::LocalFunctional::Codim1Integral< GDT::LocalEvaluation::SWIPDG::BoundaryRHS< DiffusionFactorType,
                                                                                             DirichletType,
                                                                                             DiffusionTensorType > >
        DirichletFunctionalType;
    typedef GDT::LocalAssembler::Codim1Vector< DirichletFunctionalType > DirichletVectorAssemblerType;
    std::vector< DirichletFunctionalType* > dirichlet_functionals;
    std::vector< DirichletVectorAssemblerType* > dirichlet_vector_assemblers;
    size_t component_index = this->problem().force()->num_components() + this->problem().neumann()->num_components();
    if (this->problem().diffusion_factor()->has_affine_part()) {
      for (DUNE_STUFF_SSIZE_T qq = 0; qq < this->problem().dirichlet()->num_components(); ++qq) {
        dirichlet_functionals.push_back(new DirichletFunctionalType(*(this->problem().diffusion_factor()->affine_part()),
                                                                    *(diffusion_tensor.affine_part()),
                                                                    *(this->problem().dirichlet()->component(qq))));
        dirichlet_vector_assemblers.push_back(new DirichletVectorAssemblerType(*(
            dirichlet_functionals[dirichlet_functionals.size() - 1])));
        boundary_assembler.add(*(dirichlet_vector_assemblers[dirichlet_vector_assemblers.size() - 1]),
                               *(local_vector.component(component_index)),
                               new GDT::ApplyOn::DirichletIntersections< BoundaryGridPartType >(this->boundary_info()));
        ++component_index;
      }
    }
    if (this->problem().dirichlet()->has_affine_part()) {
      for (DUNE_STUFF_SSIZE_T qq = 0; qq < this->problem().diffusion_factor()->num_components(); ++qq) {
        dirichlet_functionals.push_back(new DirichletFunctionalType(*(this->problem().diffusion_factor()->component(qq)),
                                                                    *(diffusion_tensor.affine_part()),
                                                                    *(this->problem().dirichlet()->affine_part())));
        dirichlet_vector_assemblers.push_back(new DirichletVectorAssemblerType(*(
            dirichlet_functionals[dirichlet_functionals.size() - 1])));
        boundary_assembler.add(*(dirichlet_vector_assemblers[dirichlet_vector_assemblers.size() - 1]),
                               *(local_vector.component(component_index)),
                               new GDT::ApplyOn::DirichletIntersections< BoundaryGridPartType >(this->boundary_info()));
        ++component_index;
      }
    }
    for (DUNE_STUFF_SSIZE_T pp = 0; pp < this->problem().diffusion_factor()->num_components(); ++ pp) {
      for (DUNE_STUFF_SSIZE_T qq = 0; qq < this->problem().dirichlet()->num_components(); ++qq) {
        dirichlet_functionals.push_back(new DirichletFunctionalType(*(this->problem().diffusion_factor()->component(pp)),
                                                                    *(diffusion_tensor.affine_part()),
                                                                    *(this->problem().dirichlet()->component(qq))));
        dirichlet_vector_assemblers.push_back(new DirichletVectorAssemblerType(*(
            dirichlet_functionals[dirichlet_functionals.size() - 1])));
        boundary_assembler.add(*(dirichlet_vector_assemblers[dirichlet_vector_assemblers.size() - 1]),
                               *(local_vector.component(component_index)),
                               new GDT::ApplyOn::DirichletIntersections< BoundaryGridPartType >(this->boundary_info()));
        ++component_index;
      }
    } // dirichlet boundary terms

    // do the actual work
    boundary_assembler.assemble();

    // clean up
    for (auto& element : dirichlet_vector_assemblers) delete element;
    for (auto& element : dirichlet_functionals)       delete element;
    for (auto& element : neumann_vector_assemblers)   delete element;
    for (auto& element : neumann_functionals)         delete element;
    for (auto& element : dirichlet_matrix_assemblers) delete element;
    for (auto& element : dirichlet_operators)         delete element;
  } // ... assemble_boundary_contributions(...)

  /**
   * \note  We take the matrices as input here becaus we would have to look them up in the maps otherwise. Since that
   *        has already been done above we save a little.
   */
  void assemble_coupling_contributions(const size_t subdomain,
                                       const size_t neighbour,
                                       AffinelyDecomposedMatrixType& inside_inside_matrix,
                                       AffinelyDecomposedMatrixType& inside_outside_matrix,
                                       AffinelyDecomposedMatrixType& outside_inside_matrix,
                                       AffinelyDecomposedMatrixType& outside_outside_matrix) const
  {
    typedef typename LocalDiscretizationType::TestSpaceType   LocalTestSpaceType;
    typedef typename LocalDiscretizationType::AnsatzSpaceType LocalAnsatzSpaceType;
    const LocalTestSpaceType&   inner_test_space   = *(this->local_discretizations_[subdomain]->test_space());
    const LocalAnsatzSpaceType& inner_ansatz_space = *(this->local_discretizations_[subdomain]->ansatz_space());
    const LocalTestSpaceType&   outer_test_space   = *(this->local_discretizations_[neighbour]->test_space());
    const LocalAnsatzSpaceType& outer_ansatz_space = *(this->local_discretizations_[neighbour]->ansatz_space());
    CouplingAssembler coupling_assembler(inner_test_space, inner_ansatz_space,
                                         outer_test_space, outer_ansatz_space,
                                         *(ms_grid_->couplingGridPart(subdomain, neighbour)));

    typedef typename ProblemType::DiffusionFactorType::NonparametricType DiffusionFactorType;
    typedef typename ProblemType::DiffusionTensorType::NonparametricType DiffusionTensorType;
    const auto& diffusion_tensor = *(this->problem().diffusion_tensor());
    assert(!diffusion_tensor.parametric());
    assert(diffusion_tensor.has_affine_part());
    typedef GDT::LocalOperator::Codim1CouplingIntegral< GDT::LocalEvaluation::SWIPDG::Inner< DiffusionFactorType,
                                                                                             DiffusionTensorType > >
        CouplingOperatorType;
    typedef GDT::LocalAssembler::Codim1CouplingMatrix< CouplingOperatorType > CouplingMatrixAssemblerType;
    std::vector< CouplingOperatorType* > coupling_operators;
    std::vector< CouplingMatrixAssemblerType* > coupling_matrix_assemblers;
    for (DUNE_STUFF_SSIZE_T qq = 0; qq < this->problem().diffusion_factor()->num_components(); ++qq) {
      coupling_operators.push_back(new CouplingOperatorType(*(this->problem().diffusion_factor()->component(qq)),
                                                            *(diffusion_tensor.affine_part())));
      coupling_matrix_assemblers.push_back(new CouplingMatrixAssemblerType(*(coupling_operators[qq])));
      coupling_assembler.addLocalAssembler(*(coupling_matrix_assemblers[qq]),
                                           *(inside_inside_matrix.component(qq)),
                                           *(inside_outside_matrix.component(qq)),
                                           *(outside_inside_matrix.component(qq)),
                                           *(outside_outside_matrix.component(qq)));
    }
    if (this->problem().diffusion_factor()->has_affine_part()) {
      coupling_operators.push_back(new CouplingOperatorType(*(this->problem().diffusion_factor()->affine_part()),
                                                            *(diffusion_tensor.affine_part())));
      coupling_matrix_assemblers.push_back(new CouplingMatrixAssemblerType(*(
          coupling_operators[coupling_operators.size() - 1])));
      coupling_assembler.addLocalAssembler(*(coupling_matrix_assemblers[coupling_matrix_assemblers.size() - 1]),
                                           *(inside_inside_matrix.affine_part()),
                                           *(inside_outside_matrix.affine_part()),
                                           *(outside_inside_matrix.affine_part()),
                                           *(outside_outside_matrix.affine_part()));
    }

    // do the actual work
    coupling_assembler.assemble();

    // clean up
    for (auto& element : coupling_matrix_assemblers)  delete element;
    for (auto& element : coupling_operators)          delete element;
  } // ... assemble_coupling_contributions(...)

  void build_global_containers()
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

  template< class AffinelyDecomposedContainerType >
  ssize_t find_component(const AffinelyDecomposedContainerType& container,
                         const Pymor::ParameterFunctional& coefficient) const
  {
    for (size_t qq = 0; qq < container.num_components(); ++qq)
      if (*(container.coefficient(qq)) == coefficient)
        return qq;
    return -1;
  } // ... find_component(...)

  const GridProviderType& grid_provider_;
  std::shared_ptr< const MsGridType > ms_grid_;
  PatternType pattern_;
  std::vector< std::shared_ptr< AffinelyDecomposedMatrixType > > local_matrices_;
  std::vector< std::shared_ptr< AffinelyDecomposedVectorType > > local_vectors_;
  std::vector< std::map< size_t, std::shared_ptr< PatternType > > > inside_outside_patterns_;
  std::vector< std::map< size_t, std::shared_ptr< PatternType > > > outside_inside_patterns_;
  std::vector< std::map< size_t, std::shared_ptr< AffinelyDecomposedMatrixType > > > inside_outside_matrices_;
  std::vector< std::map< size_t, std::shared_ptr< AffinelyDecomposedMatrixType > > > outside_inside_matrices_;
//  std::vector< PatternType > local_boundary_patterns_;
//  std::vector< PatternType > local_coupling_patterns_;
//  std::vector< std::map< size_t, PatternType > > inside_outside_coupling_patterns_;
//  std::vector< std::map< size_t, PatternType > > outside_inside_coupling_patterns_;
//  std::vector< std::shared_ptr< AffinelyDecomposedMatrixType > > boundary_matrices_;
//  std::vector< std::shared_ptr< AffinelyDecomposedVectorType > > boundary_vectors_;
//  std::vector< std::shared_ptr< AffinelyDecomposedMatrixType > > local_coupling_matrices_;
//  std::vector< std::map< size_t, std::shared_ptr< AffinelyDecomposedMatrixType > > > inside_outside_coupling_matrices_;
//  std::vector< std::map< size_t, std::shared_ptr< AffinelyDecomposedMatrixType > > > outside_inside_coupling_matrices_;
}; // BlockSWIPDG


} // namespace Discretizations
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_BLOCK_SWIPDG_HH
