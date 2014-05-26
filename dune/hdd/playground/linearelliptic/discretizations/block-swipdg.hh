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

#include <dune/common/timer.hh>

#include <dune/grid/multiscale/provider.hh>

#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/configtree.hh>
#include <dune/stuff/grid/layers.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/la/container.hh>
#include <dune/stuff/la/solver.hh>

#include <dune/gdt/spaces/discontinuouslagrange.hh>
#include <dune/gdt/playground/spaces/block.hh>
#include <dune/gdt/playground/localevaluation/swipdg.hh>

#include <dune/hdd/playground/linearelliptic/problems/zero-boundary.hh>

#include "../../../linearelliptic/discretizations/base.hh"
#include "swipdg.hh"

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
                , polOrder, GDT::ChooseSpaceBackend::fem_localfunctions, la_backend > DiscretizationType;
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
    , local_discretizations_(grid_provider.num_subdomains(), nullptr)
    , local_test_spaces_(grid_provider.num_subdomains(), nullptr)
    , local_ansatz_spaces_(grid_provider.num_subdomains(), nullptr)
  {
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


} // namespace internal


/**
 * \attention The given problem is replaced by a Problems::ZeroBoundary.
 * \attention The given boundary info config is replaced by a Stuff::Grid::BoundaryInfos::AllDirichlet.
 * \attention The boundary info for the local oversampled discretizations is hardwired to dirichlet zero atm!
 * \todo  The local products are assembled using a local discretization. This might not be optimal, since we also
 *        assemble a local system matric and right hand side which are never needed. But since this all happens in one
 *        grid walk the overhead seems ok.
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

  static const unsigned int dimDomain = BaseType::dimDomain;
  static const unsigned int dimRange = BaseType::dimRange;

  typedef grid::Multiscale::ProviderInterface< GridImp > GridProviderType;
  typedef typename GridProviderType::MsGridType MsGridType;

private:
  typedef typename Traits::LocalDiscretizationsContainerType::DiscretizationType LocalDiscretizationType;
  typedef typename TestSpaceType::PatternType PatternType;

  using typename BaseType::AffinelyDecomposedMatrixType;
  using typename BaseType::AffinelyDecomposedVectorType;
  typedef Pymor::LA::AffinelyDecomposedConstContainer< MatrixType > AffinelyDecomposedConstMatrixType;
  typedef Pymor::LA::AffinelyDecomposedConstContainer< VectorType > AffinelyDecomposedConstVectorType;

public:
  typedef typename LocalDiscretizationType::ProblemType LocalProblemType;
  typedef typename LocalDiscretizationType::ProductType LocalProductType;

  static std::string static_id()
  {
    return typename DiscretizationInterface< Traits >::static_id() + ".block-swipdg";
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
    if (this->problem_.diffusion_tensor().parametric())
      DUNE_THROW_COLORFULLY(NotImplemented, "The diffusion tensor must not be parametric!");
    if (!this->problem_.diffusion_tensor().has_affine_part())
      DUNE_THROW_COLORFULLY(Stuff::Exceptions::wrong_input_given, "The diffusion tensor must not be empty!");
  } // BlockSWIPDG(...)

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
        assert(!this->local_discretizations_[ss]->parametric());
//        oversampled_discretizations_[ss]->init();
//        // assemble the local products
//        // * therefore create nonparametric local discretization
//        LocalDiscretizationType local_product_discretization(grid_provider_,
//                                                             this->all_dirichlet_boundary_config_,
//                                                             this->zero_boundary_problem_,
//                                                             ss);
//        local_product_discretization.initialize();
//        // and get all of the products
//        for (auto id : local_product_discretization.available_products())
//          local_products_[ss].insert(std::make_pair(id, local_product_discretization.get_product(id)));
        // and create the local containers
        // * the matrices
        //   * just copy thos from the local discretizations
        const auto local_operator = this->local_discretizations_[ss]->get_operator();
        local_matrices_[ss] = std::make_shared< AffinelyDecomposedMatrixType >();
        //   * we take the affine part only if the diffusion has one, otherwise it contains only the dirichlet rows,
        //     thus it is empty, since the local problems are purely neumann
        if (this->problem().diffusion_factor().has_affine_part()) {
          if (!local_operator.has_affine_part())
            DUNE_THROW(Stuff::Exceptions::internal_error, "The local operator is missing the affine part!");
          local_matrices_[ss]->register_affine_part(new MatrixType(*(local_operator.affine_part().container())));
        }
//        if (local_operator.num_components() < problem_->diffusion()->num_components())
//          DUNE_PYMOR_THROW(Pymor::Exception::requirements_not_met,
//                           "The local operator should have " << problem_->diffusion()->num_components()
//                           << " components (but has only " << local_operator.num_components() << ")!");
//        for (DUNE_STUFF_SSIZE_T qq = 0; qq < problem_->diffusion()->num_components(); ++qq)
//          local_matrices_[ss]->register_component(new MatrixType(
//              local_operator.component(qq).container()->backend()),
//              problem_->diffusion()->coefficient(qq));
        // * and the vectors
        const auto local_functional = this->local_discretizations_[ss]->get_rhs();
        local_vectors_[ss] = std::make_shared< AffinelyDecomposedVectorType >();
        //   * first the affine part
        if (this->problem().force().has_affine_part() || this->problem().neumann().has_affine_part()) {
          //   * which we copy from the local discretization
          if (!local_functional.has_affine_part())
            DUNE_THROW(Stuff::Exceptions::internal_error, "The local functional is missing the affine part!");
          local_vectors_[ss]->register_affine_part(new VectorType(*(local_functional.affine_part().container())));
        } else if (this->problem().diffusion_factor().has_affine_part()
                   || this->problem().dirichlet().has_affine_part()) {
          //   * but not, if it was due to the dirichlet boundary correction
          //     In this case we just create an empty one, since the one of the local discretization has to be empty
          local_vectors_[ss]->register_affine_part(local_functional.affine_part().container()->size());
        }
//        //   * then we copy the components from the local discretizations
//        for (DUNE_STUFF_SSIZE_T qq = 0; qq < problem_->force()->num_components(); ++qq)
//          local_vectors_[ss]->register_component(new VectorType(
//              local_functional.component(qq).container()->backend()),
//              problem_->force()->coefficient(qq));
//        for (DUNE_STUFF_SSIZE_T qq = 0; qq < problem_->neumann()->num_components(); ++qq)
//          local_vectors_[ss]->register_component(new VectorType(
//              local_functional.component(problem_->force()->num_components() + qq).container()->backend()),
//              problem_->neumann()->coefficient(qq));
//        //   * and create the components due to the dirichlet boundary term
//        if (problem_->diffusion()->has_affine_part())
//          for (DUNE_STUFF_SSIZE_T qq = 0; qq < problem_->dirichlet()->num_components(); ++qq)
//            local_vectors_[ss]->register_component(new VectorType(local_discretizations_[ss]->testSpace().mapper().size()),
//                                                   problem_->dirichlet()->coefficient(qq));
//        if (problem_->dirichlet()->has_affine_part())
//          for (DUNE_STUFF_SSIZE_T qq = 0; qq < problem_->diffusion()->num_components(); ++qq)
//            local_vectors_[ss]->register_component(new VectorType(local_discretizations_[ss]->testSpace().mapper().size()),
//                                                   problem_->diffusion()->coefficient(qq));
//        if (problem_->diffusion()->num_components() > 0 && problem_->dirichlet()->num_components() > 0) {
//          Pymor::ParameterType diffusion_dirichlet_mu;
//          for (auto key : problem_->diffusion()->parameter_type().keys())
//            diffusion_dirichlet_mu.set(key, problem_->diffusion()->parameter_type().get(key));
//          for (auto key : problem_->dirichlet()->parameter_type().keys())
//            diffusion_dirichlet_mu.set(key, problem_->dirichlet()->parameter_type().get(key));
//          for (DUNE_STUFF_SSIZE_T pp = 0; pp < problem_->diffusion()->num_components(); ++ pp) {
//            for (DUNE_STUFF_SSIZE_T qq = 0; qq < problem_->dirichlet()->num_components(); ++qq) {
//              const std::string expression = "(" + problem_->diffusion()->coefficient(pp)->expression()
//                                             + ")*(" + problem_->dirichlet()->coefficient(qq)->expression() + ")";
//              local_vectors_[ss]->register_component(new VectorType(local_discretizations_[ss]->testSpace().mapper().size()),
//                                                     new Pymor::ParameterFunctional(diffusion_dirichlet_mu,
//                                                                                    expression));
//            }
//          }
//        } // create the local containers

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
            if (this->problem().diffusion_factor().has_affine_part()) {
              inside_outside_matrix->register_affine_part(new MatrixType(inner_test_mapper.size(),
                                                                         outer_ansatz_mapper.size(),
                                                                         inside_outside_pattern));
              outside_inside_matrix->register_affine_part(new MatrixType(outer_test_mapper.size(),
                                                                         inner_ansatz_mapper.size(),
                                                                         outside_inside_pattern));
            }
            for (DUNE_STUFF_SSIZE_T qq = 0; qq < this->problem().diffusion_factor().num_components(); ++qq) {
              inside_outside_matrix->register_component(new MatrixType(inner_test_mapper.size(),
                                                                       outer_ansatz_mapper.size(),
                                                                       inside_outside_pattern),
                                                        this->problem().diffusion_factor().coefficient(qq));
              outside_inside_matrix->register_component(new MatrixType(outer_test_mapper.size(),
                                                                       inner_ansatz_mapper.size(),
                                                                       outside_inside_pattern),
                                                        this->problem().diffusion_factor().coefficient(qq));
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

//      // parameter
//      this->inherit_parameter_type(*(local_matrices_[0]), "lhs");
//      this->inherit_parameter_type(*(local_vectors_[0]), "rhs");

      this->container_based_initialized_ = true;
    } // if (!this->container_based_initialized_)
  } // ... init(...)

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

  template< class M >
  void copy_local_to_global_matrix(const M& local_matrix,
                                   const PatternType& local_pattern,
                                   const size_t test_subdomain,
                                   const size_t ansatz_subdomain,
                                   M& global_matrix) const
  {
    for (size_t local_ii = 0; local_ii < local_pattern.size(); ++local_ii) {
      const size_t global_ii = this->test_space()->mapper().mapToGlobal(test_subdomain, local_ii);
      for (const size_t& local_jj : local_pattern.inner(local_ii)) {
        const size_t global_jj = this->ansatz_space()->mapper().mapToGlobal(ansatz_subdomain, local_jj);
        global_matrix.add_to_entry(global_ii, global_jj, local_matrix.get_entry(local_ii, local_jj));
      }
    }
  } // ... copy_local_to_global_matrix(...)

  template< class V >
  void copy_local_to_global_vector(const V& local_vector,
                                   const size_t subdomain,
                                   V& global_vector) const
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
    typedef typename ProblemType::DiffusionFactorType::NonparametricType  DiffusionType;
    typedef GDT::LocalOperator::Codim1BoundaryIntegral< GDT::LocalEvaluation::SWIPDG::BoundaryLHS< DiffusionType > >
        DirichletOperatorType;
    typedef GDT::LocalAssembler::Codim1BoundaryMatrix< DirichletOperatorType > DirichletMatrixAssemblerType;
    std::vector< DirichletOperatorType* > dirichlet_operators;
    std::vector< DirichletMatrixAssemblerType* > dirichlet_matrix_assemblers;
    for (DUNE_STUFF_SSIZE_T qq = 0; qq < this->problem().diffusion_factor().num_components(); ++qq) {
      dirichlet_operators.push_back(new DirichletOperatorType(*(this->problem().diffusion_factor().component(qq))));
      dirichlet_matrix_assemblers.push_back(new DirichletMatrixAssemblerType(*(dirichlet_operators[qq])));
      boundary_assembler.add(*(dirichlet_matrix_assemblers[qq]),
                             *(local_matrix.component(qq)),
                             new GDT::ApplyOn::DirichletIntersections< BoundaryGridPartType >(this->boundary_info()));
    }
    if (this->problem().diffusion_factor().has_affine_part()) {
      dirichlet_operators.push_back(new DirichletOperatorType(*(this->problem().diffusion_factor().affine_part())));
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
    for (DUNE_STUFF_SSIZE_T qq = 0; qq < this->problem().neumann().num_components(); ++qq) {
      neumann_functionals.push_back(new NeumannFunctionalType(*(this->problem().neumann().component(qq))));
      neumann_vector_assemblers.push_back(new NeumannVectorAssemblerType(*(neumann_functionals[qq])));
      boundary_assembler.add(*(neumann_vector_assemblers[qq]),
                             *(local_vector.component(this->problem().force().num_components() + qq)),
                             new GDT::ApplyOn::NeumannIntersections< BoundaryGridPartType >(this->boundary_info()));
    }
    if (this->problem().neumann().has_affine_part()) {
      neumann_functionals.push_back(new NeumannFunctionalType(*(this->problem().neumann().affine_part())));
      neumann_vector_assemblers.push_back(new NeumannVectorAssemblerType(*(
          neumann_functionals[neumann_functionals.size() - 1])));
      boundary_assembler.add(*(neumann_vector_assemblers[neumann_vector_assemblers.size() - 1]),
                             *(local_vector.affine_part()),
                             new GDT::ApplyOn::NeumannIntersections< BoundaryGridPartType >(this->boundary_info()));
    }

    // * dirichlet boundary terms
    typedef typename ProblemType::FunctionType::NonparametricType  DirichletType;
    typedef GDT::LocalFunctional::Codim1Integral< GDT::LocalEvaluation::SWIPDG::BoundaryRHS< DiffusionType,
                                                                                            DirichletType > >
        DirichletFunctionalType;
    typedef GDT::LocalAssembler::Codim1Vector< DirichletFunctionalType > DirichletVectorAssemblerType;
    std::vector< DirichletFunctionalType* > dirichlet_functionals;
    std::vector< DirichletVectorAssemblerType* > dirichlet_vector_assemblers;
    size_t component_index = this->problem().force().num_components() + this->problem().neumann().num_components();
    if (this->problem().diffusion_factor().has_affine_part()) {
      for (DUNE_STUFF_SSIZE_T qq = 0; qq < this->problem().dirichlet().num_components(); ++qq) {
        dirichlet_functionals.push_back(new DirichletFunctionalType(*(this->problem().diffusion_factor().affine_part()),
                                                                    *(this->problem().dirichlet().component(qq))));
        dirichlet_vector_assemblers.push_back(new DirichletVectorAssemblerType(*(
            dirichlet_functionals[dirichlet_functionals.size() - 1])));
        boundary_assembler.add(*(dirichlet_vector_assemblers[dirichlet_vector_assemblers.size() - 1]),
                               *(local_vector.component(component_index)),
                               new GDT::ApplyOn::DirichletIntersections< BoundaryGridPartType >(this->boundary_info()));
        ++component_index;
      }
    }
    if (this->problem().dirichlet().has_affine_part()) {
      for (DUNE_STUFF_SSIZE_T qq = 0; qq < this->problem().diffusion_factor().num_components(); ++qq) {
        dirichlet_functionals.push_back(new DirichletFunctionalType(*(this->problem().diffusion_factor().component(qq)),
                                                                    *(this->problem().dirichlet().affine_part())));
        dirichlet_vector_assemblers.push_back(new DirichletVectorAssemblerType(*(
            dirichlet_functionals[dirichlet_functionals.size() - 1])));
        boundary_assembler.add(*(dirichlet_vector_assemblers[dirichlet_vector_assemblers.size() - 1]),
                               *(local_vector.component(component_index)),
                               new GDT::ApplyOn::DirichletIntersections< BoundaryGridPartType >(this->boundary_info()));
        ++component_index;
      }
    }
    for (DUNE_STUFF_SSIZE_T pp = 0; pp < this->problem().diffusion_factor().num_components(); ++ pp) {
      for (DUNE_STUFF_SSIZE_T qq = 0; qq < this->problem().dirichlet().num_components(); ++qq) {
        dirichlet_functionals.push_back(new DirichletFunctionalType(*(this->problem().diffusion_factor().component(pp)),
                                                                    *(this->problem().dirichlet().component(qq))));
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

    typedef typename ProblemType::DiffusionFactorType::NonparametricType DiffusionType;
    typedef GDT::LocalOperator::Codim1CouplingIntegral< GDT::LocalEvaluation::SWIPDG::Inner< DiffusionType > >
        CouplingOperatorType;
    typedef GDT::LocalAssembler::Codim1CouplingMatrix< CouplingOperatorType > CouplingMatrixAssemblerType;
    std::vector< CouplingOperatorType* > coupling_operators;
    std::vector< CouplingMatrixAssemblerType* > coupling_matrix_assemblers;
    for (DUNE_STUFF_SSIZE_T qq = 0; qq < this->problem().diffusion_factor().num_components(); ++qq) {
      coupling_operators.push_back(new CouplingOperatorType(*(this->problem().diffusion_factor().component(qq))));
      coupling_matrix_assemblers.push_back(new CouplingMatrixAssemblerType(*(coupling_operators[qq])));
      coupling_assembler.addLocalAssembler(*(coupling_matrix_assemblers[qq]),
                                           *(inside_inside_matrix.component(qq)),
                                           *(inside_outside_matrix.component(qq)),
                                           *(outside_inside_matrix.component(qq)),
                                           *(outside_outside_matrix.component(qq)));
    }
    if (this->problem().diffusion_factor().has_affine_part()) {
      coupling_operators.push_back(new CouplingOperatorType(*(this->problem().diffusion_factor().affine_part())));
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
    // initialize global matrix
    this->matrix_ = std::make_shared< AffinelyDecomposedMatrixType >();
    auto& system_matrix = *(this->matrix_);
    if (local_matrices_[0]->has_affine_part())
      system_matrix.register_affine_part(new MatrixType(this->test_space()->mapper().size(),
                                                        this->ansatz_space()->mapper().size(),
                                                        pattern_));
    for (DUNE_STUFF_SSIZE_T qq = 0; qq < local_matrices_[0]->num_components(); ++qq)
      system_matrix.register_component(new MatrixType(this->test_space()->mapper().size(),
                                                      this->ansatz_space()->mapper().size(),
                                                      pattern_),
                                       local_matrices_[0]->coefficient(qq));
    // walk the subdomains
    for (size_t ss = 0; ss < ms_grid_->size(); ++ss) {
      const auto& local_matrix = *(local_matrices_[ss]);
      if (local_matrix.has_affine_part())
        copy_local_to_global_matrix(*(local_matrix.affine_part()),
                                    this->local_discretizations_[ss]->pattern(),
                                    ss, ss,
                                    *(system_matrix.affine_part()));
      for (DUNE_STUFF_SSIZE_T qq = 0; qq < local_matrix.num_components(); ++qq) {
        copy_local_to_global_matrix(*(local_matrix.component(qq)),
                                    this->local_discretizations_[ss]->pattern(),
                                    ss, ss,
                                    *(system_matrix.component(qq)));
      }

      // walk the neighbours
      for (const size_t& nn : ms_grid_->neighborsOf(ss)) {
        if (ss < nn) {
          // get the coupling patterns
          const auto result_inside_outside_pattern = inside_outside_patterns_[ss].find(nn);
          if (result_inside_outside_pattern == inside_outside_patterns_[ss].end())
            DUNE_PYMOR_THROW(Pymor::Exception::this_does_not_make_any_sense,
                             "The coupling pattern for subdomain " << ss << " and neighbour " << nn << "is missing!");
          const auto& inside_outside_pattern = *(result_inside_outside_pattern->second);
          const auto result_outside_inside_pattern = outside_inside_patterns_[nn].find(ss);
          if (result_outside_inside_pattern == outside_inside_patterns_[nn].end())
            DUNE_PYMOR_THROW(Pymor::Exception::this_does_not_make_any_sense,
                             "The coupling pattern for neighbour " << nn << " and subdomain " << ss << "is missing!");
          const auto& outside_inside_pattern = *(result_outside_inside_pattern->second);
          // and the coupling matrices
          auto result_inside_outside_matrix = inside_outside_matrices_[ss].find(nn);
          if (result_inside_outside_matrix == inside_outside_matrices_[ss].end())
            DUNE_PYMOR_THROW(Pymor::Exception::this_does_not_make_any_sense,
                             "The coupling matrix for subdomain " << ss << " and neighbour " << nn << "is missing!");
          auto& inside_outside_matrix = *(result_inside_outside_matrix->second);
          auto result_outside_inside_matrix = outside_inside_matrices_[nn].find(ss);
          if (result_outside_inside_matrix == outside_inside_matrices_[nn].end())
            DUNE_PYMOR_THROW(Pymor::Exception::this_does_not_make_any_sense,
                             "The coupling matrix for neighbour " << nn << " and subdomain " << ss << "is missing!");
          auto& outside_inside_matrix = *(result_outside_inside_matrix->second);
          // and copy them into the global matrix
          if (inside_outside_matrix.has_affine_part())
            copy_local_to_global_matrix(*(inside_outside_matrix.affine_part()),
                                        inside_outside_pattern,
                                        ss, nn,
                                        *(system_matrix.affine_part()));
          for (DUNE_STUFF_SSIZE_T qq = 0; qq < inside_outside_matrix.num_components(); ++qq) {
            copy_local_to_global_matrix(*(inside_outside_matrix.component(qq)),
                                        inside_outside_pattern,
                                        ss, nn,
                                        *(system_matrix.component(qq)));
          }
          if (outside_inside_matrix.has_affine_part())
            copy_local_to_global_matrix(*(outside_inside_matrix.affine_part()),
                                        outside_inside_pattern,
                                        nn, ss,
                                        *(system_matrix.affine_part()));
          for (DUNE_STUFF_SSIZE_T qq = 0; qq < outside_inside_matrix.num_components(); ++qq) {
            copy_local_to_global_matrix(*(outside_inside_matrix.component(qq)),
                                        outside_inside_pattern,
                                        nn, ss,
                                        *(system_matrix.component(qq)));
          }
        }
      } // walk the neighbours
    } // walk the subdomains

    // initialize global vector
    this->rhs_ = std::make_shared< AffinelyDecomposedVectorType >();
    auto& rhs_vector = *(this->rhs_);
    if (local_vectors_[0]->has_affine_part())
      rhs_vector.register_affine_part(new VectorType(this->test_space()->mapper().size()));
    for (DUNE_STUFF_SSIZE_T qq = 0; qq < local_vectors_[0]->num_components(); ++qq)
      rhs_vector.register_component(new VectorType(this->test_space()->mapper().size()),
                                    local_vectors_[0]->coefficient(qq));
    // walk the subdomains
    for (size_t ss = 0; ss < ms_grid_->size(); ++ss) {
      auto local_vector = *(local_vectors_[ss]);
      if (local_vector.has_affine_part())
        copy_local_to_global_vector(*(local_vector.affine_part()),
                                    ss,
                                    *(rhs_vector.affine_part()));
      for (DUNE_STUFF_SSIZE_T qq = 0; qq < local_vector.num_components(); ++qq)
        copy_local_to_global_vector(*(local_vector.component(qq)),
                                    ss,
                                    *(rhs_vector.component(qq)));
    } // walk the subdomains
  } // ... build_global_containers(...)

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
