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
    , local_boundary_patterns_(ms_grid_->size())
    , local_coupling_patterns_(ms_grid_->size())
    , inside_outside_coupling_patterns_(ms_grid_->size())
    , outside_inside_coupling_patterns_(ms_grid_->size())
    , boundary_matrices_(ms_grid_->size(), nullptr)
    , boundary_vectors_(ms_grid_->size(), nullptr)
    , local_coupling_matrices_(ms_grid_->size(), nullptr)
    , inside_outside_coupling_matrices_(ms_grid_->size())
    , outside_inside_coupling_matrices_(ms_grid_->size())
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
      out << prefix << "assembling locally on " << subdomains << " subdomains... " << std::flush;
      Dune::Timer timer;
      for (size_t ss = 0; ss < subdomains; ++ss) {
        // init the local discretizations (assembles matrices and patterns)
        this->local_discretizations_[ss]->init();

        // create and copy the local patterns
        add_local_to_global_pattern(this->local_discretizations_[ss]->pattern(), ss, ss);

        const auto inner_test_space = this->local_discretizations_[ss]->test_space();
        const auto inner_ansatz_space = this->local_discretizations_[ss]->ansatz_space();

        // assemble boundary
        if (ms_grid_->boundary(ss)) {
          // create local boundary patterns
          local_boundary_patterns_[ss] = inner_test_space->compute_volume_pattern(ms_grid_->boundaryGridPart(ss)->gridView(),
                                                                                  *inner_ansatz_space);
          assemble_boundary_contributions(ss);
        }

        // walk the neighbors
        for (const size_t& nn : ms_grid_->neighborsOf(ss)) {
          // visit each coupling only once (assemble primally)
          if (ss < nn) {
            const auto outer_test_space = this->local_discretizations_[nn]->test_space();
            const auto outer_ansatz_space = this->local_discretizations_[nn]->ansatz_space();
            const auto inside_outside_grid_view = ms_grid_->couplingGridPart(ss, nn);
            const auto outside_inside_grid_view = ms_grid_->couplingGridPart(nn, ss);
            // create the coupling patterns
            // * those are just a subset of the patterns from the local discretizations
            //   so the copy later on is cheaper
            //   NOTE: it is untested, if the additional grid walk that is done in 'compute_volume_pattern' is really
            //         cheaper than unneccessarily copying the whole local matrix
            local_coupling_patterns_[ss]
                = inner_test_space->compute_volume_pattern(inside_outside_grid_view->gridView(),
                                                           *inner_ansatz_space);
            local_coupling_patterns_[nn]
                = outer_test_space->compute_volume_pattern(outside_inside_grid_view->gridView(),
                                                           *outer_ansatz_space);
            // * those two are really needed
            inside_outside_coupling_patterns_[ss][nn]
                = inner_test_space->compute_face_pattern(*inside_outside_grid_view,
                                                         *outer_ansatz_space);
            outside_inside_coupling_patterns_[nn][ss]
                = outer_test_space->compute_face_pattern(*outside_inside_grid_view,
                                                         *inner_ansatz_space);
            // and copy whats needed
            add_local_to_global_pattern(inside_outside_coupling_patterns_[ss][nn], ss, nn);
            add_local_to_global_pattern(outside_inside_coupling_patterns_[nn][ss], nn, ss);

            // assemble the coupling
            assemble_coupling_contributions(ss, nn);

          } // visit each coupling only once (assemble primaly)
        } // walk the neighbors
      } // walk the subdomains for the first time
      out<< "done (took " << timer.elapsed() << "s)" << std::endl;

      // build the global matrix and vector
      // NOTE: we could just clear the local containers after copy to save memory
      out << prefix << "assembling global containers... " << std::flush;
      timer.reset();
      this->matrix_ = std::make_shared< AffinelyDecomposedMatrixType >();
      this->rhs_= std::make_shared< AffinelyDecomposedVectorType >();
      // walk the subdomains
      for (size_t ss = 0; ss < ms_grid_->size(); ++ss) {
        // copy the local containers
        copy_local_to_global_matrix(*(this->local_discretizations_[ss]->matrix_),
                                    this->local_discretizations_[ss]->pattern(),
                                    ss,
                                    ss,
                                    *(this->matrix_));
        copy_local_to_global_vector(*(this->local_discretizations_[ss]->rhs_), ss, *(this->rhs_));

        // copy the boundary containers
        if (ms_grid_->boundary(ss)) {
          copy_local_to_global_matrix(*(boundary_matrices_[ss]),
                                      local_boundary_patterns_[ss],
                                      ss,
                                      ss,
                                      *(this->matrix_));
          // no need as long as we are dirichlet zero
//          copy_local_to_global_vector(*(boundary_vectors_[ss]), ss, *(this->rhs_));
        }

//        // walk the neighbours
//        for (const size_t& nn : ms_grid_->neighborsOf(ss)) {
//          if (ss < nn) {
//            // copy the coupling containers
//            copy_local_to_global_matrix(*(local_coupling_matrices_[ss]),
//                                        local_coupling_patterns_[ss],
//                                        ss,
//                                        ss,
//                                        *(this->matrix_));
//            copy_local_to_global_matrix(*(local_coupling_matrices_[nn]),
//                                        local_coupling_patterns_[nn],
//                                        nn,
//                                        nn,
//                                        *(this->matrix_));
//            copy_local_to_global_matrix(*(inside_outside_coupling_matrices_[ss][nn]),
//                                        inside_outside_coupling_patterns_[ss][nn],
//                                        ss,
//                                        nn,
//                                        *(this->matrix_));
//            copy_local_to_global_matrix(*(outside_inside_coupling_matrices_[nn][ss]),
//                                        outside_inside_coupling_patterns_[nn][ss],
//                                        nn,
//                                        ss,
//                                        *(this->matrix_));
//          }
//        }
      }
      // walk the subdomains
      out << "done (took " << timer.elapsed() << "s)" << std::endl;

//      if (!this->matrix_->parametric() && this->matrix_->has_affine_part())
//        std::cout << "\n" << *(this->matrix_->affine_part()) << std::endl;

//      if (!this->rhs_->parametric() && this->rhs_->has_affine_part())
//        std::cout << "\n" << *(this->rhs_->affine_part()) << std::endl;

      // build parameter type
      this->inherit_parameter_type(*(this->matrix_), "lhs");
      this->inherit_parameter_type(*(this->rhs_), "rhs");

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

  void add_local_to_global_pattern(const PatternType& local, const size_t test_subdomain, const size_t ansatz_subdomain)
  {
    for (size_t local_ii = 0; local_ii < local.size(); ++local_ii) {
      const size_t global_ii = this->test_space_->mapper().mapToGlobal(test_subdomain, local_ii);
      auto& global_rows = pattern_.inner(global_ii);
      for (const size_t& local_jj : local.inner(local_ii)) {
        const size_t global_jj = this->ansatz_space_->mapper().mapToGlobal(ansatz_subdomain, local_jj);
        global_rows.insert(global_jj);
      }
    }
  } // ... add_local_to_global_pattern(...)

  void assemble_boundary_contributions(const size_t subdomain)
  {
    using namespace GDT;

    const auto local_test_space = this->local_discretizations_[subdomain]->test_space();
    const auto local_ansatz_space = this->local_discretizations_[subdomain]->ansatz_space();
    const auto& local_pattern = local_boundary_patterns_[subdomain];
    assert(local_pattern.size() > 0);

    const auto& boundary_grid_view = ms_grid_->boundaryGridPart(subdomain)->gridView();
    typedef typename MsGridType::BoundaryGridPartType::GridViewType BoundaryGridViewType;
    SystemAssembler< typename LocalDiscretizationType::TestSpaceType,
                     BoundaryGridViewType,
                     typename LocalDiscretizationType::AnsatzSpaceType >
        boundary_assembler(*local_test_space, *local_ansatz_space, boundary_grid_view);
    Stuff::Grid::BoundaryInfos::AllDirichlet< typename BoundaryGridViewType::Intersection > boundary_info;

    boundary_matrices_[subdomain] = std::make_shared< AffinelyDecomposedMatrixType >();
    boundary_vectors_[subdomain] = std::make_shared< AffinelyDecomposedVectorType >();
    auto& boundary_matrix = *(boundary_matrices_[subdomain]);
    auto& boundary_vector = *(boundary_vectors_[subdomain]);

    // boundary operator
    typedef typename ProblemType::DiffusionFactorType::NonparametricType DiffusionFactorType;
    typedef typename ProblemType::DiffusionTensorType::NonparametricType DiffusionTensorType;
    const auto& diffusion_factor = this->problem_.diffusion_factor();
    const auto& diffusion_tensor = this->problem_.diffusion_tensor();
    typedef LocalOperator::Codim1BoundaryIntegral< LocalEvaluation::SWIPDG::BoundaryLHS< DiffusionFactorType
                                                                                       , DiffusionTensorType > >
        DirichletBoundaryOperatorType;
    std::vector< std::unique_ptr< DirichletBoundaryOperatorType > > dirichlet_boundary_operators;
    typedef LocalAssembler::Codim1BoundaryMatrix< DirichletBoundaryOperatorType >
        DirichletMatrixAssemblerType;
    std::vector< std::unique_ptr< DirichletMatrixAssemblerType > > dirichlet_matrix_assemblers;
    for (size_t qq = 0; qq < diffusion_factor.num_components(); ++qq) {
      const size_t id = boundary_matrix.register_component(diffusion_factor.coefficient(qq),
                                                           local_test_space->mapper().size(),
                                                           local_ansatz_space->mapper().size(),
                                                           local_pattern);
      dirichlet_boundary_operators.emplace_back(new DirichletBoundaryOperatorType(*(diffusion_factor.component(qq)),
                                                                                  *(diffusion_tensor.affine_part())));
      dirichlet_matrix_assemblers.emplace_back(new DirichletMatrixAssemblerType(*(dirichlet_boundary_operators.back())));
      boundary_assembler.add(*(dirichlet_matrix_assemblers.back()),
                             *(boundary_matrix.component(id)),
                             new ApplyOn::DirichletIntersections< BoundaryGridViewType >(boundary_info));
    }
    if (diffusion_factor.has_affine_part()) {
      if (!boundary_matrix.has_affine_part())
        boundary_matrix.register_affine_part(local_test_space->mapper().size(),
                                             local_ansatz_space->mapper().size(),
                                             local_pattern);
      dirichlet_boundary_operators.emplace_back(new DirichletBoundaryOperatorType(*(diffusion_factor.affine_part()),
                                                                                  *(diffusion_tensor.affine_part())));
      dirichlet_matrix_assemblers.emplace_back(new DirichletMatrixAssemblerType(*(dirichlet_boundary_operators.back())));
      boundary_assembler.add(*(dirichlet_matrix_assemblers.back()),
                             *(boundary_matrix.affine_part()),
                             new ApplyOn::DirichletIntersections< BoundaryGridViewType >(boundary_info));
    }

    // no need as long as we are dirichlet zero
//    // dirichlet boundary functional
//    typedef typename ProblemType::FunctionType::NonparametricType FunctionType;
//    const auto& dirichlet = this->zero_boundary_problem_.dirichlet();
//    typedef LocalFunctional::Codim1Integral< LocalEvaluation::SWIPDG::BoundaryRHS< DiffusionFactorType, FunctionType
//                                                                                 , DiffusionTensorType > >
//        DirichletFunctionalType;
//    std::vector< std::unique_ptr< DirichletFunctionalType > > dirichlet_functionals;
//    typedef LocalAssembler::Codim1Vector< DirichletFunctionalType >
//        DirichletVectorAssemblerType;
//    std::vector< std::unique_ptr< DirichletVectorAssemblerType > > dirichlet_vector_assemblers;
//    if (diffusion_factor.has_affine_part() && dirichlet.has_affine_part()) {
//      if (!boundary_vector.has_affine_part())
//        boundary_vector.register_affine_part(local_test_space->mapper().size());
//      dirichlet_functionals.emplace_back(new DirichletFunctionalType(*(diffusion_factor.affine_part()),
//                                                                     *(diffusion_tensor.affine_part()),
//                                                                     *(dirichlet.affine_part())));
//      dirichlet_vector_assemblers.emplace_back(new DirichletVectorAssemblerType(*(dirichlet_functionals.back())));
//      boundary_assembler.add(*(dirichlet_vector_assemblers.back()),
//                             *(boundary_vector.affine_part()),
//                             new ApplyOn::DirichletIntersections< BoundaryGridViewType >(boundary_info));
//    }
//    if (diffusion_factor.has_affine_part()) {
//      for (size_t qq = 0; qq < dirichlet.num_components(); ++qq) {
//        const size_t id = boundary_vector.register_component(dirichlet.coefficient(qq),
//                                                             local_test_space->mapper().size());
//        dirichlet_functionals.emplace_back(new DirichletFunctionalType(*(diffusion_factor.affine_part()),
//                                                                       *(diffusion_tensor.affine_part()),
//                                                                       *(dirichlet.component(qq))));
//        dirichlet_vector_assemblers.emplace_back(new DirichletVectorAssemblerType(*(dirichlet_functionals.back())));
//        boundary_assembler.add(*(dirichlet_vector_assemblers.back()),
//                               *(boundary_vector.component(id)),
//                               new ApplyOn::DirichletIntersections< BoundaryGridViewType >(boundary_info));
//      }
//    }
//    if (dirichlet.has_affine_part()) {
//      for (size_t qq = 0; qq < diffusion_factor.num_components(); ++qq) {
//        const size_t id = boundary_vector.register_component(diffusion_factor.coefficient(qq),
//                                                             local_test_space->mapper().size());
//        dirichlet_functionals.emplace_back(new DirichletFunctionalType(*(diffusion_factor.component(qq)),
//                                                                       *(diffusion_tensor.affine_part()),
//                                                                       *(dirichlet.affine_part())));
//        dirichlet_vector_assemblers.emplace_back(new DirichletVectorAssemblerType(*(dirichlet_functionals.back())));
//        boundary_assembler.add(*(dirichlet_vector_assemblers.back()),
//                               *(boundary_vector.component(id)),
//                               new ApplyOn::DirichletIntersections< BoundaryGridViewType >(boundary_info));
//      }
//    }
//    Pymor::ParameterType param;
//    for (const auto& key : diffusion_factor.parameter_type().keys())
//      param.set(key, diffusion_factor.parameter_type().get(key));
//    for (const auto& key : dirichlet.parameter_type().keys())
//      param.set(key, dirichlet.parameter_type().get(key));
//    for (size_t pp = 0; pp < diffusion_factor.num_components(); ++ pp) {
//      for (size_t qq = 0; qq < dirichlet.num_components(); ++qq) {
//        const std::string expression = "(" + diffusion_factor.coefficient(pp)->expression()
//                                           + ")*(" + dirichlet.coefficient(qq)->expression() + ")";
//        const size_t id = boundary_vector.register_component(param, expression, local_test_space->mapper().size());
//        dirichlet_functionals.emplace_back(new DirichletFunctionalType(*(diffusion_factor.component(pp)),
//                                                                       *(diffusion_tensor.affine_part()),
//                                                                       *(dirichlet.component(qq))));
//        dirichlet_vector_assemblers.emplace_back(new DirichletVectorAssemblerType(*(dirichlet_functionals.back())));
//        boundary_assembler.add(*(dirichlet_vector_assemblers.back()),
//                               *(boundary_vector.component(id)),
//                               new ApplyOn::DirichletIntersections< BoundaryGridViewType >(boundary_info));
//      }
//    }

    // this is still the old variant
//    // rhs
//    // * neumann boundary terms
//    typedef typename ProblemType::NeumannType::NonparametricType NeumannType;
//    typedef GDT::LocalFunctional::Codim1Integral< GDT::LocalEvaluation::Product< NeumannType > > NeumannFunctionalType;
//    typedef GDT::LocalAssembler::Codim1Vector< NeumannFunctionalType > NeumannVectorAssemblerType;
//    std::vector< NeumannFunctionalType* > neumann_functionals;
//    std::vector< NeumannVectorAssemblerType* > neumann_vector_assemblers;
//    for (DUNE_STUFF_SSIZE_T qq = 0; qq < problem_.neumann()->num_components(); ++qq) {
//      neumann_functionals.push_back(new NeumannFunctionalType(*(problem_.neumann()->component(qq))));
//      neumann_vector_assemblers.push_back(new NeumannVectorAssemblerType(*(neumann_functionals[qq])));
//      boundary_assembler.addLocalAssembler(*(neumann_vector_assemblers[qq]),
//                                           typename BoundaryAssemblerType::AssembleOnNeumann(*(global_boundary_info_)),
//                                           *(local_vector.component(problem_.force()->num_components() + qq)));
//    }
//    if (problem_.neumann()->has_affine_part()) {
//      neumann_functionals.push_back(new NeumannFunctionalType(*(problem_.neumann()->affine_part())));
//      neumann_vector_assemblers.push_back(new NeumannVectorAssemblerType(*(
//          neumann_functionals[neumann_functionals.size() - 1])));
//      boundary_assembler.addLocalAssembler(*(neumann_vector_assemblers[neumann_vector_assemblers.size() - 1]),
//                                           typename BoundaryAssemblerType::AssembleOnNeumann(*(global_boundary_info_)),
//                                           *(local_vector.affine_part()));
//    }

    // do the actual work
    boundary_assembler.assemble();
  } // ... assemble_boundary_contributions(...)

  void assemble_coupling_contributions(const size_t subdomain, const size_t neighbour)
  {
    using namespace GDT;





    {
//    const auto inner_test_space   = this->local_discretizations_[subdomain]->test_space();
//    const auto inner_ansatz_space = this->local_discretizations_[subdomain]->ansatz_space();
//    const auto outer_test_space   = this->local_discretizations_[neighbour]->test_space();
//    const auto outer_ansatz_space = this->local_discretizations_[neighbour]->ansatz_space();
//    CouplingAssembler coupling_assembler(*inner_test_space, *inner_ansatz_space,
//                                         *outer_test_space, *outer_ansatz_space,
//                                         *(ms_grid_->couplingGridPart(subdomain, neighbour)));

//    const auto& inside_inside_pattern = local_coupling_patterns_[subdomain];
//    const auto& outside_outside_pattern = local_coupling_patterns_[neighbour];
//    const auto& inside_outside_pattern = inside_outside_coupling_patterns_[subdomain][neighbour];
//    assert(inside_outside_pattern.size() > 0);
//    const auto& outside_inside_pattern = outside_inside_coupling_patterns_[neighbour][subdomain];
//    assert(outside_inside_pattern.size() > 0);

//    local_coupling_matrices_[subdomain] = std::make_shared< AffinelyDecomposedMatrixType >();
//    local_coupling_matrices_[neighbour] = std::make_shared< AffinelyDecomposedMatrixType >();
//    inside_outside_coupling_matrices_[subdomain][neighbour] = std::make_shared< AffinelyDecomposedMatrixType >();
//    outside_inside_coupling_matrices_[neighbour][subdomain] = std::make_shared< AffinelyDecomposedMatrixType >();
//    auto& inside_inside_matrix = *(local_coupling_matrices_[subdomain]);
//    auto& outside_outside_matrix = *(local_coupling_matrices_[neighbour]);
//    auto& inside_outside_matrix = *(inside_outside_coupling_matrices_[subdomain][neighbour]);
//    auto& outside_inside_matrix = *(outside_inside_coupling_matrices_[neighbour][subdomain]);

//    typedef typename ProblemType::DiffusionFactorType::NonparametricType DiffusionFactorType;
//    typedef typename ProblemType::DiffusionTensorType::NonparametricType DiffusionTensorType;
//    const auto& diffusion_factor = this->problem_.diffusion_factor();
//    const auto& diffusion_tensor = this->problem_.diffusion_tensor();
//    typedef LocalOperator::Codim1CouplingIntegral< LocalEvaluation::SWIPDG::Inner< DiffusionFactorType
//                                                                                 , DiffusionTensorType > >
//        CouplingOperatorType;
//    std::vector< std::unique_ptr< CouplingOperatorType > > coupling_operators;
//    typedef LocalAssembler::Codim1CouplingMatrix< CouplingOperatorType > CouplingAssemblerType;
//    std::vector< std::unique_ptr< CouplingAssemblerType > > coupling_assemblers;
//    for (size_t qq = 0; qq < diffusion_factor.num_components(); ++qq) {
//      const size_t in_in_id = inside_inside_matrix.register_component(diffusion_factor.coefficient(qq),
//                                                                      inner_test_space->mapper().size(),
//                                                                      inner_test_space->mapper().size(),
//                                                                      inside_inside_pattern);
//      const size_t out_out_id = outside_outside_matrix.register_component(diffusion_factor.coefficient(qq),
//                                                                          outer_test_space->mapper().size(),
//                                                                          outer_test_space->mapper().size(),
//                                                                          outside_outside_pattern);
//      const size_t in_out_id = inside_outside_matrix.register_component(diffusion_factor.coefficient(qq),
//                                                                        inner_test_space->mapper().size(),
//                                                                        outer_test_space->mapper().size(),
//                                                                        inside_outside_pattern);
//      const size_t out_in_id = outside_inside_matrix.register_component(diffusion_factor.coefficient(qq),
//                                                                        outer_test_space->mapper().size(),
//                                                                        inner_test_space->mapper().size(),
//                                                                        outside_inside_pattern);
//      coupling_operators.emplace_back(new CouplingOperatorType(*(diffusion_factor.component(qq)),
//                                                               *(diffusion_tensor.affine_part())));
//      coupling_assemblers.emplace_back(new CouplingAssemblerType(*(coupling_operators.back())));
//      coupling_assembler.addLocalAssembler(*(coupling_assemblers.back()),
//                                           *(inside_inside_matrix.component(in_in_id)),
//                                           *(inside_outside_matrix.component(in_out_id)),
//                                           *(outside_inside_matrix.component(out_in_id)),
//                                           *(outside_outside_matrix.component(out_out_id)));
//    }
//    if (diffusion_factor.has_affine_part()) {
//      if (!inside_inside_matrix.has_affine_part())
//        inside_inside_matrix.register_affine_part(inner_test_space->mapper().size(),
//                                                  inner_test_space->mapper().size(),
//                                                  inside_inside_pattern);
//      if (!outside_outside_matrix.has_affine_part())
//        outside_outside_matrix.register_affine_part(outer_test_space->mapper().size(),
//                                                    outer_test_space->mapper().size(),
//                                                    outside_outside_pattern);
//      if(!inside_outside_matrix.has_affine_part())
//        inside_outside_matrix.register_affine_part(inner_test_space->mapper().size(),
//                                                   outer_test_space->mapper().size(),
//                                                   inside_outside_pattern);
//      if (!outside_inside_matrix.has_affine_part())
//        outside_inside_matrix.register_affine_part(outer_test_space->mapper().size(),
//                                                   inner_test_space->mapper().size(),
//                                                   outside_inside_pattern);
//      coupling_operators.emplace_back(new CouplingOperatorType(*(diffusion_factor.affine_part()),
//                                                               *(diffusion_tensor.affine_part())));
//      coupling_assemblers.emplace_back(new CouplingAssemblerType(*(coupling_operators.back())));
//      coupling_assembler.addLocalAssembler(*(coupling_assemblers.back()),
//                                           *(inside_inside_matrix.affine_part()),
//                                           *(inside_outside_matrix.affine_part()),
//                                           *(outside_inside_matrix.affine_part()),
//                                           *(outside_outside_matrix.affine_part()));
//    }

//    // do the actual work
//    coupling_assembler.assemble();
    }
  } // ... assemble_coupling_contributions(...)

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
                                                this->test_space_->mapper().size(),
                                                this->ansatz_space_->mapper().size(),
                                                pattern_);
      assert(comp >= 0);
      copy_nonparametric_local_to_global_matrix(*(local_matrix.component(qq)),
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
      copy_nonparametric_local_to_global_matrix(*(local_matrix.affine_part()),
                                                local_pattern,
                                                subdomain,
                                                neighbor,
                                                *(global_matrix.affine_part()));
    }
  } // copy_local_to_global_matrix(...)

  void copy_local_to_global_vector(const AffinelyDecomposedConstVectorType& local_vector,
                                   const size_t subdomain,
                                   AffinelyDecomposedVectorType& global_vector) const
  {
    for (size_t qq = 0; qq < local_vector.num_components(); ++qq) {
      const auto coefficient = local_vector.coefficient(qq);
      ssize_t comp = find_component(global_vector, *coefficient);
      if (comp < 0)
        comp = global_vector.register_component(coefficient,
                                                this->test_space_->mapper().size());
      assert(comp >= 0);
      copy_nonparametric_local_to_global_vector(*(local_vector.component(qq)),
                                                subdomain,
                                                *(global_vector.component(comp)));
    }
    if (local_vector.has_affine_part()) {
      if (!global_vector.has_affine_part())
        global_vector.register_affine_part(this->test_space_->mapper().size());
      copy_nonparametric_local_to_global_vector(*(local_vector.affine_part()),
                                                subdomain,
                                                *(global_vector.affine_part()));
    }
  } // copy_local_to_global_vector(...)

  template< class AffinelyDecomposedContainerType >
  ssize_t find_component(const AffinelyDecomposedContainerType& container,
                         const Pymor::ParameterFunctional& coefficient) const
  {
    for (size_t qq = 0; qq < container.num_components(); ++qq)
      if (*(container.coefficient(qq)) == coefficient)
        return qq;
    return -1;
  } // ... find_component(...)

  template< class M >
  void copy_nonparametric_local_to_global_matrix(const Stuff::LA::MatrixInterface< M >& local_matrix,
                                                 const PatternType& local_pattern,
                                                 const size_t test_subdomain,
                                                 const size_t ansatz_subdomain,
                                                 Stuff::LA::MatrixInterface< M >& global_matrix) const
  {
    for (size_t local_ii = 0; local_ii < local_pattern.size(); ++local_ii) {
      const size_t global_ii = this->test_space_->mapper().mapToGlobal(test_subdomain, local_ii);
      for (const size_t& local_jj : local_pattern.inner(local_ii)) {
        const size_t global_jj = this->ansatz_space_->mapper().mapToGlobal(ansatz_subdomain, local_jj);
        global_matrix.add_to_entry(global_ii, global_jj, local_matrix.get_entry(local_ii, local_jj));
      }
    }
  } // ... copy_nonparametric_local_to_global_matrix(...)

  template< class V >
  void copy_nonparametric_local_to_global_vector(const Stuff::LA::VectorInterface< V >& local_vector,
                                                 const size_t subdomain,
                                                 Stuff::LA::VectorInterface< V >& global_vector) const
  {
    for (size_t local_ii = 0; local_ii < local_vector.size(); ++local_ii) {
      const size_t global_ii = this->test_space_->mapper().mapToGlobal(subdomain, local_ii);
      global_vector.add_to_entry(global_ii, local_vector.get_entry(local_ii));
    }
  } // ... copy_nonparametric_local_to_global_vector(...)

  const GridProviderType& grid_provider_;
  std::shared_ptr< const MsGridType > ms_grid_;
  PatternType pattern_;
//  std::vector< std::map< std::string, LocalProductType > > local_products_;
  std::vector< PatternType > local_boundary_patterns_;
  std::vector< PatternType > local_coupling_patterns_;
  std::vector< std::map< size_t, PatternType > > inside_outside_coupling_patterns_;
  std::vector< std::map< size_t, PatternType > > outside_inside_coupling_patterns_;
  std::vector< std::shared_ptr< AffinelyDecomposedMatrixType > > boundary_matrices_;
  std::vector< std::shared_ptr< AffinelyDecomposedVectorType > > boundary_vectors_;
  std::vector< std::shared_ptr< AffinelyDecomposedMatrixType > > local_coupling_matrices_;
  std::vector< std::map< size_t, std::shared_ptr< AffinelyDecomposedMatrixType > > > inside_outside_coupling_matrices_;
  std::vector< std::map< size_t, std::shared_ptr< AffinelyDecomposedMatrixType > > > outside_inside_coupling_matrices_;
}; // BlockSWIPDG

//  void visualize_local(const DUNE_STUFF_SSIZE_T ss,
//                       const VectorType& vector,
//                       const std::string filename,
//                       const std::string name) const
//  {
//    assert(ss < (DUNE_STUFF_SSIZE_T)(ms_grid_->size()));
//    local_discretizations_[ss]->visualize(vector, filename, name);
//  } // ... visualize_local(...)

//  DUNE_STUFF_SSIZE_T num_subdomains() const
//  {
//    return ms_grid_->size();
//  }

//  VectorType localize_vector(const VectorType& global_vector, const DUNE_STUFF_SSIZE_T ss) const
//  {
//    if (ss < 0 || ss >= num_subdomains())
//      DUNE_PYMOR_THROW(Pymor::Exception::index_out_of_range,
//                       "0 <= ss < num_subdomains() = " << num_subdomains() << " is not true for ss = " << ss << "!");
//    if (global_vector.size() != ansatz_mapper_->size())
//      DUNE_PYMOR_THROW(Pymor::Exception::sizes_do_not_match,
//                       "The size() of global_vector (" << global_vector.dim()
//                       << ") does not match the size() of the ansatz space (" << ansatz_mapper_->size() << ")!");
//    assert(ss < (DUNE_STUFF_SSIZE_T)(local_discretizations_.size()));
//    VectorType local_vector = local_discretizations_[ss]->create_vector();
//    for (size_t ii = 0; ii < local_vector.size(); ++ii)
//      local_vector.set_entry(ii, global_vector.get_entry(ansatz_mapper_->mapToGlobal(ss, ii)));
//    return local_vector;
//  } // ... localize_vetor(...)

//  VectorType* localize_vector_and_return_ptr(const VectorType& global_vector, const DUNE_STUFF_SSIZE_T ss) const
//  {
//    return new VectorType(localize_vector(global_vector, ss));
//  }

//  std::vector< DUNE_STUFF_SSIZE_T > neighbouring_subdomains(const DUNE_STUFF_SSIZE_T ss) const
//  {
//    if (ss < 0 || ss >= num_subdomains())
//      DUNE_PYMOR_THROW(Pymor::Exception::index_out_of_range,
//                       "0 <= ss < num_subdomains() = " << num_subdomains() << " is not true for ss = " << ss << "!");
//    const auto set_of_neighbours = ms_grid_->neighborsOf(ss);
//    return std::vector< DUNE_STUFF_SSIZE_T >(set_of_neighbours.begin(), set_of_neighbours.end());
//  }

//public:
//  void visualize(const VectorType& vector, const std::string filename, const std::string name) const
//  {
//    if (vector.size() != ansatz_mapper_->size())
//      DUNE_PYMOR_THROW(Pymor::Exception::sizes_do_not_match,
//                       "the dim of the vector (" << vector.dim() << ") does not match the size of the ansatz space ("
//                       << ansatz_mapper_->size() << ")!");

////    std::vector< std::unique_ptr< VectorType > > local_vectors(ms_grid_->size());
////    for (size_t ss = 0; ss < ms_grid_->size(); ++ss) {
////      VectorType* local_vector_ptr = new VectorType(localize_vector(vector, ss));
////      local_vectors[ss] = std::unique_ptr< VectorType >(local_vector_ptr);
////    }

//    const MultiscaleFunction function(*this, local_discretizations_, vector, name);
//    function.visualize(filename);
////    VTKWriter< GlobalGridViewType > vtk_writer(*(ms_grid_->globalGridView()), VTK::nonconforming);
////    vtk_writer.addVertexData(function);
////    vtk_writer.write(filename);
//  } // ... visualize(...)

//  void initialize(std::ostream& out = Dune::Stuff::Common::Logger().devnull(),
//                  const std::string prefix = "")
//  {
//    if (!initialized_) {
//      const size_t subdomains = ms_grid_->size();
//      // initialize the global pattern
//      global_pattern_ = std::make_shared< PatternType >(test_mapper_->size());
//      // walk the subdomains for the first time
//      //   * to initialize the coupling pattern,
//      //   * to finalize the global sparsity pattern
//      out << prefix << "walking subdomains for the first time... " << std::flush;
//      Dune::Timer timer;
//      // this local problem and boundary info is needed for the local products
//      const auto local_product_boundary_info
//          = std::make_shared< Stuff::GridboundaryAllDirichlet< typename LocalGridPartType::IntersectionType > >();
//      typedef Problem::Default< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, scalarDiffusion >
//          LocalProblemType;
//      typedef Stuff::Function::Constant< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >
//          ConstantFunctionType;
//      typedef Pymor::Function::NonparametricDefault< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >
//          WrappedFunctionType;
//      const auto local_product_problem
//          = std::make_shared< LocalProblemType >(std::make_shared< WrappedFunctionType >(new ConstantFunctionType(RangeFieldType(0))),
//                                                 std::make_shared< WrappedFunctionType >(new ConstantFunctionType(RangeFieldType(0))),
//                                                 std::make_shared< WrappedFunctionType >(new ConstantFunctionType(RangeFieldType(0))),
//                                                 std::make_shared< WrappedFunctionType >(new ConstantFunctionType(RangeFieldType(0))));
//      for (size_t ss = 0; ss < subdomains; ++ss) {
//        // init the local discretizations (assembles matrices and patterns)
//        local_discretizations_[ss]->initialize();
//        oversampled_discretizations_[ss]->initialize();
//        // assemble the local products
//        // * therefore create nonparametric local discretization
//        LocalDiscretizationType local_product_discretization(ms_grid_->localGridPart(ss),
//                                                             local_product_boundary_info,
//                                                             local_product_problem);
//        local_product_discretization.initialize();
//        // and get all of the products
//        for (auto id : local_product_discretization.available_products())
//          local_products_[ss].insert(std::make_pair(id, local_product_discretization.get_product(id)));
//        // and create the local containers
//        // * the matrices
//        //   * just copy thos from the local discretizations
//        const auto local_operator = local_discretizations_[ss]->get_operator();
//        local_matrices_[ss] = std::make_shared< AffinelyDecomposedMatrixType >();
//        //   * we take the affine part only if the diffusion has one, otherwise it contains only the dirichlet rows,
//        //     thus it is empty, since the local problems are purely neumann
//        if (problem_.diffusion()->has_affine_part()) {
//          if (!local_operator.has_affine_part())
//            DUNE_PYMOR_THROW(Pymor::Exception::requirements_not_met, "The local operator is missing the affine part!");
//          local_matrices_[ss]->register_affine_part(new MatrixType(
//                                                      local_operator.affine_part().container()->backend()));
//        }
//        if (local_operator.num_components() < problem_.diffusion()->num_components())
//          DUNE_PYMOR_THROW(Pymor::Exception::requirements_not_met,
//                           "The local operator should have " << problem_.diffusion()->num_components()
//                           << " components (but has only " << local_operator.num_components() << ")!");
//        for (DUNE_STUFF_SSIZE_T qq = 0; qq < problem_.diffusion()->num_components(); ++qq)
//          local_matrices_[ss]->register_component(new MatrixType(
//              local_operator.component(qq).container()->backend()),
//              problem_.diffusion()->coefficient(qq));
//        // * and the vectors
//        const auto local_functional = local_discretizations_[ss]->get_rhs();
//        local_vectors_[ss] = std::make_shared< AffinelyDecomposedVectorType >();
//        //   * first the affine part
//        if (problem_.force()->has_affine_part() || problem_.neumann()->has_affine_part()) {
//          //   * which we copy from the local discretization
//          if (!local_functional.has_affine_part())
//            DUNE_PYMOR_THROW(Pymor::Exception::requirements_not_met,
//                             "The local functional is missing the affine part!");
//          local_vectors_[ss]->register_affine_part(new VectorType(
//              local_functional.affine_part().container()->backend()));
//        } else if (problem_.diffusion()->has_affine_part()
//                   || problem_.dirichlet()->has_affine_part()) {
//          //   * but not, if it was due to the dirichlet boundary correction
//          //     In this case we just create an empty one, since the one of the local discretization has to be empty
//          local_vectors_[ss]->register_affine_part(new VectorType(test_mapper_->size()));
//        }
//        //   * then we copy the components from the local discretizations
//        for (DUNE_STUFF_SSIZE_T qq = 0; qq < problem_.force()->num_components(); ++qq)
//          local_vectors_[ss]->register_component(new VectorType(
//              local_functional.component(qq).container()->backend()),
//              problem_.force()->coefficient(qq));
//        for (DUNE_STUFF_SSIZE_T qq = 0; qq < problem_.neumann()->num_components(); ++qq)
//          local_vectors_[ss]->register_component(new VectorType(
//              local_functional.component(problem_.force()->num_components() + qq).container()->backend()),
//              problem_.neumann()->coefficient(qq));
//        //   * and create the components due to the dirichlet boundary term
//        if (problem_.diffusion()->has_affine_part())
//          for (DUNE_STUFF_SSIZE_T qq = 0; qq < problem_.dirichlet()->num_components(); ++qq)
//            local_vectors_[ss]->register_component(new VectorType(local_discretizations_[ss]->testSpace().mapper().size()),
//                                                   problem_.dirichlet()->coefficient(qq));
//        if (problem_.dirichlet()->has_affine_part())
//          for (DUNE_STUFF_SSIZE_T qq = 0; qq < problem_.diffusion()->num_components(); ++qq)
//            local_vectors_[ss]->register_component(new VectorType(local_discretizations_[ss]->testSpace().mapper().size()),
//                                                   problem_.diffusion()->coefficient(qq));
//        if (problem_.diffusion()->num_components() > 0 && problem_.dirichlet()->num_components() > 0) {
//          Pymor::ParameterType diffusion_dirichlet_mu;
//          for (auto key : problem_.diffusion()->parameter_type().keys())
//            diffusion_dirichlet_mu.set(key, problem_.diffusion()->parameter_type().get(key));
//          for (auto key : problem_.dirichlet()->parameter_type().keys())
//            diffusion_dirichlet_mu.set(key, problem_.dirichlet()->parameter_type().get(key));
//          for (DUNE_STUFF_SSIZE_T pp = 0; pp < problem_.diffusion()->num_components(); ++ pp) {
//            for (DUNE_STUFF_SSIZE_T qq = 0; qq < problem_.dirichlet()->num_components(); ++qq) {
//              const std::string expression = "(" + problem_.diffusion()->coefficient(pp)->expression()
//                                             + ")*(" + problem_.dirichlet()->coefficient(qq)->expression() + ")";
//              local_vectors_[ss]->register_component(new VectorType(local_discretizations_[ss]->testSpace().mapper().size()),
//                                                     new Pymor::ParameterFunctional(diffusion_dirichlet_mu,
//                                                                                    expression));
//            }
//          }
//        } // create the local containers

//        // create and copy the local patterns
//        add_local_to_global_pattern(*(local_discretizations_[ss]->pattern()), ss, ss, *global_pattern_);
//        const auto& inner_test_space = local_discretizations_[ss]->testSpace();
//        const auto& inner_ansatz_space = local_discretizations_[ss]->ansatzSpace();
//        // walk the neighbors
//        for (const size_t& nn : ms_grid_->neighborsOf(ss)) {
//          // visit each coupling only once (assemble primally)
//          if (ss < nn) {
//            const auto& outer_test_space = local_discretizations_[nn]->testSpace();
//            const auto& outer_ansatz_space = local_discretizations_[nn]->ansatzSpace();
//            const auto& inside_outside_grid_part = *(ms_grid_->couplingGridPart(ss, nn));
//            const auto& outside_inside_grid_part = *(ms_grid_->couplingGridPart(nn, ss));
//            // create the coupling patterns
//            std::shared_ptr< PatternType > inside_outside_pattern(
//                  inner_test_space.computeCodim1Pattern(inside_outside_grid_part, outer_ansatz_space));
//            inside_outside_coupling_patterns_[ss].insert(std::make_pair(nn, inside_outside_pattern));
//            std::shared_ptr< PatternType > outside_inside_pattern(
//                  outer_test_space.computeCodim1Pattern(outside_inside_grid_part, inner_ansatz_space));
//            outside_inside_coupling_patterns_[nn].insert(std::make_pair(ss, outside_inside_pattern));
//            // and copy them
//            add_local_to_global_pattern(*inside_outside_pattern,  ss, nn, *global_pattern_);
//            add_local_to_global_pattern(*outside_inside_pattern,  nn, ss, *global_pattern_);
//          } // visit each coupling only once (assemble primaly)
//        } // walk the neighbors
//      } // walk the subdomains for the first time
//      out<< "done (took " << timer.elapsed() << " sek)" << std::endl;

//      // walk the subdomains for the second time
//      //   * to assemble the boundary matrices and vectors and
//      //   * to assemble the coupling matrices
//      out << prefix << "walking subdomains for the second time... " << std::flush;
//      for (size_t ss = 0; ss < ms_grid_->size(); ++ss) {
//        const auto& inner_test_mapper = local_discretizations_[ss]->testSpace().mapper();
//        const auto& inner_ansatz_mapper = local_discretizations_[ss]->ansatzSpace().mapper();
//        if (ms_grid_->boundary(ss))
//          assemble_boundary_contributions(ss);
//        // walk the neighbors
//        for (const size_t& nn : ms_grid_->neighborsOf(ss)) {
//          // visit each coupling only once (assemble primaly)
//          if (ss < nn) {
//            const auto& outer_test_mapper = local_discretizations_[nn]->testSpace().mapper();
//            const auto& outer_ansatz_mapper = local_discretizations_[nn]->ansatzSpace().mapper();
//            // get the patterns
//            const auto in_out_result = inside_outside_coupling_patterns_[ss].find(nn);
//            if (in_out_result == inside_outside_coupling_patterns_[ss].end())
//              DUNE_PYMOR_THROW(Pymor::Exception::this_does_not_make_any_sense,
//                               "subdomain " << ss << ", neighbour " << nn);
//            const auto& inside_outside_pattern = *(in_out_result->second);
//            const auto out_in_result = outside_inside_coupling_patterns_[nn].find(ss);
//            if (out_in_result == outside_inside_coupling_patterns_[nn].end())
//              DUNE_PYMOR_THROW(Pymor::Exception::this_does_not_make_any_sense,
//                               "subdomain " << ss << ", neighbour " << nn);
//            const auto& outside_inside_pattern = *(out_in_result->second);
//            // create the coupling matrices
//            auto inside_outside_matrix = std::make_shared< AffinelyDecomposedMatrixType >();
//            auto outside_inside_matrix = std::make_shared< AffinelyDecomposedMatrixType >();
//            if (problem_.diffusion()->has_affine_part()) {
//              inside_outside_matrix->register_affine_part(new MatrixType(inner_test_mapper.size(),
//                                                                         outer_ansatz_mapper.size(),
//                                                                         inside_outside_pattern));
//              outside_inside_matrix->register_affine_part(new MatrixType(outer_test_mapper.size(),
//                                                                         inner_ansatz_mapper.size(),
//                                                                         outside_inside_pattern));
//            }
//            for (DUNE_STUFF_SSIZE_T qq = 0; qq < problem_.diffusion()->num_components(); ++qq) {
//              inside_outside_matrix->register_component(new MatrixType(inner_test_mapper.size(),
//                                                                       outer_ansatz_mapper.size(),
//                                                                       inside_outside_pattern),
//                                                        problem_.diffusion()->coefficient(qq));
//              outside_inside_matrix->register_component(new MatrixType(outer_test_mapper.size(),
//                                                                       inner_ansatz_mapper.size(),
//                                                                       outside_inside_pattern),
//                                                        problem_.diffusion()->coefficient(qq));
//            }
//            // and assemble them
//            assemble_coupling_contributions(ss, nn,
//                                            *(local_matrices_[ss]),
//                                            *(inside_outside_matrix),
//                                            *(outside_inside_matrix),
//                                            *(local_matrices_[nn]));
//            inside_outside_coupling_matrices_[ss].insert(std::make_pair(nn, inside_outside_matrix));
//            outside_inside_coupling_matrices_[nn].insert(std::make_pair(ss, outside_inside_matrix));
//          } // visit each coupling only once
//        } // walk the neighbors
//      } // walk the subdomains for the second time
//      out<< "done (took " << timer.elapsed() << " sek)" << std::endl;
//      // parameter
//      this->inherit_parameter_type(*(local_matrices_[0]), "lhs");
//      this->inherit_parameter_type(*(local_vectors_[0]), "rhs");
//      // done
//      initialized_ = true;
//    } // if !(initialized_)
//  } // void initialize(...)

//  OperatorType get_operator() const
//  {
//    // initialize global matrix
//    AffinelyDecomposedMatrixType system_matrix;
//    if (local_matrices_[0]->has_affine_part())
//      system_matrix.register_affine_part(new MatrixType(test_mapper_->size(),
//                                                        ansatz_mapper_->size(),
//                                                        *global_pattern_));
//    for (DUNE_STUFF_SSIZE_T qq = 0; qq < local_matrices_[0]->num_components(); ++qq)
//      system_matrix.register_component(new MatrixType(test_mapper_->size(),
//                                                      ansatz_mapper_->size(),
//                                                      *global_pattern_),
//                                       local_matrices_[0]->coefficient(qq));
//    // walk the subdomains
//    for (size_t ss = 0; ss < ms_grid_->size(); ++ss) {
//      auto local_matrix = *(local_matrices_[ss]);
//      if (local_matrix.has_affine_part())
//        copy_local_to_global_matrix(local_matrix.affine_part()->backend(),
//                                    *(local_discretizations_[ss]->pattern()),
//                                    ss, ss,
//                                    system_matrix.affine_part()->backend());
//      for (DUNE_STUFF_SSIZE_T qq = 0; qq < local_matrix.num_components(); ++qq) {
//        copy_local_to_global_matrix(local_matrix.component(qq)->backend(),
//                                    *(local_discretizations_[ss]->pattern()),
//                                    ss, ss,
//                                    system_matrix.component(qq)->backend());
//      }

//      // walk the neighbours
//      for (const size_t& nn : ms_grid_->neighborsOf(ss)) {
//        if (ss < nn) {
//          // get the coupling patterns
//          const auto result_inside_outside_pattern = inside_outside_coupling_patterns_[ss].find(nn);
//          if (result_inside_outside_pattern == inside_outside_coupling_patterns_[ss].end())
//            DUNE_PYMOR_THROW(Pymor::Exception::this_does_not_make_any_sense,
//                             "The coupling pattern for subdomain " << ss << " and neighbour " << nn << "is missing!");
//          const auto& inside_outside_pattern = *(result_inside_outside_pattern->second);
//          const auto result_outside_inside_pattern = outside_inside_coupling_patterns_[nn].find(ss);
//          if (result_outside_inside_pattern == outside_inside_coupling_patterns_[nn].end())
//            DUNE_PYMOR_THROW(Pymor::Exception::this_does_not_make_any_sense,
//                             "The coupling pattern for neighbour " << nn << " and subdomain " << ss << "is missing!");
//          const auto& outside_inside_pattern = *(result_outside_inside_pattern->second);
//          // and the coupling matrices
//          auto result_inside_outside_matrix = inside_outside_coupling_matrices_[ss].find(nn);
//          if (result_inside_outside_matrix == inside_outside_coupling_matrices_[ss].end())
//            DUNE_PYMOR_THROW(Pymor::Exception::this_does_not_make_any_sense,
//                             "The coupling matrix for subdomain " << ss << " and neighbour " << nn << "is missing!");
//          auto& inside_outside_matrix = *(result_inside_outside_matrix->second);
//          auto result_outside_inside_matrix = outside_inside_coupling_matrices_[nn].find(ss);
//          if (result_outside_inside_matrix == outside_inside_coupling_matrices_[nn].end())
//            DUNE_PYMOR_THROW(Pymor::Exception::this_does_not_make_any_sense,
//                             "The coupling matrix for neighbour " << nn << " and subdomain " << ss << "is missing!");
//          auto& outside_inside_matrix = *(result_outside_inside_matrix->second);
//          // and copy them into the global matrix
//          if (inside_outside_matrix.has_affine_part())
//            copy_local_to_global_matrix(inside_outside_matrix.affine_part()->backend(),
//                                        inside_outside_pattern,
//                                        ss, nn,
//                                        system_matrix.affine_part()->backend());
//          for (DUNE_STUFF_SSIZE_T qq = 0; qq < inside_outside_matrix.num_components(); ++qq) {
//            copy_local_to_global_matrix(inside_outside_matrix.component(qq)->backend(),
//                                        inside_outside_pattern,
//                                        ss, nn,
//                                        system_matrix.component(qq)->backend());
//          }
//          if (outside_inside_matrix.has_affine_part())
//            copy_local_to_global_matrix(outside_inside_matrix.affine_part()->backend(),
//                                        outside_inside_pattern,
//                                        nn, ss,
//                                        system_matrix.affine_part()->backend());
//          for (DUNE_STUFF_SSIZE_T qq = 0; qq < outside_inside_matrix.num_components(); ++qq) {
//            copy_local_to_global_matrix(outside_inside_matrix.component(qq)->backend(),
//                                        outside_inside_pattern,
//                                        nn, ss,
//                                        system_matrix.component(qq)->backend());
//          }
//        }
//      } // walk the neighbours
//    }
//    return OperatorType(system_matrix);
//  }

//  OperatorType* get_operator_and_return_ptr() const
//  {
//    return new OperatorType(get_operator());
//  }

//  FunctionalType get_rhs() const
//  {
//    // initialize global vector
//    AffinelyDecomposedVectorType rhs_vector;
//    if (local_vectors_[0]->has_affine_part())
//      rhs_vector.register_affine_part(new VectorType(test_mapper_->size()));
//    for (DUNE_STUFF_SSIZE_T qq = 0; qq < local_vectors_[0]->num_components(); ++qq)
//      rhs_vector.register_component(new VectorType(test_mapper_->size()),
//                                    local_vectors_[0]->coefficient(qq));
//    // walk the subdomains
//    for (size_t ss = 0; ss < ms_grid_->size(); ++ss) {
//      auto local_vector = *(local_vectors_[ss]);
//      if (local_vector.has_affine_part())
//        copy_local_to_global_vector(local_vector.affine_part()->backend(),
//                                    ss,
//                                    rhs_vector.affine_part()->backend());
//      for (DUNE_STUFF_SSIZE_T qq = 0; qq < local_vector.num_components(); ++qq)
//        copy_local_to_global_vector(local_vector.component(qq)->backend(),
//                                    ss,
//                                    rhs_vector.component(qq)->backend());

//    }
//    return FunctionalType(rhs_vector);
//  }

//  FunctionalType* get_rhs_and_return_ptr() const
//  {
//    return new FunctionalType(get_rhs());
//  }

//  std::vector< std::string > available_products() const
//  {
//    return {"l2", "h1"};
//  }

//  ProductType get_product(const std::string id) const
//  {
//    if (id == "l2") {
//      auto product_matrix = new MatrixType(test_mapper_->size(),
//                                                  ansatz_mapper_->size(),
//                                                  *global_pattern_);
//      for (size_t ss = 0; ss < ms_grid_->size(); ++ss)
//        copy_local_to_global_matrix(local_discretizations_[ss]->get_product("l2").container()->backend(),
//                                    *(local_discretizations_[ss]->pattern()),
//                                    ss, ss,
//                                    *product_matrix);
//      return ProductType(product_matrix);
//    } else if (id == "h1") {
//      auto product_matrix = new MatrixType(test_mapper_->size(),
//                                                  ansatz_mapper_->size(),
//                                                  *global_pattern_);
//      for (size_t ss = 0; ss < ms_grid_->size(); ++ss)
//        copy_local_to_global_matrix(local_discretizations_[ss]->get_product("h1").container()->backend(),
//                                    *(local_discretizations_[ss]->pattern()),
//                                    ss, ss,
//                                    *product_matrix);
//      return ProductType(product_matrix);
//    } else
//      DUNE_PYMOR_THROW(Pymor::Exception::this_does_not_make_any_sense,
//                       "Product id has to be one of 'l2' or 'h1' (is " << id << ")!");
//  }

//  ProductType* get_product_and_return_ptr(const std::string id) const
//  {
//    return new ProductType(get_product(id));
//  }

//  ProductType get_local_product(const DUNE_STUFF_SSIZE_T ss, const std::string id) const
//  {
//    if (ss < 0 || ss >= num_subdomains())
//      DUNE_PYMOR_THROW(Pymor::Exception::index_out_of_range,
//                       "0 <= ss < num_subdomains() = " << num_subdomains() << " is not true for ss = " << ss << "!");
//    const auto result = local_products_[ss].find(id);
//    if (result == local_products_[ss].end())
//      DUNE_PYMOR_THROW(Pymor::Exception::key_is_not_valid, "There is no local product with id '" << id << "'!");
//    return result->second;
//  }

//  ProductType* get_local_product_and_return_ptr(const DUNE_STUFF_SSIZE_T ss, const std::string id) const
//  {
//    return new ProductType(get_local_product(ss, id));
//  }

//  OperatorType get_local_operator(const DUNE_STUFF_SSIZE_T ss) const
//  {
//    if (ss < 0 || ss >= num_subdomains())
//      DUNE_PYMOR_THROW(Pymor::Exception::index_out_of_range,
//                       "0 <= ss < num_subdomains() = " << num_subdomains() << " is not true for ss = " << ss << "!");
//    assert(ss < int(local_matrices_.size()));
//    return OperatorType(*(local_matrices_[ss]));
//  }

//  OperatorType* get_local_operator_and_return_ptr(const DUNE_STUFF_SSIZE_T ss) const
//  {
//    return new OperatorType(get_local_operator(ss));
//  }

//  OperatorType get_coupling_operator(const DUNE_STUFF_SSIZE_T ss, const DUNE_STUFF_SSIZE_T nn) const
//  {
//    if (ss < 0 || ss >= num_subdomains())
//      DUNE_PYMOR_THROW(Pymor::Exception::index_out_of_range,
//                       "0 <= ss < num_subdomains() = " << num_subdomains() << " is not true for ss = " << ss << "!");
//    const auto neighbours = ms_grid_->neighborsOf(ss);
//    if (neighbours.count(nn) == 0)
//      DUNE_PYMOR_THROW(Pymor::Exception::index_out_of_range,
//                       "Subdomain " << nn << " is not a neighbour of subdomain " << ss
//                       << " (call neighbouring_subdomains(" << ss << ") to find out)!");
//    if (ss < nn) {
//      // we need to look for this coupling operator in the inside/outside context
//      const auto result_inside_outside_matrix = inside_outside_coupling_matrices_[ss].find(nn);
//      if (result_inside_outside_matrix == inside_outside_coupling_matrices_[ss].end())
//        DUNE_PYMOR_THROW(Pymor::Exception::this_does_not_make_any_sense,
//                         "The coupling matrix for subdomain " << ss << " and neighbour " << nn << " is missing!");
//      const auto inside_outside_matrix = result_inside_outside_matrix->second;
//      return OperatorType(*inside_outside_matrix);
//    } else if (nn < ss) {
//      // we need to look for this coupling operator in the outside/inside context
//      const auto result_outside_inside_matrix = outside_inside_coupling_matrices_[ss].find(nn);
//      if (result_outside_inside_matrix == outside_inside_coupling_matrices_[ss].end())
//        DUNE_PYMOR_THROW(Pymor::Exception::this_does_not_make_any_sense,
//                         "The coupling matrix for neighbour " << nn << " and subdomain " << ss << " is missing!");
//      const auto outside_inside_matrix = result_outside_inside_matrix->second;
//      return OperatorType(*outside_inside_matrix);
//    } else {
//      // the above exception should have cought this
//      DUNE_PYMOR_THROW(Pymor::Exception::this_does_not_make_any_sense,
//                       "The multiscale grid is corrupted! Subdomain " << ss << " must not be its own neighbour!");
//    }
//  } // ... get_coupling_operator(...)

//  OperatorType* get_coupling_operator_and_return_ptr(const DUNE_STUFF_SSIZE_T ss, const DUNE_STUFF_SSIZE_T nn) const
//  {
//    return new OperatorType(get_coupling_operator(ss, nn));
//  }

//  FunctionalType get_local_functional(const DUNE_STUFF_SSIZE_T ss) const
//  {
//    if (ss < 0 || ss >= num_subdomains())
//      DUNE_PYMOR_THROW(Pymor::Exception::index_out_of_range,
//                       "0 <= ss < num_subdomains() = " << num_subdomains() << " is not true for ss = " << ss << "!");
//    assert(ss < (DUNE_STUFF_SSIZE_T)(local_vectors_.size()));
//    return FunctionalType(*(local_vectors_[ss]));
//  }

//  FunctionalType* get_local_functional_and_return_ptr(const DUNE_STUFF_SSIZE_T ss) const
//  {
//    return new FunctionalType(get_local_functional(ss));
//  }

//  void uncached_solve(VectorType& vector,
//                      const Pymor::Parameter mu = Pymor::Parameter(),
//                      std::ostream& out = Dune::Stuff::Common::Logger().devnull(),
//                      const std::string prefix = "") const
//  {
//    if (!initialized_)
//      DUNE_PYMOR_THROW(Pymor::Exception::requirements_not_met,
//                       "Do not call solve() if initialized() == false (call initialize() first)!");
//    // check input
//    if (mu.type() != Pymor::Parametric::parameter_type())
//      DUNE_PYMOR_THROW(Pymor::Exception::wrong_parameter_type,
//                       "the type of mu (" << mu.type() << ") does not match the parameter_type of this ("
//                       << Pymor::Parametric::parameter_type() << ")!");
//    if (vector.size() != ansatz_mapper_->size())
//      DUNE_PYMOR_THROW(Pymor::Exception::sizes_do_not_match,
//                       "the dim of the vector (" << vector.dim() << ") does not match the size of the ansatz space ("
//                       << ansatz_mapper_->size() << ")!");
//    out << prefix << "computing global lhs and rhs... " << std::flush;
//    Dune::Timer timer;
//    auto system_matrix = std::shared_ptr< MatrixType >(new MatrixType(test_mapper_->size(),
//                                                                      ansatz_mapper_->size(),
//                                                                      *global_pattern_));
//    VectorType rhs_vector(test_mapper_->size());
//    Pymor::Parameter mu_rhs;
//    if (local_vectors_[0]->parametric())
//      mu_rhs = Pymor::Parametric::map_parameter(mu, "rhs");
//    Pymor::Parameter mu_lhs;
//    if (local_matrices_[0]->parametric())
//      mu_lhs = Pymor::Parametric::map_parameter(mu, "lhs");
//    // walk the subdomains
//    for (size_t ss = 0; ss < ms_grid_->size(); ++ss) {
//      // get the local container and patterns
//      const auto& local_pattern = *(local_discretizations_[ss]->pattern());
//      const auto& local_matrix = *(local_matrices_[ss]);
//      const auto& local_vector = *(local_vectors_[ss]);
//      // and freeze them
//      std::shared_ptr< const VectorType > local_frozen_vector;
//      if (local_vector.parametric())
//        local_frozen_vector = std::make_shared< const VectorType >(local_vector.freeze_parameter(mu_rhs));
//      else
//        local_frozen_vector = local_vector.affine_part();
//      std::shared_ptr< const MatrixType > local_frozen_matrix;
//      if (local_matrix.parametric())
//        local_frozen_matrix = std::make_shared< const MatrixType >(local_matrix.freeze_parameter(mu_lhs));
//      else
//        local_frozen_matrix = local_matrix.affine_part();
//      // and copy them into the global containers
//      copy_local_to_global_vector(*local_frozen_vector, ss, rhs_vector);
//      copy_local_to_global_matrix(*local_frozen_matrix, local_pattern, ss, ss, *system_matrix);

//      // walk the neighbours
//      for (const size_t& nn : ms_grid_->neighborsOf(ss)) {
//        if (ss < nn) {
//          // get the coupling matrices and patterns
//          const auto result_inside_outside_pattern = inside_outside_coupling_patterns_[ss].find(nn);
//          if (result_inside_outside_pattern == inside_outside_coupling_patterns_[ss].end())
//            DUNE_PYMOR_THROW(Pymor::Exception::this_does_not_make_any_sense,
//                             "The coupling pattern for subdomain " << ss << " and neighbour " << nn << "is missing!");
//          const auto& inside_outside_pattern = *(result_inside_outside_pattern->second);
//          const auto result_outside_inside_pattern = outside_inside_coupling_patterns_[nn].find(ss);
//          if (result_outside_inside_pattern == outside_inside_coupling_patterns_[nn].end())
//            DUNE_PYMOR_THROW(Pymor::Exception::this_does_not_make_any_sense,
//                             "The coupling pattern for neighbour " << nn << " and subdomain " << ss << "is missing!");
//          const auto& outside_inside_pattern = *(result_outside_inside_pattern->second);

//          const auto result_inside_outside_matrix = inside_outside_coupling_matrices_[ss].find(nn);
//          if (result_inside_outside_matrix == inside_outside_coupling_matrices_[ss].end())
//            DUNE_PYMOR_THROW(Pymor::Exception::this_does_not_make_any_sense,
//                             "The coupling matrix for subdomain " << ss << " and neighbour " << nn << "is missing!");
//          const auto& inside_outside_matrix = *(result_inside_outside_matrix->second);
//          const auto result_outside_inside_matrix = outside_inside_coupling_matrices_[nn].find(ss);
//          if (result_outside_inside_matrix == outside_inside_coupling_matrices_[nn].end())
//            DUNE_PYMOR_THROW(Pymor::Exception::this_does_not_make_any_sense,
//                             "The coupling matrix for neighbour " << nn << " and subdomain " << ss << "is missing!");
//          const auto& outside_inside_matrix = *(result_outside_inside_matrix->second);
//          // and freeze them
//          std::shared_ptr< MatrixType > inside_outside_frozen_matrix;
//          if (inside_outside_matrix.parametric())
//            inside_outside_frozen_matrix = std::shared_ptr< MatrixType >(new MatrixType(inside_outside_matrix.freeze_parameter(mu_lhs)));
//          else
//            inside_outside_frozen_matrix = inside_outside_matrix.affine_part();
//          std::shared_ptr< MatrixType > outside_inside_frozen_matrix;
//          if (outside_inside_matrix.parametric())
//            outside_inside_frozen_matrix = std::shared_ptr< MatrixType >(new MatrixType(outside_inside_matrix.freeze_parameter(mu_lhs)));
//          else
//            outside_inside_frozen_matrix = outside_inside_matrix.affine_part();

//          // and copy them into the global matrix
//          copy_local_to_global_matrix(*inside_outside_frozen_matrix,
//                                      inside_outside_pattern,
//                                      ss, nn,
//                                      *system_matrix);
//          copy_local_to_global_matrix(*outside_inside_frozen_matrix,
//                                      outside_inside_pattern,
//                                      nn, ss,
//                                      *system_matrix);
//        }
//      } // walk the neighbours
//    } // walk the subdomains
//    out << "done (took " << timer.elapsed() << "s)" << std::endl;

//    const typename OperatorType::ComponentType op(system_matrix);
//    const std::string option = op.invert_options()[0];
//    out << prefix << "solving with '" << option << "' option... " << std::flush;
//    timer.reset();
//    op.apply_inverse(rhs_vector, vector, option);
//    out << "done (took " << timer.elapsed() << "s)" << std::endl;
//  } // ... uncached_solve(...)

//  VectorType* solve_and_return_ptr(const Pymor::Parameter mu = Pymor::Parameter()/*,
//                                   std::ostream& out = Dune::Stuff::Common::Logger().devnull(),
//                                   const std::string prefix = ""*/) const
//  {
//    VectorType* ret = create_vector_and_return_ptr();
//    this->solve(*ret, mu/*, out, prefix*/);
//    return ret;
//  }

//  /**
//   *  \note We still have an unnecessarry copy here!
//   */
//  void solve_for_local_correction(const VectorType& global_vector,
//                                  const DUNE_STUFF_SSIZE_T ss,
//                                  VectorType& local_vector,
//                                  const Pymor::Parameter mu = Pymor::Parameter(),
//                                  std::ostream&
//#ifndef NDEBUG
//                                                out
//#endif
//                                                    = Dune::Stuff::Common::Logger().devnull(),
//                                  const std::string
//#ifndef NDEBUG
//                                                    prefix
//#endif
//                                                           = "") const
//  {
//    if (!initialized_)
//      DUNE_PYMOR_THROW(Pymor::Exception::requirements_not_met,
//                       "Do not call uncached_local_solve() if initialized() == false (call initialize() first)!");
//    // check input
//    if (mu.type() != Pymor::Parametric::parameter_type())
//      DUNE_PYMOR_THROW(Pymor::Exception::wrong_parameter_type,
//                       "the type of mu (" << mu.type() << ") does not match the parameter_type of this ("
//                       << Pymor::Parametric::parameter_type() << ")!");
//    // solve locally
//#ifndef NDEBUG
////    out << prefix << "initalizing tmp solution vector on oversampled subdomain... " << std::flush;
//#endif
////    auto local_oversampled_vector = oversampled_discretizations_[ss]->create_vector();
//#ifndef NDEBUG
////    out << "done (took " << timer.elapsed() << "s)" << std::endl;
//    out << prefix << "assembling local operator and functional... " << std::flush;
//    Dune::Timer timer;
//#endif
//    // first, assemble the oversampled nonparametric system matrix and rhs
//    MatrixType oversampled_system_matrix(                 // <- this is still a copy since container() is only const atm
//          oversampled_discretizations_[ss]->get_operator().freeze_parameter(mu).container()->backend());
//    const auto oversampled_frozen_rhs_functional = oversampled_discretizations_[ss]->get_rhs().freeze_parameter(mu);
//    const VectorType& oversampled_rhs_mu = oversampled_frozen_rhs_functional.container()->backend();

//    // then we restrict the given global vector to the oversampled domain
////    std::vector< std::unique_ptr< VectorType > > tmp_local_vectors(ms_grid_->size());
////    for (size_t ii = 0; ii < ms_grid_->size(); ++ii)
////      tmp_local_vectors[ii] = std::unique_ptr< VectorType >(new VectorType(localize_vector(global_vector, ii)));
//    const MultiscaleFunction tmp_global_function(*this, local_discretizations_, global_vector, "global function");
//    auto tmp_oversampled_vector = oversampled_discretizations_[ss]->create_vector();
//    typedef GDT::DiscreteFunction< LocalAnsatzSpaceType, VectorType > LocalFunctionType;
//    LocalFunctionType tmp_oversampled_function(oversampled_discretizations_[ss]->ansatzSpace(), tmp_oversampled_vector);
//    const GDT::ProjectionOperator::Generic< LocalGridPartType >
//        oversampled_restriction_operator(*(oversampled_discretizations_[ss]->gridPart()));
//    oversampled_restriction_operator.apply(tmp_global_function, tmp_oversampled_function);
//    //   and visualize
//    tmp_oversampled_function.visualize("oversampled_mu_1_solution_subdomain_" + DSC::toString(ss));

//    // then we compute the right hand side
//    //   oversampled_rhs = oversampled_rhs_mu - oversampled_system_matrix * tmp_oversampled_vector
//    auto oversampled_rhs = oversampled_discretizations_[ss]->create_vector();
//    oversampled_system_matrix.mv(tmp_oversampled_vector, oversampled_rhs);
//    oversampled_rhs.scal(-1);
//    oversampled_rhs.iadd(oversampled_rhs_mu);
//    // and clear tmp_oversampled_vector
//    tmp_oversampled_vector.scal(0);

//    // and apply the correct boundary values
//    // so if we are on the domain boundary
////    if (ms_grid_->boundary(ss)) {
//      // we apply the dirichlet constraints
//      // * therefore we create a fake local bounaryinfo
////      //   the interior boundary has id 7
////      auto fake_boundary_id_map = std::make_shared< std::map< std::string, std::set< int > > >();
////      fake_boundary_id_map->operator[]("dirichlet") = {1, 2, 3, 4, 5, 6};
////      fake_boundary_id_map->operator[]("neumann") = {7};
//      const Stuff::GridboundaryAllDirichlet/*IdBased*/< typename LocalGridPartType::IntersectionType >
//          fake_boundary_info/*(fake_boundary_id_map)*/;
//      GDT::Constraints::Dirichlet<  typename LocalGridPartType::IntersectionType, RangeFieldType, true >
//          clear_and_set_rows(fake_boundary_info,
//                          oversampled_discretizations_[ss]->testSpace().mapper().maxNumDofs(),
//                          oversampled_discretizations_[ss]->ansatzSpace().mapper().maxNumDofs());
//      // * and a system assembler to apply those
//      GDT::SystemAssembler< LocalTestSpaceType > oversampled_assembler(oversampled_discretizations_[ss]->testSpace());
//      oversampled_assembler.addLocalConstraints(clear_and_set_rows, oversampled_system_matrix);
//      oversampled_assembler.addLocalConstraints(clear_and_set_rows, oversampled_rhs);
//      // * and do the actual work
//      oversampled_assembler.applyConstraints();
////    } else {
////      // we are in the full neumann setting and set the first DoF to zero
////      local_matrix.unit_row(0);
////      local_rhs.set_entry(0, 0);
////    }
//#ifndef NDEBUG
//    out << "done (took " << timer.elapsed() << "s)" << std::endl;
//    out << prefix << "solving " << local_matrix.rows() << "x" << local_matrix.cols() << " problem on ";
//    if (ms_grid_->oversampling())
//      out << "oversampled ";
//    out << "subdomain " << ss << "... " << std::flush;
//#endif
//    typedef Stuff::LA::BicgstabILUTSolver< MatrixType, VectorType > LASolverType;
//    const LASolverType la_solver;
//    auto settings = la_solver.defaultSettings();
//    settings["precision"] = "1e-10";
//    const auto failure = la_solver.apply(oversampled_system_matrix,
//                                         oversampled_rhs,
//                                         tmp_oversampled_vector, settings);
//    if (failure)
//      DUNE_THROW(MathError, "Linear solver failed on subdomain " << ss << "!");
//#ifndef NDEBUG
//    out << "done (took " << timer.elapsed() << "s)" << std::endl;
//#endif
//    tmp_oversampled_function.visualize("oversampled_correction_" + mu.report_for_filename() + "_subdomain_"
//                                       + DSC::toString(ss));

//    // restrict the oversampled solution to the subdomain
//    // * therefore wrap the oversampled solution into a discrete function on the oversampled subdomain
////    const ConstLocalFunctionType local_oversampled_solution(oversampled_discretizations_[ss]->ansatzSpace(),
////                                                            local_oversampled_vector);
//    // * and the target vector into a discrete function on the subdomain
//    LocalFunctionType local_solution(local_discretizations_[ss]->ansatzSpace(), local_vector);
//    // * and apply the restriction
//    const GDT::ProjectionOperator::Generic< LocalGridPartType > restriction_operator(*(ms_grid_->localGridPart(ss)));
//    restriction_operator.apply(tmp_oversampled_function, local_solution);
//    local_solution.visualize("local_correction_" + mu.report_for_filename() + "_subdomain_" + DSC::toString(ss));

////    // in case we are in the full neumann setting, move the solution to zero mean
////    if (!ms_grid_->boundary(ss)) {
////      const RangeFieldType mean = local_vector.backend().mean();
////      if (Stuff::Common::FloatCmp::ne(mean, RangeFieldType(0)))
////        local_vector.backend().rowwise() -= VectorType(1, mean).backend().transpose();
////    }

//    // now to compute errors
//    //   restrict the original global solution to the subdomain and add it
//    local_vector.iadd(localize_vector(global_vector, ss));
//    local_solution.visualize("local_solution_" + mu.report_for_filename() + "_subdomain_" + DSC::toString(ss));

//    // compute error
//    auto tmp_global_vector = create_vector();
//    this->solve(tmp_global_vector, mu);
//    const auto reference_soution_mu = localize_vector(tmp_global_vector, ss);
//    typedef GDT::ConstDiscreteFunction< LocalAnsatzSpaceType, VectorType > ConstLocalFunctionType;
//    const ConstLocalFunctionType reference_function_mu(local_discretizations_[ss]->ansatzSpace(),
//                                                       reference_soution_mu);
//    reference_function_mu.visualize("reference_solution_" + mu.report_for_filename() + "_subdomain_" + DSC::toString(ss));
//    const GDT::ProductOperator::H1Semi< LocalGridPartType > semi_h1_product(*(ms_grid_->localGridPart(ss)));
//    const Stuff::Function::Difference< ConstLocalFunctionType, ConstLocalFunctionType > difference(reference_function_mu,
//                                                                                                   local_solution);
//    const auto semi_h1_error = std::sqrt(semi_h1_product.apply2(difference, difference));
//    const auto semi_h1_norm = std::sqrt(semi_h1_product.apply2(reference_function_mu, reference_function_mu));
//    std::cout << "  relative semi h1 error on subdomain " << ss << ": " << semi_h1_error / semi_h1_norm << std::endl;
//  } // ... solve_for_local_correction(...)

//private:

//  void copy_local_to_global_matrix(const MatrixType& local_matrix,
//                                   const PatternType& local_pattern,
//                                   const size_t test_subdomain,
//                                   const size_t ansatz_subdomain,
//                                   MatrixType& global_matrix) const
//  {
//    for (size_t local_ii = 0; local_ii < local_pattern.size(); ++local_ii) {
//      const size_t global_ii = test_mapper_->mapToGlobal(test_subdomain, local_ii);
//      for (const size_t& local_jj : local_pattern.inner(local_ii)) {
//        const size_t global_jj = ansatz_mapper_->mapToGlobal(ansatz_subdomain, local_jj);
//        global_matrix.add_to_entry(global_ii, global_jj, local_matrix.get_entry(local_ii, local_jj));
//      }
//    }
//  } // ... copy_local_to_global_matrix(...)

//  void copy_local_to_global_vector(const VectorType& local_vector,
//                                   const size_t subdomain,
//                                   VectorType& global_vector) const
//  {
//    for (size_t local_ii = 0; local_ii < local_vector.size(); ++local_ii) {
//      const size_t global_ii = test_mapper_->mapToGlobal(subdomain, local_ii);
//      global_vector.add_to_entry(global_ii, local_vector.get_entry(local_ii));
//    }
//  } // ... copy_local_to_global_vector(...)

//  /**
//   * \note  We take the matrices as input here becaus we would have to look them up in the maps otherwise. Since that
//   *        has already been done above we save a little.
//   */
//  void assemble_coupling_contributions(const size_t subdomain,
//                                       const size_t neighbour,
//                                       AffinelyDecomposedMatrixType& inside_inside_matrix,
//                                       AffinelyDecomposedMatrixType& inside_outside_matrix,
//                                       AffinelyDecomposedMatrixType& outside_inside_matrix,
//                                       AffinelyDecomposedMatrixType& outside_outside_matrix) const
//  {
//    typedef typename LocalDiscretizationType::TestSpaceType   LocalTestSpaceType;
//    typedef typename LocalDiscretizationType::AnsatzSpaceType LocalAnsatzSpaceType;
//    const LocalTestSpaceType&   inner_test_space   = local_discretizations_[subdomain]->testSpace();
//    const LocalAnsatzSpaceType& inner_ansatz_space = local_discretizations_[subdomain]->ansatzSpace();
//    const LocalTestSpaceType&   outer_test_space   = local_discretizations_[neighbour]->testSpace();
//    const LocalAnsatzSpaceType& outer_ansatz_space = local_discretizations_[neighbour]->ansatzSpace();
//    CouplingAssembler coupling_assembler(inner_test_space, inner_ansatz_space,
//                                         outer_test_space, outer_ansatz_space,
//                                         *(ms_grid_->couplingGridPart(subdomain, neighbour)));

//    typedef typename ProblemType::DiffusionType::NonparametricType DiffusionType;
//    typedef GDT::LocalOperator::Codim1CouplingIntegral< GDT::LocalEvaluation::SIPDG::Inner< DiffusionType > >
//        CouplingOperatorType;
//    typedef GDT::LocalAssembler::Codim1CouplingMatrix< CouplingOperatorType > CouplingMatrixAssemblerType;
//    std::vector< CouplingOperatorType* > coupling_operators;
//    std::vector< CouplingMatrixAssemblerType* > coupling_matrix_assemblers;
//    for (DUNE_STUFF_SSIZE_T qq = 0; qq < problem_.diffusion()->num_components(); ++qq) {
//      coupling_operators.push_back(new CouplingOperatorType(*(problem_.diffusion()->component(qq)), beta_));
//      coupling_matrix_assemblers.push_back(new CouplingMatrixAssemblerType(*(coupling_operators[qq])));
//      coupling_assembler.addLocalAssembler(*(coupling_matrix_assemblers[qq]),
//                                           *(inside_inside_matrix.component(qq)),
//                                           *(inside_outside_matrix.component(qq)),
//                                           *(outside_inside_matrix.component(qq)),
//                                           *(outside_outside_matrix.component(qq)));
//    }
//    if (problem_.diffusion()->has_affine_part()) {
//      coupling_operators.push_back(new CouplingOperatorType(*(problem_.diffusion()->affine_part()), beta_));
//      coupling_matrix_assemblers.push_back(new CouplingMatrixAssemblerType(*(
//          coupling_operators[coupling_operators.size() - 1])));
//      coupling_assembler.addLocalAssembler(*(coupling_matrix_assemblers[coupling_matrix_assemblers.size() - 1]),
//                                           *(inside_inside_matrix.affine_part()),
//                                           *(inside_outside_matrix.affine_part()),
//                                           *(outside_inside_matrix.affine_part()),
//                                           *(outside_outside_matrix.affine_part()));
//    }

//    // do the actual work
//    coupling_assembler.assemble();

//    // clean up
//    for (auto& element : coupling_matrix_assemblers)  delete element;
//    for (auto& element : coupling_operators)          delete element;
//  } // ... assemble_coupling_contributions(...)


} // namespace Discretizations
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_BLOCK_SWIPDG_HH
