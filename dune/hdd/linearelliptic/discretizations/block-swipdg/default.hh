// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_BLOCK_SWIPDG_HH
#define DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_BLOCK_SWIPDG_HH

#include <boost/numeric/conversion/cast.hpp>

#include <dune/stuff/common/memory.hh>

#include <dune/grid/multiscale/provider.hh>

#include <dune/gdt/playground/spaces/block.hh>

#include "base.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Discretizations {
namespace BlockSwipdg {


// forward, needed in the Traits
template <class G, class R = double, int r = 1, int p = 1, Stuff::LA::ChooseBackend la = Stuff::LA::default_sparse_backend>
class DefaultMultiscaleGrid;


namespace internal {


template <class G>
struct DefaultMultiscaleGridLocalGridsProvider
{
  typedef grid::Multiscale::ProviderInterface<G> MultiscaleGridProviderType;
  typedef G LocalGridType;

  DefaultMultiscaleGridLocalGridsProvider(MultiscaleGridProviderType& ms_grd)
    : ms_grid_provider_(ms_grd)
  {}

  const std::shared_ptr<const typename MultiscaleGridProviderType::MsGridType>& ms_grid() const
  {
    return ms_grid_provider_.ms_grid();
  }

  size_t num_subdomains() const
  {
    return ms_grid_provider_.ms_grid()->size();
  }

  MultiscaleGridProviderType& local_grid(const size_t /*subdomain*/)
  {
    return ms_grid_provider_;
  }

  MultiscaleGridProviderType& local_grid(const size_t /*subdomain*/) const
  {
    return ms_grid_provider_;
  }

  int local_level(const size_t subdomain) const
  {
    return boost::numeric_cast<int>(subdomain);
  }

  MultiscaleGridProviderType& ms_grid_provider_;
}; // struct DefaultMultiscaleGridLocalGridsProvider


template <class G, class R, int r, int p, Stuff::LA::ChooseBackend la>
class DefaultMultiscaleGridTraits
{
public:
  typedef DefaultMultiscaleGridLocalGridsProvider<G> LocalGridsProviderType;
  static const constexpr Stuff::Grid::ChooseLayer local_layer   = Stuff::Grid::ChooseLayer::local;
  static const constexpr Stuff::Grid::ChooseLayer overlap_layer = Stuff::Grid::ChooseLayer::local_oversampled;
  typedef DefaultMultiscaleGrid<G, R, r, p, la> derived_type;

  template <class T>
  struct BlockSpace
  {
    typedef GDT::Spaces::Block<T> type;
  };
}; // class DefaultMultiscaleGridTraits


} // namespace internal


/**
 * \attention The given problem is replaced by a Problems::ZeroBoundary.
 * \attention The given boundary info config is replaced by a Stuff::Grid::BoundaryInfos::AllDirichlet.
 * \attention The boundary info for the local oversampled discretizations is hardwired to dirichlet zero atm!
 */
template <class G, class R, int r, int p, Stuff::LA::ChooseBackend la>
class DefaultMultiscaleGrid
  : DSC::StorageProvider<internal::DefaultMultiscaleGridLocalGridsProvider<G>>
  , public internal::Base<internal::DefaultMultiscaleGridTraits<G, R, r, p, la>, R, r, p, la>

{
  typedef DSC::StorageProvider<internal::DefaultMultiscaleGridLocalGridsProvider<G>> LocalGrids;
  typedef internal::Base<internal::DefaultMultiscaleGridTraits<G, R, r, p, la>, R, r, p, la> BaseType;
public:
  typedef internal::BaseTraits<internal::DefaultMultiscaleGridTraits<G, R, r, p, la>, R, r, p, la> Traits;
  typedef typename Traits::LocalGridsProviderType LocalGridsProviderType;

  using typename BaseType::PatternType;
  using typename BaseType::ProblemType;

private:
  using typename BaseType::LocalTestSpaceType;
  using typename BaseType::LocalAnsatzSpaceType;
  typedef typename LocalGridsProviderType::MultiscaleGridProviderType::MsGridType::BoundaryGridPartType BoundaryGridPartType;

public:
  static std::string static_id()
  {
    return BaseType::static_id() + ".default-multiscalegrid";
  }

  DefaultMultiscaleGrid(grid::Multiscale::ProviderInterface<G>& ms_grid_provider_provider,
                        const Stuff::Common::Configuration& bound_inf_cfg,
                        const ProblemType& prob,
                        const std::vector< std::string >& only_these_products = {})
    : LocalGrids(new LocalGridsProviderType(ms_grid_provider_provider))
    , BaseType(LocalGrids::access(), bound_inf_cfg, prob, only_these_products)
  {}

  void init(const bool prune = false)
  {
    if (this->container_based_initialized_)
      return;

    auto logger = Stuff::Common::TimedLogger().get(static_id());
    const size_t num_subdomains = local_grids_provider_.num_subdomains();

    logger.info() << "discretizing on " << num_subdomains << " subdomains..." << std::endl;
    // has to be done first, initializes the local matrices and patterns required for the coupling
    for (size_t subdomain = 0; subdomain < num_subdomains; ++subdomain)
      this->assemble_local_contributions(subdomain);

    for (size_t subdomain = 0; subdomain < num_subdomains; ++subdomain) {
      if (local_grids_provider_.ms_grid()->boundary(subdomain)) {
        GDT::SystemAssembler<LocalTestSpaceType, BoundaryGridPartType, LocalAnsatzSpaceType>
            boundary_assembler(this->local_discretizations_[subdomain]->test_space(),
                               this->local_discretizations_[subdomain]->ansatz_space(),
                               local_grids_provider_.ms_grid()->boundaryGridPart(subdomain));
        this->template assemble_boundary_contributions<BoundaryGridPartType>(subdomain, boundary_assembler);
      }

      for (const size_t& neighboring_subomain : local_grids_provider_.ms_grid()->neighborsOf(subdomain)) {
        if (subdomain < neighboring_subomain) {
          const auto& inner_test_space   = this->local_discretizations_[subdomain]->test_space();
          const auto& inner_ansatz_space = this->local_discretizations_[subdomain]->ansatz_space();
          const auto& outer_test_space   = this->local_discretizations_[neighboring_subomain]->test_space();
          const auto& outer_ansatz_space = this->local_discretizations_[neighboring_subomain]->ansatz_space();
          // compute patterns
          const auto inside_outside_grid_part = local_grids_provider_.ms_grid()->couplingGridPart(subdomain,
                                                                                                  neighboring_subomain);
          const auto outside_inside_grid_part = local_grids_provider_.ms_grid()->couplingGridPart(neighboring_subomain,
                                                                                                  subdomain);
          auto inside_outside_pattern = std::make_shared<PatternType>(
                inner_test_space.compute_face_pattern(inside_outside_grid_part, outer_ansatz_space));
          this->inside_outside_patterns_[subdomain].insert(std::make_pair(neighboring_subomain, inside_outside_pattern));
          auto outside_inside_pattern = std::make_shared<PatternType>(
                outer_test_space.compute_face_pattern(outside_inside_grid_part, inner_ansatz_space));
          this->outside_inside_patterns_[neighboring_subomain].insert(std::make_pair(subdomain, outside_inside_pattern));
          // assemble
          CouplingAssembler coupling_assembler(inside_outside_grid_part,
                                               inner_test_space, inner_ansatz_space,
                                               outer_test_space, outer_ansatz_space);
          this->assemble_coupling_contributions(subdomain, neighboring_subomain, coupling_assembler);
        } // if (subdomain < neighboring_subomain)
      } // for (const size_t& neighboring_subomain : local_grids_provider_.ms_grid()->neighborsOf(subdomain))
    } // walk the subdomains

    // copy local and coupling patterns and container
    this->build_global_containers();

    logger.info() << "assembling products... " << std::endl;
    this->assemble_products(this->only_these_products_, 2);

    // finalize
    this->finalize_init(prune);

    logger.info() << "finished!" << std::endl;
  } // ... init(...)

private:
  class CouplingAssembler
    : public BaseType::template CouplingAssemblerBase<typename LocalGridsProviderType::MultiscaleGridProviderType::MsGridType::CouplingGridPartType::IntersectionType>
  {
    typedef typename BaseType::template CouplingAssemblerBase
      <typename LocalGridsProviderType::MultiscaleGridProviderType::MsGridType::CouplingGridPartType::IntersectionType>
        CouplingAssemblerBaseType;
    typedef typename LocalGridsProviderType::MultiscaleGridProviderType::MsGridType::CouplingGridPartType
        CouplingGridPartType;

  public:
    using typename CouplingAssemblerBaseType::LocalMatrixType;

    template <class... Args>
    CouplingAssembler(const CouplingGridPartType& coupling_grid_part, Args&& ...args)
      : CouplingAssemblerBaseType(std::forward< Args >(args)...)
      , coupling_grid_part_(coupling_grid_part)
    {}

    void assemble() const
    {
      // only do something, if there are local assemblers
      if (this->localCodim1MatrixAssemblers_.size() > 0) {
        // common tmp storage for all entities
        // * for the matrix assemblers
        std::vector< size_t > numberOfTmpMatricesNeeded(2, 0);
        for (auto& localCodim1MatrixAssembler : this->localCodim1MatrixAssemblers_) {
          const auto tmp = localCodim1MatrixAssembler->numTmpObjectsRequired();
          assert(tmp.size() == 2);
          numberOfTmpMatricesNeeded[0] = std::max(numberOfTmpMatricesNeeded[0], tmp[0]);
          numberOfTmpMatricesNeeded[1] = std::max(numberOfTmpMatricesNeeded[1], tmp[1]);
        }
        const size_t maxLocalSize = std::max(this->innerTestSpace_.mapper().maxNumDofs(),
                                             std::max(this->innerAnsatzSpace_.mapper().maxNumDofs(),
                                                      std::max(this->outerTestSpace_.mapper().maxNumDofs(),
                                                               this->outerAnsatzSpace_.mapper().maxNumDofs())));
        std::vector< LocalMatrixType > tmpLocalAssemblerMatrices(numberOfTmpMatricesNeeded[0],
                                                                  LocalMatrixType(maxLocalSize,
                                                                                  maxLocalSize,
                                                                                  R(0)));
        std::vector< LocalMatrixType > tmpLocalOperatorMatrices(numberOfTmpMatricesNeeded[1],
                                                                LocalMatrixType(maxLocalSize,
                                                                                maxLocalSize,
                                                                                R(0)));
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
        const auto entityEndIt = coupling_grid_part_.template end< 0 >();
        for(auto entityIt = coupling_grid_part_.template begin< 0 >(); entityIt != entityEndIt; ++entityIt ) {
          const auto& entity = *entityIt;
          // walk the intersections
          const auto intersectionEndIt = coupling_grid_part_.iend(entity);
          for (auto intersectionIt = coupling_grid_part_.ibegin(entity);
               intersectionIt != intersectionEndIt;
               ++intersectionIt) {
            const auto& intersection = *intersectionIt;
            // for a coupling grid part, we can be sure to only get the inner coupling intersetcions
            // so no further check neccesarry then
            assert(intersection.neighbor() && !intersection.boundary());
            // call local matrix assemblers
            for (auto& localCodim1MatrixAssembler : this->localCodim1MatrixAssemblers_) {
              localCodim1MatrixAssembler->apply(this->innerTestSpace_, this->innerAnsatzSpace_,
                                                this->outerTestSpace_, this->outerAnsatzSpace_,
                                                intersection,
                                                tmpLocalMatricesContainer, tmpIndices);
            }
          } // walk the intersections
        } // walk the grid
      } // only do something, if there are local assemblers
    } // void assemble() const

  private:
    const CouplingGridPartType& coupling_grid_part_;
  }; // class CouplingAssembler

  using BaseType::local_grids_provider_;
}; // DefaultMultiscaleGrid


//#if HAVE_ALUGRID && HAVE_DUNE_FEM
//# if HAVE_DUNE_ISTL

//extern template class DefaultMultiscaleGrid< ALUGrid< 2, 2, simplex, conforming, No_Comm >,
//                                             double,
//                                             1,
//                                             1,
//                                             Stuff::LA::ChooseBackend::istl_sparse >;

//#   if HAVE_MPI

//extern template class DefaultMultiscaleGrid< ALUGrid< 2, 2, simplex, conforming, MPI_Comm >,
//                                             double,
//                                             1,
//                                             1,
//                                             Stuff::LA::ChooseBackend::istl_sparse >;

//#   endif // HAVE_MPI
//# endif // HAVE_DUNE_ISTL
//# if HAVE_EIGEN

//extern template class DefaultMultiscaleGrid< ALUGrid< 2, 2, simplex, conforming, No_Comm >,
//                                             double,
//                                             1,
//                                             1,
//                                             Stuff::LA::ChooseBackend::eigen_sparse >;

//#   if HAVE_MPI

//extern template class DefaultMultiscaleGrid< ALUGrid< 2, 2, simplex, conforming, MPI_Comm >,
//                                             double,
//                                             1,
//                                             1,
//                                             Stuff::LA::ChooseBackend::eigen_sparse >;

//#   endif // HAVE_MPI
//# endif // HAVE_EIGEN
//#endif // HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_ISTL


template< Stuff::LA::ChooseBackend la, class G, class R = double, int r = 1, int p = 1 >
DefaultMultiscaleGrid< G, R, r, p, la > make_default(grid::Multiscale::ProviderInterface< G >& grid_provider,
                                                     const DSC::Configuration& boundary_info,
                                                     const ProblemInterface< typename G::template Codim< 0 >::Entity,
                                                                             typename G::ctype, G::dimension,
                                                                             R, r >& problem,
                                                     const std::vector< std::string >& only_these_products = {})
{
  return DefaultMultiscaleGrid< G, R, r, p, la >(grid_provider, boundary_info, problem, only_these_products);
}


} // namespace BlockSwipdg
} // namespace Discretizations
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_BLOCK_SWIPDG_HH
