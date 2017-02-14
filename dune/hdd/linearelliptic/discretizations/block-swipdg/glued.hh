// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_BLOCK_SWIPDG_GLUED_HH
#define DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_BLOCK_SWIPDG_GLUED_HH

#include <dune/stuff/common/memory.hh>

#include <dune/grid/multiscale/glued.hh>

#include <dune/gdt/playground/spaces/glued-block.hh>

#include "base.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Discretizations {
namespace BlockSwipdg {


// forward, needed in the Traits
template <class MacroG, class LocalG, class R = double, int r = 1, int p = 1,
          Stuff::LA::ChooseBackend la = Stuff::LA::default_sparse_backend>
class GluedMultiscaleGrid;


namespace internal {


template <class MacroG, class LocalG>
struct GluedMultiscaleGridLocalGridsProvider
{
  typedef grid::Multiscale::Glued<MacroG, LocalG> GluedGridType;
  typedef LocalG LocalGridType;

  GluedMultiscaleGridLocalGridsProvider(GluedGridType& glued_grid)
    : glued_grid_(glued_grid)
  {}

  GluedGridType& ms_grid()
  {
    return glued_grid_;
  }

  const GluedGridType& ms_grid() const
  {
    return glued_grid_;
  }

  size_t num_subdomains() const
  {
    return glued_grid_.num_subdomains();
  }

  typename GluedGridType::LocalGridProviderType& local_grid(const size_t subdomain)
  {
    return glued_grid_.local_grid(subdomain);
  }

  const typename GluedGridType::LocalGridProviderType& local_grid(const size_t subdomain) const
  {
    return glued_grid_.local_grid(subdomain);
  }

  int local_level(const size_t subdomain) const
  {
    return glued_grid_.max_local_level(subdomain);
  }

  grid::Multiscale::Glued<MacroG, LocalG>& glued_grid_;
}; // struct GluedMultiscaleGridLocalGridsProvider


template <class MacroG, class LocalG, class R, int r, int p, Stuff::LA::ChooseBackend la>
class GluedMultiscaleGridTraits
{
public:
  typedef GluedMultiscaleGridLocalGridsProvider<MacroG, LocalG> LocalGridsProviderType;
  static const constexpr Stuff::Grid::ChooseLayer local_layer   = Stuff::Grid::ChooseLayer::level;
  static const constexpr Stuff::Grid::ChooseLayer overlap_layer = Stuff::Grid::ChooseLayer::level;
  typedef GluedMultiscaleGrid<MacroG, LocalG, R, r, p, la> derived_type;

  template <class T>
  struct BlockSpace
  {
    typedef GDT::Spaces::GluedBlock<MacroG, LocalG, T> type;
  };
}; // class GluedMultiscaleGridTraits


} // namespace internal


template <class MacroG, class LocalG, class R, int r, int p, Stuff::LA::ChooseBackend la>
class GluedMultiscaleGrid
  : DSC::StorageProvider<internal::GluedMultiscaleGridLocalGridsProvider<MacroG, LocalG>>
  , public internal::Base<internal::GluedMultiscaleGridTraits<MacroG, LocalG, R, r, p, la>, R, r, p, la>
{
  typedef DSC::StorageProvider<internal::GluedMultiscaleGridLocalGridsProvider<MacroG, LocalG>> LocalGrids;
  typedef internal::Base
    <internal::GluedMultiscaleGridTraits<MacroG, LocalG, R, r, p, la>, R, r, p, la>             BaseType;
public:
  typedef internal::BaseTraits<internal::GluedMultiscaleGridTraits<MacroG, LocalG, R, r, p, la>, R, r, p, la> Traits;
  typedef typename Traits::LocalGridsProviderType LocalGridsProviderType;
  typedef MacroG MacroGridType;
  typedef LocalG LocalGridType;
  using typename BaseType::ProblemType;
  using typename BaseType::LocalTestSpaceType;
  using typename BaseType::LocalAnsatzSpaceType;
  using typename BaseType::PatternType;
  using typename BaseType::VectorType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".glued-multiscalegrid";
  }

  GluedMultiscaleGrid(grid::Multiscale::Glued<MacroG, LocalG>& glued_grid,
                      const Stuff::Common::Configuration& bound_inf_cfg,
                      const ProblemType& prob,
                      const std::vector< std::string >& only_these_products = {})
    : LocalGrids(new LocalGridsProviderType(glued_grid))
    , BaseType(LocalGrids::access(), bound_inf_cfg, prob, only_these_products)
    , glued_grid_(glued_grid)
  {}

  void init(const bool prune = false)
  {
    if (this->container_based_initialized_)
      return;

    auto logger = Stuff::Common::TimedLogger().get(static_id());
    const size_t num_subdomains = local_grids_provider_.num_subdomains();

    logger.info() << "discretizing on " << num_subdomains << " subdomains..." << std::endl;
    // has to be done first, initializes the local matrices and patterns required for the coupling
    for (auto&& macro_entity : DSC::entityRange(glued_grid_.macro_grid_view())) {
      const size_t subdomain = glued_grid_.subdomain(macro_entity);
      this->assemble_local_contributions(subdomain);
    }

    for (auto&& macro_entity : DSC::entityRange(glued_grid_.macro_grid_view())) {
      const size_t subdomain = glued_grid_.subdomain(macro_entity);
      if (glued_grid_.boundary(macro_entity)) {
        auto boundary_assembler = make_boundary_assembler(this->local_discretizations_[subdomain]->test_space(),
                                                          this->local_discretizations_[subdomain]->ansatz_space(),
                                                          glued_grid_.local_boundary_entities(macro_entity,
                                                                                              glued_grid_.max_local_level(macro_entity)),
                                                          this->local_discretizations_[subdomain]->test_space().grid_view());
        this->template assemble_boundary_contributions<typename LocalTestSpaceType::GridViewType>(subdomain, boundary_assembler);
      }

      for (auto&& macro_intersection : DSC::intersectionRange(glued_grid_.macro_grid_view(), macro_entity)) {
        if (!macro_intersection.boundary() && macro_intersection.neighbor()) {
          const auto macro_neighbor_ptr = macro_intersection.outside();
          const auto& macro_neighbor = *macro_neighbor_ptr;
          const size_t neighboring_subomain = glued_grid_.subdomain(macro_neighbor);
          if (subdomain < neighboring_subomain) {
            const auto& inner_test_space   = this->local_discretizations_[subdomain]->test_space();
            const auto& inner_ansatz_space = this->local_discretizations_[subdomain]->ansatz_space();
            const auto& outer_test_space   = this->local_discretizations_[neighboring_subomain]->test_space();
            const auto& outer_ansatz_space = this->local_discretizations_[neighboring_subomain]->ansatz_space();
            const auto& inside_outside_coupling = glued_grid_.coupling(macro_entity, glued_grid_.max_local_level(macro_entity),
                                                                       macro_neighbor, glued_grid_.max_local_level(macro_neighbor));
            const auto& outside_inside_coupling = glued_grid_.coupling(macro_neighbor, glued_grid_.max_local_level(macro_neighbor),
                                                                       macro_entity, glued_grid_.max_local_level(macro_entity));
            // compute patterns
            auto inside_outside_pattern = create_coupling_pattern(inside_outside_coupling, inner_test_space, outer_ansatz_space);
            this->inside_outside_patterns_[subdomain].insert(std::make_pair(neighboring_subomain, inside_outside_pattern));
            auto outside_inside_pattern = create_coupling_pattern(outside_inside_coupling, outer_test_space, inner_ansatz_space);
            this->outside_inside_patterns_[neighboring_subomain].insert(std::make_pair(subdomain, outside_inside_pattern));
            // assemble
            auto coupling_assembler = make_coupling_assembler(inner_test_space, inner_ansatz_space,
                                                              outer_test_space, outer_ansatz_space,
                                                              inside_outside_coupling);
            this->assemble_coupling_contributions(subdomain, neighboring_subomain, coupling_assembler);
          } // if (subdomain < neighboring_subomain)
        } // if (!intersection.boundary() && intersection.neighbor())
      } // for (auto&& macro_intersection : DSC::IntersectionRange(glued_grid_.macro_grid_view(), macro_entity))
    } // walk the subdomains

    // copy local and coupling patterns and container
    this->build_global_containers();

//    logger.info() << "assembling products... " << std::endl;
//    this->assemble_products(this->only_these_products_, 2);

    // finalize
    this->finalize_init(prune);

    logger.info() << "finished!" << std::endl;
  } // ... init(...)

  using BaseType::visualize;

  virtual void visualize(const VectorType& vector,
                         const std::string filename,
                         const std::string name,
                         const bool add_dirichlet = true,
                         Pymor::Parameter mu = Pymor::Parameter()) const
  {
    VectorType tmp = vector.copy();
    const auto vectors = this->available_vectors();
    if (add_dirichlet && std::find(vectors.begin(), vectors.end(), "dirichlet") != vectors.end()) {
      const auto dirichlet_vector = this->get_vector("dirichlet");
      if (dirichlet_vector.parametric()) {
        const Pymor::Parameter mu_dirichlet = this->map_parameter(mu, "dirichlet");
        if (mu_dirichlet.type() != dirichlet_vector.parameter_type())
          DUNE_THROW(Pymor::Exceptions::wrong_parameter_type,
                     mu_dirichlet.type() << " vs. " << dirichlet_vector.parameter_type());
        tmp += dirichlet_vector.freeze_parameter(mu);
      } else
        tmp += *(dirichlet_vector.affine_part());
    }
    auto local_vectors = localize_vector(tmp);
    typedef GDT::ConstDiscreteFunction<LocalAnsatzSpaceType, VectorType>
        LocalDiscreteFunctionType;
    std::vector<std::shared_ptr<LocalDiscreteFunctionType>> local_discrete_functions;
    // those are required if the local spaces live on grid parts
    typedef Stuff::Functions::VisualizationAdapter<typename LocalGridType::LevelGridView, r>
        LocalVisualizationAdapterType;
    std::vector<std::shared_ptr<LocalVisualizationAdapterType>> local_adapters;
    for (size_t ss = 0; ss < glued_grid_.num_subdomains(); ++ss) {
      local_discrete_functions.emplace_back(new LocalDiscreteFunctionType(*this->local_ansatz_spaces_[ss],
                                                                          local_vectors[ss],
                                                                          name));
      local_adapters.emplace_back(new LocalVisualizationAdapterType(*local_discrete_functions.back()));
    }
    grid::Multiscale::GluedVTKWriter<MacroGridType, LocalGridType> vtk_writer(glued_grid_);
    vtk_writer.addVertexData(local_adapters);
    vtk_writer.write(filename, VTK::appendedraw);
  } // ... visualize(...)

  VectorType localize_vector(const VectorType& global_vector, const size_t ss) const
  {
    if (ss >= glued_grid_.num_subdomains())
      DUNE_THROW(Stuff::Exceptions::index_out_of_range,
                 "0 <= ss < num_subdomains() = " << glued_grid_.num_subdomains() << " is not true for ss = " << ss << "!");
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

  std::vector<VectorType> localize_vector(const VectorType& global_vector) const
  {
    std::vector<VectorType> ret(glued_grid_.num_subdomains());
    for (size_t ss = 0; ss < glued_grid_.num_subdomains(); ++ss)
      ret[ss] = localize_vector(global_vector, ss);
    return ret;
  }

private:
  // We only use the GDT::SystemAssembler to collect the local functors via the the add(...) methods, so the spaces and
  // grid views do not matter, apart from their types.
  template <class EntityPointerType, class GridViewType>
  class BoundaryAssembler
    : public GDT::SystemAssembler<LocalTestSpaceType, GridViewType, LocalAnsatzSpaceType>
  {
    typedef GDT::SystemAssembler<LocalTestSpaceType, GridViewType, LocalAnsatzSpaceType> BaseType;
  public:
    BoundaryAssembler(const LocalTestSpaceType& tst_spc,
                      const LocalTestSpaceType& anstz_spc,
                      const std::vector<std::pair<EntityPointerType, std::vector<int>>>& boundary_entity_information,
                      const GridViewType& unused_grid_view)
      : BaseType(tst_spc, anstz_spc, unused_grid_view)
      , boundary_entity_information_(boundary_entity_information)
    {}

    void assemble()
    {
      prepare();
      if ((codim0_functors_.size() + codim1_functors_.size()) > 0) {
        // walk boundary entities
        for (const auto& element : boundary_entity_information_) {
          const auto& entity_pointer = element.first;
          const auto& entity = *entity_pointer;
          // apply codim 0 functors
          apply_local(entity);
          // walk those intersections which lie on the domain boundary
          const auto& local_boundary_intersections = element.second;
          for (auto&& intersection : DSC::intersectionRange(grid_view_, entity)) {
            const int local_intersection_index = intersection.indexInInside();
            if (std::find(local_boundary_intersections.begin(),
                          local_boundary_intersections.end(),
                          local_intersection_index) != local_boundary_intersections.end())
              apply_local(intersection, entity, entity);
          }
        }
      }
      finalize();
      clear();
    } // ... assemble(...)

  private:
    using BaseType::prepare;
    using BaseType::apply_local;
    using BaseType::finalize;
    using BaseType::clear;
    using BaseType::grid_view_;
    using BaseType::codim0_functors_;
    using BaseType::codim1_functors_;
    const std::vector<std::pair<EntityPointerType, std::vector<int>>>& boundary_entity_information_;
  }; // class BoundaryAssembler

  template <class EntityPointerType, class GV>
  BoundaryAssembler<EntityPointerType, GV> make_boundary_assembler(const LocalTestSpaceType& tst_spc,
                                                                   const LocalAnsatzSpaceType& anstz_spc,
                                                                   const std::vector<std::pair<EntityPointerType, std::vector<int>>>& boundary_entity_information,
                                                                   const GV& unused_grid_view) const
  {
    return BoundaryAssembler<EntityPointerType, GV>(tst_spc, anstz_spc, boundary_entity_information, unused_grid_view);
  }

  template <class CouplingGlueType>
  class CouplingAssembler
  {
    typedef Dune::DynamicMatrix< R > LocalMatrixType;
    typedef Dune::DynamicVector< R > LocalVectorType;
    typedef std::vector< std::vector< LocalMatrixType > > LocalMatricesContainerType;
    typedef std::vector< std::vector< LocalVectorType > > LocalVectorsContainerType;
    typedef std::vector< Dune::DynamicVector< size_t > > IndicesContainer;

    typedef typename CouplingGlueType::Intersection IntersectionType;

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
    CouplingAssembler(const LocalTestSpaceType& inner_test_space,
                      const LocalAnsatzSpaceType& inner_ansatz_space,
                      const LocalTestSpaceType& outer_test_space,
                      const LocalAnsatzSpaceType& outer_ansatz_space,
                      const CouplingGlueType& coupling_glue)
      : innerTestSpace_(inner_test_space)
      , innerAnsatzSpace_(inner_ansatz_space)
      , outerTestSpace_(outer_test_space)
      , outerAnsatzSpace_(outer_ansatz_space)
      , coupling_glue_(coupling_glue)
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
                                                                                  0));
        std::vector< LocalMatrixType > tmpLocalOperatorMatrices(numberOfTmpMatricesNeeded[1],
                                                                LocalMatrixType(maxLocalSize,
                                                                                maxLocalSize,
                                                                                0));
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

        const auto coupling_intersection_it_end = coupling_glue_.template iend<0>();
        for (auto coupling_intersection_it = coupling_glue_.template ibegin<0>();
             coupling_intersection_it != coupling_intersection_it_end;
             ++coupling_intersection_it) {
          const auto& intersection = *coupling_intersection_it;
          // call local matrix assemblers
          for (auto& localCodim1MatrixAssembler : localCodim1MatrixAssemblers_) {
            localCodim1MatrixAssembler->apply(innerTestSpace_, innerAnsatzSpace_,
                                              outerTestSpace_, outerAnsatzSpace_,
                                              intersection,
                                              tmpLocalMatricesContainer, tmpIndices);
          }
        } // walk the grid
      } // only do something, if there are local assemblers
    } // void assemble() const

  private:
    const LocalTestSpaceType& innerTestSpace_;
    const LocalAnsatzSpaceType& innerAnsatzSpace_;
    const LocalTestSpaceType& outerTestSpace_;
    const LocalAnsatzSpaceType& outerAnsatzSpace_;
    const CouplingGlueType coupling_glue_;
    std::vector< LocalCodim1MatrixAssemblerApplication* > localCodim1MatrixAssemblers_;
  }; // class CouplingAssembler

  template <class CouplingGlueType>
  CouplingAssembler<CouplingGlueType> make_coupling_assembler(const LocalTestSpaceType& innr_tst_spc,
                                                              const LocalAnsatzSpaceType& innr_anstz_spc,
                                                              const LocalTestSpaceType& outr_tst_spc,
                                                              const LocalAnsatzSpaceType& outr_anstz_spc,
                                                              const CouplingGlueType& coupling_glue) const
  {
    return CouplingAssembler<CouplingGlueType>(innr_tst_spc, innr_anstz_spc, outr_tst_spc, outr_anstz_spc, coupling_glue);
  }

  template <class CouplingGlueType, class InnerTestSpaceType, class OuterAnsatzSpaceType>
  std::shared_ptr<PatternType> create_coupling_pattern(const CouplingGlueType& coupling_glue,
                                                       const InnerTestSpaceType& inner_test_space,
                                                       const OuterAnsatzSpaceType& outer_ansatz_space) const
  {
    PatternType pattern(inner_test_space.mapper().size());
    DynamicVector<size_t> global_rows(inner_test_space.mapper().maxNumDofs(), 0);
    DynamicVector<size_t> global_cols(outer_ansatz_space.mapper().maxNumDofs(), 0);
    // walk the coupling
    const auto coupling_intersection_it_end = coupling_glue.template iend<0>();
    for (auto coupling_intersection_it = coupling_glue.template ibegin<0>();
         coupling_intersection_it != coupling_intersection_it_end;
         ++coupling_intersection_it) {
      const auto& coupling_intersection = *coupling_intersection_it;
      const auto inner_entity_ptr = coupling_intersection.inside();
      const auto outer_neighbor_ptr = coupling_intersection.outside();
      const auto& inner_entity = *inner_entity_ptr;
      const auto& outer_neighbor = *outer_neighbor_ptr;
      inner_test_space.mapper().globalIndices(inner_entity, global_rows);
      outer_ansatz_space.mapper().globalIndices(outer_neighbor, global_cols);
      const auto test_base_entity = inner_test_space.base_function_set(inner_entity);
      const auto ansatz_base_neighbour = outer_ansatz_space.base_function_set(outer_neighbor);
      for (size_t ii = 0; ii < test_base_entity.size(); ++ii)
        for (size_t jj = 0; jj < ansatz_base_neighbour.size(); ++jj)
          pattern.insert(global_rows[ii], global_cols[jj]);

    } // walk the coupling
    pattern.sort();
    return std::make_shared<PatternType>(std::move(pattern));
  } // ... create_coupling_pattern(...)

  using BaseType::local_grids_provider_;
  using BaseType::only_these_products_;

  grid::Multiscale::Glued<MacroG, LocalG>& glued_grid_;
}; // GluedMultiscaleGrid


template< Stuff::LA::ChooseBackend la, class M, class L, class R, int r = 1, int p = 1 >
GluedMultiscaleGrid< M, L, R, r, p, la > make_glued(grid::Multiscale::Glued<M, L>& glued_grid,
                                                              const DSC::Configuration& boundary_info,
                                                              const ProblemInterface< typename L::template Codim< 0 >::Entity,
                                                              typename L::ctype, L::dimension,
                                                              R, r >& problem,
                                                              const std::vector< std::string >& only_these_products = {})
{
  return GluedMultiscaleGrid< M, L, R, r, p, la >(glued_grid, boundary_info, problem, only_these_products);
}


} // namespace BlockSwipdg
} // namespace Discretizations
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_BLOCK_SWIPDG_GLUED_HH
