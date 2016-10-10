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

//  typedef GluedMultiscaleGrid< MacroGridImp, LocalGridImp, RangeFieldImp, rangeDim, polynomialOrder, la_backend > derived_type;
//  typedef RangeFieldImp     RangeFieldType;
//  static const unsigned int dimRange = rangeDim;
//  static const unsigned int polOrder = polynomialOrder;
//private:
//  friend class GluedMultiscaleGrid< MacroGridImp, LocalGridImp, RangeFieldImp, rangeDim, polynomialOrder, la_backend >;
//  typedef GluedMultiscaleGridLocalDiscretizationsContainer
//      < MacroGridImp, LocalGridImp, RangeFieldType, dimRange, polOrder, la_backend >
//      GluedMultiscaleGridLocalDiscretizationsContainerType;
//  typedef typename GluedMultiscaleGridLocalDiscretizationsContainerType::TestSpaceType   LocalTestSpaceType;
//  typedef typename GluedMultiscaleGridLocalDiscretizationsContainerType::AnsatzSpaceType LocalAnsatzSpaceType;
//public:
//  typedef GDT::Spaces::GluedBlock<MacroGridImp, LocalGridImp, LocalTestSpaceType>   TestSpaceType;
//  typedef GDT::Spaces::GluedBlock<MacroGridImp, LocalGridImp, LocalAnsatzSpaceType> AnsatzSpaceType;
//  typedef typename AnsatzSpaceType::GridViewType        GridViewType;
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
//  typedef internal::GluedMultiscaleGridTraits
//      < MacroGridImp, LocalGridImp, RangeFieldImp, rangeDim, polynomialOrder, la_backend >                     Traits;
//  typedef MacroGridImp MacroGridType;
//  typedef LocalGridImp LocalGridType;
  using typename BaseType::ProblemType;
//  using typename BaseType::GridViewType    GridViewType;
//  using typename BaseType::TestSpaceType   TestSpaceType;
//  using typename BaseType::AnsatzSpaceType AnsatzSpaceType;
//  using typename BaseType::EntityType      EntityType;
//  using typename BaseType::DomainFieldType DomainFieldType;
//  using typename BaseType::RangeFieldType  RangeFieldType;
//  using typename BaseType::MatrixType      MatrixType;
//  using typename BaseType::VectorType      VectorType;
//  using typename BaseType::OperatorType    OperatorType;
//  using typename BaseType::ProductType     ProductType;
//  using typename BaseType::FunctionalType  FunctionalType;

//  static const unsigned int dimDomain = BaseType::dimDomain;
//  static const unsigned int dimRange  = BaseType::dimRange;

////  typedef grid::Multiscale::ProviderInterface< GridImp > GridProviderType;
////  typedef typename GridProviderType::GridType   GridType;
////  typedef typename GridProviderType::MsGridType MsGridType;

//  typedef typename LocalDiscretizationsBaseType::DiscretizationType LocalDiscretizationType;
////  typedef typename Traits::LocalDiscretizationsContainerType::OversampledDiscretizationType OversampledDiscretizationType;
//  typedef typename TestSpaceType::PatternType                       PatternType;

//private:
//  using typename BaseType::AffinelyDecomposedMatrixType;
//  using typename BaseType::AffinelyDecomposedVectorType;
//  typedef Pymor::LA::AffinelyDecomposedConstContainer< MatrixType > AffinelyDecomposedConstMatrixType;
//  typedef Pymor::LA::AffinelyDecomposedConstContainer< VectorType > AffinelyDecomposedConstVectorType;

//public:
//  typedef typename LocalDiscretizationType::ProblemType LocalProblemType;
//  typedef typename LocalDiscretizationType::ProductType LocalProductType;

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
  {}

////  const std::vector< std::shared_ptr< LocalDiscretizationType > >& local_discretizations() const;

//  void init(const bool prune = false)
//  {
//    if (this->container_based_initialized_)
//      return;

//    auto logger = Stuff::Common::TimedLogger().get("hdd.linearelliptic.discretizations.glued-block-swipdg.init");
//    logger.info() << "discretizing on " << glued_grid_.num_subdomains() << " subdomains..." << std::endl;

//    pattern_ = std::make_shared< PatternType >(BaseType::test_space().mapper().size());

//    // walk the subdomains for the first time
//    //   * to initialize the coupling pattern,
//    //   * to finalize the global sparsity pattern
//    logger.info() << "computing patterns and local contributions... " << std::endl;
//    for (auto&& macro_entity : DSC::entityRange(glued_grid_.macro_grid_view())) {
//      const size_t ss = glued_grid_.macro_grid_view().indexSet().index(macro_entity);
//      // init the local discretizations (assembles matrices and patterns)
//      this->local_discretizations_[ss]->init(false);
//      // and the local boundary containers
//      local_boundary_matrices_[ss] = std::make_shared<AffinelyDecomposedMatrixType>();
//      local_boundary_vectors_[ss]  = std::make_shared<AffinelyDecomposedVectorType>();
////      // * the matrices
////      //   * just copy those from the local discretizations
////      const auto local_operator = this->local_discretizations_[ss]->get_operator();
////      local_boundary_matrices_[ss] = std::make_shared< AffinelyDecomposedMatrixType >();
////      //   * we take the affine part only if the diffusion has one, otherwise it contains only the dirichlet rows,
////      //     thus it is empty, since the local problems are purely neumann
////      if (this->problem().diffusion_factor()->has_affine_part()) {
////        if (!local_operator.has_affine_part())
////          DUNE_THROW(Stuff::Exceptions::internal_error, "The local operator is missing the affine part!");
////        local_boundary_matrices_[ss]->register_affine_part(new MatrixType(*(local_operator.affine_part().container())));
////      }
////      if (local_operator.num_components() < this->problem().diffusion_factor()->num_components())
////        DUNE_THROW(Stuff::Exceptions::requirements_not_met,
////                   "The local operator should have " << this->problem().diffusion_factor()->num_components()
////                   << " components (but has only " << local_operator.num_components() << ")!");
////      for (ssize_t qq = 0; qq < this->problem().diffusion_factor()->num_components(); ++qq) {
////        local_boundary_matrices_[ss]->register_component(new MatrixType(
////            local_operator.component(qq).container()->backend()),
////            this->problem().diffusion_factor()->coefficient(qq));
////      }
////      // * and the vectors
////      const auto local_functional = this->local_discretizations_[ss]->get_rhs();
////      local_boundary_vectors_[ss] = std::make_shared< AffinelyDecomposedVectorType >();
////      for (size_t qq = 0; qq < boost::numeric_cast< size_t >(local_functional.num_components()); ++qq)
////        local_boundary_vectors_[ss]->register_component(new VectorType(*(local_functional.component(qq).container())),
////                                               new Pymor::ParameterFunctional(local_functional.coefficient(qq)));
////      if (local_functional.has_affine_part())
////        local_boundary_vectors_[ss]->register_affine_part(new VectorType(*(local_functional.affine_part().container())));

//      // create and copy the local patterns
//      add_local_to_global_pattern(this->local_discretizations_[ss]->pattern(), ss, ss, *pattern_);
//      const auto& inner_test_space = this->local_discretizations_[ss]->test_space();
//      const auto& inner_ansatz_space = this->local_discretizations_[ss]->ansatz_space();
//      // walk the neighbors
//      for (auto&& intersection : DSC::intersectionRange(glued_grid_.macro_grid_view(), macro_entity)) {
//        if (intersection.neighbor() && !intersection.boundary()) {
//          const auto macro_neighbor_ptr = intersection.outside();
//          const auto& macro_neighbor = *macro_neighbor_ptr;
//          const size_t nn = glued_grid_.macro_grid_view().indexSet().index(macro_neighbor);
//          // visit each coupling only once (assemble primally)
//          if (ss < nn) {
//            const auto& outer_test_space = this->local_discretizations_[nn]->test_space();
//            const auto& outer_ansatz_space = this->local_discretizations_[nn]->ansatz_space();
//            const auto& inside_outside_coupling = glued_grid_.coupling(macro_entity, glued_grid_.max_local_level(macro_entity),
//                                                                       macro_neighbor, glued_grid_.max_local_level(macro_neighbor));
//            const auto& outside_inside_coupling = glued_grid_.coupling(macro_neighbor, glued_grid_.max_local_level(macro_neighbor),
//                                                                       macro_entity, glued_grid_.max_local_level(macro_entity));
//            // create the coupling patterns
//            auto inside_outside_pattern = create_coupling_pattern(inside_outside_coupling, inner_test_space, outer_ansatz_space);
//            inside_outside_patterns_[ss].insert(std::make_pair(nn, inside_outside_pattern));
//            auto outside_inside_pattern = create_coupling_pattern(outside_inside_coupling, outer_test_space, inner_ansatz_space);
//            outside_inside_patterns_[nn].insert(std::make_pair(ss, outside_inside_pattern));
//            // and copy them
//            add_local_to_global_pattern(*inside_outside_pattern,  ss, nn, *pattern_);
//            add_local_to_global_pattern(*outside_inside_pattern,  nn, ss, *pattern_);
//          } // visit each coupling only once (assemble primaly)
//        }
//      } // walk the neighbors
//    } // walk the subdomains for the first time

//    // walk the subdomains for the second time
//    //   * to assemble the boundary matrices and vectors and
//    //   * to assemble the coupling matrices
//    logger.info() << "computing coupling and boundary contributions... " << std::endl;
//    for (auto&& macro_entity : DSC::entityRange(glued_grid_.macro_grid_view())) {
//      const size_t ss = glued_grid_.macro_grid_view().indexSet().index(macro_entity);
//      const auto& inner_test_mapper = this->local_discretizations_[ss]->test_space().mapper();
//      const auto& inner_ansatz_mapper = this->local_discretizations_[ss]->ansatz_space().mapper();
//      for (auto&& macro_intersection : DSC::intersectionRange(glued_grid_.macro_grid_view(), macro_entity)) {
//        if (macro_intersection.boundary() && !macro_intersection.neighbor())
//          assemble_boundary_contributions(macro_entity);
//        else if (!macro_intersection.boundary() && macro_intersection.neighbor()) {
//          const auto macro_neighbor_ptr = macro_intersection.outside();
//          const auto& macro_neighbor = *macro_neighbor_ptr;
//          const size_t nn = glued_grid_.macro_grid_view().indexSet().index(macro_neighbor);
//          // visit each coupling only once (assemble primaly)
//          if (ss < nn) {
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
//          assemble_coupling_contributions(macro_entity, macro_neighbor,
//                                          *(this->local_discretizations_[ss]->system_matrix()),
//                                          *(inside_outside_matrix),
//                                          *(outside_inside_matrix),
//                                          *(this->local_discretizations_[nn]->system_matrix()));
//          inside_outside_matrices_[ss].insert(std::make_pair(nn, inside_outside_matrix));
//          outside_inside_matrices_[nn].insert(std::make_pair(ss, outside_inside_matrix));
//          } // visit each coupling only once
//        } else
//          DUNE_THROW(Stuff::Exceptions::internal_error, "Unknown intersection type encountered!");
//      } // walk the intersection of the macro entity
//    } // walk the subdomains for the second time

//    // build global containers
//    pattern_->sort();
//    build_global_containers();

//    logger.info() << "assembling products... " << std::endl;
//    build_products();

//    // finalize
//    this->finalize_init(prune);

//    logger.info() << "finished!" << std::endl;
//  } // ... init(...)

//  using BaseType::visualize;

//  virtual void visualize(const VectorType& vector,
//                         const std::string filename,
//                         const std::string name,
//                         const bool add_dirichlet = true,
//                         Pymor::Parameter mu = Pymor::Parameter()) const
//  {
//    VectorType tmp = vector.copy();
//    const auto vectors = this->available_vectors();
//    if (add_dirichlet && std::find(vectors.begin(), vectors.end(), "dirichlet") != vectors.end()) {
//      const auto dirichlet_vector = this->get_vector("dirichlet");
//      if (dirichlet_vector.parametric()) {
//        const Pymor::Parameter mu_dirichlet = this->map_parameter(mu, "dirichlet");
//        if (mu_dirichlet.type() != dirichlet_vector.parameter_type())
//          DUNE_THROW(Pymor::Exceptions::wrong_parameter_type,
//                     mu_dirichlet.type() << " vs. " << dirichlet_vector.parameter_type());
//        tmp += dirichlet_vector.freeze_parameter(mu);
//      } else
//        tmp += *(dirichlet_vector.affine_part());
//    }
//    auto local_vectors = localize_vector(tmp);
//    typedef GDT::ConstDiscreteFunction<typename LocalDiscretizationType::AnsatzSpaceType, VectorType>
//        LocalDiscreteFunctionType;
//    std::vector<std::shared_ptr<LocalDiscreteFunctionType>> local_discrete_functions;
//    // those are required if the local spaces live on grid parts
//    typedef Stuff::Functions::VisualizationAdapter<typename LocalGridType::LevelGridView, dimRange>
//        LocalVisualizationAdapterType;
//    std::vector<std::shared_ptr<LocalVisualizationAdapterType>> local_adapters;
//    for (size_t ss = 0; ss < glued_grid_.num_subdomains(); ++ss) {
//      local_discrete_functions.emplace_back(new LocalDiscreteFunctionType(*this->local_ansatz_spaces_[ss],
//                                                                          local_vectors[ss],
//                                                                          name));
//      local_adapters.emplace_back(new LocalVisualizationAdapterType(*local_discrete_functions.back()));
//    }
//    grid::Multiscale::GluedVTKWriter<MacroGridType, LocalGridType> vtk_writer(glued_grid_);
//    vtk_writer.addVertexData(local_adapters);
//    vtk_writer.write(filename, VTK::appendedraw);
//  } // ... visualize(...)

////  ssize_t num_subdomains() const;

////  std::vector< ssize_t > neighbouring_subdomains(const ssize_t ss) const;

//  VectorType localize_vector(const VectorType& global_vector, const size_t ss) const
//  {
//    if (ss >= glued_grid_.num_subdomains())
//      DUNE_THROW(Stuff::Exceptions::index_out_of_range,
//                 "0 <= ss < num_subdomains() = " << glued_grid_.num_subdomains() << " is not true for ss = " << ss << "!");
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

//  std::vector<VectorType> localize_vector(const VectorType& global_vector) const
//  {
//    std::vector<VectorType> ret(glued_grid_.num_subdomains());
//    for (size_t ss = 0; ss < glued_grid_.num_subdomains(); ++ss)
//      ret[ss] = localize_vector(global_vector, ss);
//    return ret;
//  }

////  VectorType globalize_vectors(const std::vector< VectorType >& local_vectors) const;

////  VectorType* globalize_vectors_and_return_ptr(const std::vector< VectorType >& local_vectors) const;

////  VectorType* localize_vector_and_return_ptr(const VectorType& global_vector, const ssize_t ss) const;

////  ProductType get_local_product(const size_t ss, const std::string id) const;

////  ProductType* get_local_product_and_return_ptr(const ssize_t ss, const std::string id) const;

////  OperatorType get_local_operator(const size_t ss) const;

////  OperatorType* get_local_operator_and_return_ptr(const ssize_t ss) const;

////  OperatorType get_coupling_operator(const size_t ss, const size_t nn) const;

////  OperatorType* get_coupling_operator_and_return_ptr(const ssize_t ss, const ssize_t nn) const;

////  FunctionalType get_local_functional(const size_t ss) const;

////  FunctionalType* get_local_functional_and_return_ptr(const ssize_t ss) const;

////  VectorType solve_for_local_correction(const std::vector< VectorType >& local_vectors,
////                                        const size_t subdomain,
////                                        const Pymor::Parameter mu = Pymor::Parameter()) const;

////  LocalDiscretizationType get_local_discretization(const size_t subdomain) const;

////  LocalDiscretizationType* pb_get_local_discretization(const ssize_t subdomain) const;

////  LocalDiscretizationType get_oversampled_discretization(const size_t subdomain,
////                                                               const std::string boundary_value_type) const;

////  LocalDiscretizationType* pb_get_oversampled_discretization(const ssize_t subdomain,
////                                                                   const std::string boundary_value_type) const;

//private:
//  // we only use the GDT::SystemAssembler to obtain the add(...) methods, so the spaces and grid views do not matter,
//  // apart from their types
//  template <class EntityPointerType>
//  class BoundaryAssembler
//    : public GDT::SystemAssembler<typename LocalDiscretizationType::TestSpaceType,
//                                  GridViewType,
//                                  typename LocalDiscretizationType::AnsatzSpaceType>
//  {
//    typedef GDT::SystemAssembler<typename LocalDiscretizationType::TestSpaceType,
//                                 GridViewType,
//                                 typename LocalDiscretizationType::AnsatzSpaceType> BaseType;
//  public:
//    BoundaryAssembler(const typename LocalDiscretizationType::TestSpaceType& tst_spc,
//                      const typename LocalDiscretizationType::AnsatzSpaceType& anstz_spc,
//                      const std::vector<std::pair<EntityPointerType, std::vector<int>>>& boundary_entity_information,
//                      const GridViewType& unused_grid_view)
//      : BaseType(tst_spc, anstz_spc, unused_grid_view)
//      , boundary_entity_information_(boundary_entity_information)
//    {}

//    void assemble()
//    {
//      prepare();
//      if ((codim0_functors_.size() + codim1_functors_.size()) > 0) {
//        // walk boundary entities
//        for (const auto& element : boundary_entity_information_) {
//          const auto& entity_pointer = element.first;
//          const auto& entity = *entity_pointer;
//          // apply codim 0 functors
//          apply_local(entity);
//          // walk those intersections which lie on the domain boundary
//          const auto& local_boundary_intersections = element.second;
//          for (auto&& intersection : DSC::intersectionRange(grid_view_, entity)) {
//            const int local_intersection_index = intersection.indexInInside();
//            if (std::find(local_boundary_intersections.begin(),
//                          local_boundary_intersections.end(),
//                          local_intersection_index) != local_boundary_intersections.end())
//              apply_local(intersection, entity, entity);
//          }
//        }
//      }
//      finalize();
//      clear();
//    } // ... assemble(...)

//  private:
//    using BaseType::prepare;
//    using BaseType::apply_local;
//    using BaseType::finalize;
//    using BaseType::clear;
//    using BaseType::grid_view_;
//    using BaseType::codim0_functors_;
//    using BaseType::codim1_functors_;
//    const std::vector<std::pair<EntityPointerType, std::vector<int>>>& boundary_entity_information_;
//  }; // class BoundaryAssembler

//  template <class EntityPointerType>
//  BoundaryAssembler<EntityPointerType> make_boundary_assembler(const typename LocalDiscretizationType::TestSpaceType& tst_spc,
//                                                               const typename LocalDiscretizationType::AnsatzSpaceType& anstz_spc,
//                                                               const std::vector<std::pair<EntityPointerType, std::vector<int>>>& boundary_entity_information,
//                                                               const GridViewType& unused_grid_view) const
//  {
//    return BoundaryAssembler<EntityPointerType>(tst_spc, anstz_spc, boundary_entity_information, unused_grid_view);
//  }

//  template <class CouplingGlueType>
//  class CouplingAssembler
//  {
//    typedef Dune::DynamicMatrix< RangeFieldType > LocalMatrixType;
//    typedef Dune::DynamicVector< RangeFieldType > LocalVectorType;
//    typedef std::vector< std::vector< LocalMatrixType > > LocalMatricesContainerType;
//    typedef std::vector< std::vector< LocalVectorType > > LocalVectorsContainerType;
//    typedef std::vector< Dune::DynamicVector< size_t > > IndicesContainer;

//    typedef typename LocalDiscretizationType::TestSpaceType LocalTestSpaceType;
//    typedef typename LocalDiscretizationType::AnsatzSpaceType LocalAnsatzSpaceType;

//    typedef typename CouplingGlueType::Intersection IntersectionType;

//    class LocalCodim1MatrixAssemblerApplication
//    {
//    public:
//      virtual ~LocalCodim1MatrixAssemblerApplication(){}

//      virtual void apply(const LocalTestSpaceType& /*inner_test_space*/,
//                         const LocalAnsatzSpaceType& /*inner_ansatz_space*/,
//                         const LocalTestSpaceType& /*outer_test_space*/,
//                         const LocalAnsatzSpaceType& /*outer_ansatz_space*/,
//                         const IntersectionType& /*_intersection*/,
//                         LocalMatricesContainerType& /*_localMatricesContainer*/,
//                         IndicesContainer& /*indicesContainer*/) const = 0;

//      virtual std::vector< size_t > numTmpObjectsRequired() const = 0;
//    };

//    template< class LocalAssemblerType, class M >
//    class LocalCodim1MatrixAssemblerWrapper
//      : public LocalCodim1MatrixAssemblerApplication
//    {
//    public:
//      LocalCodim1MatrixAssemblerWrapper(const LocalAssemblerType& localAssembler,
//                                        Dune::Stuff::LA::MatrixInterface< M >& in_in_matrix,
//                                        Dune::Stuff::LA::MatrixInterface< M >& in_out_matrix,
//                                        Dune::Stuff::LA::MatrixInterface< M >& out_in_matrix,
//                                        Dune::Stuff::LA::MatrixInterface< M >& out_out_matrix)
//        : localMatrixAssembler_(localAssembler)
//        , in_in_matrix_(in_in_matrix)
//        , out_out_matrix_(out_out_matrix)
//        , in_out_matrix_(in_out_matrix)
//        , out_in_matrix_(out_in_matrix)
//      {}

//      virtual void apply(const LocalTestSpaceType& inner_test_space,
//                         const LocalAnsatzSpaceType& inner_ansatz_space,
//                         const LocalTestSpaceType& outer_test_space,
//                         const LocalAnsatzSpaceType& outer_ansatz_space,
//                         const IntersectionType& intersection,
//                         LocalMatricesContainerType& localMatricesContainer,
//                         IndicesContainer& indicesContainer) const
//      {
//        localMatrixAssembler_.assembleLocal(inner_test_space, inner_ansatz_space,
//                                            outer_test_space, outer_ansatz_space,
//                                            intersection,
//                                            in_in_matrix_, out_out_matrix_, in_out_matrix_, out_in_matrix_,
//                                            localMatricesContainer, indicesContainer);
//      }

//      virtual std::vector< size_t > numTmpObjectsRequired() const
//      {
//        return localMatrixAssembler_.numTmpObjectsRequired();
//      }

//    private:
//      const LocalAssemblerType& localMatrixAssembler_;
//      Dune::Stuff::LA::MatrixInterface< M >& in_in_matrix_;
//      Dune::Stuff::LA::MatrixInterface< M >& out_out_matrix_;
//      Dune::Stuff::LA::MatrixInterface< M >& in_out_matrix_;
//      Dune::Stuff::LA::MatrixInterface< M >& out_in_matrix_;
//    }; // class LocalCodim1MatrixAssemblerWrapper

//  public:
//    CouplingAssembler(const LocalTestSpaceType& inner_test_space,
//                      const LocalAnsatzSpaceType& inner_ansatz_space,
//                      const LocalTestSpaceType& outer_test_space,
//                      const LocalAnsatzSpaceType& outer_ansatz_space,
//                      const CouplingGlueType& coupling_glue)
//      : innerTestSpace_(inner_test_space)
//      , innerAnsatzSpace_(inner_ansatz_space)
//      , outerTestSpace_(outer_test_space)
//      , outerAnsatzSpace_(outer_ansatz_space)
//      , coupling_glue_(coupling_glue)
//    {}

//    ~CouplingAssembler()
//    {
//      clearLocalAssemblers();
//    }

//    void clearLocalAssemblers()
//    {
//      for (auto& element: localCodim1MatrixAssemblers_)
//        delete element;
//    }

//    template< class L, class M >
//    void addLocalAssembler(const GDT::LocalAssembler::Codim1CouplingMatrix< L >& localAssembler,
//                           Dune::Stuff::LA::MatrixInterface< M >& in_in_matrix,
//                           Dune::Stuff::LA::MatrixInterface< M >& in_out_matrix,
//                           Dune::Stuff::LA::MatrixInterface< M >& out_in_matrix,
//                           Dune::Stuff::LA::MatrixInterface< M >& out_out_matrix)
//    {
//      assert(in_in_matrix.rows() == innerTestSpace_.mapper().size());
//      assert(in_in_matrix.cols() == innerAnsatzSpace_.mapper().size());
//      assert(in_out_matrix.rows() == innerTestSpace_.mapper().size());
//      assert(in_out_matrix.cols() == outerAnsatzSpace_.mapper().size());
//      assert(out_in_matrix.rows() == outerTestSpace_.mapper().size());
//      assert(out_in_matrix.cols() == innerAnsatzSpace_.mapper().size());
//      assert(out_out_matrix.rows() == outerTestSpace_.mapper().size());
//      assert(out_out_matrix.cols() == outerAnsatzSpace_.mapper().size());
//      localCodim1MatrixAssemblers_.push_back(
//            new LocalCodim1MatrixAssemblerWrapper< GDT::LocalAssembler::Codim1CouplingMatrix< L >, M >(
//              localAssembler, in_in_matrix, in_out_matrix, out_in_matrix, out_out_matrix));
//    }

//    void assemble() const
//    {
//      // only do something, if there are local assemblers
//      if (localCodim1MatrixAssemblers_.size() > 0) {
//        // common tmp storage for all entities
//        // * for the matrix assemblers
//        std::vector< size_t > numberOfTmpMatricesNeeded(2, 0);
//        for (auto& localCodim1MatrixAssembler : localCodim1MatrixAssemblers_) {
//          const auto tmp = localCodim1MatrixAssembler->numTmpObjectsRequired();
//          assert(tmp.size() == 2);
//          numberOfTmpMatricesNeeded[0] = std::max(numberOfTmpMatricesNeeded[0], tmp[0]);
//          numberOfTmpMatricesNeeded[1] = std::max(numberOfTmpMatricesNeeded[1], tmp[1]);
//        }
//        const size_t maxLocalSize = std::max(innerTestSpace_.mapper().maxNumDofs(),
//                                             std::max(innerAnsatzSpace_.mapper().maxNumDofs(),
//                                                      std::max(outerTestSpace_.mapper().maxNumDofs(),
//                                                               outerAnsatzSpace_.mapper().maxNumDofs())));
//        std::vector< LocalMatrixType > tmpLocalAssemblerMatrices( numberOfTmpMatricesNeeded[0],
//                                                                  LocalMatrixType(maxLocalSize,
//                                                                                  maxLocalSize,
//                                                                                  RangeFieldType(0)));
//        std::vector< LocalMatrixType > tmpLocalOperatorMatrices(numberOfTmpMatricesNeeded[1],
//                                                                LocalMatrixType(maxLocalSize,
//                                                                                maxLocalSize,
//                                                                                RangeFieldType(0)));
//        std::vector< std::vector< LocalMatrixType > > tmpLocalMatricesContainer;
//        tmpLocalMatricesContainer.push_back(tmpLocalAssemblerMatrices);
//        tmpLocalMatricesContainer.push_back(tmpLocalOperatorMatrices);
//        // * for the global indices
//        std::vector< Dune::DynamicVector< size_t > > tmpIndices = {
//            Dune::DynamicVector< size_t >(maxLocalSize)
//          , Dune::DynamicVector< size_t >(maxLocalSize)
//          , Dune::DynamicVector< size_t >(maxLocalSize)
//          , Dune::DynamicVector< size_t >(maxLocalSize)
//        };

//        const auto coupling_intersection_it_end = coupling_glue_.template iend<0>();
//        for (auto coupling_intersection_it = coupling_glue_.template ibegin<0>();
//             coupling_intersection_it != coupling_intersection_it_end;
//             ++coupling_intersection_it) {
//          const auto& intersection = *coupling_intersection_it;
//          // call local matrix assemblers
//          for (auto& localCodim1MatrixAssembler : localCodim1MatrixAssemblers_) {
//            localCodim1MatrixAssembler->apply(innerTestSpace_, innerAnsatzSpace_,
//                                              outerTestSpace_, outerAnsatzSpace_,
//                                              intersection,
//                                              tmpLocalMatricesContainer, tmpIndices);
//          }
//        } // walk the grid
//      } // only do something, if there are local assemblers
//    } // void assemble() const

//  private:
//    const LocalTestSpaceType& innerTestSpace_;
//    const LocalAnsatzSpaceType& innerAnsatzSpace_;
//    const LocalTestSpaceType& outerTestSpace_;
//    const LocalAnsatzSpaceType& outerAnsatzSpace_;
//    const CouplingGlueType coupling_glue_;
//    std::vector< LocalCodim1MatrixAssemblerApplication* > localCodim1MatrixAssemblers_;
//  }; // class CouplingAssembler

//  template <class CouplingGlueType>
//  CouplingAssembler<CouplingGlueType> make_coupling_assembler(const typename LocalDiscretizationType::TestSpaceType& innr_tst_spc,
//                                                              const typename LocalDiscretizationType::AnsatzSpaceType& innr_anstz_spc,
//                                                              const typename LocalDiscretizationType::TestSpaceType& outr_tst_spc,
//                                                              const typename LocalDiscretizationType::AnsatzSpaceType& outr_anstz_spc,
//                                                              const CouplingGlueType& coupling_glue) const
//  {
//    return CouplingAssembler<CouplingGlueType>(innr_tst_spc, innr_anstz_spc, outr_tst_spc, outr_anstz_spc, coupling_glue);
//  }

//  void add_local_to_global_pattern(const PatternType& local,
//                                   const size_t test_subdomain,
//                                   const size_t ansatz_subdomain,
//                                   PatternType& global) const
//  {
//    for (size_t local_ii = 0; local_ii < local.size(); ++local_ii) {
//      const size_t global_ii = this->test_space().mapper().mapToGlobal(test_subdomain, local_ii);
//      const auto& local_rows = local.inner(local_ii);
//      for (const auto& local_jj : local_rows) {
//        const size_t global_jj = this->ansatz_space().mapper().mapToGlobal(ansatz_subdomain, local_jj);
//        global.insert(global_ii, global_jj);
//      }
//    }
//  } // ... add_local_to_global_pattern(...)

//  template <class CouplingGlueType, class InnerTestSpaceType, class OuterAnsatzSpaceType>
//  std::shared_ptr<PatternType> create_coupling_pattern(const CouplingGlueType& coupling_glue,
//                                                       const InnerTestSpaceType& inner_test_space,
//                                                       const OuterAnsatzSpaceType& outer_ansatz_space) const
//  {
//    PatternType pattern(inner_test_space.mapper().size());
//    DynamicVector<size_t> global_rows(inner_test_space.mapper().maxNumDofs(), 0);
//    DynamicVector<size_t> global_cols(outer_ansatz_space.mapper().maxNumDofs(), 0);
//    // walk the coupling
//    const auto coupling_intersection_it_end = coupling_glue.template iend<0>();
//    for (auto coupling_intersection_it = coupling_glue.template ibegin<0>();
//         coupling_intersection_it != coupling_intersection_it_end;
//         ++coupling_intersection_it) {
//      const auto& coupling_intersection = *coupling_intersection_it;
//      const auto inner_entity_ptr = coupling_intersection.inside();
//      const auto outer_neighbor_ptr = coupling_intersection.outside();
//      const auto& inner_entity = *inner_entity_ptr;
//      const auto& outer_neighbor = *outer_neighbor_ptr;
//      inner_test_space.mapper().globalIndices(inner_entity, global_rows);
//      outer_ansatz_space.mapper().globalIndices(outer_neighbor, global_cols);
//      const auto test_base_entity = inner_test_space.base_function_set(inner_entity);
//      const auto ansatz_base_neighbour = outer_ansatz_space.base_function_set(outer_neighbor);
//      for (size_t ii = 0; ii < test_base_entity.size(); ++ii)
//        for (size_t jj = 0; jj < ansatz_base_neighbour.size(); ++jj)
//          pattern.insert(global_rows[ii], global_cols[jj]);

//    } // walk the coupling
//    pattern.sort();
//    return std::make_shared<PatternType>(std::move(pattern));
//  } // ... create_coupling_pattern(...)

//  void copy_local_to_global_matrix(const AffinelyDecomposedConstMatrixType& local_matrix,
//                                   const PatternType& local_pattern,
//                                   const size_t subdomain,
//                                   const size_t neighbor,
//                                   AffinelyDecomposedMatrixType& global_matrix) const
//  {
//    for (size_t qq = 0; qq < boost::numeric_cast< size_t >(local_matrix.num_components()); ++qq) {
//      const auto coefficient = local_matrix.coefficient(qq);
//      ssize_t comp = find_component(global_matrix, *coefficient);
//      if (comp < 0)
//        comp = global_matrix.register_component(coefficient,
//                                                this->test_space().mapper().size(),
//                                                this->ansatz_space().mapper().size(),
//                                                *pattern_);
//      assert(comp >= 0);
//      copy_local_to_global_matrix(*(local_matrix.component(qq)),
//                                  local_pattern,
//                                  subdomain,
//                                  neighbor,
//                                  *(global_matrix.component(comp)));
//    }
//    if (local_matrix.has_affine_part()) {
//      if (!global_matrix.has_affine_part())
//        global_matrix.register_affine_part(this->test_space_.mapper().size(),
//                                           this->ansatz_space_.mapper().size(),
//                                           *pattern_);
//      copy_local_to_global_matrix(*(local_matrix.affine_part()),
//                                  local_pattern,
//                                  subdomain,
//                                  neighbor,
//                                  *(global_matrix.affine_part()));
//    }
//  } // copy_local_to_global_matrix(...)

//  template< class ML, class MG >
//  void copy_local_to_global_matrix(const Stuff::LA::MatrixInterface< ML >& local_matrix,
//                                   const PatternType& local_pattern,
//                                   const size_t test_subdomain,
//                                   const size_t ansatz_subdomain,
//                                   Stuff::LA::MatrixInterface< MG >& global_matrix) const
//  {
//    for (size_t local_ii = 0; local_ii < local_pattern.size(); ++local_ii) {
//      const size_t global_ii = this->test_space().mapper().mapToGlobal(test_subdomain, local_ii);
//      for (const size_t& local_jj : local_pattern.inner(local_ii)) {
//        const size_t global_jj = this->ansatz_space().mapper().mapToGlobal(ansatz_subdomain, local_jj);
//        global_matrix.add_to_entry(global_ii, global_jj, local_matrix.get_entry(local_ii, local_jj));
//      }
//    }
//  } // ... copy_local_to_global_matrix(...)

//  void copy_local_to_global_vector(const AffinelyDecomposedConstVectorType& local_vector,
//                                   const size_t subdomain,
//                                   AffinelyDecomposedVectorType& global_vector) const
//  {
//    for (size_t qq = 0; qq < boost::numeric_cast< size_t >(local_vector.num_components()); ++qq) {
//      const auto coefficient = local_vector.coefficient(qq);
//      ssize_t comp = find_component(global_vector, *coefficient);
//      if (comp < 0)
//        comp = global_vector.register_component(coefficient,
//                                                this->test_space().mapper().size());
//      assert(comp >= 0);
//      copy_local_to_global_vector(*(local_vector.component(qq)),
//                                  subdomain,
//                                  *(global_vector.component(comp)));
//    }
//    if (local_vector.has_affine_part()) {
//      if (!global_vector.has_affine_part())
//        global_vector.register_affine_part(this->test_space().mapper().size());
//      copy_local_to_global_vector(*(local_vector.affine_part()),
//                                  subdomain,
//                                  *(global_vector.affine_part()));
//    }
//  } // copy_local_to_global_vector(...)

//  template< class VL, class VG >
//  void copy_local_to_global_vector(const Stuff::LA::VectorInterface< VL >& local_vector,
//                                   const size_t subdomain,
//                                   Stuff::LA::VectorInterface< VG >& global_vector) const
//  {
//    for (size_t local_ii = 0; local_ii < local_vector.size(); ++local_ii) {
//      const size_t global_ii = this->test_space().mapper().mapToGlobal(subdomain, local_ii);
//      global_vector.add_to_entry(global_ii, local_vector.get_entry(local_ii));
//    }
//  } // ... copy_local_to_global_vector(...)

//  template <class MacroEntityType>
//  void assemble_boundary_contributions(const MacroEntityType& macro_entity) const
//  {
////    auto logger = Stuff::Common::TimedLogger().get("hdd.linearelliptic.discretizations.block-swipdg.assemble_boundary_contributions");
//    const size_t subdomain = glued_grid_.macro_grid_view().indexSet().index(macro_entity);
//    const auto& local_test_space   = this->local_discretizations_[subdomain]->test_space();
//    const auto& local_ansatz_space = this->local_discretizations_[subdomain]->ansatz_space();
//    const auto& local_pattern = this->local_discretizations_[subdomain]->pattern();
//    auto& local_matrix = *this->local_boundary_matrices_[subdomain];
//    auto& local_vector = *this->local_boundary_vectors_[subdomain];
//    const auto& boundary_entities = glued_grid_.local_boundary_entities(macro_entity,
//                                                                        glued_grid_.max_local_level(macro_entity));
//    auto unused_local_view = glued_grid_.local_grid(subdomain).level_view(0);
//    auto boundary_assembler = make_boundary_assembler(local_test_space,
//                                                      local_ansatz_space,
//                                                      boundary_entities,
//                                                      unused_local_view); // the last arguments is required but not used
//    // lhs
//    // * dirichlet boundary terms
//    typedef typename ProblemType::DiffusionFactorType::NonparametricType DiffusionFactorType;
//    typedef typename ProblemType::DiffusionTensorType::NonparametricType DiffusionTensorType;
//    const auto& diffusion_tensor = *(this->problem().diffusion_tensor());
//    assert(!diffusion_tensor.parametric());
//    assert(diffusion_tensor.has_affine_part());
//    typedef GDT::LocalOperator::Codim1BoundaryIntegral<GDT::LocalEvaluation::SWIPDG::BoundaryLHS
//        <DiffusionFactorType, DiffusionTensorType>> DirichletOperatorType;
//    typedef GDT::LocalAssembler::Codim1BoundaryMatrix<DirichletOperatorType> DirichletMatrixAssemblerType;
//    std::vector<std::unique_ptr<DirichletOperatorType>> dirichlet_operators;
//    std::vector<std::unique_ptr<DirichletMatrixAssemblerType>> dirichlet_matrix_assemblers;
//    for (ssize_t qq = 0; qq < this->problem().diffusion_factor()->num_components(); ++qq) {
//      dirichlet_operators.emplace_back(new DirichletOperatorType(*(this->problem().diffusion_factor()->component(qq)),
//                                                                 *(diffusion_tensor.affine_part())));
//      dirichlet_matrix_assemblers.emplace_back(new DirichletMatrixAssemblerType(*dirichlet_operators.back()));
//      auto comp = local_matrix.register_component(*this->problem().diffusion_factor()->coefficient(qq),
//                                                  local_test_space.mapper().size(),
//                                                  local_ansatz_space.mapper().size(),
//                                                  local_pattern);
//      boundary_assembler.add(*dirichlet_matrix_assemblers.back(),
//                             *(local_matrix.component(comp)),
//                             new Stuff::Grid::ApplyOn::DirichletIntersections<GridViewType>(this->boundary_info()));
//    }
//    if (this->problem().diffusion_factor()->has_affine_part()) {
//      dirichlet_operators.emplace_back(new DirichletOperatorType(*(this->problem().diffusion_factor()->affine_part()),
//                                                                 *(diffusion_tensor.affine_part())));
//      dirichlet_matrix_assemblers.emplace_back(new DirichletMatrixAssemblerType(*dirichlet_operators.back()));
//      if (!local_matrix.has_affine_part())
//        local_matrix.register_affine_part(new MatrixType(local_test_space.mapper().size(),
//                                                         local_ansatz_space.mapper().size(),
//                                                         local_pattern));
//      boundary_assembler.add(*dirichlet_matrix_assemblers.back(),
//                             *(local_matrix.affine_part()),
//                             new Stuff::Grid::ApplyOn::DirichletIntersections<GridViewType>(this->boundary_info()));
//    }

//    // rhs
//    // * neumann boundary terms
//    typedef typename ProblemType::FunctionType::NonparametricType NeumannType;
//    typedef GDT::LocalFunctional::Codim1Integral<GDT::LocalEvaluation::Product<NeumannType>> NeumannFunctionalType;
//    typedef GDT::LocalAssembler::Codim1Vector<NeumannFunctionalType> NeumannVectorAssemblerType;
//    std::vector<std::unique_ptr<NeumannFunctionalType>> neumann_functionals;
//    std::vector<std::unique_ptr<NeumannVectorAssemblerType>> neumann_vector_assemblers;
//    for (ssize_t qq = 0; qq < this->problem().neumann()->num_components(); ++qq) {
//      neumann_functionals.emplace_back(new NeumannFunctionalType(*(this->problem().neumann()->component(qq))));
//      neumann_vector_assemblers.emplace_back(new NeumannVectorAssemblerType(*(neumann_functionals[qq])));
//      auto comp = local_vector.register_component(*this->problem().neumann()->coefficient(qq),
//                                                  local_test_space.mapper().size());
//      boundary_assembler.add(*(neumann_vector_assemblers[qq]),
//                             *(local_vector.component(comp)),
//                             new Stuff::Grid::ApplyOn::NeumannIntersections<GridViewType>(this->boundary_info()));
//    }
//    if (this->problem().neumann()->has_affine_part()) {
//      neumann_functionals.emplace_back(new NeumannFunctionalType(*(this->problem().neumann()->affine_part())));
//      neumann_vector_assemblers.emplace_back(new NeumannVectorAssemblerType(*(neumann_functionals.back())));
//      if (!local_vector.has_affine_part())
//        local_vector.register_affine_part(new VectorType(local_test_space.mapper().size()));
//      boundary_assembler.add(*(neumann_vector_assemblers.back()),
//                             *(local_vector.affine_part()),
//                             new Stuff::Grid::ApplyOn::NeumannIntersections<GridViewType>(this->boundary_info()));
//    }

//    // * dirichlet boundary terms
//    typedef typename ProblemType::FunctionType::NonparametricType DirichletType;
//    typedef GDT::LocalFunctional::Codim1Integral< GDT::LocalEvaluation::SWIPDG::BoundaryRHS
//        <DiffusionFactorType, DirichletType, DiffusionTensorType>> DirichletFunctionalType;
//    typedef GDT::LocalAssembler::Codim1Vector<DirichletFunctionalType> DirichletVectorAssemblerType;
//    std::vector<std::unique_ptr<DirichletFunctionalType>> dirichlet_functionals;
//    std::vector<std::unique_ptr<DirichletVectorAssemblerType>> dirichlet_vector_assemblers;
//    if (this->problem().diffusion_factor()->has_affine_part()) {
//      for (ssize_t qq = 0; qq < this->problem().dirichlet()->num_components(); ++qq) {
//        dirichlet_functionals.emplace_back(new DirichletFunctionalType(*(this->problem().diffusion_factor()->affine_part()),
//                                                                       *(diffusion_tensor.affine_part()),
//                                                                       *(this->problem().dirichlet()->component(qq))));
//        dirichlet_vector_assemblers.emplace_back(new DirichletVectorAssemblerType(*(dirichlet_functionals.back())));
//        auto comp = local_vector.register_component(*this->problem().dirichlet()->coefficient(qq),
//                                                    local_test_space.mapper().size());
//        boundary_assembler.add(*(dirichlet_vector_assemblers.back()),
//                               *(local_vector.component(comp)),
//                               new Stuff::Grid::ApplyOn::DirichletIntersections<GridViewType>(this->boundary_info()));
//      }
//    }
//    if (this->problem().dirichlet()->has_affine_part()) {
//      for (ssize_t qq = 0; qq < this->problem().diffusion_factor()->num_components(); ++qq) {
//        dirichlet_functionals.emplace_back(new DirichletFunctionalType(*(this->problem().diffusion_factor()->component(qq)),
//                                                                       *(diffusion_tensor.affine_part()),
//                                                                       *(this->problem().dirichlet()->affine_part())));
//        dirichlet_vector_assemblers.emplace_back(new DirichletVectorAssemblerType(*(
//            dirichlet_functionals[dirichlet_functionals.size() - 1])));
//        auto comp = local_vector.register_component(*this->problem().diffusion_factor()->coefficient(qq),
//                                                    local_test_space.mapper().size());
//        boundary_assembler.add(*(dirichlet_vector_assemblers.back()),
//                               *(local_vector.component(comp)),
//                               new Stuff::Grid::ApplyOn::DirichletIntersections<GridViewType>(this->boundary_info()));
//      }
//    }
//    Pymor::ParameterType param;
//    for (auto key : this->problem().diffusion_factor()->parameter_type().keys())
//      param.set(key, this->problem().diffusion_factor()->parameter_type().get(key));
//    for (auto key : this->problem().dirichlet()->parameter_type().keys())
//      param.set(key, this->problem().dirichlet()->parameter_type().get(key));
//    for (ssize_t pp = 0; pp < this->problem().diffusion_factor()->num_components(); ++ pp) {
//      for (ssize_t qq = 0; qq < this->problem().dirichlet()->num_components(); ++qq) {
//        dirichlet_functionals.emplace_back(new DirichletFunctionalType(*(this->problem().diffusion_factor()->component(pp)),
//                                                                       *(diffusion_tensor.affine_part()),
//                                                                       *(this->problem().dirichlet()->component(qq))));
//        dirichlet_vector_assemblers.emplace_back(new DirichletVectorAssemblerType(*(
//            dirichlet_functionals[dirichlet_functionals.size() - 1])));
//        const std::string expression = "(" + this->problem().diffusion_factor()->coefficient(pp)->expression()
//                                       + ")*(" + this->problem().dirichlet()->coefficient(qq)->expression() + ")";
//        auto comp = local_vector.register_component(new Pymor::ParameterFunctional(param, expression),
//                                                    local_test_space.mapper().size());
//        boundary_assembler.add(*(dirichlet_vector_assemblers.back()),
//                               *(local_vector.component(comp)),
//                               new Stuff::Grid::ApplyOn::DirichletIntersections<GridViewType>(this->boundary_info()));
//      }
//    }
//    if (this->problem().dirichlet()->has_affine_part() && this->problem().diffusion_factor()->has_affine_part()) {
//      if (!local_vector.has_affine_part())
//        local_vector.register_affine_part(new VectorType(local_test_space.mapper().size()));
//      dirichlet_functionals.emplace_back(new DirichletFunctionalType(*(this->problem().diffusion_factor()->affine_part()),
//                                                                     *(diffusion_tensor.affine_part()),
//                                                                     *(this->problem().dirichlet()->affine_part())));
//      dirichlet_vector_assemblers.emplace_back(new DirichletVectorAssemblerType(*(dirichlet_functionals.back())));
//      boundary_assembler.add(*(dirichlet_vector_assemblers.back()),
//                             *(local_vector.affine_part()),
//                             new Stuff::Grid::ApplyOn::DirichletIntersections<GridViewType>(this->boundary_info()));
//    } // dirichlet boundary terms

//    // do the actual work
//    boundary_assembler.assemble();
//  } // ... assemble_boundary_contributions(...)

//  template <class MacroEntityType>
//  void assemble_coupling_contributions(const MacroEntityType& macro_entity,
//                                       const MacroEntityType& macro_neighbor,
//                                       AffinelyDecomposedMatrixType& inside_inside_matrix,
//                                       AffinelyDecomposedMatrixType& inside_outside_matrix,
//                                       AffinelyDecomposedMatrixType& outside_inside_matrix,
//                                       AffinelyDecomposedMatrixType& outside_outside_matrix) const
//  {
//    const size_t subdomain = glued_grid_.macro_grid_view().indexSet().index(macro_entity);
//    const size_t neighbour = glued_grid_.macro_grid_view().indexSet().index(macro_neighbor);
//    typedef typename LocalDiscretizationType::TestSpaceType   LocalTestSpaceType;
//    typedef typename LocalDiscretizationType::AnsatzSpaceType LocalAnsatzSpaceType;
//    const LocalTestSpaceType&   inner_test_space   = this->local_discretizations_[subdomain]->test_space();
//    const LocalAnsatzSpaceType& inner_ansatz_space = this->local_discretizations_[subdomain]->ansatz_space();
//    const LocalTestSpaceType&   outer_test_space   = this->local_discretizations_[neighbour]->test_space();
//    const LocalAnsatzSpaceType& outer_ansatz_space = this->local_discretizations_[neighbour]->ansatz_space();
//    auto coupling_assembler = make_coupling_assembler(inner_test_space, inner_ansatz_space,
//                                                      outer_test_space, outer_ansatz_space,
//                                                      glued_grid_.coupling(macro_entity,
//                                                                           glued_grid_.max_local_level(macro_entity),
//                                                                           macro_neighbor,
//                                                                           glued_grid_.max_local_level(macro_neighbor)));

//    typedef typename ProblemType::DiffusionFactorType::NonparametricType DiffusionFactorType;
//    typedef typename ProblemType::DiffusionTensorType::NonparametricType DiffusionTensorType;
//    const auto& diffusion_tensor = *(this->problem().diffusion_tensor());
//    assert(!diffusion_tensor.parametric());
//    assert(diffusion_tensor.has_affine_part());
//    typedef GDT::LocalOperator::Codim1CouplingIntegral< GDT::LocalEvaluation::SWIPDG::Inner< DiffusionFactorType,
//                                                                                             DiffusionTensorType > >
//        CouplingOperatorType;
//    typedef GDT::LocalAssembler::Codim1CouplingMatrix< CouplingOperatorType > CouplingMatrixAssemblerType;
//    std::vector< std::unique_ptr< CouplingOperatorType > > coupling_operators;
//    std::vector< std::unique_ptr< CouplingMatrixAssemblerType > > coupling_matrix_assemblers;
//    for (ssize_t qq = 0; qq < this->problem().diffusion_factor()->num_components(); ++qq) {
//      coupling_operators.emplace_back(new CouplingOperatorType(*(this->problem().diffusion_factor()->component(qq)),
//                                                               *(diffusion_tensor.affine_part())));
//      coupling_matrix_assemblers.emplace_back(new CouplingMatrixAssemblerType(*(coupling_operators[qq])));
//      coupling_assembler.addLocalAssembler(*(coupling_matrix_assemblers[qq]),
//                                           *(inside_inside_matrix.component(qq)),
//                                           *(inside_outside_matrix.component(qq)),
//                                           *(outside_inside_matrix.component(qq)),
//                                           *(outside_outside_matrix.component(qq)));
//    }
//    if (this->problem().diffusion_factor()->has_affine_part()) {
//      coupling_operators.emplace_back(new CouplingOperatorType(*(this->problem().diffusion_factor()->affine_part()),
//                                                               *(diffusion_tensor.affine_part())));
//      coupling_matrix_assemblers.emplace_back(new CouplingMatrixAssemblerType(*(
//          coupling_operators[coupling_operators.size() - 1])));
//      coupling_assembler.addLocalAssembler(*(coupling_matrix_assemblers[coupling_matrix_assemblers.size() - 1]),
//                                           *(inside_inside_matrix.affine_part()),
//                                           *(inside_outside_matrix.affine_part()),
//                                           *(outside_inside_matrix.affine_part()),
//                                           *(outside_outside_matrix.affine_part()));
//    }

//    // do the actual work
//    coupling_assembler.assemble();
//  } // ... assemble_coupling_contributions(...)

//  void build_products()
//  {
//    for (const std::string& product_id : this->local_discretizations_.at(0)->available_products()) {
//      auto product_matrix = std::make_shared<AffinelyDecomposedMatrixType>();
//      for (size_t ss = 0; ss < glued_grid_.num_subdomains(); ++ss)
//        copy_local_to_global_matrix(this->local_discretizations_[ss]->get_product(product_id).container(),
//                                    this->local_discretizations_[ss]->pattern(),
//                                    ss,
//                                    ss,
//                                    *product_matrix);
//      this->products_.insert(std::make_pair(product_id, product_matrix));
//    }
//  } // ... build_products(...)

//  void build_global_containers()
//  {
//    // walk the subdomains
//    for (auto&& macro_entity : DSC::entityRange(glued_grid_.macro_grid_view())) {
//      const size_t ss = glued_grid_.macro_grid_view().indexSet().index(macro_entity);

//      copy_local_to_global_matrix(*this->local_discretizations_[ss]->system_matrix(),
//                                  this->local_discretizations_[ss]->pattern(),
//                                  ss,
//                                  ss,
//                                  *(this->matrix_));
//      copy_local_to_global_vector(*this->local_discretizations_[ss]->rhs(),
//                                  ss,
//                                  *(this->rhs_));

//      copy_local_to_global_matrix(*this->local_boundary_matrices_[ss],
//                                  this->local_discretizations_[ss]->pattern(),
//                                  ss,
//                                  ss,
//                                  *(this->matrix_));
//      copy_local_to_global_vector(*this->local_boundary_vectors_[ss],
//                                  ss,
//                                  *(this->rhs_));

//      // walk the neighbours
//      for (auto&& intersection : DSC::intersectionRange(glued_grid_.macro_grid_view(), macro_entity)) {
//        if (intersection.neighbor() && !intersection.boundary()) {
//          const auto macro_neighbor_ptr = intersection.outside();
//          const auto& macro_neighbor = *macro_neighbor_ptr;
//          const size_t nn = glued_grid_.macro_grid_view().indexSet().index(macro_neighbor);
//          if (ss < nn) {
//            // get the coupling patterns
//            const auto result_inside_outside_pattern = inside_outside_patterns_[ss].find(nn);
//            if (result_inside_outside_pattern == inside_outside_patterns_[ss].end())
//              DUNE_THROW(Stuff::Exceptions::internal_error,
//                         "The coupling pattern for subdomain " << ss << " and neighbour " << nn << "is missing!");
//            const auto& inside_outside_pattern = *(result_inside_outside_pattern->second);
//            const auto result_outside_inside_pattern = outside_inside_patterns_[nn].find(ss);
//            if (result_outside_inside_pattern == outside_inside_patterns_[nn].end())
//              DUNE_THROW(Stuff::Exceptions::internal_error,
//                         "The coupling pattern for neighbour " << nn << " and subdomain " << ss << "is missing!");
//            const auto& outside_inside_pattern = *(result_outside_inside_pattern->second);
//            // and the coupling matrices
//            auto result_inside_outside_matrix = inside_outside_matrices_[ss].find(nn);
//            if (result_inside_outside_matrix == inside_outside_matrices_[ss].end())
//              DUNE_THROW(Stuff::Exceptions::internal_error,
//                         "The coupling matrix for subdomain " << ss << " and neighbour " << nn << "is missing!");
//            auto& inside_outside_matrix = *(result_inside_outside_matrix->second);
//            auto result_outside_inside_matrix = outside_inside_matrices_[nn].find(ss);
//            if (result_outside_inside_matrix == outside_inside_matrices_[nn].end())
//              DUNE_THROW(Stuff::Exceptions::internal_error,
//                         "The coupling matrix for neighbour " << nn << " and subdomain " << ss << "is missing!");
//            auto& outside_inside_matrix = *(result_outside_inside_matrix->second);
//            // and copy them into the global matrix
//            copy_local_to_global_matrix(inside_outside_matrix,
//                                        inside_outside_pattern,
//                                        ss, nn,
//                                        *(this->matrix_));
//            copy_local_to_global_matrix(outside_inside_matrix,
//                                        outside_inside_pattern,
//                                        nn, ss,
//                                        *(this->matrix_));
//          }
//        }
//      } // walk the neighbours
//    } // walk the subdomains
//  } // ... build_global_containers(...)

//  template< class AffinelyDecomposedContainerType >
//  ssize_t find_component(const AffinelyDecomposedContainerType& container,
//                         const Pymor::ParameterFunctional& coefficient) const
//  {
//    for (size_t qq = 0; qq < boost::numeric_cast< size_t >(container.num_components()); ++qq)
//      if (*(container.coefficient(qq)) == coefficient)
//        return boost::numeric_cast<ssize_t>(qq);
//    return -1;
//  } // ... find_component(...)

//  using BaseType::pattern_;

//  grid::Multiscale::Glued<MacroGridType, LocalGridType>& glued_grid_;
//  const std::vector< std::string > only_these_products_;
//  std::vector< std::shared_ptr< AffinelyDecomposedMatrixType > > local_boundary_matrices_;
//  std::vector< std::shared_ptr< AffinelyDecomposedVectorType > > local_boundary_vectors_;
//  std::vector< std::map< size_t, std::shared_ptr< PatternType > > > inside_outside_patterns_;
//  std::vector< std::map< size_t, std::shared_ptr< PatternType > > > outside_inside_patterns_;
//  std::vector< std::map< size_t, std::shared_ptr< AffinelyDecomposedMatrixType > > > inside_outside_matrices_;
//  std::vector< std::map< size_t, std::shared_ptr< AffinelyDecomposedMatrixType > > > outside_inside_matrices_;
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
