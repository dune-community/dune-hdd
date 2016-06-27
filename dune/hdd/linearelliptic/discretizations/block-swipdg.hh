// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_BLOCK_SWIPDG_HH
#define DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_BLOCK_SWIPDG_HH

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

# include <dune/common/timer.hh>
# include <dune/common/dynmatrix.hh>

# if HAVE_ALUGRID
#   include <dune/grid/alugrid.hh>
# endif

# include <dune/geometry/quadraturerules.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/grid/multiscale/provider.hh>

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
#include <dune/gdt/playground/spaces/block.hh>
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

#include "base.hh"
#include "swipdg.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Discretizations {


// forward, needed in the Traits
template< class GridImp, class RangeFieldImp, int rangeDim, int polynomialOrder = 1
        , Stuff::LA::ChooseBackend la_backend = Stuff::LA::default_sparse_backend >
class BlockSWIPDG;


namespace internal {


template< class GridType, class RangeFieldType, int dimRange, int polOrder, Stuff::LA::ChooseBackend la_backend >
class LocalDiscretizationsContainer
{
  typedef grid::Multiscale::ProviderInterface< GridType > GridProviderType;
public:
  typedef SWIPDG< GridType, Stuff::Grid::ChooseLayer::local, RangeFieldType, dimRange
                , polOrder, GDT::ChooseSpaceBackend::fem, la_backend > DiscretizationType;
  typedef SWIPDG< GridType, Stuff::Grid::ChooseLayer::local_oversampled, RangeFieldType, dimRange
                , polOrder, GDT::ChooseSpaceBackend::fem, la_backend > OversampledDiscretizationType;
  typedef typename DiscretizationType::ProblemType     ProblemType;
  typedef typename DiscretizationType::TestSpaceType   TestSpaceType;
  typedef typename DiscretizationType::AnsatzSpaceType AnsatzSpaceType;

private:
  typedef Problems::ZeroBoundary< ProblemType > FakeProblemType;
  typedef typename DiscretizationType::GridViewType::Intersection IntersectionType;

public:
  LocalDiscretizationsContainer(const GridProviderType& grid_provider,
                                const ProblemType& prob,
                                const std::vector< std::string >& only_these_products)
    : zero_boundary_problem_(prob)
    , all_dirichlet_boundary_config_(Stuff::Grid::BoundaryInfos::AllDirichlet< IntersectionType >::default_config())
    , all_neumann_boundary_config_(Stuff::Grid::BoundaryInfos::AllNeumann< IntersectionType >::default_config())
    , multiscale_boundary_config_(Stuff::Grid::BoundaryInfoConfigs::IdBased::default_config())
    , local_discretizations_(grid_provider.num_subdomains(), nullptr)
    , oversampled_discretizations_dirichlet_(grid_provider.num_subdomains(), nullptr)
    , oversampled_discretizations_neumann_(grid_provider.num_subdomains(), nullptr)
    , local_test_spaces_(grid_provider.num_subdomains(), nullptr)
    , local_ansatz_spaces_(grid_provider.num_subdomains(), nullptr)
  {
    multiscale_boundary_config_["neumann"] = "7";
    for (size_t ss = 0; ss < grid_provider.num_subdomains(); ++ss) {
      local_discretizations_[ss] = std::make_shared< DiscretizationType >(grid_provider,
                                                                          all_neumann_boundary_config_,
                                                                          zero_boundary_problem_,
                                                                          ss,
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
  typedef BlockSWIPDG< GridImp, RangeFieldImp, rangeDim, polynomialOrder, la_backend >        ThisType;
public:
  typedef internal::BlockSWIPDGTraits< GridImp, RangeFieldImp, rangeDim, polynomialOrder, la_backend > Traits;
  typedef typename BaseType::ProblemType     ProblemType;
  typedef typename BaseType::GridViewType    GridViewType;
  typedef typename BaseType::TestSpaceType   TestSpaceType;
  typedef typename BaseType::AnsatzSpaceType AnsatzSpaceType;
  typedef typename BaseType::EntityType      EntityType;
  typedef typename BaseType::DomainFieldType DomainFieldType;
  typedef typename BaseType::RangeFieldType  RangeFieldType;
  typedef typename BaseType::MatrixType      MatrixType;
  typedef typename BaseType::VectorType      VectorType;
  typedef typename BaseType::OperatorType    OperatorType;
  typedef typename BaseType::ProductType     ProductType;
  typedef typename BaseType::FunctionalType  FunctionalType;

  static const unsigned int dimDomain = BaseType::dimDomain;
  static const unsigned int dimRange  = BaseType::dimRange;

  typedef grid::Multiscale::ProviderInterface< GridImp > GridProviderType;
  typedef typename GridProviderType::GridType   GridType;
  typedef typename GridProviderType::MsGridType MsGridType;

  typedef typename Traits::LocalDiscretizationsContainerType::DiscretizationType            LocalDiscretizationType;
  typedef typename Traits::LocalDiscretizationsContainerType::OversampledDiscretizationType OversampledDiscretizationType;
  typedef typename TestSpaceType::PatternType PatternType;

private:
  using typename BaseType::AffinelyDecomposedMatrixType;
  using typename BaseType::AffinelyDecomposedVectorType;
  typedef Pymor::LA::AffinelyDecomposedConstContainer< MatrixType > AffinelyDecomposedConstMatrixType;
  typedef Pymor::LA::AffinelyDecomposedConstContainer< VectorType > AffinelyDecomposedConstVectorType;

public:
  typedef typename LocalDiscretizationType::ProblemType LocalProblemType;
  typedef typename LocalDiscretizationType::ProductType LocalProductType;

  static std::string static_id();

  BlockSWIPDG(const GridProviderType& grid_provider,
              const Stuff::Common::Configuration& bound_inf_cfg,
              const ProblemType& prob,
              const std::vector< std::string >& only_these_products = {});

  const std::vector< std::shared_ptr< LocalDiscretizationType > >& local_discretizations() const;

  void init(const bool prune = false);

  ssize_t num_subdomains() const;

  std::vector< ssize_t > neighbouring_subdomains(const ssize_t ss) const;

  VectorType localize_vector(const VectorType& global_vector, const size_t ss) const;

  VectorType globalize_vectors(const std::vector< VectorType >& local_vectors) const;

  VectorType* globalize_vectors_and_return_ptr(const std::vector< VectorType >& local_vectors) const;

  VectorType* localize_vector_and_return_ptr(const VectorType& global_vector, const ssize_t ss) const;

  ProductType get_local_product(const size_t ss, const std::string id) const;

  ProductType* get_local_product_and_return_ptr(const ssize_t ss, const std::string id) const;

  OperatorType get_local_operator(const size_t ss) const;

  OperatorType* get_local_operator_and_return_ptr(const ssize_t ss) const;

  OperatorType get_coupling_operator(const size_t ss, const size_t nn) const;

  OperatorType* get_coupling_operator_and_return_ptr(const ssize_t ss, const ssize_t nn) const;

  FunctionalType get_local_functional(const size_t ss) const;

  FunctionalType* get_local_functional_and_return_ptr(const ssize_t ss) const;

  VectorType solve_for_local_correction(const std::vector< VectorType >& local_vectors,
                                        const size_t subdomain,
                                        const Pymor::Parameter mu = Pymor::Parameter()) const;

  LocalDiscretizationType get_local_discretization(const size_t subdomain) const;

  LocalDiscretizationType* pb_get_local_discretization(const ssize_t subdomain) const;

  OversampledDiscretizationType get_oversampled_discretization(const size_t subdomain,
                                                               const std::string boundary_value_type) const;

  OversampledDiscretizationType* pb_get_oversampled_discretization(const ssize_t subdomain,
                                                                   const std::string boundary_value_type) const;

//  OversampledDiscretizationType* pb_get_oversampled_discretization(const ssize_t subdomain,
//                                                                   const std::string boundary_value_type,
//                                                                   const VectorType& boundary_values) const;

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
            // so no further check neccesarry then
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
    const CouplingGridPartType grid_part_;
    std::vector< LocalCodim1MatrixAssemblerApplication* > localCodim1MatrixAssemblers_;
  }; // class CouplingAssembler

  void add_local_to_global_pattern(const PatternType& local,
                                   const size_t test_subdomain,
                                   const size_t ansatz_subdomain,
                                   PatternType& global) const;

  void copy_local_to_global_matrix(const AffinelyDecomposedConstMatrixType& local_matrix,
                                   const PatternType& local_pattern,
                                   const size_t subdomain,
                                   const size_t neighbor,
                                   AffinelyDecomposedMatrixType& global_matrix) const;

  template< class ML, class MG >
  void copy_local_to_global_matrix(const Stuff::LA::MatrixInterface< ML >& local_matrix,
                                   const PatternType& local_pattern,
                                   const size_t test_subdomain,
                                   const size_t ansatz_subdomain,
                                   Stuff::LA::MatrixInterface< MG >& global_matrix) const
  {
    for (size_t local_ii = 0; local_ii < local_pattern.size(); ++local_ii) {
      const size_t global_ii = this->test_space().mapper().mapToGlobal(test_subdomain, local_ii);
      for (const size_t& local_jj : local_pattern.inner(local_ii)) {
        const size_t global_jj = this->ansatz_space().mapper().mapToGlobal(ansatz_subdomain, local_jj);
        global_matrix.add_to_entry(global_ii, global_jj, local_matrix.get_entry(local_ii, local_jj));
      }
    }
  } // ... copy_local_to_global_matrix(...)

  void copy_local_to_global_vector(const AffinelyDecomposedConstVectorType& local_vector,
                                   const size_t subdomain,
                                   AffinelyDecomposedVectorType& global_vector) const;

  template< class VL, class VG >
  void copy_local_to_global_vector(const Stuff::LA::VectorInterface< VL >& local_vector,
                                   const size_t subdomain,
                                   Stuff::LA::VectorInterface< VG >& global_vector) const
  {
    for (size_t local_ii = 0; local_ii < local_vector.size(); ++local_ii) {
      const size_t global_ii = this->test_space().mapper().mapToGlobal(subdomain, local_ii);
      global_vector.add_to_entry(global_ii, local_vector.get_entry(local_ii));
    }
  } // ... copy_local_to_global_vector(...)

  void assemble_boundary_contributions(const size_t subdomain) const;

  /**
   * \note  We take the matrices as input here becaus we would have to look them up in the maps otherwise. Since that
   *        has already been done above we save a little.
   */
  void assemble_coupling_contributions(const size_t subdomain,
                                       const size_t neighbour,
                                       AffinelyDecomposedMatrixType& inside_inside_matrix,
                                       AffinelyDecomposedMatrixType& inside_outside_matrix,
                                       AffinelyDecomposedMatrixType& outside_inside_matrix,
                                       AffinelyDecomposedMatrixType& outside_outside_matrix) const;

  void build_global_containers();

  template <class T, class A>
  ssize_t find_add_component(AffinelyDecomposedMatrixType& matrix,
                             const Pymor::ParameterFunctional& coefficient,
                             const T& t_space,
                             const A& a_space,
                             const PatternType& ptrn) const
  {
    for (size_t qq = 0; qq < boost::numeric_cast<size_t>(matrix.num_components()); ++qq)
      if (*(matrix.coefficient(qq)) == coefficient)
        return qq;
    return matrix.register_component(coefficient, t_space.mapper().size(), a_space.mapper().size(), ptrn);
  } // ... find_add_component(...)

  template <class S>
  ssize_t find_add_component(AffinelyDecomposedVectorType& vector,
                             const Pymor::ParameterFunctional& coefficient,
                             const S& space) const
  {
    for (size_t qq = 0; qq < boost::numeric_cast<size_t>(vector.num_components()); ++qq)
      if (*(vector.coefficient(qq)) == coefficient)
        return qq;
    return vector.register_component(coefficient, space.mapper().size());
  } // ... find_add_component(...)

  const GridProviderType& grid_provider_;
  std::shared_ptr< const MsGridType > ms_grid_;
  const std::vector< std::string > only_these_products_;
  using BaseType::pattern_;
  std::vector< std::shared_ptr< AffinelyDecomposedMatrixType > > local_matrices_;
  std::vector< std::shared_ptr< AffinelyDecomposedVectorType > > local_vectors_;
  std::vector< std::map< size_t, std::shared_ptr< PatternType > > > inside_outside_patterns_;
  std::vector< std::map< size_t, std::shared_ptr< PatternType > > > outside_inside_patterns_;
  std::vector< std::map< size_t, std::shared_ptr< AffinelyDecomposedMatrixType > > > inside_outside_matrices_;
  std::vector< std::map< size_t, std::shared_ptr< AffinelyDecomposedMatrixType > > > outside_inside_matrices_;
}; // BlockSWIPDG


#if HAVE_ALUGRID && HAVE_DUNE_FEM
# if HAVE_DUNE_ISTL

extern template class BlockSWIPDG< ALUGrid< 2, 2, simplex, conforming, No_Comm >,
                                   double,
                                   1,
                                   1,
                                   Stuff::LA::ChooseBackend::istl_sparse >;

#   if HAVE_MPI

extern template class BlockSWIPDG< ALUGrid< 2, 2, simplex, conforming, MPI_Comm >,
                                   double,
                                   1,
                                   1,
                                   Stuff::LA::ChooseBackend::istl_sparse >;

#   endif // HAVE_MPI
# endif // HAVE_DUNE_ISTL
# if HAVE_EIGEN

extern template class BlockSWIPDG< ALUGrid< 2, 2, simplex, conforming, No_Comm >,
                                   double,
                                   1,
                                   1,
                                   Stuff::LA::ChooseBackend::eigen_sparse >;

#   if HAVE_MPI

extern template class BlockSWIPDG< ALUGrid< 2, 2, simplex, conforming, MPI_Comm >,
                                   double,
                                   1,
                                   1,
                                   Stuff::LA::ChooseBackend::eigen_sparse >;

#   endif // HAVE_MPI
# endif // HAVE_EIGEN
#endif // HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_ISTL


template< Stuff::LA::ChooseBackend la, class G, class R, int r = 1, int p = 1 >
BlockSWIPDG< G, R, r, p, la > make_block_swipdg(const grid::Multiscale::ProviderInterface< G >& grid_provider,
                                                const DSC::Configuration& boundary_info,
                                                const ProblemInterface< typename G::template Codim< 0 >::Entity,
                                                                        typename G::ctype, G::dimension,
                                                                        R, r >& problem,
                                                const std::vector< std::string >& only_these_products = {})
{
  return BlockSWIPDG< G, R, r, p, la >(grid_provider, boundary_info, problem, only_these_products);
}


} // namespace Discretizations
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_BLOCK_SWIPDG_HH
