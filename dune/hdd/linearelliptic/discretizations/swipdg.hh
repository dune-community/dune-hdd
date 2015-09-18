// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_SWIPDG_HH
#define DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_SWIPDG_HH

#include <memory>
#include <vector>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/timer.hh>

#include <dune/stuff/common/disable_warnings.hh>
# if HAVE_ALUGRID
#   include <dune/grid/alugrid.hh>
# endif
#include <dune/stuff/common/reenable_warnings.hh>

#if HAVE_DUNE_GRID_MULTISCALE
# include <dune/grid/multiscale/provider/interface.hh>
#endif

#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/grid/layers.hh>
#include <dune/stuff/grid/provider.hh>
#include <dune/stuff/la/container.hh>
#include <dune/stuff/la/solver.hh>
#include <dune/stuff/grid/walker/functors.hh>

#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/assembler/tmp-storage.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/functionals/l2.hh>
#include <dune/gdt/operators/oswaldinterpolation.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/playground/functionals/swipdg.hh>
#include <dune/gdt/playground/operators/elliptic-swipdg.hh>
#include <dune/gdt/playground/operators/fluxreconstruction.hh>
#include <dune/gdt/playground/products/swipdgpenalty.hh>
#include <dune/gdt/spaces/fv/default.hh>
#include <dune/gdt/spaces/rt/pdelab.hh>
#include <dune/gdt/products/boundaryl2.hh>
#include <dune/gdt/products/elliptic.hh>
#include <dune/gdt/products/h1.hh>
#include <dune/gdt/products/l2.hh>
#include <dune/gdt/spaces/dg.hh>

#include "base.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Discretizations {


// forward, for friendlyness
template< class GridImp, class RangeFieldImp, int rangeDim, int polynomialOrder, Stuff::LA::ChooseBackend la_backend >
class BlockSWIPDG;


// forward, needed in the Traits
template< class GridImp, Stuff::Grid::ChooseLayer layer, class RangeFieldImp, int rangeDim, int polynomialOrder = 1,
#if HAVE_DUNE_FEM
          GDT::ChooseSpaceBackend space_backend = GDT::ChooseSpaceBackend::fem,
#else
# error No suitable space backend available!
#endif
          Stuff::LA::ChooseBackend la_backend  = Stuff::LA::default_sparse_backend >
class SWIPDG;


namespace internal {


template< class GridImp, Stuff::Grid::ChooseLayer layer, class RangeFieldImp, int rangeDim, int polynomialOrder,
          GDT::ChooseSpaceBackend space_backend,
          Stuff::LA::ChooseBackend la_backend >
class SWIPDGTraits
  : public internal::ContainerBasedDefaultTraits< typename Stuff::LA::Container< RangeFieldImp, la_backend >::MatrixType,
                                                  typename Stuff::LA::Container< RangeFieldImp, la_backend >::VectorType>
{
public:
  typedef SWIPDG< GridImp, layer, RangeFieldImp, rangeDim, polynomialOrder, space_backend, la_backend > derived_type;
  typedef GridImp GridType;
  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRange = rangeDim;
  static const unsigned int polOrder = polynomialOrder;

private:
  typedef GDT::Spaces::DGProvider< GridType, layer, space_backend, polOrder, RangeFieldType, dimRange > SpaceProvider;

  friend class SWIPDG< GridImp, layer, RangeFieldImp, rangeDim, polynomialOrder, space_backend, la_backend >;

public:
  typedef typename SpaceProvider::Type TestSpaceType;
  typedef TestSpaceType AnsatzSpaceType;
  typedef typename TestSpaceType::GridViewType GridViewType;
}; // class SWIPDGTraits


} // namespace internal


template< class GridImp, Stuff::Grid::ChooseLayer layer, class RangeFieldImp, int rangeDim, int polynomialOrder,
          GDT::ChooseSpaceBackend space_backend,
          Stuff::LA::ChooseBackend la_backend >
class SWIPDG
  : public ContainerBasedDefault< internal::SWIPDGTraits< GridImp, layer, RangeFieldImp, rangeDim, polynomialOrder,
                                                          space_backend, la_backend > >
{
  typedef ContainerBasedDefault< internal::SWIPDGTraits< GridImp, layer, RangeFieldImp, rangeDim, polynomialOrder,
                                                         space_backend, la_backend > > BaseType;
  typedef SWIPDG< GridImp, layer, RangeFieldImp, rangeDim, polynomialOrder, space_backend, la_backend > ThisType;
public:
  typedef internal::SWIPDGTraits< GridImp, layer, RangeFieldImp, rangeDim, polynomialOrder, space_backend, la_backend >
      Traits;
  typedef GridImp GridType;
  using typename BaseType::ProblemType;
  using typename BaseType::BoundaryInfoType;
  using typename BaseType::TestSpaceType;
  using typename BaseType::AnsatzSpaceType;
  using typename BaseType::MatrixType;
  using typename BaseType::VectorType;
  using typename BaseType::GridViewType;
  using typename BaseType::RangeFieldType;

  typedef typename TestSpaceType::PatternType PatternType;

  static const unsigned int dimDomain = BaseType::dimDomain;
  static const unsigned int dimRange = BaseType::dimRange;

private:
  typedef typename Traits::SpaceProvider SpaceProvider;

  typedef Stuff::Grid::ProviderInterface< GridType >      GridProviderType;
#if HAVE_DUNE_GRID_MULTISCALE
  typedef grid::Multiscale::ProviderInterface< GridType > MsGridProviderType;
#endif
  using typename BaseType::AffinelyDecomposedMatrixType;
  using typename BaseType::AffinelyDecomposedVectorType;

  typedef typename ProblemType::DiffusionFactorType::NonparametricType DiffusionFactorType;
  typedef typename ProblemType::DiffusionTensorType::NonparametricType DiffusionTensorType;
  typedef GDT::Operators::EllipticSWIPDG< DiffusionFactorType, MatrixType
                                        , TestSpaceType, AnsatzSpaceType
                                        , GridViewType, DiffusionTensorType > EllipticOperatorType;

public:
  static std::string static_id()
  {
    return DiscretizationInterface< Traits >::static_id() + ".swipdg";
  }

  SWIPDG(GridProviderType& grid_provider,
         const Stuff::Common::Configuration& bound_inf_cfg,
         const ProblemType& prob,
         const int level_or_subdomain = 0,
         const std::vector< std::string >& only_these_products = {})
    : BaseType(SpaceProvider::create(grid_provider, level_or_subdomain),
               SpaceProvider::create(grid_provider, level_or_subdomain),
               bound_inf_cfg,
               prob)
    , beta_(GDT::LocalEvaluation::SIPDG::internal::default_beta(dimDomain))
    , pattern_(EllipticOperatorType::pattern(BaseType::test_space(), BaseType::ansatz_space()))
    , only_these_products_(only_these_products)
  {
    // in case of parametric diffusion tensor this discretization is not affinely decomposable any more
    if (this->problem_.diffusion_tensor()->parametric())
      DUNE_THROW(NotImplemented, "The diffusion tensor must not be parametric!");
    if (!this->problem_.diffusion_tensor()->has_affine_part())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given, "The diffusion tensor must not be empty!");
  } // SWIPDG(...)

#if HAVE_DUNE_GRID_MULTISCALE
  SWIPDG(const MsGridProviderType& grid_provider,
         const Stuff::Common::Configuration& bound_inf_cfg,
         const ProblemType& prob,
         const int level_or_subdomain = 0,
         const std::vector< std::string >& only_these_products = {})
    : BaseType(SpaceProvider::create(grid_provider, level_or_subdomain),
               SpaceProvider::create(grid_provider, level_or_subdomain),
               bound_inf_cfg,
               prob)
    , beta_(GDT::LocalEvaluation::SIPDG::internal::default_beta(dimDomain))
    , pattern_(EllipticOperatorType::pattern(BaseType::test_space(), BaseType::ansatz_space()))
    , only_these_products_(only_these_products)
  {
    // in case of parametric diffusion tensor this discretization is not affinely decomposable any more
    if (this->problem_.diffusion_tensor()->parametric())
      DUNE_THROW(NotImplemented, "The diffusion tensor must not be parametric!");
    if (!this->problem_.diffusion_tensor()->has_affine_part())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given, "The diffusion tensor must not be empty!");
  } // SWIPDG(...)
#endif // HAVE_DUNE_GRID_MULTISCALE

  const PatternType& pattern() const
  {
    return pattern_;
  }

  void init(std::ostream& out = Stuff::Common::Logger().devnull(), const std::string prefix = "")
  {
    if (!this->container_based_initialized_) {
      using namespace GDT;

      auto& matrix = *(this->matrix_);
      auto& rhs = *(this->rhs_);
      const auto& space = this->test_space_;
      const auto& boundary_info = *(this->boundary_info_);

      out << prefix << "assembling... " << std::flush;
      Dune::Timer timer;
      SystemAssembler< TestSpaceType > system_assembler(space);
      Stuff::Grid::Functor::DirichletDetector< GridViewType > dirichlet_detector(this->boundary_info());
      system_assembler.add(dirichlet_detector);

      // lhs operator
      const auto& diffusion_factor = *(this->problem_.diffusion_factor());
      const auto& diffusion_tensor = *(this->problem_.diffusion_tensor());
      assert(!diffusion_tensor.parametric());
      assert(diffusion_tensor.has_affine_part());
      std::vector< std::unique_ptr< EllipticOperatorType > > elliptic_operators;
      for (size_t qq = 0; qq < boost::numeric_cast< size_t >(diffusion_factor.num_components()); ++qq) {
        const size_t id = matrix.register_component(diffusion_factor.coefficient(qq),
                                                    space.mapper().size(), space.mapper().size(), pattern_);
        elliptic_operators.emplace_back(new EllipticOperatorType(
            *(diffusion_factor.component(qq)),
            *(diffusion_tensor.affine_part()),
            boundary_info,
            *(matrix.component(id)),
            space));
      }
      if (diffusion_factor.has_affine_part()) {
        if (!matrix.has_affine_part())
          matrix.register_affine_part(space.mapper().size(), space.mapper().size(), pattern_);
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

      // products
      const size_t over_integrate = 2;
      // * L2
      typedef Products::L2Assemblable< MatrixType, TestSpaceType, GridViewType, AnsatzSpaceType > L2ProductType;
      std::unique_ptr< L2ProductType > l2_product;
      auto l2_product_matrix = std::make_shared< AffinelyDecomposedMatrixType >();
      if (std::find(only_these_products_.begin(), only_these_products_.end(), "l2") != only_these_products_.end()) {
        l2_product_matrix->register_affine_part(this->test_space().mapper().size(),
                                                this->ansatz_space().mapper().size(),
                                                L2ProductType::pattern(this->test_space(),
                                                                       this->ansatz_space()));
        l2_product = DSC::make_unique< L2ProductType >(*(l2_product_matrix->affine_part()),
                                                       this->test_space(),
                                                       this->grid_view(),
                                                       this->ansatz_space(),
                                                       over_integrate);
        system_assembler.add(*l2_product);
      }
      // * H1 semi
      typedef Products::H1SemiAssemblable< MatrixType, TestSpaceType, GridViewType, AnsatzSpaceType >
          H1ProductType;
      std::unique_ptr< H1ProductType > h1_product;
      auto h1_product_matrix = std::make_shared< AffinelyDecomposedMatrixType >();
      if (std::find(only_these_products_.begin(), only_these_products_.end(), "h1_semi") != only_these_products_.end()) {
        h1_product_matrix->register_affine_part(this->test_space().mapper().size(),
                                                this->ansatz_space().mapper().size(),
                                                H1ProductType::pattern(this->test_space(),
                                                                       this->ansatz_space()));
        h1_product = DSC::make_unique< H1ProductType >(*(h1_product_matrix->affine_part()),
                                                       this->test_space(),
                                                       this->grid_view(),
                                                       this->ansatz_space(),
                                                       over_integrate);
        system_assembler.add(*h1_product);
      }
      // * elliptic
      assert(!diffusion_tensor.parametric());
      typedef Products::EllipticAssemblable< MatrixType, DiffusionFactorType, TestSpaceType, GridViewType,
                                                  AnsatzSpaceType, RangeFieldType, DiffusionTensorType >
          EllipticProductType;
      std::vector< std::unique_ptr< EllipticProductType > > elliptic_products;
      auto elliptic_product_matrix = std::make_shared< AffinelyDecomposedMatrixType >();
      if (std::find(only_these_products_.begin(), only_these_products_.end(), "elliptic") != only_these_products_.end()) {
        for (DUNE_STUFF_SSIZE_T qq = 0; qq < diffusion_factor.num_components(); ++qq) {
          const auto id = elliptic_product_matrix->register_component(diffusion_factor.coefficient(qq),
                                                                      this->test_space().mapper().size(),
                                                                      this->ansatz_space().mapper().size(),
                                                                      EllipticProductType::pattern(this->test_space(),
                                                                                                   this->ansatz_space()));
          elliptic_products.emplace_back(new EllipticProductType(*elliptic_product_matrix->component(id),
                                                                 this->test_space(),
                                                                 this->grid_view(),
                                                                 this->ansatz_space(),
                                                                 *diffusion_factor.component(qq),
                                                                 *diffusion_tensor.affine_part(),
                                                                 over_integrate));
        }
        if (diffusion_factor.has_affine_part()) {
          elliptic_product_matrix->register_affine_part(this->test_space().mapper().size(),
                                                        this->ansatz_space().mapper().size(),
                                                        EllipticProductType::pattern(this->test_space(),
                                                                                     this->ansatz_space()));
          elliptic_products.emplace_back(new EllipticProductType(*elliptic_product_matrix->affine_part(),
                                                                 this->test_space(),
                                                                 this->grid_view(),
                                                                 this->ansatz_space(),
                                                                 *diffusion_factor.affine_part(),
                                                                 *diffusion_tensor.affine_part(),
                                                                 over_integrate));
        }
        for (auto& product : elliptic_products)
          system_assembler.add(*product);
      }
      // * boundary L2
      typedef Products::BoundaryL2Assemblable< MatrixType, TestSpaceType, GridViewType, AnsatzSpaceType >
          BoundaryL2ProductType;
      std::unique_ptr< BoundaryL2ProductType > boundary_l2_product;
      auto boundary_l2_product_matrix = std::make_shared< AffinelyDecomposedMatrixType >();
      if (std::find(only_these_products_.begin(), only_these_products_.end(), "boundary_l2") != only_these_products_.end()) {
        boundary_l2_product_matrix->register_affine_part(this->test_space().mapper().size(),
                                                         this->ansatz_space().mapper().size(),
                                                         BoundaryL2ProductType::pattern(this->test_space(),
                                                                                        this->ansatz_space()));
        boundary_l2_product = DSC::make_unique< BoundaryL2ProductType >(*(boundary_l2_product_matrix->affine_part()),
                                                                        this->test_space(),
                                                                        this->grid_view(),
                                                                        this->ansatz_space(),
                                                                        over_integrate);
        system_assembler.add(*boundary_l2_product);
      }
      // * penalty term
      typedef Products::SwipdgPenaltyAssemblable
          < MatrixType, DiffusionFactorType, DiffusionTensorType, TestSpaceType > PenaltyProductType;
      std::vector< std::unique_ptr< PenaltyProductType > > penalty_products;
      auto penalty_product_matrix = std::make_shared< AffinelyDecomposedMatrixType >();
      if (std::find(only_these_products_.begin(), only_these_products_.end(), "penalty") != only_these_products_.end()) {
        for (DUNE_STUFF_SSIZE_T qq = 0; qq < diffusion_factor.num_components(); ++qq) {
          const auto id = penalty_product_matrix->register_component(diffusion_factor.coefficient(qq),
                                                                     this->test_space().mapper().size(),
                                                                     this->ansatz_space().mapper().size(),
                                                                     PenaltyProductType::pattern(this->test_space(),
                                                                                                 this->ansatz_space()));
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
                                                       PenaltyProductType::pattern(this->test_space(),
                                                                                   this->ansatz_space()));
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
      out << "done (took " << timer.elapsed() << "s)" << std::endl;

      if (!dirichlet_detector.found())
        this->purely_neumann_ = true;

      // build parameter type
      this->inherit_parameter_type(matrix.parameter_type(), "lhs");
      this->inherit_parameter_type(rhs.parameter_type(),    "rhs");

      // finalize
      if (std::find(only_these_products_.begin(), only_these_products_.end(), "l2") != only_these_products_.end())
        this->products_.insert(std::make_pair("l2", l2_product_matrix));
      if (std::find(only_these_products_.begin(), only_these_products_.end(), "h1_semi") != only_these_products_.end())
        this->products_.insert(std::make_pair("h1_semi", h1_product_matrix));
      if (std::find(only_these_products_.begin(), only_these_products_.end(), "elliptic") != only_these_products_.end())
        this->products_.insert(std::make_pair("elliptic", elliptic_product_matrix));
      if (std::find(only_these_products_.begin(), only_these_products_.end(), "boundary_l2") != only_these_products_.end())
        this->products_.insert(std::make_pair("boundary_l2", boundary_l2_product_matrix));
      if (std::find(only_these_products_.begin(), only_these_products_.end(), "penalty") != only_these_products_.end())
        this->products_.insert(std::make_pair("penalty", penalty_product_matrix));
      if (std::find(only_these_products_.begin(), only_these_products_.end(), "energy") != only_these_products_.end())
        this->products_.insert(std::make_pair("energy",
                                              std::make_shared< AffinelyDecomposedMatrixType >(this->matrix_->copy())));

      this->finalize_init();
    } // if (!this->container_based_initialized_)
  } // ... init(...)

private:
  friend class BlockSWIPDG< GridImp, RangeFieldImp, rangeDim, polynomialOrder, la_backend >;

  const RangeFieldType beta_;
  PatternType pattern_;
  const std::vector< std::string > only_these_products_;
}; // class SWIPDG


} // namespace Discretizations
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_SWIPDG_HH
