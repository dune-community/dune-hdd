// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_CG_HH
#define DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_CG_HH

#include <memory>
#include <vector>
//#include <sstream>

//#include <dune/common/exceptions.hh>
//#include <dune/common/parametertree.hh>
#include <dune/common/timer.hh>

//#include <dune/grid/part/interface.hh>

//#include <dune/fem/misc/mpimanager.hh>

//#include <dune/stuff/common/logging.hh>
//#include <dune/stuff/discretefunction/projection/dirichlet.hh>
//#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/configtree.hh>
#include <dune/stuff/grid/partview.hh>
#include <dune/stuff/grid/provider.hh>
#include <dune/stuff/la/container.hh>
#include <dune/stuff/la/solver.hh>

#include <dune/gdt/spaces/continuouslagrange/pdelab.hh>
#include <dune/gdt/spaces/continuouslagrange/fem.hh>
#include <dune/gdt/spaces/continuouslagrange/fem-localfunctions.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/operators/elliptic.hh>
#include <dune/gdt/functionals/l2.hh>
#include <dune/gdt/products/l2.hh>
#include <dune/gdt/products/h1.hh>
#include <dune/gdt/products/elliptic.hh>

//#include <dune/pymor/common/exceptions.hh>
//#include <dune/pymor/la/container/eigen.hh>
//#include <dune/pymor/la/container/affine.hh>
//#include <dune/pymor/functionals/default.hh>
//#include <dune/pymor/functionals/affine.hh>
//#include <dune/pymor/operators/eigen.hh>
//#include <dune/pymor/operators/affine.hh>

//#include "../problems/interfaces.hh"
#include "base.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Discretizations {


// forward, needed in the Traits
template< class GridImp, class RangeFieldImp, int rangeDim, int polynomialOrder = 1,
          GDT::ChooseSpaceBackend space_backend = GDT::ChooseSpaceBackend::pdelab,
          Stuff::LA::ChooseBackend la_backend = Stuff::LA::ChooseBackend::eigen_sparse >
class CG;


namespace internal {


template< class GridImp, class RangeFieldImp, int rangeDim, int polynomialOrder,
          GDT::ChooseSpaceBackend space_backend,
          Stuff::LA::ChooseBackend la_backend >
class CGTraits
{
public:
  typedef CG< GridImp, RangeFieldImp, rangeDim, polynomialOrder, space_backend, la_backend > derived_type;
  typedef GridImp GridType;
  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRange = rangeDim;
  static const unsigned int polOrder = polynomialOrder;

private:
  template< class G, class R, int r, int p, GDT::ChooseSpaceBackend b >
  class SpaceChooser
  {
    static_assert(AlwaysFalse< G >::value, "No space available for this backend!");
  };

  template< class G, class R, int r, int p >
  class SpaceChooser< G, R, r, p, GDT::ChooseSpaceBackend::fem >
  {
    typedef typename Stuff::Grid::LeafPartView< G, Stuff::Grid::ChoosePartView::part >::Type GridPartViewType;
  public:
    typedef GDT::Spaces::ContinuousLagrange::FemBased< GridPartViewType, p, R, r > Type;
  };

  template< class G, class R, int r, int p >
  class SpaceChooser< G, R, r, p, GDT::ChooseSpaceBackend::fem_localfunction >
  {
    typedef typename Stuff::Grid::LeafPartView< G, Stuff::Grid::ChoosePartView::part >::Type GridPartViewType;
  public:
    typedef GDT::Spaces::ContinuousLagrange::FemLocalfunctionsBased< GridPartViewType, p, R, r > Type;
  };

  template< class G, class R, int r, int p >
  class SpaceChooser< G, R, r, p, GDT::ChooseSpaceBackend::pdelab >
  {
    typedef typename Stuff::Grid::LeafPartView< G, Stuff::Grid::ChoosePartView::view >::Type GridPartViewType;
  public:
    typedef GDT::Spaces::ContinuousLagrange::PdelabBased< GridPartViewType, p, R, r > Type;
  };

  template< class R, Stuff::LA::ChooseBackend b >
  struct ContainerChooser
  {
    static_assert(AlwaysFalse< R >::value, "No container available for this backend!");
  };

  template< class R >
  struct ContainerChooser< R, Stuff::LA::ChooseBackend::common_dense >
  {
    typedef Stuff::LA::CommonDenseMatrix< R > MatrixType;
    typedef Stuff::LA::CommonDenseVector< R > VectorType;
  };

  template< class R >
  struct ContainerChooser< R, Stuff::LA::ChooseBackend::istl_sparse >
  {
    typedef Stuff::LA::IstlRowMajorSparseMatrix< R > MatrixType;
    typedef Stuff::LA::IstlDenseVector< R > VectorType;
  };

  template< class R >
  struct ContainerChooser< R, Stuff::LA::ChooseBackend::eigen_dense >
  {
    typedef Stuff::LA::EigenDenseMatrix< R > MatrixType;
    typedef Stuff::LA::EigenDenseVector< R > VectorType;
  };

  template< class R >
  struct ContainerChooser< R, Stuff::LA::ChooseBackend::eigen_sparse >
  {
    typedef Stuff::LA::EigenRowMajorSparseMatrix< R > MatrixType;
    typedef Stuff::LA::EigenDenseVector< R > VectorType;
  };

public:
  typedef typename ContainerChooser< RangeFieldType, la_backend >::MatrixType MatrixType;
  typedef typename ContainerChooser< RangeFieldType, la_backend >::VectorType VectorType;
  typedef typename SpaceChooser< GridType, RangeFieldType, dimRange, polOrder, space_backend >::Type TestSpaceType;
  typedef TestSpaceType AnsatzSpaceType;
  typedef typename TestSpaceType::GridViewType GridViewType;
}; // class CGTraits


} // namespace internal


template< class GridImp, class RangeFieldImp, int rangeDim, int polynomialOrder,
          GDT::ChooseSpaceBackend space_backend,
          Stuff::LA::ChooseBackend la_backend >
class CG
  : public ContainerBasedDefault< internal::CGTraits< GridImp, RangeFieldImp, rangeDim, polynomialOrder, space_backend, la_backend > >
{
  typedef ContainerBasedDefault< internal::CGTraits< GridImp, RangeFieldImp, rangeDim, polynomialOrder, space_backend, la_backend > > BaseType;
public:
  typedef internal::CGTraits< GridImp, RangeFieldImp, rangeDim, polynomialOrder, space_backend, la_backend > Traits;
  typedef GridImp GridType;
  using typename BaseType::ProblemType;
  using typename BaseType::BoundaryInfoType;
  using typename BaseType::TestSpaceType;
  using typename BaseType::AnsatzSpaceType;
  using typename BaseType::MatrixType;
  using typename BaseType::VectorType;
  using typename BaseType::GridViewType;
  using typename BaseType::RangeFieldType;

private:
  typedef Stuff::Grid::ProviderInterface< GridType > GridProviderType;
  using typename BaseType::AffinelyDecomposedMatrixType;
  using typename BaseType::AffinelyDecomposedVectorType;

public:
  static std::string static_id()
  {
    return typename DiscretizationInterface< Traits >::static_id() + ".cg";
  }

  CG(const GridProviderType& grid_provider,
     const Stuff::Common::ConfigTree& bound_inf_cfg,
     const ProblemType& prob)
    : BaseType(std::make_shared< TestSpaceType >(grid_provider.template leaf< TestSpaceType::part_view_type >()),
               std::make_shared< AnsatzSpaceType >(grid_provider.template leaf< AnsatzSpaceType::part_view_type >()),
               bound_inf_cfg,
               prob)
  {}

  void init(std::ostream& out = Stuff::Common::Logger().devnull(), const std::string prefix = "")
  {
    if (!this->matrix_ || !this->rhs_) {
      Dune::Timer timer;
      this->matrix_ = std::make_shared< AffinelyDecomposedMatrixType >();
      this->rhs_ = std::make_shared< AffinelyDecomposedVectorType >();
      auto dirichlet_vector = std::make_shared< AffinelyDecomposedVectorType >();

      auto& matrix = *(this->matrix_);
      auto& rhs = *(this->rhs_);
      const auto& space = *(this->test_space_);
      const auto& grid_view = *(space.grid_view());
      const auto& boundary_info = *(this->boundary_info_);

      out << prefix << "assembling... " << std::flush;
      timer.reset();
      GDT::SystemAssembler< TestSpaceType > system_assembler(space);

      // project dirichlet boundary values
      typedef typename ProblemType::FunctionType::NonparametricType DirichletType;
      const auto& dirichlet = this->problem_.dirichlet();
      typedef GDT::DiscreteFunction< AnsatzSpaceType, VectorType > DiscreteFunctionType;
      typedef GDT::Operators::DirichletProjectionLocalizable< GridViewType, DirichletType, DiscreteFunctionType >
          DirichletProjectionOperator;
      std::vector< std::unique_ptr< DiscreteFunctionType > > dirichlet_projections;
      std::vector< std::unique_ptr< DirichletProjectionOperator > > dirichlet_projection_operators;
      for (size_t qq = 0; qq < dirichlet.num_components(); ++qq) {
        const size_t id = dirichlet_vector->register_component(new VectorType(space.mapper().size()),
                                                               dirichlet.coefficient(qq));
        dirichlet_projections.emplace_back(new DiscreteFunctionType(space, *(dirichlet_vector->component(id))));
        dirichlet_projection_operators.emplace_back(
              new DirichletProjectionOperator(grid_view,
                                              boundary_info,
                                              *(dirichlet.component(qq)),
                                              *(dirichlet_projections[dirichlet_projections.size() - 1])));
      }
      if (dirichlet.has_affine_part()) {
        dirichlet_vector->register_affine_part(new VectorType(space.mapper().size()));
        dirichlet_projections.emplace_back(new DiscreteFunctionType(space, *(dirichlet_vector->affine_part())));
        dirichlet_projection_operators.emplace_back(
              new DirichletProjectionOperator(grid_view,
                                              boundary_info,
                                              *(dirichlet.affine_part()),
                                              *(dirichlet_projections[dirichlet_projections.size() - 1])));
      }
      for (auto& projection_operator : dirichlet_projection_operators)
        system_assembler.add(*projection_operator, new GDT::ApplyOn::BoundaryEntities< GridViewType >());

      // lhs operator
      typedef typename ProblemType::DiffusionFactorType::NonparametricType DiffusionFactorType;
      const auto& diffusion_factor = this->problem_.diffusion_factor();
      typedef GDT::Operators::EllipticCG< DiffusionFactorType, MatrixType, TestSpaceType > EllipticOperatorType;
      const auto pattern = EllipticOperatorType::pattern(space);
      std::vector< std::unique_ptr< EllipticOperatorType > > elliptic_operators;
      for (size_t qq = 0; qq < diffusion_factor.num_components(); ++qq) {
        const size_t id = matrix.register_component(new MatrixType(space.mapper().size(),
                                                                   space.mapper().size(),
                                                                   pattern),
                                                    diffusion_factor.coefficient(qq));
        elliptic_operators.emplace_back(new EllipticOperatorType(*(diffusion_factor.component(qq)),
                                                                 *(matrix.component(id)),
                                                                 space));
      }
      if (diffusion_factor.has_affine_part()) {
        matrix.register_affine_part(new MatrixType(space.mapper().size(), space.mapper().size(), pattern));
        elliptic_operators.emplace_back(new EllipticOperatorType(*(diffusion_factor.affine_part()),
                                                                 *(matrix.affine_part()),
                                                                 space));
      }
      for (auto& elliptic_operator : elliptic_operators)
        system_assembler.add(*elliptic_operator);

      // rhs functional
      // * force
      typedef typename ProblemType::FunctionType::NonparametricType ForceType;
      const auto& force = this->problem_.force();
      typedef GDT::Functionals::L2Volume< ForceType, VectorType, TestSpaceType > L2VolumeFunctionalType;
      std::vector< std::unique_ptr< L2VolumeFunctionalType > > force_functionals;
      for (size_t qq = 0; qq < force.num_components(); ++qq) {
        const size_t id = rhs.register_component(new VectorType(space.mapper().size()), force.coefficient(qq));
        force_functionals.emplace_back(new L2VolumeFunctionalType(*(force.component(qq)),
                                                                  *(rhs.component(id)),
                                                                  space));
      }
      if (force.has_affine_part()) {
        rhs.register_affine_part(new VectorType(space.mapper().size()));
        force_functionals.emplace_back(new L2VolumeFunctionalType(*(force.affine_part()),
                                                                  *(rhs.affine_part()),
                                                                  space));
      }
      for (auto& force_functional : force_functionals)
        system_assembler.add(*force_functional);
      // * neumann
      typedef typename ProblemType::FunctionType::NonparametricType NeumannType;
      const auto& neumann = this->problem_.neumann();
      typedef GDT::Functionals::L2Face< NeumannType, VectorType, TestSpaceType > L2FaceFunctionalType;
      std::vector< std::unique_ptr< L2FaceFunctionalType > > neumann_functionals;
      for (size_t qq = 0; qq < neumann.num_components(); ++qq) {
        const size_t id = rhs.register_component(new VectorType(space.mapper().size()), neumann.coefficient(qq));
        neumann_functionals.emplace_back(new L2FaceFunctionalType(*(neumann.component(qq)),
                                                                  *(rhs.component(id)),
                                                                  space));
      }
      if (neumann.has_affine_part()) {
        if (!rhs.has_affine_part())
          rhs.register_affine_part(new VectorType(space.mapper().size()));
        neumann_functionals.emplace_back(new L2FaceFunctionalType(*(neumann.affine_part()),
                                                                  *(rhs.affine_part()),
                                                                  space));
      }
      for (auto& neumann_functional : neumann_functionals)
        system_assembler.add(*neumann_functional,
                             new GDT::ApplyOn::NeumannIntersections< GridViewType >(boundary_info));

      // products
      // * L2
      typedef GDT::Products::L2Assemblable< MatrixType, TestSpaceType > L2ProductType;
      auto l2_product_matrix = std::make_shared< AffinelyDecomposedMatrixType >();
      l2_product_matrix->register_affine_part(new MatrixType(space.mapper().size(),
                                                             space.mapper().size(),
                                                             L2ProductType::pattern(space)));
      L2ProductType l2_product(*(l2_product_matrix->affine_part()), space);
      system_assembler.add(l2_product);
      // * H1 semi
      typedef GDT::Products::H1SemiAssemblable< MatrixType, TestSpaceType > H1ProductType;
      auto h1_product_matrix = std::make_shared< AffinelyDecomposedMatrixType >();
      h1_product_matrix->register_affine_part(new MatrixType(space.mapper().size(),
                                                             space.mapper().size(),
                                                             H1ProductType::pattern(space)));
      H1ProductType h1_product(*(h1_product_matrix->affine_part()), space);
      system_assembler.add(h1_product);
      // * energy
      typedef GDT::Products::EllipticAssemblable< DiffusionFactorType, MatrixType, TestSpaceType > EnergyProductType;
      auto energy_product_matrix = std::make_shared< AffinelyDecomposedMatrixType >();
      std::vector< std::unique_ptr< EnergyProductType > > energy_products;
      for (size_t qq = 0; qq < diffusion_factor.num_components(); ++qq) {
        const size_t id = energy_product_matrix->register_component(new MatrixType(space.mapper().size(),
                                                                                   space.mapper().size(),
                                                                                   EnergyProductType::pattern(space)),
                                                                    diffusion_factor.coefficient(qq));
        energy_products.emplace_back(new EnergyProductType(*(diffusion_factor.component(qq)),
                                                           *(energy_product_matrix->component(id)),
                                                           space));
      }
      if (diffusion_factor.has_affine_part()) {
        energy_product_matrix->register_affine_part(new MatrixType(space.mapper().size(),
                                                                   space.mapper().size(),
                                                                   EnergyProductType::pattern(space)));
        energy_products.emplace_back(new EnergyProductType(*(diffusion_factor.affine_part()),
                                                           *(energy_product_matrix->affine_part()),
                                                           space));
      }
      for (auto& energy_product : energy_products)
        system_assembler.add(*energy_product);

      // do the actual assembling
      system_assembler.walk();
      out << "done (took " << timer.elapsed() << "s)" << std::endl;

      out << prefix << "computing dirichlet shift... " << std::flush;
      VectorType tmp(space.mapper().size());
      if (matrix.has_affine_part() && dirichlet_vector->has_affine_part()) {
        if (!rhs.has_affine_part())
          rhs.register_affine_part(new VectorType(space.mapper().size()));
        matrix.affine_part()->mv(*(dirichlet_vector->affine_part()), tmp);
        *(rhs.affine_part()) -= tmp;
      }
      if (matrix.has_affine_part()) {
        for (size_t qq = 0; qq < dirichlet_vector->num_components(); ++qq) {
          const size_t ind = rhs.register_component(new VectorType(space.mapper().size()),
                                                    dirichlet_vector->coefficient(qq));
          matrix.affine_part()->mv(*(dirichlet_vector->component(qq)), tmp);
          *(rhs.component(ind)) -= tmp;
        }
      }
      if (dirichlet_vector->has_affine_part()) {
        for (size_t qq = 0; qq < matrix.num_components(); ++qq) {
          const size_t ind = rhs.register_component(new VectorType(space.mapper().size()), matrix.coefficient(qq));
          matrix.component(qq)->mv(*(dirichlet_vector->affine_part()), tmp);
          *(rhs.component(ind)) -= tmp;
        }
      }
      Pymor::ParameterType param;
      for (auto key : matrix.parameter_type().keys())
        param.set(key, matrix.parameter_type().get(key));
      for (auto key : dirichlet_vector->parameter_type().keys())
        param.set(key, dirichlet_vector->parameter_type().get(key));
      for (size_t pp = 0; pp < matrix.num_components(); ++ pp) {
        for (size_t qq = 0; qq < dirichlet_vector->num_components(); ++qq) {
          const std::string expression = "(" + matrix.coefficient(pp)->expression()
                                         + ")*(" + dirichlet_vector->coefficient(qq)->expression() + ")";
          const size_t ind = rhs.register_component(new VectorType(space.mapper().size()),
                                                    new Pymor::ParameterFunctional(param,
                                                                                   expression));
          matrix.component(pp)->mv(*(dirichlet_vector->component(qq)), tmp);
          *(rhs.component(ind)) -= tmp;
        }
      }
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      out << prefix << "applying constraints... " << std::flush;
      GDT::Constraints::Dirichlet< typename GridViewType::Intersection, RangeFieldType, true >
        clear_and_set_dirichlet_rows(boundary_info, space.mapper().maxNumDofs(), space.mapper().maxNumDofs());
      GDT::Constraints::Dirichlet< typename GridViewType::Intersection, RangeFieldType, false >
        clear_dirichlet_rows(boundary_info, space.mapper().maxNumDofs(), space.mapper().maxNumDofs());
      // we always need an affine part in the system matrix for the dirichlet rows
      if (!matrix.has_affine_part())
        matrix.register_affine_part(new MatrixType(space.mapper().size(), space.mapper().size(), pattern));
      system_assembler.add(clear_and_set_dirichlet_rows,
                           *(matrix.affine_part()),
                           new GDT::ApplyOn::BoundaryEntities< GridViewType >());
      for (size_t qq = 0; qq < matrix.num_components(); ++qq)
        system_assembler.add(clear_dirichlet_rows, *(matrix.component(qq)),
                             new GDT::ApplyOn::BoundaryEntities< GridViewType >());
      if (rhs.has_affine_part())
        system_assembler.add(clear_dirichlet_rows, *(rhs.affine_part()),
                             new GDT::ApplyOn::BoundaryEntities< GridViewType >());
      for (size_t qq = 0; qq < matrix.num_components(); ++qq)
        system_assembler.add(clear_dirichlet_rows, *(rhs.component(qq)),
                             new GDT::ApplyOn::BoundaryEntities< GridViewType >());
      system_assembler.walk();
      out << "done (took " << timer.elapsed() << " sec)" << std::endl;

      // build parameter type
      this->inherit_parameter_type(matrix.parameter_type(), "lhs");
      this->inherit_parameter_type(rhs.parameter_type(), "rhs");
      this->inherit_parameter_type(dirichlet_vector->parameter_type(), "dirichlet");

      // finalize
      this->products_.insert(std::make_pair("l2", l2_product_matrix));
      this->products_.insert(std::make_pair("h1_semi", h1_product_matrix));
      this->products_.insert(std::make_pair("energy", energy_product_matrix));
      this->vectors_.insert(std::make_pair("dirichlet", dirichlet_vector));
    } // if !(!this->matrix_ || !this->rhs_)
  } // ... init(...)
}; // class CG


} // namespace Discretizations
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_CG_HH
