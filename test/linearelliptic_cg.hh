// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_TEST_LINEARELLIPTIC_CG_HH
#define DUNE_HDD_TEST_LINEARELLIPTIC_CG_HH

#include <algorithm>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/localization-study.hh>

#include <dune/gdt/products/l2.hh>
#include <dune/gdt/products/h1.hh>
#include <dune/gdt/products/elliptic.hh>
#include <dune/gdt/spaces/cg/pdelab.hh>
#include <dune/gdt/discretefunction/default.hh>

#include <dune/hdd/linearelliptic/discretizations/cg.hh>
#include <dune/hdd/linearelliptic/estimators/swipdg.hh>
#include <dune/hdd/linearelliptic/testcases/ESV2007.hh>
#include <dune/hdd/linearelliptic/testcases/OS2014.hh>

#include "linearelliptic.hh"
#include "linearelliptic-swipdg-expectations.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Tests {


template< class TestCaseType, int polOrder, bool  = true >
class CGStudyExpectations
{
public:
  static size_t rate(const TestCaseType& /*test_case*/, const std::string type)
  {
    return type == "L2" ? polOrder + 1 : polOrder;
  } // ... rate(...)

  static std::vector< double > results(const TestCaseType& /*test_case*/, const std::string type)
  {
    return {};
  } // ... results(...)
}; // CGStudyExpectations
namespace internal {
template< class TestCaseType, int polOrder, GDT::ChooseSpaceBackend space_backend, Stuff::LA::ChooseBackend la_backend >
class DiscretizationCG
{
  typedef typename TestCaseType::GridType GridType;
  typedef typename TestCaseType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange = TestCaseType::dimRange;
public:
  typedef Discretizations::CG
      < GridType, Stuff::Grid::ChooseLayer::level, RangeFieldType, dimRange, polOrder, space_backend, la_backend > Type;
}; // class DiscretizationCG


} // namespace internal


template< class TestCaseType,
          int polOrder = 1,
          GDT::ChooseSpaceBackend space_backend = GDT::ChooseSpaceBackend::pdelab,
          Stuff::LA::ChooseBackend la_backend = Stuff::LA::default_sparse_backend >
class CGStudy
  : public DHLT::EocStudyBase< TestCaseType,
                         typename internal::DiscretizationCG< TestCaseType,
                                                                  polOrder,
                                                                  space_backend,
                                                                  la_backend >::Type >
  , public Stuff::Common::LocalizationStudy
{
  typedef DHLT::EocStudyBase< TestCaseType,
                        typename internal::DiscretizationCG< TestCaseType,
                                                                 polOrder,
                                                                 space_backend,
                                                                 la_backend >::Type > StudyBaseType;
  typedef Stuff::Common::LocalizationStudy                                            LocalizationBaseType;
  typedef CGStudy< TestCaseType, polOrder, space_backend, la_backend >            ThisType;
public:
  typedef typename StudyBaseType::DiscretizationType DiscretizationType;
  typedef typename StudyBaseType::GridViewType GridViewType;
  typedef typename StudyBaseType::FunctionType FunctionType;
  typedef typename StudyBaseType::VectorType   VectorType;

public:
  CGStudy(TestCaseType& test_case,
              const std::vector< std::string > only_these_norms = {},
              const std::vector< std::string > only_these_local_norms = {},
              const std::string visualize_prefix = "")
    : StudyBaseType(test_case, only_these_norms, visualize_prefix)
    , LocalizationBaseType(only_these_local_norms)
  {}

  virtual ~CGStudy() {}

  virtual std::string identifier() const override final
  {
    return DiscretizationType::static_id() + " (polorder " + Stuff::Common::toString(polOrder) + ")";
  }

  virtual DS::LA::CommonDenseVector< double > compute_indicators(const std::string type) const {
  return DS::LA::CommonDenseVector< double >();
}

  virtual size_t expected_rate(const std::string type) const override final
  {
    // If you get an undefined reference here from the linker you are missing the appropriate
    // specialization of BlockCGStudyExpectations!
    // For a new TestCaseType you have to add a specialization in a separate object file
    // (see linearelliptic-block-CG-expectations_os2014_2daluconform.cxx for example) and adjust the
    // CMakeLists.txt accordingly. For a new polOrder add
    //     template class BlockCGStudyExpectations< TestCasesType, polOrder >;
    // in the appropriate (existing) object file and implement a specialization for this polOrder, if needed!
    //
    // Oh: and do not forget to add
    //   'extern template class BlockCGStudyExpectations< ... >'
    // to each test source using these results!
    return CGStudyExpectations< TestCaseType, polOrder >::rate(this->test_case_, type);
  } // ... expected_rate(...)

  virtual std::vector< double > expected_results(const std::string type) const override final
  {
    // If you get an undefined reference here from the linker you are missing the appropriate
    // specialization of BlockCGStudyExpectations!
    // For a new TestCaseType you have to add a specialization in a separate object file
    // (see linearelliptic-block-CG-expectations_os2014_2daluconform.cxx for example) and adjust the
    // CMakeLists.txt accordingly. For a new polOrder add
    //     template class BlockCGStudyExpectations< TestCasesType, polOrder >;
    // in the appropriate (existing) object file and implement a specialization for this polOrder, if needed!
    //
    // Oh: and do not forget to add
    //   'extern template class BlockCGStudyExpectations< ... >'
    // to each test source using these results!
    return CGStudyExpectations< TestCaseType, polOrder >::results(this->test_case_, type);
  } // ... expected_results(...)

  virtual Stuff::LA::CommonDenseVector< double > compute_reference_indicators() const override final
  {
    typedef typename TestCaseType::ProblemType ProblemType;
    typedef typename ProblemType::DiffusionFactorType::NonparametricType DiffusionFactorType;
    typedef typename ProblemType::DiffusionTensorType::NonparametricType DiffusionTensorType;
    typedef typename DiffusionFactorType::RangeFieldType RangeFieldType;
    assert(!this->test_case_.problem().diffusion_factor()->parametric());
    assert(!this->test_case_.problem().diffusion_tensor()->parametric());
    const auto diffusion_factor = this->test_case_.problem().diffusion_factor()->affine_part();
    const auto diffusion_tensor = this->test_case_.problem().diffusion_tensor()->affine_part();

    // get current solution
    assert(this->current_refinement_ <= this->num_refinements());
    if (this->last_computed_refinement_ != this->current_refinement_)
      const_cast< ThisType& >(*this).compute_on_current_refinement();
    assert(this->last_computed_refinement_ == this->current_refinement_);
    assert(this->current_solution_vector_);
    typedef typename StudyBaseType::ConstDiscreteFunctionType ConstDiscreteFunctionType;
    const ConstDiscreteFunctionType current_solution(this->reference_discretization_->ansatz_space(),
                                                     *this->current_solution_vector_,
                                                     "current solution");
    // compute error
    if (this->test_case_.provides_exact_solution()) {
      DUNE_THROW(Stuff::Exceptions::you_have_to_implement_this, "Felix!");
    } else {
      if (!this->reference_solution_computed_)
        const_cast< ThisType& >(*this).compute_reference_solution();
      assert(this->reference_discretization_);
      assert(this->reference_solution_vector_);
      // get reference solution
      const ConstDiscreteFunctionType reference_solution(this->reference_discretization_->ansatz_space(),
                                                         *this->reference_solution_vector_,
                                                         "reference solution");
      // define error norm
      const auto difference = reference_solution - current_solution;
      typedef typename ConstDiscreteFunctionType::DifferenceType DifferenceType;
      GDT::Products::EllipticLocalizable< GridViewType, DiffusionFactorType,
                                          DifferenceType, DifferenceType, RangeFieldType,
                                          DiffusionTensorType >
          local_energy_norm(this->test_case_.reference_grid_view(),
                            difference,
                            difference,
                            *diffusion_factor,
                            *diffusion_tensor,
                            false);
      local_energy_norm.prepare();
      // prepare
      const auto& current_grid_view = this->current_discretization_->grid_view();
      const int current_level = current_grid_view.template begin< 0 >()->level();
      Stuff::LA::CommonDenseVector< double > error_indicators(boost::numeric_cast< size_t >(current_grid_view.indexSet().size(0)),
                                                              0.0);
      std::vector< size_t > fine_entities_per_coarse_entity(error_indicators.size(), 0);
      RangeFieldType energy_error_squared = 0.0;
      // walk the reference grid
      auto reference_grid_view = this->test_case_.reference_grid_view();
      const auto entity_it_end = reference_grid_view.template end< 0 >();
      for (auto entity_it = reference_grid_view.template begin< 0 >();
           entity_it != entity_it_end;
           ++entity_it) {
        const auto& entity = *entity_it;
        int fine_level = entity.level();
        typename GridViewType::template Codim< 0 >::EntityPointer father_entity_ptr(entity);
        // find father entity in coarse grid view
        assert(fine_level >= current_level);
        while (fine_level > current_level) {
          assert(father_entity_ptr->hasFather());
          father_entity_ptr = father_entity_ptr->father();
          fine_level = father_entity_ptr.level();
        }
        assert(fine_level == current_level);
        const auto& father_entity = *father_entity_ptr;
        assert(current_grid_view.indexSet().contains(father_entity));
        const size_t index = current_grid_view.indexSet().index(father_entity);
        ++fine_entities_per_coarse_entity[index];
        const RangeFieldType local_energy_error_squared = local_energy_norm.compute_locally(entity);
        energy_error_squared += local_energy_error_squared;
        error_indicators[index] += local_energy_error_squared;
      } // walk the reference grid
      // average
      for (size_t ii = 0; ii < error_indicators.size(); ++ii)
        error_indicators[ii] /= (energy_error_squared * fine_entities_per_coarse_entity[ii]);
      if (!this->visualize_prefix_.empty())
        visualize_indicators(current_grid_view,
                             error_indicators,
                             "energy",
                             this->visualize_prefix_ + "_indicators_energy_error_"
                                + DSC::toString(this->current_refinement_));
      return error_indicators;
    }
  } // ... compute_reference_indicators(...)

  virtual std::vector< std::string > provided_indicators() const override final
  {
    return {};
  }

  std::map< std::string, std::vector< double > > run_eoc(std::ostream& out)
  {
    return StudyBaseType::run(out);
  }

  void run_localization(std::ostream& out)
  {
    LocalizationBaseType::run(out);
  }

private:
  virtual std::vector< std::string > available_norms() const override final
  {
    return {"L2", "H1_semi", "energy"};
  }

  virtual double compute_norm(const GridViewType& grid_view,
                              const FunctionType& function,
                              const std::string type) const override final
  {
    using namespace GDT;
    typedef typename TestCaseType::ProblemType::DiffusionFactorType::NonparametricType DiffusionFactorType;
    typedef typename TestCaseType::ProblemType::DiffusionTensorType::NonparametricType DiffusionTensorType;
    if (type == "L2") {
      return Products::L2< GridViewType >(grid_view).induced_norm(function);
    } else if (type == "H1_semi") {
      return Products::H1Semi< GridViewType >(grid_view).induced_norm(function);
    } else if (type == "energy") {
      const auto& diffusion_factor = *(this->test_case_.problem().diffusion_factor());
      assert(!diffusion_factor.parametric());
      assert(diffusion_factor.has_affine_part());
      const auto& diffusion_tensor = *(this->test_case_.problem().diffusion_tensor());
      assert(!diffusion_tensor.parametric());
      assert(diffusion_tensor.has_affine_part());
      Products::Elliptic< GridViewType, DiffusionFactorType, double, DiffusionTensorType >
          elliptic_product(grid_view, *diffusion_factor.affine_part(), *diffusion_tensor.affine_part());
      return elliptic_product.induced_norm(function);
    } else
      DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Wrong type '" << type << "' requested!");
  } // ... compute_norm(...)

  virtual std::vector< std::string > available_estimators() const override final
  {
    return {};
  }

  virtual double estimate(const VectorType& vector, const std::string type) const override final
  {
    return 1;
  } // ... estimate(...)

  template< class GV, class VV >
  void visualize_indicators(const GV& grid_view,
                            const VV& vector,
                            const std::string name,
                            const std::string filename) const
  {
    typedef GDT::Spaces::FV::Default< GV, typename VV::ScalarType, 1 > FVSpaceType;
    const FVSpaceType fv_space(grid_view);
    GDT::ConstDiscreteFunction< FVSpaceType, VV > discrete_function(fv_space, vector, name);
    discrete_function.visualize(filename);
  } // ... visualize_indicators(...)
}; // class CGStudy


} // namespace Tests
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_TEST_LINEARELLIPTIC_CG_HH
