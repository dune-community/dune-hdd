// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_TEST_LINEARELLIPTIC_SWIPDG_HH
#define DUNE_HDD_TEST_LINEARELLIPTIC_SWIPDG_HH

#include <algorithm>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/localization-study.hh>

#include <dune/gdt/products/l2.hh>
#include <dune/gdt/products/h1.hh>
#include <dune/gdt/playground/products/elliptic.hh>
#include <dune/gdt/playground/spaces/finitevolume/default.hh>
#include <dune/gdt/discretefunction/default.hh>

#include <dune/hdd/linearelliptic/discretizations/swipdg.hh>
#include <dune/hdd/linearelliptic/discretizations/swipdg-estimator.hh>
#include <dune/hdd/playground/linearelliptic/testcases/ESV2007.hh>
#include <dune/hdd/playground/linearelliptic/testcases/OS2014.hh>

#include "linearelliptic.hh"
#include "linearelliptic-swipdg-expectations.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Tests {
namespace internal {


template< class TestCaseType, int polOrder, GDT::ChooseSpaceBackend space_backend, Stuff::LA::ChooseBackend la_backend >
class DiscretizationSWIPDG
{
  typedef typename TestCaseType::GridType GridType;
  typedef typename TestCaseType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange = TestCaseType::dimRange;
public:
  typedef Discretizations::SWIPDG
      < GridType, Stuff::Grid::ChooseLayer::level, RangeFieldType, dimRange, polOrder, space_backend, la_backend > Type;
  typedef Discretizations::SWIPDGEstimator< typename Type::AnsatzSpaceType,
                                            typename Type::VectorType,
                                            typename Type::ProblemType,
                                            GridType >
      EstimatorType;
}; // class DiscretizationSWIPDG


} // namespace internal


template< class TestCaseType, int polOrder, GDT::ChooseSpaceBackend space_backend, Stuff::LA::ChooseBackend la_backend >
class SWIPDGStudy
  : public EocStudyBase< TestCaseType,
                         typename internal::DiscretizationSWIPDG< TestCaseType,
                                                                  polOrder,
                                                                  space_backend,
                                                                  la_backend >::Type >
  , public Stuff::Common::LocalizationStudy
{
  typedef EocStudyBase< TestCaseType,
                        typename internal::DiscretizationSWIPDG< TestCaseType,
                                                                 polOrder,
                                                                 space_backend,
                                                                 la_backend >::Type > StudyBaseType;
  typedef Stuff::Common::LocalizationStudy                                            LocalizationBaseType;
  typedef SWIPDGStudy< TestCaseType, polOrder, space_backend, la_backend >            ThisType;

  typedef typename StudyBaseType::DiscretizationType DiscretizationType;
  typedef typename internal::DiscretizationSWIPDG< TestCaseType, polOrder, space_backend, la_backend >::EstimatorType
                                                     EstimatorType;
  typedef typename StudyBaseType::GridViewType GridViewType;
  typedef typename StudyBaseType::FunctionType FunctionType;
  typedef typename StudyBaseType::VectorType   VectorType;

public:
  SWIPDGStudy(const TestCaseType& test_case,
              const std::vector< std::string > only_these_norms = {},
              const std::vector< std::string > only_these_local_norms = {},
              const std::string visualize_prefix = "")
    : StudyBaseType(test_case, only_these_norms, visualize_prefix)
    , LocalizationBaseType(only_these_local_norms)
  {}

  virtual ~SWIPDGStudy() {}

  virtual std::string identifier() const DS_OVERRIDE DS_FINAL
  {
    return DiscretizationType::static_id() + " (polorder " + Stuff::Common::toString(polOrder) + ")";
  }

  virtual size_t expected_rate(const std::string type) const DS_OVERRIDE DS_FINAL
  {
    // If you get an undefined reference here from the linker you are missing the appropriate
    // specialization of BlockSWIPDGStudyExpectations!
    // For a new TestCaseType you have to add a specialization in a separate object file
    // (see linearelliptic-block-swipdg-expectations_os2014_2daluconform.cxx for example) and adjust the
    // CMakeLists.txt accordingly. For a new polOrder add
    //     template class BlockSWIPDGStudyExpectations< TestCasesType, polOrder >;
    // in the appropriate (existing) object file and implement a specialization for this polOrder, if needed!
    //
    // Oh: and do not forget to add
    //   'extern template class BlockSWIPDGStudyExpectations< ... >'
    // to each test source using these results!
    return SWIPDGStudyExpectations< TestCaseType, polOrder >::rate(this->test_case_, type);
  } // ... expected_rate(...)

  virtual std::vector< double > expected_results(const std::string type) const DS_OVERRIDE DS_FINAL
  {
    // If you get an undefined reference here from the linker you are missing the appropriate
    // specialization of BlockSWIPDGStudyExpectations!
    // For a new TestCaseType you have to add a specialization in a separate object file
    // (see linearelliptic-block-swipdg-expectations_os2014_2daluconform.cxx for example) and adjust the
    // CMakeLists.txt accordingly. For a new polOrder add
    //     template class BlockSWIPDGStudyExpectations< TestCasesType, polOrder >;
    // in the appropriate (existing) object file and implement a specialization for this polOrder, if needed!
    //
    // Oh: and do not forget to add
    //   'extern template class BlockSWIPDGStudyExpectations< ... >'
    // to each test source using these results!
    return SWIPDGStudyExpectations< TestCaseType, polOrder >::results(this->test_case_, type);
  } // ... expected_results(...)

  virtual Stuff::LA::CommonDenseVector< double > compute_reference_indicators() const
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
    typedef typename StudyBaseType::DiscreteFunctionType      DiscreteFunctionType;
    typedef typename StudyBaseType::ConstDiscreteFunctionType ConstDiscreteFunctionType;
    const ConstDiscreteFunctionType current_solution(*(this->reference_discretization_->ansatz_space()),
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
      const ConstDiscreteFunctionType reference_solution(*(this->reference_discretization_->ansatz_space()),
                                                         *this->reference_solution_vector_,
                                                         "reference solution");
      // define error norm
      const auto difference = reference_solution - current_solution;
      typedef typename ConstDiscreteFunctionType::DifferenceType DifferenceType;
      GDT::Products::EllipticLocalizable< DiffusionFactorType, GridViewType,
                                          DifferenceType, DifferenceType, RangeFieldType,
                                          DiffusionTensorType >
          local_energy_norm(*diffusion_factor,
                            *diffusion_tensor,
                            *(this->test_case_.reference_grid_view()),
                            difference,
                            difference,
                            Discretizations::internal::SWIPDGEstimators::over_integrate);
      local_energy_norm.prepare();
      // prepare
      const auto current_grid_view = this->current_discretization_->grid_view();
      const int current_level = current_grid_view->template begin< 0 >()->level();
      Stuff::LA::CommonDenseVector< double > error_indicators(boost::numeric_cast< size_t >(current_grid_view->indexSet().size(0)),
                                                              0.0);
      std::vector< size_t > fine_entities_per_coarse_entity(error_indicators.size(), 0);
      RangeFieldType energy_error_squared = 0.0;
      // walk the reference grid
      const auto entity_it_end = this->test_case_.reference_grid_view()->template end< 0 >();
      for (auto entity_it = this->test_case_.reference_grid_view()->template begin< 0 >();
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
        assert(current_grid_view->indexSet().contains(father_entity));
        const size_t index = current_grid_view->indexSet().index(father_entity);
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

  virtual std::vector< std::string > provided_indicators() const
  {
    return EstimatorType::available_local();
  }

  virtual Stuff::LA::CommonDenseVector< double > compute_indicators(const std::string type) const
  {
    // get current solution
    assert(this->current_refinement_ <= this->num_refinements());
    if (this->last_computed_refinement_ != this->current_refinement_) {
      const_cast< ThisType& >(*this).compute_on_current_refinement();
    }
    assert(this->current_solution_vector_on_level_);
    auto indicators = EstimatorType::estimate_local(*this->current_discretization_->ansatz_space(),
                                                    *this->current_solution_vector_on_level_,
                                                    this->test_case_.problem(),
                                                    type);
    if (!this->visualize_prefix_.empty())
      visualize_indicators(this->current_discretization_->grid_view(),
                           indicators,
                           type,
                           this->visualize_prefix_ + "_indicators_" + type + "_"
                           + DSC::toString(this->current_refinement_));
    return indicators;
  } // ... compute_indicators(...)

  std::map< std::string, std::vector< double > > run_eoc(std::ostream& out)
  {
    return StudyBaseType::run(out);
  }

  void run_localization(std::ostream& out)
  {
    LocalizationBaseType::run(out);
  }

private:
  virtual std::vector< std::string > available_norms() const DS_OVERRIDE DS_FINAL
  {
    return {"L2", "H1_semi", "energy"};
  }

  virtual double compute_norm(const GridViewType& grid_view,
                              const FunctionType& function,
                              const std::string type) const DS_OVERRIDE DS_FINAL
  {
    using namespace GDT;
    typedef typename TestCaseType::ProblemType::DiffusionFactorType::NonparametricType DiffusionFactorType;
    typedef typename TestCaseType::ProblemType::DiffusionTensorType::NonparametricType DiffusionTensorType;
    if (type == "L2") {
      return Products::L2< GridViewType >(grid_view).induced_norm(function);
    } else if (type == "H1_semi") {
      return Products::H1SemiGeneric< GridViewType >(grid_view).induced_norm(function);
    } else if (type == "energy") {
      const auto& diffusion_factor = *(this->test_case_.problem().diffusion_factor());
      assert(!diffusion_factor.parametric());
      assert(diffusion_factor.has_affine_part());
      const auto& diffusion_tensor = *(this->test_case_.problem().diffusion_tensor());
      assert(!diffusion_tensor.parametric());
      assert(diffusion_tensor.has_affine_part());
      Products::Elliptic< DiffusionFactorType, GridViewType, double, DiffusionTensorType >
          elliptic_product(*diffusion_factor.affine_part(), *diffusion_tensor.affine_part(), grid_view);
      return elliptic_product.induced_norm(function);
    } else
      DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Wrong type '" << type << "' requested!");
  } // ... compute_norm(...)

  virtual std::vector< std::string > available_estimators() const DS_OVERRIDE DS_FINAL
  {
    auto ret = EstimatorType::available();
    if (std::find(ret.begin(), ret.end(), "eta_ESV2007") != ret.end())
      ret.push_back("eff_ESV2007");
    if (std::find(ret.begin(), ret.end(), "eta_ESV2007_alt") != ret.end())
      ret.push_back("eff_ESV2007_alt");
    return ret;
  }

  virtual double estimate(const VectorType& vector, const std::string type) const DS_OVERRIDE DS_FINAL
  {
    if (type == "eff_ESV2007")
      return estimate(vector, "eta_ESV2007") / const_cast< ThisType& >(*this).current_error_norm("energy");
    else if (type == "eff_ESV2007_alt")
      return estimate(vector, "eta_ESV2007_alt") / const_cast< ThisType& >(*this).current_error_norm("energy");
    else {
      assert(this->current_discretization_);
      return EstimatorType::estimate(*this->current_discretization_->ansatz_space(),
                                     vector,
                                     this->test_case_.problem(),
                                     type);
    }
  } // ... estimate(...)

  template< class GV, class VV >
  void visualize_indicators(const std::shared_ptr< const GV >& grid_view_ptr,
                            const VV& vector,
                            const std::string name,
                            const std::string filename) const
  {
    typedef GDT::Spaces::FiniteVolume::Default< GV, typename VV::ScalarType, 1 > FVSpaceType;
    const FVSpaceType fv_space(grid_view_ptr);
    GDT::ConstDiscreteFunction< FVSpaceType, VV > discrete_function(fv_space, vector, name);
    discrete_function.visualize(filename);
  } // ... visualize_indicators(...)
}; // class SWIPDGStudy


} // namespace Tests
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_TEST_LINEARELLIPTIC_SWIPDG_HH
