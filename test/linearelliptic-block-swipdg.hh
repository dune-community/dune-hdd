// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_TEST_LINEARELLIPTIC_BLOCK_SWIPDG_HH
#define DUNE_HDD_TEST_LINEARELLIPTIC_BLOCK_SWIPDG_HH

#include <boost/numeric/conversion/cast.hpp>

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/common/localization-study.hh>
#include <dune/stuff/grid/search.hh>
#include <dune/stuff/common/ranges.hh>

#include <dune/gdt/products/l2.hh>
#include <dune/gdt/products/h1.hh>
#include <dune/gdt/products/elliptic.hh>

#include <dune/hdd/linearelliptic/discretizations/block-swipdg.hh>
#include <dune/hdd/linearelliptic/estimators/block-swipdg.hh>
#include <dune/hdd/linearelliptic/testcases/ESV2007.hh>
#include <dune/hdd/linearelliptic/testcases/OS2014.hh>

#include "linearelliptic.hh"
#include "linearelliptic-block-swipdg-expectations.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Tests {
namespace internal {


template< class TestCaseType, int polOrder, Stuff::LA::ChooseBackend la_backend >
class DiscretizationBlockSWIPDG
{
  typedef typename TestCaseType::GridType GridType;
  typedef typename TestCaseType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange = TestCaseType::dimRange;
public:
  typedef Discretizations::BlockSWIPDG< GridType, RangeFieldType, dimRange, polOrder, la_backend >            Type;
  typedef Estimators::BlockSWIPDG
      < typename Type::AnsatzSpaceType, typename Type::VectorType, typename Type::ProblemType, GridType > EstimatorType;
}; // class DiscretizationBlockSWIPDG


} // namespace internal


template< class TestCaseType, int polOrder = 1, Stuff::LA::ChooseBackend la_backend = Stuff::LA::default_sparse_backend >
class BlockSWIPDGStudy
  : public MultiscaleEocStudyBase< TestCaseType,
                                   typename internal::DiscretizationBlockSWIPDG< TestCaseType,
                                                                                 polOrder,
                                                                                 la_backend >::Type >
  , public Stuff::Common::LocalizationStudy
{
  typedef BlockSWIPDGStudy< TestCaseType, polOrder, la_backend > ThisType;
  typedef MultiscaleEocStudyBase
      < TestCaseType,
        typename internal::DiscretizationBlockSWIPDG< TestCaseType, polOrder, la_backend >::Type > StudyBaseType;
  typedef Stuff::Common::LocalizationStudy                                                         LocalizationBaseType;

  typedef typename StudyBaseType::DiscretizationType DiscretizationType;
  typedef typename internal::DiscretizationBlockSWIPDG< TestCaseType, polOrder, la_backend >::EstimatorType
                                                     EstimatorType;
  typedef typename DiscretizationType::GridViewType GridViewType;
  typedef typename StudyBaseType::FunctionType      FunctionType;
  typedef typename StudyBaseType::VectorType        VectorType;

public:
  BlockSWIPDGStudy(const TestCaseType& test_case,
                   const std::vector< std::string > only_these_norms = std::vector< std::string >(),
                   const std::vector< std::string > only_these_local_indicators = std::vector< std::string >(),
                   const std::string visualize_prefix = "")
    : StudyBaseType(test_case, only_these_norms, visualize_prefix)
    , LocalizationBaseType(only_these_local_indicators)
  {}

  virtual ~BlockSWIPDGStudy() {}

  virtual std::string identifier() const override final
  {
    return DiscretizationType::static_id()
        + " (polorder " + Stuff::Common::toString(polOrder)
        + ", " + this->test_case_.partitioning() + " partitioning)";
  } // ... identifier(...)

  virtual size_t expected_rate(const std::string type) const override final
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
    return BlockSWIPDGStudyExpectations< TestCaseType, polOrder >::rate(this->test_case_, type);
  } // ... expected_rate(...)

  virtual std::vector< double > expected_results(const std::string type) const override final
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
    return BlockSWIPDGStudyExpectations< TestCaseType, polOrder >::results(this->test_case_, type);
  } // ... expected_results(...)

  virtual Stuff::LA::CommonDenseVector< double > compute_reference_indicators() const override final
  {
    typedef typename TestCaseType::ProblemType ProblemType;
    typedef typename ProblemType::DiffusionFactorType::NonparametricType DiffusionFactorType;
    typedef typename ProblemType::DiffusionTensorType::NonparametricType DiffusionTensorType;
    typedef typename DiffusionFactorType::DomainType     DomainType;
    typedef typename DiffusionFactorType::RangeFieldType RangeFieldType;
    const auto mu_bar = this->get_mu(this->test_case_);
    const auto problem_mu_bar = this->test_case_.problem().with_mu(mu_bar);
    const auto diffusion_factor_mu_bar = problem_mu_bar->diffusion_factor()->affine_part();
    const auto diffusion_tensor_mu_bar = problem_mu_bar->diffusion_tensor()->affine_part();
    // get current solution
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
      GDT::Products::EllipticLocalizable< typename DiscretizationType::GridViewType, DiffusionFactorType,
                                          DifferenceType, DifferenceType, RangeFieldType,
                                          DiffusionTensorType >
          local_energy_norm(this->reference_discretization_->grid_view(),
                            difference,
                            difference,
                            *diffusion_factor_mu_bar,
                            *diffusion_tensor_mu_bar,
                            Estimators::internal::SWIPDG::over_integrate);
      // prepare
      local_energy_norm.prepare();
      const auto ms_grid = this->current_discretization_->ansatz_space()->ms_grid();
      const auto current_grid_view = this->current_discretization_->grid_view();
      Stuff::Grid::EntityInlevelSearch< typename DiscretizationType::GridViewType > entity_search(current_grid_view);
      Stuff::LA::CommonDenseVector< double > error_indicators(ms_grid->size(), 0.0);
      std::vector< size_t > fine_entities_per_subdomain(error_indicators.size(), 0);
      RangeFieldType energy_error_squared = 0.0;
      // walk the reference grid
      const auto entity_it_end = this->reference_discretization_->grid_view().template end< 0 >();
      for (auto entity_it = this->reference_discretization_->grid_view().template begin< 0 >();
           entity_it != entity_it_end;
           ++entity_it) {
        const auto& entity = *entity_it;
        // search for the father entity in the current (coarser) grid
        std::vector< DomainType > center = {entity.geometry().center()};
        const auto father_entity_ptr_ptrs = entity_search(center);
        assert(father_entity_ptr_ptrs.size() == 1);
        const auto& father_entity_ptr_ptr = father_entity_ptr_ptrs[0];
        assert(father_entity_ptr_ptr);
        const auto father_entity_ptr = *father_entity_ptr_ptr;
        const auto& father_entity = *father_entity_ptr;
        assert(current_grid_view.contains(father_entity));
        const size_t current_subdomain = ms_grid->subdomainOf(father_entity);
        ++fine_entities_per_subdomain[current_subdomain];
        const RangeFieldType local_energy_error_squared = local_energy_norm.compute_locally(entity);
        energy_error_squared += local_energy_error_squared;
        error_indicators[current_subdomain] += local_energy_error_squared;
      } // walk the reference grid
      // average
      for (size_t ii = 0; ii < error_indicators.size(); ++ii)
        error_indicators[ii] /= (energy_error_squared * fine_entities_per_subdomain[ii]);
      if (!this->visualize_prefix_.empty())
      visualize_indicators(*ms_grid,
                           error_indicators,
                           "energy",
                           this->visualize_prefix_ + "_indicators_energy_error_"
                              + DSC::toString(this->current_refinement_));
      return error_indicators;
    } // this->test_case_.provides_exact_solution()
  } // ... compute_reference_indicators(...)

  virtual std::vector< std::string > provided_indicators() const
  {
    return EstimatorType::available_local();
  }

  virtual Stuff::LA::CommonDenseVector< double > compute_indicators(const std::string type) const
  {
    // get current solution
    const_cast< ThisType& >(*this).compute_on_current_refinement();
    assert(this->current_solution_vector_on_level_);
    auto indicators = EstimatorType::estimate_local(*this->current_discretization_->ansatz_space(),
                                                    *this->current_solution_vector_on_level_,
                                                    this->test_case_.problem(),
                                                    type,
                                                    this->test_case_.parameters());
    if (!this->visualize_prefix_.empty())
      visualize_indicators(*this->current_discretization_->ansatz_space()->ms_grid(),
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
  virtual std::vector< std::string > available_norms() const override final
  {
    std::vector< std::string > norms = {"L2", "H1_semi"};
    if (this->test_case_.parametric()) {
      for (auto parameter : this->test_case_.parameters())
        norms.push_back("energy_" + parameter.first);
    } else
      norms.push_back("energy");
    return norms;
  } // ... available_norms(...)

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
    } else if (type.substr(0, 7) == "energy_" && type.size() > 7) {
      const auto parameter_id = type.substr(7);
      assert(this->test_case_.parametric());
      const auto mu = this->test_case_.parameters().at(parameter_id);
      const auto nonparametric_problem = this->test_case_.problem().with_mu(mu);
      Products::Elliptic< GridViewType, DiffusionFactorType, double, DiffusionTensorType >
          elliptic_product(grid_view,
                           *nonparametric_problem->diffusion_factor()->affine_part(),
                           *nonparametric_problem->diffusion_tensor()->affine_part());
      return elliptic_product.induced_norm(function);
    } else
      DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Wrong type '" << type << "' requested!");
  } // ... compute_norm(...)

  /**
   * \attention Have to be returned in the correct order: if one type string is contained in another, the larger one
   *            has to appear first!
   */
  std::vector< std::string > effectivities() const
  {
    return { "OS2014_alt", "OS2014_*", "OS2014" };
  }

  virtual std::vector< std::string > available_estimators() const override final
  {
    auto ret = EstimatorType::available();
    for (auto id : effectivities()) {
      if (std::find(ret.begin(), ret.end(), "eta_" + id) != ret.end()) {
        if (this->test_case_.parametric()) {
          for (auto parameter : this->test_case_.parameters())
            ret.push_back("eff_" + id + "_" + parameter.first);
        } else
          ret.push_back("eff_" + id);
      }
    }
    return ret;
  } // ... available_estimators(..)

  virtual double estimate(const VectorType& vector, const std::string type) const override final
  {
    // process all effectivities
    for (auto id : effectivities()) {
      if (type.substr(0, ("eff_" + id).size()) == "eff_" + id) {
        if (this->test_case_.parametric()
            && type.substr(0, ("eff_" + id + "_").size()) == "eff_" + id + "_"
            && type.size() > ("eff_" + id + "_").size()) {
          const auto parameter_id = type.substr(("eff_" + id + "_").size());
          return estimate(vector, "eta_" + id)
              / const_cast< ThisType& >(*this).current_error_norm("energy_" + parameter_id);
        } else {
          assert(type == "eff_" + id);
          assert(!this->test_case_.parametric());
          return estimate(vector, "eta_" + id) / const_cast< ThisType& >(*this).current_error_norm("energy");
        }
      }
    }
    // else we have a normal estimator
    assert(this->current_discretization_);
    return EstimatorType::estimate(*this->current_discretization_->ansatz_space(),
                                   vector,
                                   this->test_case_.problem(),
                                   type,
                                   this->test_case_.parameters());
  } // ... estimate(...)

  template< class MSG, class VV >
  void visualize_indicators(const MSG& ms_grid,
                            const VV& vector,
                            const std::string name,
                            const std::string filename) const
  {
    assert(vector.size() == ms_grid.size());
    const auto grid_view = ms_grid.globalGridView();
    VV fine_vector(boost::numeric_cast< size_t >(grid_view.indexSet().size(0)), 0.0);
    for (const auto& entity : Stuff::Common::entityRange(grid_view)) {
      const size_t index = grid_view.indexSet().index(entity);
      const size_t subdomain = ms_grid.subdomainOf(entity);
      fine_vector[index] = vector[subdomain];
    }
    typedef GDT::Spaces::FiniteVolume::Default< typename MSG::GlobalGridViewType, typename VV::ScalarType, 1 >
        FVSpaceType;
    const FVSpaceType fv_space(grid_view);
    GDT::ConstDiscreteFunction< FVSpaceType, VV > discrete_function(fv_space, fine_vector, name);
    discrete_function.visualize(filename);
  } // ... visualize_indicators(...)
}; // class BlockSWIPDGStudy


} // namespace Tests
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_TEST_LINEARELLIPTIC_BLOCK_SWIPDG_HH
