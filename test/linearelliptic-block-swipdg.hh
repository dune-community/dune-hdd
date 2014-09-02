// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_TEST_LINEARELLIPTIC_BLOCK_SWIPDG_HH
#define DUNE_HDD_TEST_LINEARELLIPTIC_BLOCK_SWIPDG_HH

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/type_utils.hh>

#include <dune/gdt/products/l2.hh>
#include <dune/gdt/products/h1.hh>
#include <dune/gdt/playground/products/elliptic.hh>

#include <dune/hdd/playground/linearelliptic/discretizations/block-swipdg.hh>
#include <dune/hdd/playground/linearelliptic/discretizations/block-swipdg-estimator.hh>
#include <dune/hdd/playground/linearelliptic/testcases/ESV2007.hh>
#include <dune/hdd/playground/linearelliptic/testcases/OS2014.hh>

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
  typedef Discretizations::BlockSWIPDG< GridType, RangeFieldType, dimRange, polOrder, la_backend > Type;
  typedef Discretizations::BlockSWIPDGEstimator< typename Type::AnsatzSpaceType,
                                                 typename Type::VectorType,
                                                 typename Type::ProblemType,
                                                 GridType >
      EstimatorType;
}; // class DiscretizationBlockSWIPDG


} // namespace internal


template< class TestCaseType, int polOrder, Stuff::LA::ChooseBackend la_backend >
class EocStudyBlockSWIPDG
  : public MultiscaleEocStudyBase< TestCaseType,
                                   typename internal::DiscretizationBlockSWIPDG< TestCaseType,
                                                                                 polOrder,
                                                                                 la_backend >::Type >
{
  typedef EocStudyBlockSWIPDG< TestCaseType, polOrder, la_backend > ThisType;
  typedef MultiscaleEocStudyBase
      < TestCaseType,
        typename internal::DiscretizationBlockSWIPDG< TestCaseType, polOrder, la_backend >::Type > BaseType;

  typedef typename BaseType::DiscretizationType      DiscretizationType;
  typedef typename internal::DiscretizationBlockSWIPDG< TestCaseType, polOrder, la_backend >::EstimatorType
      EstimatorType;
  typedef typename DiscretizationType::GridViewType GridViewType;
  typedef typename BaseType::FunctionType FunctionType;
  typedef typename BaseType::VectorType   VectorType;

public:
  EocStudyBlockSWIPDG(const TestCaseType& test_case,
                      const std::vector< std::string > only_these_norms = std::vector< std::string >())
    : BaseType(test_case, only_these_norms)
  {}

  virtual ~EocStudyBlockSWIPDG() {}

  virtual std::string identifier() const DS_OVERRIDE DS_FINAL
  {
    return DiscretizationType::static_id()
        + " (polorder " + Stuff::Common::toString(polOrder)
        + ", " + this->test_case_.partitioning() + " partitioning)";
  } // ... identifier(...)

  virtual size_t expected_rate(const std::string type) const DS_OVERRIDE DS_FINAL
  {
    // If you get an undefined reference here from the linker you are missing the appropriate
    // specialization of BlockSWIPDGStudyExpectations!
    // For a new TestCaseType you have to add a specialization in a separate object file
    // (see linearelliptic-block-swipdg-expectations_os2014_2daluconform.cxx for example) and adjust the
    // CMakeLists.txt accordingly. For a new polOrder add
    //     template class BlockSWIPDGStudyExpectations< TestCasesType, polOrder, true >;
    // in the appropriate (existing) object file and implement a specialization for this polOrder, if needed!
    return BlockSWIPDGStudyExpectations< TestCaseType, polOrder >::rate(this->test_case_, type);
  } // ... expected_rate(...)

  virtual std::vector< double > expected_results(const std::string type) const DS_OVERRIDE DS_FINAL
  {
    // If you get an undefined reference here from the linker you are missing the appropriate
    // specialization of BlockSWIPDGStudyExpectations!
    // For a new TestCaseType you have to add a specialization in a separate object file
    // (see linearelliptic-block-swipdg-expectations_os2014_2daluconform.cxx for example) and adjust the
    // CMakeLists.txt accordingly. For a new polOrder add
    //     template class BlockSWIPDGStudyExpectations< TestCasesType, polOrder, true >;
    // in the appropriate (existing) object file and implement a specialization for this polOrder, if needed!
    return BlockSWIPDGStudyExpectations< TestCaseType, polOrder >::results(this->test_case_, type);
  } // ... expected_results(...)

private:
  virtual std::vector< std::string > available_norms() const DS_OVERRIDE DS_FINAL
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
    } else if (type.substr(0, 7) == "energy_" && type.size() > 7) {
      const auto parameter_id = type.substr(7);
      assert(this->test_case_.parametric());
      const auto mu = this->test_case_.parameters().at(parameter_id);
      const auto nonparametric_problem = this->test_case_.problem().with_mu(mu);
      Products::Elliptic< DiffusionFactorType, GridViewType, double, DiffusionTensorType >
          elliptic_product(*nonparametric_problem->diffusion_factor()->affine_part(),
                           *nonparametric_problem->diffusion_tensor()->affine_part(), grid_view);
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

  virtual std::vector< std::string > available_estimators() const DS_OVERRIDE DS_FINAL
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

  virtual double estimate(const VectorType& vector, const std::string type) const DS_OVERRIDE DS_FINAL
  {
    // process all efficitvities
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
}; // class EocStudyBlockSWIPDG


//#if HAVE_ALUGRID && HAVE_DUNE_GRID_MULTISCALE && HAVE_DUNE_FEM && HAVE_EIGEN


//extern template class EocStudyBlockSWIPDG< TestCases::ESV2007Multiscale< ALUGrid< 2, 2, simplex, conforming > >,
//                                           1,
//                                           Stuff::LA::ChooseBackend::eigen_sparse >;

//extern template class EocStudyBlockSWIPDG< TestCases::OS2014Multiscale< ALUGrid< 2, 2, simplex, conforming > >,
//                                           1,
//                                           Stuff::LA::ChooseBackend::eigen_sparse >;


//#endif // HAVE_ALUGRID && HAVE_DUNE_GRID_MULTISCALE && HAVE_DUNE_FEM && HAVE_EIGEN


} // namespace Tests
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_TEST_LINEARELLIPTIC_BLOCK_SWIPDG_HH
