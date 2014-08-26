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

#include "linearelliptic.hh"

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
    if (type == "L2")
      return polOrder + 1;
    else if (type == "H1_semi")
      return polOrder;
    else if (type.substr(0, 6) == "energy")
      return polOrder;
    else if (type == "eta_NC_OS2014")
      return polOrder;
    else if (type == "eta_R_OS2014")
      return polOrder + 1;
    else if (type == "eta_DF_OS2014")
      return polOrder;
    else if (type == "eta_OS2014")
      return polOrder;
    else if (type == "eff_OS2014")
      return 0;
    else
      DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Wrong type '" << type << "' requested!");
  } // ... expected_rate(...)

  virtual std::vector< double > expected_results(const std::string type) const DS_OVERRIDE DS_FINAL
  {
    if (std::is_same< TestCaseType, TestCases::ESV2007Multiscale< ALUConformGrid< 2, 2 > > >::value
        || std::is_same< TestCaseType, TestCases::ESV2007Multiscale< ALUGrid< 2, 2, simplex, conforming > > >::value) {
      if (this->test_case_.partitioning() == "[1 1 1]") {
        if (polOrder == 1) {
          if (type == "energy")
            return {3.29e-01, 1.63e-01, 8.05e-02, 4.02e-02};
          else if (type == "eta_NC_OS2014")
            return {1.67e-01, 7.90e-02, 3.92e-02, 1.96e-02};
          else if (type == "eta_R_OS2014")
            return {5.80e-01, 2.91e-01, 1.46e-01, 7.28e-02};
          else if (type == "eta_DF_OS2014")
            return {3.56e-01, 1.77e-01, 8.74e-02, 4.36e-02};
          else if (type == "eta_OS2014")
            return {1.11e+00, 5.46e-01, 2.73e-01, 1.37e-01};
          else if (type == "eff_OS2014")
            return {3.36, 3.38, 3.39, 3.40};
          else
            DUNE_THROW(Stuff::Exceptions::test_results_missing, type);
        } else
          DUNE_THROW(Stuff::Exceptions::test_results_missing, polOrder);
      } else if (this->test_case_.partitioning() == "[2 2 1]") {
        if (polOrder == 1) {
          if (type == "energy")
            return {3.29e-01, 1.63e-01, 8.05e-02, 4.02e-02};
          else if (type == "eta_NC_OS2014")
            return {1.67e-01, 7.90e-02, 3.92e-02, 1.96e-02};
          else if (type == "eta_R_OS2014")
            return {2.90e-01, 1.46e-01, 7.28e-02, 3.64e-02};
          else if (type == "eta_DF_OS2014")
            return {3.56e-01, 1.77e-01, 8.74e-02, 4.36e-02};
          else if (type == "eta_OS2014")
            return {1.11e+00, 5.46e-01, 2.73e-01, 1.37e-01};
          else if (type == "eff_OS2014")
            return {2.48, 2.48, 2.49, 2.49};
          else
            DUNE_THROW(Stuff::Exceptions::test_results_missing, type);
        } else
          DUNE_THROW(Stuff::Exceptions::test_results_missing, polOrder);
      } else if (this->test_case_.partitioning() == "[4 4 1]") {
        if (polOrder == 1) {
          if (type == "energy")
            return {3.29e-01, 1.63e-01, 8.05e-02, 4.02e-02};
          else if (type == "eta_NC_OS2014")
            return {1.67e-01, 7.90e-02, 3.92e-02, 1.96e-02};
          else if (type == "eta_R_OS2014")
            return {1.46e-01, 7.27e-02, 3.64e-02, 1.82e-02};
          else if (type == "eta_DF_OS2014")
            return {3.56e-01, 1.77e-01, 8.74e-02, 4.36e-02};
          else if (type == "eta_OS2014")
            return {1.11e+00, 5.46e-01, 2.73e-01, 1.37e-01};
          else if (type == "eff_OS2014")
            return {2.04, 2.03, 2.03, 2.04};
          else
            DUNE_THROW(Stuff::Exceptions::test_results_missing, type);
        } else
          DUNE_THROW(Stuff::Exceptions::test_results_missing, polOrder);
      } else if (this->test_case_.partitioning() == "[8 8 1]") {
        if (polOrder == 1) {
          if (type == "energy")
            return {3.29e-01, 1.63e-01, 8.05e-02, 4.02e-02};
          else if (type == "eta_NC_OS2014")
            return {1.67e-01, 7.90e-02, 3.92e-02, 1.96e-02};
          else if (type == "eta_R_OS2014")
            return {7.24e-02, 3.64e-02, 1.83e-02, 9.10e-03};
          else if (type == "eta_DF_OS2014")
            return {3.56e-01, 1.77e-01, 8.74e-02, 4.36e-02};
          else if (type == "eta_OS2014")
            return {1.11e+00, 5.46e-01, 2.73e-01, 1.37e-01};
          else if (type == "eff_OS2014")
            return {1.82, 1.81, 1.81, 1.81};
          else
            DUNE_THROW(Stuff::Exceptions::test_results_missing, type);
        } else
          DUNE_THROW(Stuff::Exceptions::test_results_missing, polOrder);
      } else
        DUNE_THROW(Stuff::Exceptions::test_results_missing, this->test_case_.partitioning());
    } else
      DUNE_THROW(Stuff::Exceptions::test_results_missing, Stuff::Common::Typename< TestCaseType >::value());
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

  virtual std::vector< std::string > available_estimators() const DS_OVERRIDE DS_FINAL
  {
    auto ret = EstimatorType::available();
    if (std::find(ret.begin(), ret.end(), "eta_OS2014") != ret.end())
      ret.push_back("eff_OS2014");
    if (std::find(ret.begin(), ret.end(), "eta_OS2014_alt") != ret.end())
      ret.push_back("eff_OS2014_alt");
    return ret;
  } // ... available_estimators(..)

  virtual double estimate(const VectorType& vector, const std::string type) const DS_OVERRIDE DS_FINAL
  {
    if (type == "eff_OS2014")
      return estimate(vector, "eta_OS2014") / const_cast< ThisType& >(*this).current_error_norm("energy");
    else if (type == "eff_OS2014_alt")
      return estimate(vector, "eta_OS2014_alt") / const_cast< ThisType& >(*this).current_error_norm("energy");
    else {
      assert(this->current_discretization_);
      return EstimatorType::estimate(*this->current_discretization_->ansatz_space(),
                                     vector,
                                     this->test_case_.problem(),
                                     type);
    }
  } // ... estimate(...)
}; // class EocStudyBlockSWIPDG


#if HAVE_ALUGRID && HAVE_DUNE_GRID_MULTISCALE && HAVE_DUNE_FEM && HAVE_EIGEN


extern template class EocStudyBlockSWIPDG< TestCases::ESV2007Multiscale< ALUGrid< 2, 2, simplex, conforming > >,
                                           1,
                                           Stuff::LA::ChooseBackend::eigen_sparse >;


#endif // HAVE_ALUGRID && HAVE_DUNE_GRID_MULTISCALE && HAVE_DUNE_FEM && HAVE_EIGEN


} // namespace Tests
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_TEST_LINEARELLIPTIC_BLOCK_SWIPDG_HH
