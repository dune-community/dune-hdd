// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_TEST_LINEARELLIPTIC_BLOCK_SWIPDG_HH
#define DUNE_HDD_TEST_LINEARELLIPTIC_BLOCK_SWIPDG_HH

#include <dune/stuff/common/exceptions.hh>

#include <dune/gdt/products/l2.hh>
#include <dune/gdt/products/h1.hh>
#include <dune/gdt/playground/products/elliptic.hh>

#include <dune/hdd/playground/linearelliptic/discretizations/block-swipdg.hh>
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
}; // class DiscretizationBlockSWIPDG


} // namespace internal


template< class TestCaseType, int polOrder, Stuff::LA::ChooseBackend la_backend >
class EocStudyBlockSWIPDG
  : public MultiscaleEocStudyBase< TestCaseType, typename internal::DiscretizationBlockSWIPDG< TestCaseType, polOrder, la_backend >::Type >
{
  typedef MultiscaleEocStudyBase
      < TestCaseType, typename internal::DiscretizationBlockSWIPDG< TestCaseType, polOrder, la_backend >::Type > BaseType;

  typedef typename BaseType::DiscretizationType DiscretizationType;
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
    else if (type == "energy")
      return polOrder;
    else if (type == "eta_NC")
      return polOrder;
    else if (type == "eta_R")
      return polOrder + 1;
    else if (type == "eta_DF")
      return polOrder;
    else if (type == "eta")
      return polOrder;
    else if (type == "efficiency")
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
          else if (type == "eta_NC")
            return {1.67e-01, 7.90e-02, 3.92e-02, 1.96e-02};
          else if (type == "eta_R")
            return {5.80e-01, 2.91e-01, 1.46e-01, 7.28e-02};
          else if (type == "eta_DF") {
            return {3.56e-01, 1.77e-01, 8.74e-02, 4.36e-02};
          } else if (type == "eta")
            return {1.11e+00, 5.46e-01, 2.73e-01, 1.37e-01};
          else
            DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Wrong type '" << type << "' requested!");
        } else
          DUNE_THROW(NotImplemented, "Please record the expected results for this polOrder!");
      } else if (this->test_case_.partitioning() == "[2 2 1]") {
        if (polOrder == 1) {
          if (type == "energy")
            return {3.29e-01, 1.63e-01, 8.05e-02, 4.02e-02};
          else if (type == "eta_NC")
            return {1.67e-01, 7.90e-02, 3.92e-02, 1.96e-02};
          else if (type == "eta_R")
            return {2.90e-01, 1.46e-01, 7.28e-02, 3.64e-02};
          else if (type == "eta_DF") {
            return {3.56e-01, 1.77e-01, 8.74e-02, 4.36e-02};
          } else if (type == "eta")
            return {1.11e+00, 5.46e-01, 2.73e-01, 1.37e-01};
          else
            DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Wrong type '" << type << "' requested!");
        } else
          DUNE_THROW(NotImplemented, "Please record the expected results for this polOrder!");
      } else if (this->test_case_.partitioning() == "[4 4 1]") {
        if (polOrder == 1) {
          if (type == "energy")
            return {3.29e-01, 1.63e-01, 8.05e-02, 4.02e-02};
          else if (type == "eta_NC")
            return {1.67e-01, 7.90e-02, 3.92e-02, 1.96e-02};
          else if (type == "eta_R")
            return {1.46e-01, 7.27e-02, 3.64e-02, 1.82e-02};
          else if (type == "eta_DF") {
            return {3.56e-01, 1.77e-01, 8.74e-02, 4.36e-02};
          } else if (type == "eta")
            return {1.11e+00, 5.46e-01, 2.73e-01, 1.37e-01};
          else
            DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Wrong type '" << type << "' requested!");
        } else
          DUNE_THROW(NotImplemented, "Please record the expected results for this polOrder!");
      } else if (this->test_case_.partitioning() == "[8 8 1]") {
        if (polOrder == 1) {
          if (type == "energy")
            return {3.29e-01, 1.63e-01, 8.05e-02, 4.02e-02};
          else if (type == "eta_NC")
            return {1.67e-01, 7.90e-02, 3.92e-02, 1.96e-02};
          else if (type == "eta_R")
            return {7.24e-02, 3.64e-02, 1.83e-02, 9.10e-03};
          else if (type == "eta_DF") {
            return {3.56e-01, 1.77e-01, 8.74e-02, 4.36e-02};
          } else if (type == "eta")
            return {1.11e+00, 5.46e-01, 2.73e-01, 1.37e-01};
          else
            DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Wrong type '" << type << "' requested!");
        } else
          DUNE_THROW(NotImplemented, "Please record the expected results for this polOrder!");
      } else
        DUNE_THROW(NotImplemented, "Please record the expected results for this partitioning!");
    } else
      DUNE_THROW(NotImplemented, "Please record the expected results for this TestCaseType/GridType combination!");
  } // ... expected_results(...)

private:
  virtual std::vector< std::string > available_norms_() const DS_OVERRIDE DS_FINAL
  {
    return {
      "energy"
    };
  } // ... available_norms_(...)

  virtual double compute_norm_(const GridViewType& grid_view,
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
  } // ... compute_norm_(...)

  virtual std::vector< std::string > available_estimators_() const DS_OVERRIDE DS_FINAL
  {
    return DiscretizationType::available_estimators();
  }

  virtual double estimate_(const VectorType& vector, const std::string type) const DS_OVERRIDE DS_FINAL
  {
    assert(this->current_discretization_);
    return this->current_discretization_->estimate(vector, type);
  }
}; // class EocStudyBlockSWIPDG


} // namespace Tests
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_TEST_LINEARELLIPTIC_BLOCK_SWIPDG_HH
