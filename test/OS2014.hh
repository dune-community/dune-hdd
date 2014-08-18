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
#include <dune/hdd/playground/linearelliptic/testcases/OS2014.hh>

#include "linearelliptic.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Tests {
namespace internal {


template< class TestCaseType, int polOrder, Stuff::LA::ChooseBackend la_backend >
class Discretization
{
  typedef typename TestCaseType::GridType GridType;
  typedef typename TestCaseType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange = TestCaseType::dimRange;
public:
  typedef Discretizations::BlockSWIPDG< GridType, RangeFieldType, dimRange, polOrder, la_backend > Type;
}; // class Discretization


} // namespace internal


template< class TestCaseType, int polOrder, Stuff::LA::ChooseBackend la_backend >
class OS2014EstimatorConvergenceStudy
  : public MultiscaleEocStudyBase< TestCaseType, typename internal::Discretization< TestCaseType, polOrder, la_backend >::Type >
{
  typedef MultiscaleEocStudyBase
      < TestCaseType, typename internal::Discretization< TestCaseType, polOrder, la_backend >::Type > BaseType;

  typedef typename BaseType::DiscretizationType DiscretizationType;
  typedef typename DiscretizationType::GridViewType GridViewType;
  typedef typename BaseType::FunctionType FunctionType;
  typedef typename BaseType::VectorType   VectorType;

  typedef typename BaseType::ConstDiscreteFunctionType ConstDiscreteFunctionType;

public:
  OS2014EstimatorConvergenceStudy(const TestCaseType& test_case)
    : BaseType(test_case)
  {}

  virtual ~OS2014EstimatorConvergenceStudy() {}

  virtual std::string identifier() const DS_OVERRIDE DS_FINAL
  {
    return DiscretizationType::static_id() + " (polorder " + Stuff::Common::toString(polOrder) + ")";
  }

  virtual size_t expected_rate(const std::string type) const DS_OVERRIDE DS_FINAL
  {
    if (type == "energy")
      return polOrder;
//    else if (type == "eta_NC")
//      return polOrder;
//    else if (type == "eta_R")
//      return polOrder + 1;
//    else if (type == "eta_DF")
//      return polOrder;
//    else if (type == "eta")
//      return polOrder;
//    else if (type == "efficiency")
//      return 0;
    else
      DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Wrong type '" << type << "' requested!");
  } // ... expected_rate(...)

  virtual std::vector< double > expected_results(const std::string
#if HAVE_ALUGRID
                                                                   type
#endif
                                                                       ) const DS_OVERRIDE DS_FINAL
  {
#if HAVE_ALUGRID
    if (std::is_same< TestCaseType, TestCases::OS2014Multiscale< ALUConformGrid< 2, 2 > > >::value
        || std::is_same< TestCaseType, TestCases::OS2014Multiscale< ALUGrid< 2, 2, simplex, conforming > > >::value) {
      if (polOrder == 1) {
        if (type == "energy")
          return {1, 1, 1, 1};
        else
          DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Wrong type '" << type << "' requested!");
      } else
        DUNE_THROW(NotImplemented, "Please record the expected results for this polOrder!");
    } else
#endif // HAVE_ALUGRID
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
    return {
        "eta_nc gamma"
      , "eta_r"
      , "eta_df old gamma_tilde"
      , "eta_df"
      , "eta_df alpha"
      , "eta_old"
      , "eta"
    };
  }

  virtual double estimate_(const VectorType& vector, const std::string type) const DS_OVERRIDE DS_FINAL
  {
    std::abort();
    assert(this->current_discretization_);
    return this->current_discretization_->estimate(vector, type);
  }
}; // class OS2014EstimatorConvergenceStudy


} // namespace Tests
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_TEST_LINEARELLIPTIC_BLOCK_SWIPDG_HH
