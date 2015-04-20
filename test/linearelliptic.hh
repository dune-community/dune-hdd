// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_test_case_LINEARELLIPTIC_EOC_STUDY_HH
#define DUNE_HDD_test_case_LINEARELLIPTIC_EOC_STUDY_HH

#include <memory>
#include <vector>
#include <limits>
#include <algorithm>

#include <dune/common/timer.hh>

#if HAVE_DUNE_FEM
# include <dune/fem/misc/gridwidth.hh>
#else
# error "We need dune-fem until we have a reliable and fast method to compute the grid width!"
#endif

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/common/convergence-study.hh>
#include <dune/stuff/grid/layers.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/prolongations.hh>

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Tests {


/**
 *  Assumes that the DiscretizationImp is derived from Discretizations::ContainerBasedDefault and has a ctor of type
 *  DiscretizationImp(Stuff::Grid::ConstProviderInterface, BoundaryInfoConfig, ProblemType, level).
 *  Assumes that TestCaseType is derived from TestCases::Base and additionally provides
 *  * const Stuff::Common::Configuration& boundary_info() const
 *  * const ProblemType& problem() const
 *  * bool provides_exact_solution() const
 *  * const ExactSolutionType& exact_solution() const
 *  and the appropriate types.
 *
 *  Any derived class only has to provide identifier(), provided_norms(), expected_rate(), compute_norm() and
 *  expected_results().
 */
template< class TestCaseType, class DiscretizationImp >
class EocStudyBase
  : public Stuff::Common::ConvergenceStudy
{
  typedef Stuff::Common::ConvergenceStudy BaseType;
protected:

  typedef DiscretizationImp DiscretizationType;
  typedef typename DiscretizationType::VectorType VectorType;
  typedef typename DiscretizationType::AnsatzSpaceType SpaceType;
  typedef GDT::DiscreteFunction< SpaceType, VectorType >      DiscreteFunctionType;
  typedef GDT::ConstDiscreteFunction< SpaceType, VectorType > ConstDiscreteFunctionType;

  typedef typename TestCaseType::ExactSolutionType ExactSolutionType;
  typedef typename TestCaseType::template Level< Stuff::Grid::ChoosePartView::view >::Type GridViewType;
  typedef typename TestCaseType::FunctionType FunctionType;

public:
  EocStudyBase(TestCaseType& test_case,
               const std::vector< std::string > only_these_norms = {},
               const std::string visualize_prefix = "")
    : BaseType(only_these_norms)
    , test_case_(test_case)
    , current_refinement_(0)
    , last_computed_refinement_(std::numeric_limits< size_t >::max())
    , time_to_solution_(0)
    , reference_solution_computed_(false)
    , current_discretization_(nullptr)
    , current_solution_vector_on_level_(nullptr)
    , reference_discretization_(nullptr)
    , reference_solution_vector_(nullptr)
    , current_solution_vector_(nullptr)
    , visualize_prefix_(visualize_prefix)
  {
    if (test_case_.problem().parametric())
      DUNE_THROW(NotImplemented, "Parametric problems are not implemented yet!");
  }

  virtual ~EocStudyBase() {}

  virtual size_t num_refinements() const override final
  {
    return test_case_.num_refinements();
  }

  virtual std::vector< std::string > provided_norms() const override final
  {
    std::vector< std::string > ret = available_norms();
    for (auto estimator : available_estimators()) {
      if (is_norm(estimator))
        DUNE_THROW(Stuff::Exceptions::internal_error,
                   "We do not want to handle the case that norms and estimators have the same name!");
      ret.push_back(estimator);
    }
    return ret;
  } // ... provided_norms(...)

  virtual double norm_reference_solution(const std::string type) override final
  {
    if (is_norm(type)) {
      if (test_case_.provides_exact_solution()) {
        // visualize
        if (!visualize_prefix_.empty()) {
          test_case_.exact_solution().visualize(test_case_.reference_grid_view(),
                                                visualize_prefix_ + "_exact_solution");
        }
        return compute_norm(test_case_.reference_grid_view(), test_case_.exact_solution(), type);
      } else {
        compute_reference_solution();
        assert(reference_discretization_);
        assert(reference_solution_vector_);
        const ConstDiscreteFunctionType reference_solution(*(reference_discretization_->ansatz_space()),
                                                           *reference_solution_vector_,
                                                           "reference solution");
        return compute_norm(test_case_.reference_grid_view(), reference_solution, type);
      }
    } else
      return 1.0;
  } // ... norm_reference_solution(...)

  virtual size_t current_grid_size() const override final
  {
    assert(current_refinement_ <= num_refinements());
    const int level = test_case_.level_of(current_refinement_);
    return test_case_.grid().size(level, 0);
  } // ... current_grid_size(...)

  virtual double current_grid_width() const override final
  {
    assert(current_refinement_ <= num_refinements());
    const int level = test_case_.level_of(current_refinement_);
    const auto grid_part = test_case_.template level< Stuff::Grid::ChoosePartView::part >(level);
    return Fem::GridWidth::calcGridWidth(grid_part);
  } // ... current_grid_width(...)

  virtual double compute_on_current_refinement() override final
  {
    using namespace Dune;
    using namespace Dune::GDT;
    if (current_refinement_ != last_computed_refinement_) {
      assert(current_refinement_ <= num_refinements());
      // compute solution
      Timer timer;
      current_discretization_
          = Stuff::Common::make_unique< DiscretizationType >(test_case_,
                                                             test_case_.boundary_info(),
                                                             test_case_.problem(),
                                                             test_case_.level_of(current_refinement_));
      current_discretization_->init();
      current_solution_vector_on_level_
          = Stuff::Common::make_unique< VectorType >(current_discretization_->create_vector());
      current_discretization_->solve(*current_solution_vector_on_level_);
      time_to_solution_ = timer.elapsed();
      const ConstDiscreteFunctionType current_refinement_solution(*current_discretization_->ansatz_space(),
                                                                  *current_solution_vector_on_level_,
                                                                  "solution on current level");
      // prolong to reference grid part
      if (!reference_solution_computed_)
        compute_reference_solution();
      auto reference_grid_view = test_case_.reference_grid_view();
      const Operators::Prolongation< GridViewType > prolongation_operator(reference_grid_view);
      assert(reference_discretization_);
      if (!current_solution_vector_)
        current_solution_vector_ = Stuff::Common::make_unique< VectorType >(reference_discretization_->create_vector());
      DiscreteFunctionType reference_refinement_solution(*(reference_discretization_->ansatz_space()),
                                                         *current_solution_vector_,
                                                         "solution on reference grid part");
      prolongation_operator.apply(current_refinement_solution, reference_refinement_solution);
      last_computed_refinement_ = current_refinement_;
      // visualize
      if (!visualize_prefix_.empty()) {
        this->test_case_.problem().visualize(current_discretization_->grid_view(),
                                             visualize_prefix_ + "_problem_" + DSC::toString(current_refinement_));
        current_refinement_solution.visualize(visualize_prefix_ + "_solution_" + DSC::toString(current_refinement_));
      }
    }
    return time_to_solution_;
  } // ... compute_on_current_refinement(...)

  virtual double current_error_norm(const std::string type) override final
  {
    // get current solution
    assert(current_refinement_ <= num_refinements());
    if (last_computed_refinement_ != current_refinement_) {
      compute_on_current_refinement();
    }
    assert(last_computed_refinement_ == current_refinement_);
    if (is_norm(type)) {
      assert(current_solution_vector_);
      if (!reference_solution_computed_)
        compute_reference_solution();
      assert(reference_discretization_);
      const ConstDiscreteFunctionType current_solution(*(reference_discretization_->ansatz_space()),
                                                       *current_solution_vector_,
                                                       "current solution");
      // compute error
      if (test_case_.provides_exact_solution()) {
        return compute_norm(test_case_.reference_grid_view(), test_case_.exact_solution() - current_solution, type);
      } else {
        // get reference solution
        compute_reference_solution();
        assert(reference_discretization_);
        assert(reference_solution_vector_);
        const ConstDiscreteFunctionType reference_solution(*(reference_discretization_->ansatz_space()),
                                                           *reference_solution_vector_,
                                                           "reference solution");
        return compute_norm(test_case_.reference_grid_view(), reference_solution- current_solution, type);
      }
    } else {
      assert(current_solution_vector_on_level_);
      return estimate(*current_solution_vector_on_level_, type);
    }
  } // ... current_error_norm(...)

  virtual void refine() override final
  {
    if (current_refinement_ <= num_refinements())
      ++current_refinement_;
  } // ... refine()

  virtual std::vector< double > expected_results(const std::string type) const override = 0;

  std::map< std::string, std::vector< double > > run(std::ostream& out, const bool print_timings = false)
  {
    return BaseType::run(false, out, print_timings);
  }

protected:
  void compute_reference_solution()
  {
    if (!reference_solution_computed_) {
      reference_discretization_
          = Stuff::Common::make_unique< DiscretizationType >(test_case_,
                                                             test_case_.boundary_info(),
                                                             test_case_.problem(),
                                                             test_case_.reference_level());
      reference_discretization_->init();
      reference_solution_vector_ = Stuff::Common::make_unique< VectorType >(reference_discretization_->create_vector());
      reference_discretization_->solve(*reference_solution_vector_);
      reference_solution_computed_ = true;
      // visualize
      if (!visualize_prefix_.empty()) {
        this->test_case_.problem().visualize(reference_discretization_->grid_view(),
                                             visualize_prefix_ + "_problem_reference");
        ConstDiscreteFunctionType(*reference_discretization_->ansatz_space(),
                                  *reference_solution_vector_,
                                  "reference solution").visualize(visualize_prefix_ + "_reference_solution");
      }
    }
  } // ... compute_reference_solution()

  bool is_norm(const std::string type) const
  {
    const auto norms = available_norms();
    return std::find(norms.begin(), norms.end(), type) != norms.end();
  } // ... is_norm(...)

  virtual std::vector< std::string > available_norms() const = 0;

  virtual std::vector< std::string > available_estimators() const = 0;

  virtual double estimate(const VectorType& vector, const std::string type) const = 0;

  virtual double compute_norm(const GridViewType& grid_view,
                              const FunctionType& function,
                              const std::string type) const = 0;

  TestCaseType& test_case_;
  size_t current_refinement_;
  size_t last_computed_refinement_;
  double time_to_solution_;
  bool reference_solution_computed_;
  std::unique_ptr< DiscretizationType > current_discretization_;
  std::unique_ptr< VectorType > current_solution_vector_on_level_;
  std::unique_ptr< DiscretizationType > reference_discretization_;
  std::unique_ptr< VectorType > reference_solution_vector_;
  std::unique_ptr< VectorType > current_solution_vector_;
  const std::string visualize_prefix_;
}; // class EocStudyBase


#if HAVE_DUNE_GRID_MULTISCALE


template< class MultiscaleTestCaseType, class DiscretizationImp >
class MultiscaleEocStudyBase
  : public Stuff::Common::ConvergenceStudy
{
  typedef Stuff::Common::ConvergenceStudy BaseType;
protected:

  typedef DiscretizationImp DiscretizationType;
  typedef typename DiscretizationType::VectorType VectorType;
  typedef typename DiscretizationType::AnsatzSpaceType SpaceType;
  typedef GDT::DiscreteFunction< SpaceType, VectorType >      DiscreteFunctionType;
  typedef GDT::ConstDiscreteFunction< SpaceType, VectorType > ConstDiscreteFunctionType;

  typedef typename MultiscaleTestCaseType::ExactSolutionType ExactSolutionType;
  typedef typename DiscretizationType::GridViewType GridViewType;
  typedef typename MultiscaleTestCaseType::FunctionType FunctionType;
  typedef typename MultiscaleTestCaseType::ProblemType::NonparametricType NonparametricProblemType;

  static const MultiscaleTestCaseType& check_test_case(const MultiscaleTestCaseType& test_case)
  {
    if (test_case.parametric() && test_case.parameters().find("mu") == test_case.parameters().end())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given,
                 "The parameter mu has to be provided for a parametric testcase!");
    return test_case;
  } // ... check_test_case(...)

  static Pymor::Parameter get_mu(const MultiscaleTestCaseType& test_case)
  {
    if (test_case.parametric())
      return test_case.parameters().at("mu");
    else
      return Pymor::Parameter();
  } // ... get_mu(...)

public:
  MultiscaleEocStudyBase(const MultiscaleTestCaseType& test_case,
                         const std::vector< std::string > only_these_norms = std::vector< std::string >(),
                         const std::string visualize_prefix = "")
    : BaseType(only_these_norms)
    , test_case_(check_test_case(test_case))
    , current_refinement_(0)
    , last_computed_refinement_(std::numeric_limits< size_t >::max())
    , time_to_solution_(0)
    , nonparametric_problem_(test_case_.problem().with_mu(get_mu(test_case_)))
    , reference_solution_computed_(false)
    , current_discretization_(nullptr)
    , current_solution_vector_on_level_(nullptr)
    , reference_discretization_(nullptr)
    , reference_solution_vector_(nullptr)
    , current_solution_vector_(nullptr)
    , visualize_prefix_(visualize_prefix)
  {}

  virtual ~MultiscaleEocStudyBase() {}

  virtual size_t num_refinements() const override final
  {
    return test_case_.num_refinements();
  }

  virtual std::vector< std::string > provided_norms() const override final
  {
    std::vector< std::string > ret = available_norms();
    for (auto estimator : available_estimators()) {
      if (is_norm(estimator))
        DUNE_THROW(Stuff::Exceptions::internal_error,
                   "We do not want to handle the case that norms and estimators have the same name!");
      ret.push_back(estimator);
    }
    return ret;
  } // ... provided_norms(...)

  virtual double norm_reference_solution(const std::string type) override final
  {
    if (is_norm(type)) {
      const auto reference_grid_view
          = test_case_.reference_provider()->template global< Stuff::Grid::ChoosePartView::view >();
      if (test_case_.provides_exact_solution()) {
        // visualize
        if (!visualize_prefix_.empty()) {
          test_case_.exact_solution().visualize(reference_grid_view,
                                                visualize_prefix_ + "_exact_solution");
        }
        return compute_norm(reference_grid_view, test_case_.exact_solution(), type);
      } else {
        compute_reference_solution();
        assert(reference_discretization_);
        assert(reference_solution_vector_);
        const ConstDiscreteFunctionType reference_solution(*(reference_discretization_->ansatz_space()),
                                                           *reference_solution_vector_,
                                                           "reference solution");
        return compute_norm(reference_grid_view, reference_solution, type);
      }
    } else
      return 1.0;
  } // ... norm_reference_solution(...)

  virtual size_t current_grid_size() const override final
  {
    assert(current_refinement_ <= num_refinements());
    const auto grid_part
        = test_case_.level_provider(current_refinement_)->template global< Stuff::Grid::ChoosePartView::part >();
    return grid_part.gridView().size(0);
  } // ... current_grid_size(...)

  virtual double current_grid_width() const override final
  {
    assert(current_refinement_ <= num_refinements());
    const auto grid_part
        = test_case_.level_provider(current_refinement_)->template global< Stuff::Grid::ChoosePartView::part >();
    return Fem::GridWidth::calcGridWidth(grid_part);
  } // ... current_grid_width(...)

  virtual double compute_on_current_refinement() override final
  {
    using namespace Dune;
    using namespace Dune::GDT;
    if (current_refinement_ != last_computed_refinement_) {
      assert(current_refinement_ <= num_refinements());
      // compute solution
      Timer timer;
      current_discretization_
          = Stuff::Common::make_unique< DiscretizationType >(*(test_case_.level_provider(current_refinement_)),
                                                             test_case_.boundary_info(),
                                                             *nonparametric_problem_);
      current_discretization_->init();
      current_solution_vector_on_level_
          = Stuff::Common::make_unique< VectorType >(current_discretization_->create_vector());
      current_discretization_->solve(*current_solution_vector_on_level_);
      time_to_solution_ = timer.elapsed();
      const ConstDiscreteFunctionType current_refinement_solution(*current_discretization_->ansatz_space(),
                                                                  *current_solution_vector_on_level_,
                                                                  "solution on current level");
      // prolong to reference grid part
      if (!reference_solution_computed_)
        compute_reference_solution();
      assert(reference_discretization_);
      const auto reference_grid_view = test_case_.reference_provider()->template global< Stuff::Grid::ChoosePartView::view >();
      const Operators::Prolongation< GridViewType > prolongation_operator(reference_grid_view);
      if (!current_solution_vector_)
        current_solution_vector_ = Stuff::Common::make_unique< VectorType >(reference_discretization_->create_vector());
      DiscreteFunctionType reference_refinement_solution(*(reference_discretization_->ansatz_space()),
                                                         *current_solution_vector_,
                                                         "solution on reference grid part");
      prolongation_operator.apply(current_refinement_solution, reference_refinement_solution);
      last_computed_refinement_ = current_refinement_;
      // visualize
      if (!visualize_prefix_.empty()) {
        this->test_case_.problem().visualize(current_discretization_->grid_view(),
                                             visualize_prefix_ + "_problem_" + DSC::toString(current_refinement_));
        current_refinement_solution.visualize(visualize_prefix_ + "_solution_" + DSC::toString(current_refinement_));
      }
    }
    return time_to_solution_;
  } // ... compute_on_current_refinement(...)

  virtual double current_error_norm(const std::string type) override final
  {
    // get current solution
    assert(current_refinement_ <= num_refinements());
    if (last_computed_refinement_ != current_refinement_) {
      compute_on_current_refinement();
    }
    assert(last_computed_refinement_ == current_refinement_);
    if (is_norm(type)) {
      assert(current_solution_vector_);
      if (!reference_solution_computed_)
        compute_reference_solution();
      assert(reference_discretization_);
      const ConstDiscreteFunctionType current_solution(*(reference_discretization_->ansatz_space()),
                                                       *current_solution_vector_,
                                                       "current solution");
      // compute error
      const auto reference_grid_view
          = test_case_.reference_provider()-> template global< Stuff::Grid::ChoosePartView::view >();
      if (test_case_.provides_exact_solution()) {
        return compute_norm(reference_grid_view, test_case_.exact_solution() - current_solution, type);
      } else {
        // get reference solution
        compute_reference_solution();
        assert(reference_discretization_);
        assert(reference_solution_vector_);
        const ConstDiscreteFunctionType reference_solution(*(reference_discretization_->ansatz_space()),
                                                           *reference_solution_vector_,
                                                           "reference solution");
        return compute_norm(reference_grid_view, reference_solution- current_solution, type);
      }
    } else {
      assert(current_solution_vector_on_level_);
      return estimate(*current_solution_vector_on_level_, type);
    }
  } // ... current_error_norm(...)

  virtual void refine() override final
  {
    if (current_refinement_ <= num_refinements())
      ++current_refinement_;
  } // ... refine()

  virtual std::vector< double > expected_results(const std::string type) const override = 0;

  std::map< std::string, std::vector< double > > run(std::ostream& out, const bool print_timings = false)
  {
    return BaseType::run(false, out, print_timings);
  }

protected:
  void compute_reference_solution()
  {
    if (!reference_solution_computed_) {
      reference_discretization_
          = Stuff::Common::make_unique< DiscretizationType >(*(test_case_.reference_provider()),
                                                             test_case_.boundary_info(),
                                                             *nonparametric_problem_);
      reference_discretization_->init();
      reference_solution_vector_ = Stuff::Common::make_unique< VectorType >(reference_discretization_->create_vector());
      reference_discretization_->solve(*reference_solution_vector_);
      reference_solution_computed_ = true;
      // visualize
      if (!visualize_prefix_.empty()) {
        this->test_case_.problem().visualize(reference_discretization_->grid_view(),
                                             visualize_prefix_ + "_problem_reference");
        ConstDiscreteFunctionType(*reference_discretization_->ansatz_space(),
                                  *reference_solution_vector_,
                                  "reference solution").visualize(visualize_prefix_ + "_reference_solution");
      }
    }
  } // ... compute_reference_solution()

  bool is_norm(const std::string type) const
  {
    const auto norms = available_norms();
    return std::find(norms.begin(), norms.end(), type) != norms.end();
  } // ... is_norm(...)

  virtual std::vector< std::string > available_norms() const = 0;

  virtual std::vector< std::string > available_estimators() const = 0;

  virtual double estimate(const VectorType& vector, const std::string type) const = 0;

  virtual double compute_norm(const GridViewType& grid_view,
                              const FunctionType& function,
                              const std::string type) const = 0;

  const MultiscaleTestCaseType& test_case_;
  size_t current_refinement_;
  size_t last_computed_refinement_;
  double time_to_solution_;
  std::shared_ptr< NonparametricProblemType > nonparametric_problem_;
  bool reference_solution_computed_;
  std::unique_ptr< DiscretizationType > current_discretization_;
  std::unique_ptr< VectorType > current_solution_vector_on_level_;
  std::unique_ptr< DiscretizationType > reference_discretization_;
  std::unique_ptr< VectorType > reference_solution_vector_;
  std::unique_ptr< VectorType > current_solution_vector_;
  const std::string visualize_prefix_;
}; // class MultiscaleEocStudyBase


#endif // HAVE_DUNE_GRID_MULTISCALE

} // namespace Tests
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_test_case_LINEARELLIPTIC_EOC_STUDY_HH
