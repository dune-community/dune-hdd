// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_TESTCASES_OS2014_HH
#define DUNE_HDD_LINEARELLIPTIC_TESTCASES_OS2014_HH

#include <algorithm>
#include <cmath>

#if HAVE_ALUGRID
# include <dune/grid/alugrid.hh>
#endif

#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/common/fixed_map.hh>

#include <dune/pymor/functions/default.hh>

#include <dune/hdd/playground/linearelliptic/problems/OS2014.hh>

#include "base.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace TestCases {
namespace internal {


template< class GridType >
class OS2014Base
{
  static_assert(GridType::dimension == 2, "This test case is only available in 2d!");
public:
  typedef typename GridType::template Codim< 0 >::Entity EntityType;
  typedef typename GridType::ctype                       DomainFieldType;
  static const unsigned int                              dimDomain = GridType::dimension;
  typedef double            RangeFieldType;
  static const unsigned int dimRange = 1;
protected:
  typedef Problems::OS2014< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > ParametricProblemType;
public:
  typedef ProblemInterface< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > ProblemType;
  typedef Stuff::Functions::Constant< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >
      ExactSolutionType;
  typedef typename ProblemType::FunctionType::NonparametricType FunctionType;

  typedef Stuff::Common::FixedMap< std::string, Pymor::Parameter, 3 > ParametersType;

protected:
  static ParametersType default_parameters()
  {
    return ParametersType({{"current", Pymor::Parameter("mu", 0.5)},
                           {"fixed", Pymor::Parameter("mu", 1)},
                           {"minimizing", Pymor::Parameter("mu", 0.1)}});
  }

  static int initial_refinements()
  {
    int ret = 1;
#if HAVE_ALUGRID
    if (std::is_same< GridType, ALUConformGrid< 2, 2 > >::value
        || std::is_same< GridType, ALUGrid< 2, 2, simplex, conforming > >::value)
      ret += 1;
#endif // HAVE_ALUGRID
    return ret;
  } // ... initial_refinements()

public:
  OS2014Base(ParametersType parameters = default_parameters())
    : parameters_(parameters)
    , boundary_info_cfg_(Stuff::Grid::BoundaryInfoConfigs::AllDirichlet::default_config())
    , parametric_problem_(3)
    , problem_(parametric_problem_.with_mu(parameters_["current"]))
    , exact_solution_(0)
  {
    assert(parameters_.find("current") != parameters_.end());
    assert(parameters_.find("fixed") != parameters_.end());
    assert(parameters_.find("minimizing") != parameters_.end());
    for (auto parameter : parameters_)
      assert(parameter.second.type() == parametric_problem_.parameter_type());
    const auto& diffusion_factor = *parametric_problem_.diffusion_factor();
    alpha_ = diffusion_factor.alpha(parameters_["current"], parameters_["fixed"]);
    gamma_ = diffusion_factor.gamma(parameters_["current"], parameters_["fixed"]);
    gamma_tilde_sqrt_ = std::max(std::sqrt(gamma_), 1.0/std::sqrt(alpha_));
  }

  void print_header(std::ostream& out = std::cout) const
  {
    out << "+==================================================================+\n"
        << "|+================================================================+|\n"
        << "||  Testcase OS2014: (parametric) estimator convergence study     ||\n"
        << "||  (see Ohlberger, Schindler, 2014)                              ||\n"
        << "|+----------------------------------------------------------------+|\n"
        << "||  domain = [-1, 1] x [-1, 1]                                    ||\n"
        << "||  diffusion = 1 + (1 - mu) 3/4 sin(4 pi (x + 1/2 y))            ||\n"
        << "||  force     = 1/2 pi^2 cos(1/2 pi x) cos(1/2 pi y)              ||\n"
        << "||  dirichlet = 0                                                 ||\n"
        << "||  reference solution: discrete solution on finest grid          ||\n"
        << "||  parameter: mu in [0.1, 1]                                     ||\n"
        << "|+================================================================+|\n"
        << "+==================================================================+\n"
        << "    solving for mu: " << parameters_["current"] << "\n"
        << "    fixed mu:       " << parameters_["fixed"] << "\n"
        << "    sqrt(alpha(mu, mu_fixed))^-1:    "
            << std::setprecision(2) << std::scientific << 1.0/std::sqrt(alpha_) << "\n"
        << "    sqrt(gamma(mu, mu_fixed)):       "
            << std::setprecision(2) << std::scientific << std::sqrt(gamma_) << "\n"
        << "    sqrt(gamma_tilde(mu, mu_fixed)): "
            << std::setprecision(2) << std::scientific << gamma_tilde_sqrt_ << "\n"
        << "====================================================================" << std::endl;
  } // ... print_header(...)

  const Stuff::Common::Configuration& boundary_info() const
  {
    return boundary_info_cfg_;
  }

  const ProblemType& parametric_problem() const
  {
    return parametric_problem_;
  }

  const ProblemType& problem() const
  {
    return *problem_;
  }

  bool provides_exact_solution() const
  {
    return false;
  }

  const ExactSolutionType& exact_solution() const
  {
    DUNE_THROW(Stuff::Exceptions::you_are_using_this_wrong,
               "Do not call exact_solution() if provides_exact_solution() is false!");
  }

  const ParametersType& parameters() const
  {
    return parameters_;
  }

protected:
  const ParametersType parameters_;
  const Stuff::Common::Configuration boundary_info_cfg_;
  const ParametricProblemType parametric_problem_;
  const std::shared_ptr< ProblemType > problem_;
  const ExactSolutionType exact_solution_;
  double alpha_;
  double gamma_;
  double gamma_tilde_sqrt_;
}; // class OS2014Base


} // namespace internal


template< class GridType >
class OS2014
  : public internal::OS2014Base< GridType >
  , public Base< GridType >
{
  typedef internal::OS2014Base< GridType > OS2014BaseType;
  typedef Base< GridType >                 TestCaseBaseType;

  static std::shared_ptr< GridType > create_initial_grid(const int refinements)
  {
    Stuff::Grid::Providers::Cube< GridType > grid_provider(-1, 1, 4);
    auto grid = grid_provider.grid_ptr();
    grid->globalRefine(refinements);
    return grid;
  } // ... create_initial_grid(...)

public:
  typedef typename OS2014BaseType::ParametersType ParametersType;

  OS2014(ParametersType parameters = OS2014BaseType::default_parameters(), const size_t num_refinements = 3)
    : OS2014BaseType(parameters)
    , TestCaseBaseType(create_initial_grid(OS2014BaseType::initial_refinements()), num_refinements)
  {}
}; // class OS2014


//#if HAVE_DUNE_GRID_MULTISCALE


//template< class GridType >
//class OS2014Multiscale
//  : public MultiscaleCubeBase< GridType >
//{
//  typedef MultiscaleCubeBase< GridType > BaseType;
//  static_assert(BaseType::dimDomain == 2, "This test case is only available in 2d!");
//public:
//  typedef typename BaseType::EntityType EntityType;
//  typedef typename BaseType::DomainFieldType DomainFieldType;
//  static const unsigned int                  dimDomain = BaseType::dimDomain;
//  typedef double            RangeFieldType;
//  static const unsigned int dimRange = 1;

//  typedef ProblemInterface< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > ProblemType;
//  typedef Stuff::Functions::Constant< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >
//      ExactSolutionType;
//  typedef typename ProblemType::FunctionType::NonparametricType FunctionType;
//  typedef Stuff::Common::FixedMap< std::string, Pymor::Parameter, 3 > ParametersType;

//private:
//  typedef Problems::OS2014< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > RealProblemType;

//public:
//  OS2014Multiscale(ParametersType parameters = ParametersType({{"fixed", Pymor::Parameter("mu", 1)},
//                                                               {"solve", Pymor::Parameter("mu", 1)},
//                                                               {"minimizing", Pymor::Parameter("mu", 0.1)}}),
//                   const std::string num_partitions = "[1 1 1]",
//                   const size_t num_refinements = 3)
//    : BaseType(initial_grid_cfg(num_partitions), initial_refinements(), num_refinements)
//    , parameters_(parameters)
//    , boundary_info_cfg_(Stuff::Grid::BoundaryInfoConfigs::AllDirichlet::default_config())
//    , real_problem_(3)
//    , problem_(real_problem_.with_mu(parameters_["solve"]))
//    , exact_solution_(0)
//  {
//    for (auto parameter : parameters_)
//      assert(parameter.second.type() == real_problem_.parameter_type());
//    assert(parameters_.find("fixed") != parameters_.end());
//    assert(parameters_.find("solve") != parameters_.end());
//    assert(parameters_.find("minimizing") != parameters_.end());
//    const auto& diffusion_factor = *real_problem_.diffusion_factor();
//    alpha_ = diffusion_factor.alpha(parameters_["solve"], parameters_["fixed"]);
//    gamma_ = diffusion_factor.gamma(parameters_["solve"], parameters_["fixed"]);
//    gamma_tilde_sqrt_ = std::max(std::sqrt(gamma_), 1.0/std::sqrt(alpha_));
//  }

//  void print_header(std::ostream& out = std::cout) const
//  {
//    out << "+==================================================================+\n"
//        << "|+================================================================+|\n"
//        << "||  Testcase OS2014: (parametric) estimator convergence study     ||\n"
//        << "||  (see Ohlberger, Schindler, 2014)                              ||\n"
//        << "|+----------------------------------------------------------------+|\n"
//        << "||  domain = [-1, 1] x [-1, 1]                                    ||\n"
//        << "||  diffusion = 1 + (1 - mu) sin(2 pi x) sin(2 pi y)              ||\n"
//        << "||  force     = 1                                                 ||\n"
//        << "||  dirichlet = 0                                                 ||\n"
//        << "||  reference solution: discrete solution on finest grid          ||\n"
//        << "||  parameter: mu in [0.1, 1]                                     ||\n"
//        << "|+================================================================+|\n"
//        << "+==================================================================+\n"
//        << "    parameter type: " << real_problem_.parameter_type() << "\n"
//        << "    fixed mu:       " << parameters_["fixed"] << "\n"
//        << "    solving for mu: " << parameters_["solve"] << "\n"
//        << "    sqrt(alpha(mu, mu_fixed))^-1:    "
//            << std::setprecision(2) << std::scientific << 1.0/std::sqrt(alpha_) << "\n"
//        << "    sqrt(gamma(mu, mu_fixed)):       "
//            << std::setprecision(2) << std::scientific << std::sqrt(gamma_) << "\n"
//        << "    sqrt(gamma_tilde(mu, mu_fixed)): "
//            << std::setprecision(2) << std::scientific << gamma_tilde_sqrt_  << std::endl;
//  } // ... print_header(...)

//  const Stuff::Common::Configuration& boundary_info() const
//  {
//    return boundary_info_cfg_;
//  }

//  const ProblemType& parametric_problem() const
//  {
//    return real_problem_;
//  }

//  const ProblemType& problem() const
//  {
//    return *problem_;
//  }

//  bool provides_exact_solution() const
//  {
//    return false;
//  }

//  const ExactSolutionType& exact_solution() const
//  {
//    DUNE_THROW(Stuff::Exceptions::you_are_using_this_wrong,
//               "Do not call exact_solution() if provides_exact_solution() is false!");
//  }

//  const ParametersType& parameters() const
//  {
//    return parameters_;
//  }

//private:
//  static Stuff::Common::Configuration initial_grid_cfg(const std::string num_partitions)
//  {
//    Stuff::Common::Configuration grid_cfg = Stuff::Grid::Providers::Cube< GridType >::default_config();
//    grid_cfg["lower_left"] = "-1";
//    grid_cfg["upper_right"] = "1";
//    grid_cfg["num_elements"] = "4";
//    grid_cfg["num_partitions"] = num_partitions;
//    return grid_cfg;
//  } // ... initial_grid_cfg(...)

//  static int initial_refinements()
//  {
//    int ret = 1;
//#if HAVE_ALUGRID
//    if (std::is_same< GridType, ALUConformGrid< 2, 2 > >::value
//        || std::is_same< GridType, ALUGrid< 2, 2, simplex, conforming > >::value)
//      ret += 1;
//#endif // HAVE_ALUGRID
//    return ret;
//  } // ... initial_refinements()

//  const ParametersType parameters_;
//  const Stuff::Common::Configuration boundary_info_cfg_;
//  const RealProblemType real_problem_;
//  const std::shared_ptr< ProblemType > problem_;
//  const ExactSolutionType exact_solution_;
//  double alpha_;
//  double gamma_;
//  double gamma_tilde_sqrt_;
//}; // class OS2014Multiscale


//#endif // HAVE_DUNE_GRID_MULTISCALE

} // namespace TestCases
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_TESTCASES_OS2014_HH
