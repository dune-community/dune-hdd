// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_TESTCASES_OS2014_HH
#define DUNE_HDD_LINEARELLIPTIC_TESTCASES_OS2014_HH

#include <algorithm>
#include <cmath>
#include <map>

#include <dune/stuff/common/disable_warnings.hh>
# if HAVE_ALUGRID
#   include <dune/grid/alugrid.hh>
# endif
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/stuff/functions/constant.hh>

#include <dune/pymor/functions/default.hh>

#include <dune/hdd/playground/linearelliptic/problems/OS2014.hh>

#include "base.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace TestCases {
namespace OS2014 {
namespace internal {


template< class GridType >
class ParametricConvergenceBase
{
  static_assert(GridType::dimension == 2, "This test case is only available in 2d!");
public:
  typedef typename GridType::template Codim< 0 >::Entity EntityType;
  typedef typename GridType::ctype                       DomainFieldType;
  static const unsigned int                              dimDomain = GridType::dimension;
  typedef double            RangeFieldType;
  static const unsigned int dimRange = 1;
public:
  typedef Problems::OS2014::ParametricESV2007
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > ProblemType;
  typedef Stuff::Functions::Constant
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > ExactSolutionType;
  typedef typename ProblemType::FunctionType::NonparametricType            FunctionType;

  typedef std::map< std::string, Pymor::ParameterType > ParameterTypesMapType;
  typedef std::map< std::string, Pymor::Parameter >     ParametersMapType;

protected:
  static const size_t default_num_refinements_ = 3;

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
  static ParameterTypesMapType required_parameters()
  {
    return ParameterTypesMapType({{"mu",            Pymor::ParameterType("mu", 1)},
                                  {"mu_bar",        Pymor::ParameterType("mu", 1)},
                                  {"mu_hat",        Pymor::ParameterType("mu", 1)},
                                  {"mu_minimizing", Pymor::ParameterType("mu", 1)}});
  }

  ParametricConvergenceBase(const ParametersMapType parameters)
    : parameters_(parameters)
    , boundary_info_cfg_(Stuff::Grid::BoundaryInfoConfigs::AllDirichlet::default_config())
    , problem_(3)
    , exact_solution_(0)
  {}

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
        << "+==================================================================+" << std::endl;
  } // ... print_header(...)

  const Stuff::Common::Configuration& boundary_info() const
  {
    return boundary_info_cfg_;
  }

  const ProblemType& problem() const
  {
    return problem_;
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

  const ParametersMapType& parameters() const
  {
    return parameters_;
  }

protected:
  const ParametersMapType parameters_;
  const Stuff::Common::Configuration boundary_info_cfg_;
  const ProblemType problem_;
  const ExactSolutionType exact_solution_;
}; // class ParametricConvergenceBase


} // namespace internal


template< class GridType >
class ParametricConvergence
  : public internal::ParametricConvergenceBase< GridType >
  , public Base< GridType >
{
  typedef internal::ParametricConvergenceBase< GridType > ParametricConvergenceBaseType;
  typedef Base< GridType >                 TestCaseBaseType;

  static std::shared_ptr< GridType > create_initial_grid(const int refinements)
  {
    Stuff::Grid::Providers::Cube< GridType > grid_provider(-1, 1, 4);
    auto grid = grid_provider.grid_ptr();
    grid->globalRefine(refinements);
    return grid;
  } // ... create_initial_grid(...)

public:
  typedef typename TestCaseBaseType::ParametersMapType ParametersMapType;

  using ParametricConvergenceBaseType::required_parameters;
  using ParametricConvergenceBaseType::parameters;

  ParametricConvergence(const ParametersMapType parameters,
         const size_t num_refinements = ParametricConvergenceBaseType::default_num_refinements_)
    : ParametricConvergenceBaseType(parameters)
    , TestCaseBaseType(create_initial_grid(ParametricConvergenceBaseType::initial_refinements()), num_refinements)
  {
    this->check_parameters(ParametricConvergenceBaseType::required_parameters(), parameters);
    this->inherit_parameter_type(this->problem_, "problem");
  }
}; // class ParametricConvergence


#if HAVE_DUNE_GRID_MULTISCALE


template< class GridType >
class ParametricBlockConvergence
  : public internal::ParametricConvergenceBase< GridType >
  , public MultiscaleCubeBase< GridType >
{
  typedef internal::ParametricConvergenceBase< GridType > ParametricConvergenceBaseType;
  typedef MultiscaleCubeBase< GridType >   TestCaseBaseType;

  static Stuff::Common::Configuration initial_grid_cfg(const std::string num_partitions)
  {
    Stuff::Common::Configuration grid_cfg = Stuff::Grid::Providers::Cube< GridType >::default_config();
    grid_cfg["lower_left"] = "-1";
    grid_cfg["upper_right"] = "1";
    grid_cfg["num_elements"] = "4";
    grid_cfg["num_partitions"] = num_partitions;
    return grid_cfg;
  } // ... initial_grid_cfg(...)

public:
  typedef typename TestCaseBaseType::ParametersMapType ParametersMapType;

  using ParametricConvergenceBaseType::required_parameters;
  using ParametricConvergenceBaseType::parameters;

  ParametricBlockConvergence(const ParametersMapType parameters,
                   const std::string num_partitions = "[1 1 1]",
                   const size_t num_refinements = ParametricConvergenceBaseType::default_num_refinements_)
    : ParametricConvergenceBaseType(parameters)
    , TestCaseBaseType(initial_grid_cfg(num_partitions), ParametricConvergenceBaseType::initial_refinements(), num_refinements)
  {
    this->check_parameters(ParametricConvergenceBaseType::required_parameters(), parameters);
    this->inherit_parameter_type(this->problem_, "problem");
  }
}; // class ParametricBlockConvergence


#endif // HAVE_DUNE_GRID_MULTISCALE

} // namespace OS2014
} // namespace TestCases
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_TESTCASES_OS2014_HH
