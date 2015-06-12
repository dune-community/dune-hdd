// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_TESTCASES_RANDOMBLOCK_HH
#define DUNE_HDD_LINEARELLIPTIC_TESTCASES_RANDOMBLOCK_HH

#include <algorithm>
#include <cmath>
#include <map>

#include <dune/stuff/common/disable_warnings.hh>
# if HAVE_ALUGRID
#   include <dune/grid/alugrid.hh>
# endif
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/stuff/common/fvector.hh>
#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/grid/provider.hh>

#include <dune/pymor/functions/default.hh>

#include <dune/hdd/linearelliptic/problems/default.hh>
#include <dune/hdd/linearelliptic/problems/random_block_problem.hh>

#include "base.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace TestCases {
namespace internal {


template< class GridType >
class RandomBlockTestcaseBase
{
  typedef RandomBlockTestcaseBase<GridType> ThisType;
public:
  typedef typename GridType::template Codim< 0 >::Entity EntityType;
  typedef typename GridType::ctype                       DomainFieldType;
  static const unsigned int                              dimDomain = GridType::dimension;
  typedef double            RangeFieldType;
  static const unsigned int dimRange = 1;
  typedef Dune::HDD::LinearElliptic::Problems::RandomBlockProblem
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > RandomBlockProblemType;

public:
  typedef Dune::HDD::LinearElliptic::Problems::Default
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > ProblemType;
  typedef Stuff::Functions::Constant
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > ExactSolutionType;
  typedef typename ProblemType::FunctionType::NonparametricType            FunctionType;

  typedef std::map< std::string, Pymor::ParameterType > ParameterTypesMapType;
  typedef std::map< std::string, Pymor::Parameter >     ParametersMapType;

protected:

  static ParametersMapType add_parameter_range(const ParametersMapType& parameters)
  {
    ParametersMapType ret = parameters;
    ret["parameter_range_min"] = Pymor::Parameter("mu", 0.1);
    ret["parameter_range_max"] = Pymor::Parameter("mu", 1.0);
    return ret;
  } // ... add_parameter_range(...)

public:
  static ParametersMapType default_parameters()
  {
    return ParametersMapType({{"mu",     Pymor::Parameter("mu", 1.0)},
                              {"mu_bar", Pymor::Parameter("mu", 1.0)},
                              {"mu_hat", Pymor::Parameter("mu", 1.0)}});
  } // ... default_parameters(...)

  static ParameterTypesMapType required_parameters()
  {
    const size_t parameter_size = 1;
    return ParameterTypesMapType({{"mu",     Pymor::ParameterType("mu", parameter_size)},
                                  {"mu_bar", Pymor::ParameterType("mu", parameter_size)},
                                  {"mu_hat", Pymor::ParameterType("mu", parameter_size)}});
  } // ... required_parameters(...)

  RandomBlockTestcaseBase(const ParametersMapType& parameters,
                   DSC::Configuration config )
    : parameters_(add_parameter_range(parameters))
    , boundary_info_cfg_(Stuff::Grid::BoundaryInfoConfigs::AllDirichlet::default_config())
    , problem_(*RandomBlockProblemType::create(config))
    , exact_solution_(0)
  {}

  void print_header(std::ostream& out = std::cout) const
  {
    out << "+==================================================================+\n"
        << "|+================================================================+|\n"
        << "||  Testcase RandomBlock                                         ||\n"
        << "|+----------------------------------------------------------------+|\n"
        << "||  domain = [0, 0] x [-1, 1]                                     ||\n"
        << "||  diffusion = mu_i in each block                                ||\n"
        << "||  force     = 1                                                 ||\n"
        << "||  dirichlet = 0                                                 ||\n"
        << "||  reference solution: discrete solution on finest grid          ||\n"
        << "||  parameter: mu in [0.1, 1]^p                                   ||\n"
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
}; // class RandomBlockBase


} // namespace internal

template< class GridType >
class RandomBlockTestcase
  : public internal::RandomBlockTestcaseBase< GridType >
  , public DSG::Providers::Cube< GridType >
  , public  internal::ParametricBase
{
  typedef internal::RandomBlockTestcaseBase< GridType > BaseType;
  typedef DSG::Providers::Cube< GridType > GridProviderType;

public:
  typedef typename Base<GridType>::ParametersMapType ParametersMapType;
  using BaseType::required_parameters;
  using BaseType::parameters;
  using BaseType::dimDomain;

  RandomBlockTestcase(const Stuff::Common::Configuration config,
               const ParametersMapType parameters
               = BaseType::default_parameters())
    : BaseType(parameters, config)
    , GridProviderType(*GridProviderType::create(config.sub("grids")))
  {
    this->check_parameters(BaseType::required_parameters(), parameters);
    this->inherit_parameter_type(this->problem_, "problem");
  }
}; // class BlockRandomBlockTestcase

} // namespace TestCases
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_TESTCASES_RANDOMBLOCK_HH
