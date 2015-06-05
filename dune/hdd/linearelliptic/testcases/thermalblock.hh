// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_TESTCASES_THERMALBLOCK_HH
#define DUNE_HDD_LINEARELLIPTIC_TESTCASES_THERMALBLOCK_HH

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
#include <dune/hdd/linearelliptic/problems/thermalblock.hh>

#include "base.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace TestCases {
namespace internal {


template< class GridType >
class ThermalblockBase
{
  typedef ThermalblockBase<GridType> ThisType;
public:
  typedef typename GridType::template Codim< 0 >::Entity EntityType;
  typedef typename GridType::ctype                       DomainFieldType;
  static const unsigned int                              dimDomain = GridType::dimension;
  typedef double            RangeFieldType;
  static const unsigned int dimRange = 1;
  typedef Problems::Thermalblock
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > ThermalblockProblemType;

public:
  typedef Problems::Default
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > ProblemType;
  typedef Stuff::Functions::Constant
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > ExactSolutionType;
  typedef typename ProblemType::FunctionType::NonparametricType            FunctionType;

  typedef std::map< std::string, Pymor::ParameterType > ParameterTypesMapType;
  typedef std::map< std::string, Pymor::Parameter >     ParametersMapType;

  static const size_t default_num_refinements = 3;

protected:
  static int initial_refinements()
  {
    int ret = 1;
#if HAVE_ALUGRID
    if (std::is_same< GridType, ALUGrid< 2, 2, simplex, conforming > >::value)
      ret += 1;
#endif // HAVE_ALUGRID
    return ret;
  } // ... initial_refinements()

  static ParametersMapType add_parameter_range(const DSC::FieldVector< size_t, dimDomain >& num_blocks,
                                               const ParametersMapType& parameters)
  {
    size_t parameter_size = 1;
    for (const auto& element : num_blocks)
      parameter_size *= element;
    ParametersMapType ret = parameters;
    ret["parameter_range_min"] = Pymor::Parameter("mu", std::vector< double >(parameter_size, 0.1));
    ret["parameter_range_max"] = Pymor::Parameter("mu", std::vector< double >(parameter_size, 1.0));
    return ret;
  } // ... add_parameter_range(...)

  static DSC::Configuration problem_config(const DSC::FieldVector< size_t, dimDomain >& num_blocks)
  {
    DSC::Configuration ret = ThermalblockProblemType::default_config();
    assert(ret.has_key("diffusion_factor.num_elements"));
    ret.set("diffusion_factor.num_elements", num_blocks, /*overwrite=*/true);
    return ret;
  } // ... problem_config(...)

public:
  static ParametersMapType default_parameters(const DSC::FieldVector< size_t, dimDomain >& num_blocks)
  {
    size_t parameter_size = 1;
    for (const auto& element : num_blocks)
      parameter_size *= element;
    return ParametersMapType({{"mu",     Pymor::Parameter("mu", std::vector< double >(parameter_size, 1.0))},
                              {"mu_bar", Pymor::Parameter("mu", std::vector< double >(parameter_size, 1.0))},
                              {"mu_hat", Pymor::Parameter("mu", std::vector< double >(parameter_size, 1.0))}});
  } // ... default_parameters(...)

  static ParameterTypesMapType required_parameters(const DSC::FieldVector< size_t, dimDomain >& num_blocks)
  {
    size_t parameter_size = 1;
    for (const auto& element : num_blocks)
      parameter_size *= element;
    return ParameterTypesMapType({{"mu",     Pymor::ParameterType("mu", parameter_size)},
                                  {"mu_bar", Pymor::ParameterType("mu", parameter_size)},
                                  {"mu_hat", Pymor::ParameterType("mu", parameter_size)}});
  } // ... required_parameters(...)

  ThermalblockBase(const DSC::FieldVector< size_t, dimDomain >& num_blocks,
                   const ParametersMapType& parameters
                   = ThisType::default_parameters(DSC::FieldVector<size_t, dimDomain>{{2u,2u}}))
    : parameters_(add_parameter_range(num_blocks, parameters))
    , boundary_info_cfg_(Stuff::Grid::BoundaryInfoConfigs::AllDirichlet::default_config())
    , problem_(*ThermalblockProblemType::create(problem_config(num_blocks)))
    , exact_solution_(0)
  {}

  void print_header(std::ostream& out = std::cout) const
  {
    out << "+==================================================================+\n"
        << "|+================================================================+|\n"
        << "||  Testcase thermalblock                                         ||\n"
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
}; // class ThermalblockBase


} // namespace internal


#if HAVE_DUNE_GRID_MULTISCALE


template< class GridType >
class BlockThermalblock
  : public internal::ThermalblockBase< GridType >
  , public MultiscaleCubeBase< GridType >
{
  typedef internal::ThermalblockBase< GridType > ThermalblockBaseType;
  typedef MultiscaleCubeBase< GridType >         TestCaseBaseType;

  static Stuff::Common::Configuration initial_grid_cfg(const std::string num_partitions)
  {
    Stuff::Common::Configuration grid_cfg = Stuff::Grid::Providers::Cube< GridType >::default_config();
    grid_cfg["lower_left"] = "0";
    grid_cfg["upper_right"] = "1";
    grid_cfg["num_elements"] = "4";
    grid_cfg["num_partitions"] = num_partitions;
    return grid_cfg;
  } // ... initial_grid_cfg(...)

public:
  typedef typename TestCaseBaseType::ParametersMapType ParametersMapType;

  using ThermalblockBaseType::required_parameters;
  using ThermalblockBaseType::parameters;
  using ThermalblockBaseType::dimDomain;

  BlockThermalblock(const ParametersMapType parameters,
                    const DSC::FieldVector< size_t, dimDomain >& num_blocks = {2, 2},
                    const std::string num_partitions = "[1 1 1]",
                    const size_t num_refinements = ThermalblockBaseType::default_num_refinements)
    : ThermalblockBaseType(num_blocks, parameters)
    , TestCaseBaseType(initial_grid_cfg(num_partitions), ThermalblockBaseType::initial_refinements(), num_refinements)
  {
    this->check_parameters(ThermalblockBaseType::required_parameters(num_blocks), parameters);
    this->inherit_parameter_type(this->problem_, "problem");
  }
}; // class BlockThermalblock


#else // HAVE_DUNE_GRID_MULTISCALE


template< class GridType >
class BlockThermalblock
{
  static_assert(AlwaysFalse< GridType >::value, "You are missing dune-grid-multiscale!");
};


#endif // HAVE_DUNE_GRID_MULTISCALE

template< class GridType >
class Thermalblock
  : public internal::ThermalblockBase< GridType >
  , public DSG::Providers::Cube< GridType >
  , public  internal::ParametricBase
{
  typedef internal::ThermalblockBase< GridType > ThermalblockBaseType;
  typedef DSG::Providers::Cube< GridType > GridProviderType;

  static Stuff::Common::Configuration initial_grid_cfg(const std::string num_partitions)
  {
    Stuff::Common::Configuration grid_cfg = Stuff::Grid::Providers::Cube< GridType >::default_config();
    grid_cfg["lower_left"] = "0";
    grid_cfg["upper_right"] = "1";
    grid_cfg["num_elements"] = "4";
    grid_cfg["num_partitions"] = num_partitions;
    return grid_cfg;
  } // ... initial_grid_cfg(...)

public:
  typedef typename Base<GridType>::ParametersMapType ParametersMapType;

  using ThermalblockBaseType::required_parameters;
  using ThermalblockBaseType::parameters;
  using ThermalblockBaseType::dimDomain;

  Thermalblock(const ParametersMapType parameters,
                    const DSC::FieldVector< size_t, dimDomain >& num_blocks = {2, 2},
                    const std::string num_partitions = "[1 1 1]",
                    const size_t num_refinements = ThermalblockBaseType::default_num_refinements)
    : ThermalblockBaseType(num_blocks, parameters)
    , GridProviderType(*GridProviderType::create(initial_grid_cfg(num_partitions)))
  {
    this->check_parameters(ThermalblockBaseType::required_parameters(num_blocks), parameters);
    this->inherit_parameter_type(this->problem_, "problem");
  }
}; // class BlockThermalblock

} // namespace TestCases
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_TESTCASES_OS2014_HH
