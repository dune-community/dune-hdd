// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_TESTCASES_BATTERY_HH
#define DUNE_HDD_LINEARELLIPTIC_TESTCASES_BATTERY_HH

#include <map>
#include <string>
#include <memory>

#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/functions/constant.hh>

#include <dune/hdd/linearelliptic/problems/ORS2016.hh>

#include "base.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace TestCases {
namespace internal {


template< class GridType >
class ParametricBatteryBase
{
  static_assert(GridType::dimension == 3, "This test case is only available in 3d!");
public:
  typedef typename GridType::template Codim< 0 >::Entity EntityType;
  typedef typename GridType::ctype                       DomainFieldType;
  static const unsigned int                              dimDomain = GridType::dimension;
  typedef double            RangeFieldType;
  static const unsigned int dimRange = 1;
  typedef Problems::ORS2016< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > ProblemType;
  typedef Stuff::Functions::Constant
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >                            ExactSolutionType;
  typedef typename ProblemType::FunctionType::NonparametricType                                       FunctionType;

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

  static Stuff::Common::Configuration configuration(const std::string filename)
  {
    Stuff::Common::Configuration config = ProblemType::default_config();
    config["diffusion_factor.lower_left"]   = "[0      0     0]";
    config["diffusion_factor.upper_right"]  = "[0.0184 0.008 0.008]";
    config["diffusion_factor.num_elements"] = "[46     20    20]";
    config["diffusion_factor.separator"]    = "[0.0084 0.01; 0 0.008; 0 0.008]";
    config["diffusion_factor.filename"]     = filename;
    return config;
  } // ... configuration()

public:
  static ParameterTypesMapType required_parameters()
  {
    return ParameterTypesMapType({{"mu",     Pymor::ParameterType("mu", 1)},
                                  {"mu_bar", Pymor::ParameterType("mu", 1)},
                                  {"mu_hat", Pymor::ParameterType("mu", 1)}});
  }

  ParametricBatteryBase(const ParametersMapType parameters, const std::string& filename)
    : parameters_(add_parameter_range(parameters))
    , boundary_info_cfg_(Stuff::Common::Configuration({"type",                                "default", "dirichlet.0", "dirichlet.1"},
                                                      {"stuff.grid.boundaryinfo.normalbased", "neumann", "[-1 0 0]",    "[1 0 0]"}))
    , problem_(ProblemType::create(configuration(filename)))
    , exact_solution_(0)
  {}

  void print_header(std::ostream& out = std::cout) const
  {
    out << "+========================================+\n"
        << "|+======================================+|\n"
        << "||  Testcase: Battery                   ||\n"
        << "||  parameter: ELECTROLYTE in [0.1, 1]  ||\n"
        << "|+======================================+|\n"
        << "+========================================+" << std::endl;
  } // ... print_header(...)

  const Stuff::Common::Configuration& boundary_info() const
  {
    return boundary_info_cfg_;
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
    return exact_solution_;
  }

  const ParametersMapType& parameters() const
  {
    return parameters_;
  }

protected:
  const ParametersMapType parameters_;
  const Stuff::Common::Configuration boundary_info_cfg_;
  const std::unique_ptr< const ProblemType > problem_;
  const ExactSolutionType exact_solution_;
}; // class ParametricBatteryBase


} // namespace internal


#if HAVE_DUNE_GRID_MULTISCALE


template< class GridType >
class ParametricBlockBattery
  : public internal::ParametricBatteryBase< GridType >
  , public MultiscaleCubeBase< GridType >
{
  typedef internal::ParametricBatteryBase< GridType > BatteryBaseType;
  typedef MultiscaleCubeBase< GridType >              TestCaseBaseType;

  static Stuff::Common::Configuration initial_grid_cfg(const std::string num_partitions,
                                                       const size_t oversampling_layers)
  {
    Stuff::Common::Configuration grid_cfg = BatteryBaseType::configuration("");
    grid_cfg["type"] = TestCaseBaseType::MsGridProviderType::static_id();
    grid_cfg["num_partitions"] = num_partitions;
    grid_cfg.set("oversampling_layers", oversampling_layers, /*overwrite=*/true);
    return grid_cfg;
  } // ... initial_grid_cfg(...)

public:
  typedef typename TestCaseBaseType::ParametersMapType ParametersMapType;

  using BatteryBaseType::required_parameters;
  using BatteryBaseType::parameters;

  ParametricBlockBattery(const ParametersMapType parameters,
                         const std::string num_partitions = "[1 1 1]",
                         const size_t num_refinements = 0,
                         const size_t oversampling_layers = 0,
                         const bool H_with_h = false,
                         const std::string filename = "geometry__46x20x20_h4e-6m")
    : BatteryBaseType(parameters, filename)
    , TestCaseBaseType(initial_grid_cfg(num_partitions, oversampling_layers),
                       0,
                       num_refinements,
                       H_with_h)
  {
    this->check_parameters(BatteryBaseType::required_parameters(), parameters);
    this->inherit_parameter_type(*this->problem_, "problem");
  }
}; // class ParametricBlockBattery


#endif // HAVE_DUNE_GRID_MULTISCALE


} // namespace TestCases
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_TESTCASES_BATTERY_HH
