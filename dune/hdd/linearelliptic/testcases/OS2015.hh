// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_TESTCASES_OS2015_HH
#define DUNE_HDD_LINEARELLIPTIC_TESTCASES_OS2015_HH

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

#include <dune/hdd/linearelliptic/problems/OS2015.hh>

#include "base.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace TestCases {
namespace OS2015 {
namespace internal {


static inline Stuff::Common::Configuration multiscale_problem_cfg()
{
  std::istringstream ss(
      "# forces\n"
      "forces.0.domain = [0.95 1.10; 0.30 0.45]\n"
      "forces.0.value = 2000\n"
      "forces.1.domain = [3.00 3.15; 0.75 0.90]\n"
      "forces.1.value = -1000\n"
      "forces.2.domain = [4.25 4.40; 0.25 0.40]\n"
      "forces.2.value = -1000\n"
      "# top channel\n"
      "channel.0.value = -1.07763239495\n"
      "channel.1.value = -1.07699512772\n"
      "channel.2.value = -1.07356156439\n"
      "channel.3.value = -1.06602281736\n"
      "channel.4.value = -1.06503683743\n"
      "channel.5.value = -1.07974870426\n"
      "channel.6.value = -1.05665895923\n"
      "channel.7.value = -1.08310334837\n"
      "channel.8.value = -1.05865484973\n"
      "channel.9.value = -1.05871039535\n"
      "channel.10.value = -1.08136695901\n"
      "channel.11.value = -1.08490172721\n"
      "channel.12.value = -1.06641120758\n"
      "channel.13.value = -1.06812773298\n"
      "channel.14.value = -1.07695652049\n"
      "channel.15.value = -1.08630079205\n"
      "channel.16.value = -1.08273722112\n"
      "channel.17.value = -1.07500402155\n"
      "channel.18.value = -1.08607142562\n"
      "channel.19.value = -1.07268761799\n"
      "channel.20.value = -1.08537037362\n"
      "channel.21.value = -1.08466927273\n"
      "channel.22.value = -1.08444661815\n"
      "channel.23.value = -1.08957037967\n"
      "channel.24.value = -1.08047394052\n"
      "channel.25.value = -1.08221229083\n"
      "channel.26.value = -1.08568599863\n"
      "channel.27.value = -1.08428347872\n"
      "channel.28.value = -1.09104098734\n"
      "channel.29.value = -1.09492700673\n"
      "channel.30.value = -1.09760440537\n"
      "channel.31.value = -1.09644989453\n"
      "channel.32.value = -1.09441681025\n"
      "channel.33.value = -1.09533290654\n"
      "channel.34.value = -1.1001430808\n"
      "channel.35.value = -1.10065627621\n"
      "channel.36.value = -1.10125877186\n"
      "channel.37.value = -1.10057485893\n"
      "channel.38.value = -1.10002261906\n"
      "channel.39.value = -1.10219154209\n"
      "channel.40.value = -1.09994463801\n"
      "channel.41.value = -1.10265630533\n"
      "channel.42.value = -1.10448566526\n"
      "channel.43.value = -1.10735820121\n"
      "channel.44.value = -1.1070022367\n"
      "channel.45.value = -1.10777650387\n"
      "channel.46.value = -1.10892785562\n"
      "channel.0.domain = [1.7 1.75; 0.5 0.55]\n"
      "channel.1.domain = [1.75 1.8; 0.5 0.55]\n"
      "channel.2.domain = [1.8 1.85; 0.5 0.55]\n"
      "channel.3.domain = [1.85 1.9; 0.5 0.55]\n"
      "channel.4.domain = [1.9 1.95; 0.5 0.55]\n"
      "channel.5.domain = [1.95 2.0; 0.5 0.55]\n"
      "channel.6.domain = [2.0 2.05; 0.5 0.55]\n"
      "channel.7.domain = [2.05 2.1; 0.5 0.55]\n"
      "channel.8.domain = [2.1 2.15; 0.5 0.55]\n"
      "channel.9.domain = [2.15 2.2; 0.5 0.55]\n"
      "channel.10.domain = [2.2 2.25; 0.5 0.55]\n"
      "channel.11.domain = [2.25 2.3; 0.5 0.55]\n"
      "channel.12.domain = [2.3 2.35; 0.5 0.55]\n"
      "channel.13.domain = [2.35 2.4; 0.5 0.55]\n"
      "channel.14.domain = [2.4 2.45; 0.5 0.55]\n"
      "channel.15.domain = [2.45 2.5; 0.5 0.55]\n"
      "channel.16.domain = [2.5 2.55; 0.5 0.55]\n"
      "channel.17.domain = [2.55 2.6; 0.5 0.55]\n"
      "channel.18.domain = [2.6 2.65; 0.5 0.55]\n"
      "channel.19.domain = [2.65 2.7; 0.5 0.55]\n"
      "channel.20.domain = [2.7 2.75; 0.5 0.55]\n"
      "channel.21.domain = [2.75 2.8; 0.5 0.55]\n"
      "channel.22.domain = [2.8 2.85; 0.5 0.55]\n"
      "channel.23.domain = [2.85 2.9; 0.5 0.55]\n"
      "channel.24.domain = [2.9 2.95; 0.5 0.55]\n"
      "channel.25.domain = [2.95 3.0; 0.5 0.55]\n"
      "channel.26.domain = [3.0 3.05; 0.5 0.55]\n"
      "channel.27.domain = [3.05 3.1; 0.5 0.55]\n"
      "channel.28.domain = [3.1 3.15; 0.5 0.55]\n"
      "channel.29.domain = [3.15 3.2; 0.5 0.55]\n"
      "channel.30.domain = [3.2 3.25; 0.5 0.55]\n"
      "channel.31.domain = [3.25 3.3; 0.5 0.55]\n"
      "channel.32.domain = [3.3 3.35; 0.5 0.55]\n"
      "channel.33.domain = [3.35 3.4; 0.5 0.55]\n"
      "channel.34.domain = [3.4 3.45; 0.5 0.55]\n"
      "channel.35.domain = [3.45 3.5; 0.5 0.55]\n"
      "channel.36.domain = [3.5 3.55; 0.5 0.55]\n"
      "channel.37.domain = [3.55 3.6; 0.5 0.55]\n"
      "channel.38.domain = [3.6 3.65; 0.5 0.55]\n"
      "channel.39.domain = [3.65 3.7; 0.5 0.55]\n"
      "channel.40.domain = [3.7 3.75; 0.5 0.55]\n"
      "channel.41.domain = [3.75 3.8; 0.5 0.55]\n"
      "channel.42.domain = [3.8 3.85; 0.5 0.55]\n"
      "channel.43.domain = [3.85 3.9; 0.5 0.55]\n"
      "channel.44.domain = [3.9 3.95; 0.5 0.55]\n"
      "channel.45.domain = [3.95 4.0; 0.5 0.55]\n"
      "channel.46.domain = [4.0 4.05; 0.5 0.55]\n"
      "# middle/top channel\n"
      "channel.47.value = -1.10372589211\n"
      "channel.48.value = -1.1020889988\n"
      "channel.49.value = -1.09806955069\n"
      "channel.50.value = -1.10000902421\n"
      "channel.51.value = -1.08797468724\n"
      "channel.52.value = -1.08827472176\n"
      "channel.53.value = -1.08692237109\n"
      "channel.54.value = -1.07893190093\n"
      "channel.55.value = -1.08748373853\n"
      "channel.56.value = -1.07445197324\n"
      "channel.57.value = -1.08246613163\n"
      "channel.58.value = -1.06726790504\n"
      "channel.59.value = -1.07891217847\n"
      "channel.60.value = -1.07260827126\n"
      "channel.61.value = -1.07094062748\n"
      "channel.62.value = -1.0692399429\n"
      "channel.63.value = -1.00099885701\n"
      "channel.64.value = -1.00109544002\n"
      "channel.65.value = -0.966491003242\n"
      "channel.66.value = -0.802284684014\n"
      "channel.67.value = -0.980790923021\n"
      "channel.68.value = -0.614478271687\n"
      "channel.69.value = -0.288129858959\n"
      "channel.70.value = -0.929509396842\n"
      "channel.71.value = -0.992376505995\n"
      "channel.72.value = -0.968162494855\n"
      "channel.73.value = -0.397316938901\n"
      "channel.74.value = -0.970934956609\n"
      "channel.75.value = -0.784344730096\n"
      "channel.76.value = -0.539725422323\n"
      "channel.77.value = -0.915632282372\n"
      "channel.78.value = -0.275089177273\n"
      "channel.79.value = -0.949684959286\n"
      "channel.80.value = -0.936132529794\n"
      "channel.47.domain = [2.6 2.65; 0.45 0.50]\n"
      "channel.48.domain = [2.65 2.7; 0.45 0.50]\n"
      "channel.49.domain = [2.7 2.75; 0.45 0.50]\n"
      "channel.50.domain = [2.75 2.8; 0.45 0.50]\n"
      "channel.51.domain = [2.8 2.85; 0.45 0.50]\n"
      "channel.52.domain = [2.85 2.9; 0.45 0.50]\n"
      "channel.53.domain = [2.9 2.95; 0.45 0.50]\n"
      "channel.54.domain = [2.95 3.0; 0.45 0.50]\n"
      "channel.55.domain = [3.0 3.05; 0.45 0.50]\n"
      "channel.56.domain = [3.05 3.1; 0.45 0.50]\n"
      "channel.57.domain = [3.1 3.15; 0.45 0.50]\n"
      "channel.58.domain = [3.15 3.2; 0.45 0.50]\n"
      "channel.59.domain = [3.2 3.25; 0.45 0.50]\n"
      "channel.60.domain = [3.25 3.3; 0.45 0.50]\n"
      "channel.61.domain = [3.3 3.35; 0.45 0.50]\n"
      "channel.62.domain = [3.35 3.4; 0.45 0.50]\n"
      "channel.63.domain = [3.4 3.45; 0.45 0.50]\n"
      "channel.64.domain = [3.45 3.5; 0.45 0.50]\n"
      "channel.65.domain = [3.5 3.55; 0.45 0.50]\n"
      "channel.66.domain = [3.55 3.6; 0.45 0.50]\n"
      "channel.67.domain = [3.6 3.65; 0.45 0.50]\n"
      "channel.68.domain = [3.65 3.7; 0.45 0.50]\n"
      "channel.69.domain = [3.7 3.75; 0.45 0.50]\n"
      "channel.70.domain = [3.75 3.8; 0.45 0.50]\n"
      "channel.71.domain = [3.8 3.85; 0.45 0.50]\n"
      "channel.72.domain = [3.85 3.9; 0.45 0.50]\n"
      "channel.73.domain = [3.9 3.95; 0.45 0.50]\n"
      "channel.74.domain = [3.95 4.0; 0.45 0.50]\n"
      "channel.75.domain = [4.0 4.05; 0.45 0.50]\n"
      "channel.76.domain = [4.05 4.1; 0.45 0.50]\n"
      "channel.77.domain = [4.1 4.15; 0.45 0.50]\n"
      "channel.78.domain = [4.15 4.2; 0.45 0.50]\n"
      "channel.79.domain = [4.2 4.25; 0.45 0.50]\n"
      "channel.80.domain = [4.25 4.3; 0.45 0.50]\n"
      "# middle/bottom channel\n"
      "channel.81.value = -1.10923642795\n"
      "channel.82.value = -1.10685618623\n"
      "channel.83.value = -1.1057800376\n"
      "channel.84.value = -1.10187723629\n"
      "channel.85.value = -1.10351710464\n"
      "channel.86.value = -1.10037551137\n"
      "channel.87.value = -1.09724407076\n"
      "channel.88.value = -1.09604600208\n"
      "channel.89.value = -1.09354469656\n"
      "channel.90.value = -1.08934455354\n"
      "channel.91.value = -1.08155476586\n"
      "channel.92.value = -1.07815397899\n"
      "channel.93.value = -1.09174062023\n"
      "channel.94.value = -1.07433616068\n"
      "channel.95.value = -1.08030587701\n"
      "channel.81.domain = [1.95 2.0; 0.40 0.45]\n"
      "channel.82.domain = [2.0 2.05; 0.40 0.45]\n"
      "channel.83.domain = [2.05 2.1; 0.40 0.45]\n"
      "channel.84.domain = [2.1 2.15; 0.40 0.45]\n"
      "channel.85.domain = [2.15 2.2; 0.40 0.45]\n"
      "channel.86.domain = [2.2 2.25; 0.40 0.45]\n"
      "channel.87.domain = [2.25 2.3; 0.40 0.45]\n"
      "channel.88.domain = [2.3 2.35; 0.40 0.45]\n"
      "channel.89.domain = [2.35 2.4; 0.40 0.45]\n"
      "channel.90.domain = [2.4 2.45; 0.40 0.45]\n"
      "channel.91.domain = [2.45 2.5; 0.40 0.45]\n"
      "channel.92.domain = [2.5 2.55; 0.40 0.45]\n"
      "channel.93.domain = [2.55 2.6; 0.40 0.45]\n"
      "channel.94.domain = [2.6 2.65; 0.40 0.45]\n"
      "channel.95.domain = [2.65 2.7; 0.40 0.45]\n"
      "# bottom channel\n"
      "channel.96.value = -1.00032869407\n"
      "channel.97.value = -1.01175908905\n"
      "channel.98.value = -1.04954395793\n"
      "channel.99.value = -1.017967697\n"
      "channel.100.value = -1.04647184091\n"
      "channel.101.value = -1.01911894831\n"
      "channel.102.value = -1.00699340158\n"
      "channel.103.value = -0.995492960025\n"
      "channel.104.value = -1.0373059007\n"
      "channel.96.domain = [2.25 2.3; 0.35 0.40]\n"
      "channel.97.domain = [2.3 2.35; 0.35 0.40]\n"
      "channel.98.domain = [2.35 2.4; 0.35 0.40]\n"
      "channel.99.domain = [2.4 2.45; 0.35 0.40]\n"
      "channel.100.domain = [2.45 2.5; 0.35 0.40]\n"
      "channel.101.domain = [2.5 2.55; 0.35 0.40]\n"
      "channel.102.domain = [2.55 2.6; 0.35 0.40]\n"
      "channel.103.domain = [2.6 2.65; 0.35 0.40]\n"
      "channel.104.domain = [2.65 2.7; 0.35 0.40]"
  );
  Stuff::Common::Configuration config(ss);
  config["parametric_channel"] = "true";
  config["lower_left"] = "[0 0]";
  config["upper_right"] = "[5 1]";
  config["channel_boundary_layer"] = "0";
  return config;
} // ... multiscale_problem_cfg(...)


static inline Stuff::Common::Configuration multiscale_grid_cfg()
{
  Stuff::Common::Configuration config;
  config["lower_left"]   = "[0 0]";
  config["upper_right"]  = "[5 1]";
  config["num_elements"] = "[100 20]";
  return config;
} // ... multiscale_grid_cfg(...)


template< class GridType >
class AcademicBase
{
  static_assert(GridType::dimension == 2, "This test case is only available in 2d!");
public:
  typedef typename GridType::template Codim< 0 >::Entity EntityType;
  typedef typename GridType::ctype                       DomainFieldType;
  static const unsigned int                              dimDomain = GridType::dimension;
  typedef double            RangeFieldType;
  static const unsigned int dimRange = 1;
public:
  typedef Problems::OS2015::ParametricESV2007
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

  static ParametersMapType add_parameter_range(const ParametersMapType& parameters)
  {
    ParametersMapType ret = parameters;
    ret["parameter_range_min"] = Pymor::Parameter("mu", 0.1);
    ret["parameter_range_max"] = Pymor::Parameter("mu", 1.0);
    return ret;
  } // ... add_parameter_range(...)

public:
  static ParameterTypesMapType required_parameters()
  {
    return ParameterTypesMapType({{"mu",            Pymor::ParameterType("mu", 1)},
                                  {"mu_bar",        Pymor::ParameterType("mu", 1)},
                                  {"mu_hat",        Pymor::ParameterType("mu", 1)}});
  }

  AcademicBase(const ParametersMapType& parameters)
    : parameters_(add_parameter_range(parameters))
    , boundary_info_cfg_(Stuff::Grid::BoundaryInfoConfigs::AllDirichlet::default_config())
    , problem_(3)
    , exact_solution_(0)
  {}

  void print_header(std::ostream& out = std::cout) const
  {
    out << "+==========================================================+\n"
        << "|+========================================================+|\n"
        << "||  Testcase OS2015: academic                             ||\n"
        << "||  (see Ohlberger, Schindler, 2016, sec. 6)              ||\n"
        << "|+--------------------------------------------------------+|\n"
        << "||  domain = [-1, 1] x [-1, 1]                            ||\n"
        << "||  diffusion = 1 + (1 - mu) cos(1/2 pi x) cos(1/2 pi y)  ||\n"
        << "||  force     = 1/2 pi^2 cos(1/2 pi x) cos(1/2 pi y)      ||\n"
        << "||  dirichlet = 0                                         ||\n"
        << "||  reference solution: discrete solution on finest grid  ||\n"
        << "||  parameter: mu in [0.1, 1]                             ||\n"
        << "|+========================================================+|\n"
        << "+==========================================================+" << std::endl;
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
}; // class AcademicBase


template< class GridType >
class MultiscaleBase
{
  static_assert(GridType::dimension == 2, "This test case is only available in 2d!");
public:
  typedef typename GridType::template Codim< 0 >::Entity EntityType;
  typedef typename GridType::ctype                       DomainFieldType;
  static const unsigned int                              dimDomain = GridType::dimension;
  typedef double            RangeFieldType;
  static const unsigned int dimRange = 1;
  typedef Problems::OS2015::Spe10Model1< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > ProblemType;
  typedef Stuff::Functions::Constant
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > ExactSolutionType;
  typedef typename ProblemType::FunctionType::NonparametricType            FunctionType;

  typedef std::map< std::string, Pymor::ParameterType > ParameterTypesMapType;
  typedef std::map< std::string, Pymor::Parameter >     ParametersMapType;

protected:
  static const size_t default_num_refinements_ = 1;

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

  static ParametersMapType add_parameter_range(const ParametersMapType& parameters)
  {
    ParametersMapType ret = parameters;
    ret["parameter_range_min"] = Pymor::Parameter("mu", 0.1);
    ret["parameter_range_max"] = Pymor::Parameter("mu", 1.0);
    return ret;
  } // ... add_parameter_range(...)

  static Stuff::Common::Configuration configuration(const std::string filename)
  {
    auto config = multiscale_problem_cfg();
    config.add(multiscale_grid_cfg(), "", true);
    config["filename"] = filename;
    return config;
  } // ... configuration()

public:
  static ParameterTypesMapType required_parameters()
  {
    return ParameterTypesMapType({{"mu",     Pymor::ParameterType("mu", 1)},
                                  {"mu_bar", Pymor::ParameterType("mu", 1)},
                                  {"mu_hat", Pymor::ParameterType("mu", 1)}});
  }

  MultiscaleBase(const ParametersMapType parameters, const std::string filename)
    : parameters_(add_parameter_range(parameters))
    , boundary_info_cfg_(Stuff::Grid::BoundaryInfoConfigs::AllDirichlet::default_config())
    , problem_(ProblemType::create(configuration(filename)))
    , exact_solution_(0)
  {}

  void print_header(std::ostream& out = std::cout) const
  {
    out << "+==========================================================+\n"
        << "|+========================================================+|\n"
        << "||  Testcase OS2015: multiscale                           ||\n"
        << "||  (see Ohlberger, Schindler, 2016, sec. 6)              ||\n"
        << "|+--------------------------------------------------------+|\n"
        << "||  domain = [0, 5] x [0, 1]                              ||\n"
        << "||  diffusion factor: 1 + (1 - mu)*channel                ||\n"
        << "||  diffusion tensor: spe10 model 1 scalar data           ||\n"
        << "||    (http://www.spe.org/web/csp/datasets/set01.htm)     ||\n"
        << "||         |  2000 in [0.55, 0.70] x [0.70, 0.85]         ||\n"
        << "||  force: | -1000 in [3.00, 3.15] x [0.77, 0.90]         ||\n"
        << "||         | -1000 in [4.30, 4.45] x [0.50, 0.65]         ||\n"
        << "||  dirichlet = 0                                         ||\n"
        << "||  reference solution: discrete solution on finest grid  ||\n"
        << "||  parameter: mu in [0.1, 1]                             ||\n"
        << "|+========================================================+|\n"
        << "+==========================================================+" << std::endl;
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
}; // class MultiscaleBase


} // namespace internal


#if HAVE_DUNE_GRID_MULTISCALE


template< class GridType >
class Academic
  : public internal::AcademicBase< GridType >
  , public MultiscaleCubeBase< GridType >
{
  typedef internal::AcademicBase< GridType > AcademicBaseType;
  typedef MultiscaleCubeBase< GridType >   TestCaseBaseType;

  static Stuff::Common::Configuration initial_grid_cfg(const std::string num_partitions,
                                                       const size_t oversampling_layers)
  {
    Stuff::Common::Configuration grid_cfg = Stuff::Grid::Providers::Cube< GridType >::default_config();
    grid_cfg["lower_left"] = "-1";
    grid_cfg["upper_right"] = "1";
    grid_cfg["num_elements"] = "4";
    grid_cfg["num_partitions"] = num_partitions;
    grid_cfg.set("oversampling_layers", oversampling_layers, /*overwrite=*/true);
    return grid_cfg;
  } // ... initial_grid_cfg(...)

public:
  typedef typename TestCaseBaseType::ParametersMapType ParametersMapType;

  using AcademicBaseType::required_parameters;
  using AcademicBaseType::parameters;

  Academic(const ParametersMapType params,
           const std::string num_partitions = "[1 1 1]",
           const size_t num_refinements = AcademicBaseType::default_num_refinements_,
           const size_t oversampling_layers = 0,
           const bool H_with_h = false)
    : AcademicBaseType(params)
    , TestCaseBaseType(initial_grid_cfg(num_partitions, oversampling_layers),
                       AcademicBaseType::initial_refinements(),
                       num_refinements,
                       H_with_h)
  {
    this->check_parameters(AcademicBaseType::required_parameters(), params);
    this->inherit_parameter_type(this->problem_, "problem");
  }
}; // class Academic


template< class GridType >
class Multiscale
  : public internal::MultiscaleBase< GridType >
  , public MultiscaleCubeBase< GridType >
{
  typedef internal::MultiscaleBase< GridType > MultiscaleBaseType;
  typedef MultiscaleCubeBase< GridType >       TestCaseBaseType;

  static Stuff::Common::Configuration initial_grid_cfg(const std::string num_partitions,
                                                       const size_t oversampling_layers)
  {
    Stuff::Common::Configuration grid_cfg = MultiscaleBaseType::configuration("");
    grid_cfg["num_partitions"] = num_partitions;
    grid_cfg.set("oversampling_layers", oversampling_layers, /*overwrite=*/true);
    return grid_cfg;
  } // ... initial_grid_cfg(...)

public:
  typedef typename TestCaseBaseType::ParametersMapType ParametersMapType;

  using MultiscaleBaseType::required_parameters;
  using MultiscaleBaseType::parameters;

  Multiscale(const ParametersMapType params,
             const std::string num_partitions = "[1 1 1]",
             const size_t num_refinements = MultiscaleBaseType::default_num_refinements_,
             const size_t oversampling_layers = 0,
             const bool H_with_h = false,
             const std::string filename = Stuff::Functions::Spe10::internal::model1_filename)
    : MultiscaleBaseType(params, filename)
    , TestCaseBaseType(initial_grid_cfg(num_partitions, oversampling_layers),
                       MultiscaleBaseType::initial_refinements(),
                       num_refinements,
                       H_with_h)
  {
    this->check_parameters(MultiscaleBaseType::required_parameters(), params);
    this->inherit_parameter_type(*this->problem_, "problem");
  }
}; // class Multiscale


#endif // HAVE_DUNE_GRID_MULTISCALE

} // namespace OS2015
} // namespace TestCases
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_TESTCASES_OS2015_HH
