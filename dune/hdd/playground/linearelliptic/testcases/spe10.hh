// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_TESTCASES_SPE10_HH
#define DUNE_HDD_LINEARELLIPTIC_TESTCASES_SPE10_HH

#include <algorithm>
#include <cmath>
#include <map>

#include <dune/stuff/common/disable_warnings.hh>
# if HAVE_ALUGRID
#   include <dune/grid/alugrid.hh>
# endif
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/spe10.hh>
#include <dune/stuff/playground/functions/indicator.hh>

#include <dune/hdd/linearelliptic/problems/default.hh>

#include "base.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace TestCases {
namespace internal {


template< class GridType >
class Spe10Model1Base
{
  static_assert(GridType::dimension == 2, "This test case is only available in 2d!");
public:
  typedef typename GridType::template Codim< 0 >::Entity EntityType;
  typedef typename GridType::ctype                       DomainFieldType;
  static const unsigned int                              dimDomain = GridType::dimension;
  typedef double            RangeFieldType;
  static const unsigned int dimRange = 1;
  typedef Problems::Default< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > ProblemType;
  typedef Stuff::Functions::Constant< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >
      ExactSolutionType;
  typedef typename ProblemType::FunctionType::NonparametricType FunctionType;
private:
  typedef Stuff::Functions::Indicator < EntityType, DomainFieldType, dimDomain, RangeFieldType > IndicatorFunctionType;
protected:
  typedef Stuff::Functions::Spe10Model1
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimDomain, dimDomain > Spe10Model1FunctionType;

  static const size_t default_num_refinements_ = 1;

  static int initial_refinements()
  {
    int ret = 0;
#if HAVE_ALUGRID
    if (std::is_same< GridType, ALUConformGrid< 2, 2 > >::value
        || std::is_same< GridType, ALUGrid< 2, 2, simplex, conforming > >::value)
      ret += 1;
#endif // HAVE_ALUGRID
    return ret;
  } // ... initial_refinements()

  static Stuff::Common::Configuration configuration(const std::string filename)
  {
    Stuff::Common::Configuration config = Spe10Model1FunctionType::default_config();
    config["lower_left"]   = "[0 0]";
    config["upper_right"]  = "[5 1]";
    config["num_elements"] = "[100 20]";
    config["name"] = "diffusion_tensor";
    config["filename"] = filename;
    return config;
  } // ... configuration()

public:
  Spe10Model1Base(const std::string filename)
    : boundary_info_cfg_(Stuff::Grid::BoundaryInfoConfigs::AllDirichlet::default_config())
    , problem_(std::make_shared< ExactSolutionType >(1, "diffusion_factor"),
               Spe10Model1FunctionType::create(configuration(filename)),
               std::shared_ptr< IndicatorFunctionType >(new IndicatorFunctionType(
                   {{{{0.55, 0.70}, {0.70, 0.85}},  1.0},
                    {{{0.45, 0.15}, {0.60, 0.30}},  1.0},
                    {{{3.00, 0.75}, {3.15, 0.90}}, -1.0},
                    {{{4.30, 0.50}, {4.45, 0.65}}, -1.0}},
                   "force")),
               std::make_shared< ExactSolutionType >(0, "dirichlet"),
               std::make_shared< ExactSolutionType >(0, "neumann"))
    , exact_solution_(0)
  {}

  void print_header(std::ostream& out = std::cout) const
  {
    out << "+==========================================================+\n"
        << "|+========================================================+|\n"
        << "||  Testcase: SPE10, Model1                               ||\n"
        << "||  (see http://www.spe.org/web/csp/datasets/set01.htm)   ||\n"
        << "|+--------------------------------------------------------+|\n"
        << "||  domain = [0, 5] x [0, 1]                              ||\n"
        << "||  diffusion: spe10 model 1 scalar data                  ||\n"
        << "||         /  1 in [0.55, 0.70] x [0.70, 0.85]            ||\n"
        << "||         |  1 in [0.45, 0.60] x [0.15, 0.30]            ||\n"
        << "||  force: { -1 in [3.00, 3.15] x [0.77, 0.90]            ||\n"
        << "||         | -1 in [4.30, 4.45] x [0.50, 0.65]            ||\n"
        << "||         \\  0 else                                      ||\n"
        << "||  dirichlet = 0                                         ||\n"
        << "||  reference solution: discrete solution on finest grid  ||\n"
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
    return exact_solution_;
  }

protected:
  const Stuff::Common::Configuration boundary_info_cfg_;
  const ProblemType problem_;
  const ExactSolutionType exact_solution_;
}; // class Spe10Model1Base


} // namespace internal


template< class GridType >
class Spe10Model1
  : public internal::Spe10Model1Base< GridType >
  , public Base< GridType >
{
  typedef internal::Spe10Model1Base< GridType > Spe10Model1BaseType;
  typedef Base< GridType >                      TestCaseBaseType;

  typedef typename Spe10Model1BaseType::Spe10Model1FunctionType Spe10Model1FunctionType;

  static std::shared_ptr< GridType > create_initial_grid(const int refinements)
  {
    auto grid = Stuff::Grid::Providers::Cube< GridType >::create(Spe10Model1BaseType::configuration(""))->grid_ptr();
    grid->globalRefine(refinements);
    return grid;
  } // ... create_initial_grid(...)

public:
  Spe10Model1(const std::string filename
                  = Spe10Model1FunctionType::default_config().template get< std::string >("filename"),
              const size_t num_refinements = Spe10Model1BaseType::default_num_refinements_)
    : Spe10Model1BaseType(filename)
    , TestCaseBaseType(create_initial_grid(Spe10Model1BaseType::initial_refinements()), num_refinements)
  {}
}; // class Spe10Model1


} // namespace TestCases
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_TESTCASES_SPE10_HH
