// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_TESTCASES_SPE10MODEL2_HH
#define DUNE_HDD_LINEARELLIPTIC_TESTCASES_SPE10MODEL2_HH

#include <map>
#include <string>
#include <memory>

#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/functions/constant.hh>

#include <dune/hdd/linearelliptic/problems/spe10model2.hh>

#include "base.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace TestCases {
namespace Spe10 {


template< class GridType >
class Model2 :public Base< GridType >
{
  static_assert(GridType::dimension == 3, "This test case is only available in 3D!");
  typedef Base< GridType > BaseType;
public:
  typedef typename GridType::template Codim< 0 >::Entity EntityType;
  typedef typename GridType::ctype                       DomainFieldType;
  static const unsigned int                              dimDomain = GridType::dimension;
  typedef double            RangeFieldType;
  static const unsigned int dimRange = 1;
  typedef Problems::Spe10Model2< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > ProblemType;
  typedef Stuff::Functions::Constant
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >                            ExactSolutionType;
  typedef typename ProblemType::FunctionType::NonparametricType                                       FunctionType;

protected:
  static const size_t default_num_refinements_ = 1;

  static int initial_refinements()
  {
    int ret = 0;
#if HAVE_ALUGRID
    if (std::is_same< GridType, ALUGrid< 3, 3, simplex, conforming > >::value)
      ret += 1;
#endif // HAVE_ALUGRID
    return ret;
  } // ... initial_refinements()

  static Stuff::Common::Configuration configuration(const std::string filename)
  {
    Stuff::Common::Configuration config = ProblemType::default_config();
    config.set("num_elements", ProblemType::Spe10Model2::Spe10FunctionType::num_elements);
    config["filename"] = filename;
    return config;
  } // ... configuration()

  static std::shared_ptr< GridType > create_initial_grid(const int refinements)
  {
    auto grid = Stuff::Grid::Providers::Cube< GridType >::create(configuration(""))->grid_ptr();
    grid->preAdapt();
    grid->globalRefine(refinements);
    grid->postAdapt();
    grid->loadBalance();
    return grid;
  } // ... create_initial_grid(...)

public:
  Model2(std::string filename = "",  size_t num_refinements = 0)
    : BaseType(create_initial_grid(initial_refinements()), num_refinements)
    , boundary_info_cfg_()
    , problem_(ProblemType::create(configuration(filename)))
    , exact_solution_(0)
  {
    boundary_info_cfg_["type"] = Stuff::Grid::BoundaryInfoConfigs::NormalBased::static_id();
    boundary_info_cfg_["default"] = "neumann";
    boundary_info_cfg_["compare_tolerance"] = "1e-10";
    boundary_info_cfg_["dirichlet.0"] = "[0.0 -1.0 0.0]";
  }

  void print_header(std::ostream& out = std::cout) const
  {
    out << "+==========================================================+\n"
        << "|+========================================================+|\n"
        << "||  Testcase: SPE10, Model2                               ||\n"
        << "||  (see http://www.spe.org/web/csp/datasets/set01.htm)   ||\n"
        << "|+--------------------------------------------------------+|\n"
        << "||  domain = [0, 365.76] x [0, 670.56] x [0, 51.816]      ||\n"
        << "||  diffusion: spe10 model 1 scalar data                  ||\n"
        << "||         |  2000 in [0.55, 0.70] x [0.70, 0.85]         ||\n"
        << "||  force: | -1000 in [3.00, 3.15] x [0.77, 0.90]         ||\n"
        << "||         | -1000 in [4.30, 4.45] x [0.50, 0.65]         ||\n"
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

protected:
  Stuff::Common::Configuration boundary_info_cfg_;
  const std::unique_ptr< const ProblemType > problem_;
  const ExactSolutionType exact_solution_;
}; // class Spe10Model1Base

} // namespace Spe10
} // namespace TestCases
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_TESTCASES_SPE10MODEL2_HH
