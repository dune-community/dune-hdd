// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_TESTCASES_OS2014_HH
#define DUNE_HDD_LINEARELLIPTIC_TESTCASES_OS2014_HH

#include <dune/stuff/common/disable_warnings.hh>
# if HAVE_ALUGRID
#   include <dune/grid/alugrid.hh>
# endif
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/stuff/functions/ESV2007.hh>
#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/common/configuration.hh>

#include <dune/hdd/linearelliptic/problems/ESV2007.hh>
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

  typedef Problems::ESV2007< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > ProblemType;
  typedef Stuff::Functions::ESV2007::Testcase1ExactSolution
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > ExactSolutionType;

  typedef Stuff::LocalizableFunctionInterface < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >
      FunctionType;

protected:
  static const size_t default_num_refinements = 3;

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
  OS2014Base()
    : boundary_info_cfg_(Stuff::Grid::BoundaryInfoConfigs::AllDirichlet::default_config())
    , problem_(3)
    , exact_solution_(2, problem_.static_id() + ".exact_solution")
  {}

  void print_header(std::ostream& out = std::cout) const
  {
    out << "+==================================================================+\n"
        << "|+================================================================+|\n"
        << "||  Testcase OS2014: smooth data, homogeneous dirichlet          ||\n"
        << "||  (see testcase 1, page 23 in Ern, Stephansen, Vohralik, 2007)  ||\n"
        << "|+----------------------------------------------------------------+|\n"
        << "||  domain = [-1, 1] x [-1, 1]                                    ||\n"
        << "||  diffusion = 1                                                 ||\n"
        << "||  force     = 1/2 pi^2 cos(1/2 pi x) cos(1/2 pi y)              ||\n"
        << "||  dirichlet = 0                                                 ||\n"
        << "||  exact solution = cos(1/2 pi x) cos(1/2 pi y)                  ||\n"
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
    return true;
  }

  const ExactSolutionType& exact_solution() const
  {
    return exact_solution_;
  }

private:
  const Stuff::Common::Configuration boundary_info_cfg_;
  const ProblemType problem_;
  const ExactSolutionType exact_solution_;
}; // class OS2014


} // namespace internal


template< class GridType >
class OS2014
  : public internal::OS2014Base< GridType >
  , public Base< GridType >
{
  typedef internal::OS2014Base< GridType > OS2014BaseType;
  typedef Base< GridType >                  TestCaseBaseType;

private:
  static std::shared_ptr< GridType > create_initial_grid(const int refinements, const int start_cell_count)
  {
    Stuff::Grid::Providers::Cube< GridType > grid_provider(-1, 1, start_cell_count);
    auto grid = grid_provider.grid_ptr();
    grid->preAdapt();
    grid->globalRefine(refinements);
    grid->postAdapt();
    grid->loadBalance();
    return grid;
  } // ... create_initial_grid(...)

public:
  OS2014(const size_t num_refinements = OS2014BaseType::default_num_refinements, const int start_cell_count = 4)
    : TestCaseBaseType(create_initial_grid(OS2014BaseType::initial_refinements(), start_cell_count), num_refinements)
  {}
}; // class OS2014


# if HAVE_DUNE_GRID_MULTISCALE


template< class GridType >
class OS2014Multiscale
  : public internal::OS2014Base< GridType >
  , public MultiscaleCubeBase< GridType >
{
  typedef internal::OS2014Base< GridType > OS2014BaseType;
  typedef MultiscaleCubeBase< GridType >    TestCaseBaseType;

private:
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
  OS2014Multiscale(const std::string num_partitions = "[1 1 1]",
                    const size_t num_refinements = OS2014BaseType::default_num_refinements)
    : TestCaseBaseType(initial_grid_cfg(num_partitions), OS2014BaseType::initial_refinements(), num_refinements)
  {}


}; // class OS2014Multiscale


# endif // HAVE_DUNE_GRID_MULTISCALE

} // namespace TestCases
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_TESTCASES_OS2014_HH
