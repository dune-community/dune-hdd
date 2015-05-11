// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_TESTCASES_ESV2007_HH
#define DUNE_HDD_LINEARELLIPTIC_TESTCASES_ESV2007_HH

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
class ESV2007Base
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
    int ret = 3;
#if HAVE_ALUGRID
    if (std::is_same< GridType, ALUConformGrid< 2, 2 > >::value
        || std::is_same< GridType, ALUGrid< 2, 2, simplex, conforming > >::value)
      ret += 1;
#endif // HAVE_ALUGRID
    return ret;
  } // ... initial_refinements()

public:
  ESV2007Base()
    : boundary_info_cfg_(Stuff::Grid::BoundaryInfoConfigs::AllDirichlet::default_config())
    , problem_(3)
    , exact_solution_(2, problem_.static_id() + ".exact_solution")
  {}

  void print_header(std::ostream& out = DSC_LOG_INFO_0) const
  {
    out << "+==================================================================+\n"
        << "|+================================================================+|\n"
        << "||  Testcase ESV2007: smooth data, homogeneous dirichlet          ||\n"
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
}; // class ESV2007


} // namespace internal


template< class GridType >
class ESV2007
  : public internal::ESV2007Base< GridType >
  , public Base< GridType >
{
  typedef internal::ESV2007Base< GridType > ESV2007BaseType;
  typedef Base< GridType >                  TestCaseBaseType;

private:
  static std::shared_ptr< GridType > create_initial_grid(const int refinements)
  {
    Stuff::Grid::Providers::Cube< GridType > grid_provider(-1, 1, 4);
    auto grid = grid_provider.grid_ptr();
    grid->globalRefine(refinements);
    return grid;
  } // ... create_initial_grid(...)

public:
  ESV2007(const size_t num_refinements = ESV2007BaseType::default_num_refinements,
          const size_t initial_refine = ESV2007BaseType::initial_refinements())
    : TestCaseBaseType(create_initial_grid(initial_refine), num_refinements)
  {}
}; // class ESV2007


# if HAVE_DUNE_GRID_MULTISCALE


template< class GridType >
class ESV2007Multiscale
  : public internal::ESV2007Base< GridType >
  , public MultiscaleCubeBase< GridType >
{
  typedef internal::ESV2007Base< GridType > ESV2007BaseType;
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
  ESV2007Multiscale(const std::string num_partitions = "[1 1 1]",
                    const size_t num_refinements = ESV2007BaseType::default_num_refinements)
    : TestCaseBaseType(initial_grid_cfg(num_partitions), ESV2007BaseType::initial_refinements(), num_refinements)
  {}


}; // class ESV2007Multiscale


# endif // HAVE_DUNE_GRID_MULTISCALE

} // namespace TestCases
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_TESTCASES_ESV2007_HH
