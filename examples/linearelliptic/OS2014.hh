// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_EXAMPLES_LINEARELLIPTIC_OS2014_HH
#define DUNE_HDD_EXAMPLES_LINEARELLIPTIC_OS2014_HH

#define DUNE_HDD_LINEARELLIPTIC_TESTCASES_BASE_DISABLE_WARNING

#include <boost/numeric/conversion/cast.hpp>

#include <dune/fem/misc/mpimanager.hh>

#include <dune/stuff/common/logging.hh>

#include <dune/pymor/parameters/base.hh>

#include <dune/hdd/linearelliptic/testcases/spe10.hh>
#include <dune/hdd/linearelliptic/discretizations/block-swipdg.hh>

namespace internal {


class Initializer
{
public:
  Initializer(const ssize_t info_log_levels,
              const ssize_t debug_log_levels,
              const bool enable_warnings,
              const bool enable_colors,
              const std::string info_color,
              const std::string debug_color,
              const std::string warn_color)
  {
    try {
      int argc = 0;
      char** argv = new char* [0];
      Dune::Fem::MPIManager::initialize(argc, argv);
    } catch (...) {}
    DSC::TimedLogger().create(info_log_levels,
                              debug_log_levels,
                              enable_warnings,
                              enable_colors,
                              info_color,
                              debug_color,
                              warn_color);
    DSC::TimedLogger().get("OS2014.spe10model1example").info() << "creating grid and problem... " << std::endl;
  }
}; // class Initializer


} // namespace internal


template< class GridType >
class OS2014Spe10Model1Example
  : internal::Initializer
{
  static_assert(GridType::dimension == 2, "Only available in 2d!");
public:
  typedef Dune::HDD::LinearElliptic::TestCases::Spe10::ParametricBlockModel1< GridType > TestCaseType;
  typedef Dune::HDD::LinearElliptic::Discretizations::BlockSWIPDG< GridType, double, 1 > DiscretizationType;

  OS2014Spe10Model1Example(const std::string partitioning = "[1 1 1]",
                           const DUNE_STUFF_SSIZE_T num_refinements = 0,
                           const ssize_t info_log_levels  = 0,
                           const ssize_t debug_log_levels = -1,
                           const bool enable_warnings = true,
                           const bool enable_colors   = true,
                           const std::string info_color  = DSC::TimedLogging::default_info_color(),
                           const std::string debug_color = DSC::TimedLogging::default_debug_color(),
                           const std::string warn_color  = DSC::TimedLogging::default_warning_color())
    : internal::Initializer(info_log_levels,
                            debug_log_levels,
                            enable_warnings,
                            enable_colors,
                            info_color,
                            debug_color,
                            warn_color)
    , test_case_({{"mu", Dune::Pymor::Parameter("mu", 1)},     // <- it does not matter which parameters we give to the
                  {"mu_hat", Dune::Pymor::Parameter("mu", 1)}, //    test case here, since we use test_case_.problem()
                  {"mu_bar", Dune::Pymor::Parameter("mu", 1)}, //    which is the parametric problem anyway
                  {"mu_minimizing", Dune::Pymor::Parameter("mu", 1)}},
                 partitioning,
                 boost::numeric_cast< size_t >(num_refinements))
    , discretization_(*test_case_.reference_provider(),
                      test_case_.boundary_info(),
                      test_case_.problem())
  {
    auto logger = DSC::TimedLogger().get("OS2014.spe10model1example");
    logger.info() << "initializing discretization... " << std::flush;
    discretization_.init();
    logger.info() << "done (grid has " << discretization_.grid_view()->indexSet().size(0)
                  << " elements, discretization has " << discretization_.ansatz_space()->mapper().size() << " DoFs)"
                  << std::endl;
  } // ... OS2014Spe10Model1Example(...)

  const TestCaseType& test_case() const
  {
    return test_case_;
  }

  DiscretizationType& discretization()
  {
    return discretization_;
  }

  DiscretizationType* discretization_and_return_ptr() const
  {
    return new DiscretizationType(discretization_);
  }

public:
  TestCaseType test_case_;
  DiscretizationType discretization_;
}; // class OS2014Spe10Model1Example


#endif // DUNE_HDD_EXAMPLES_LINEARELLIPTIC_OS2014_HH
