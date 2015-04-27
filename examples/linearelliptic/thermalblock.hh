// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_EXAMPLES_LINEARELLIPTIC_THERMALBLOCK_HH
#define DUNE_HDD_EXAMPLES_LINEARELLIPTIC_THERMALBLOCK_HH

#include <boost/numeric/conversion/cast.hpp>

#if HAVE_DUNE_FEM
# include <dune/fem/misc/mpimanager.hh>
#endif

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/grid/provider/cube.hh>

#include <dune/hdd/linearelliptic/problems/thermalblock.hh>
#include <dune/hdd/linearelliptic/discretizations/cg.hh>

namespace internal {


class Initializer
{
public:
  Initializer(const DUNE_STUFF_SSIZE_T info_log_levels,
              const DUNE_STUFF_SSIZE_T debug_log_levels,
              const bool enable_warnings,
              const bool enable_colors,
              const std::string info_color,
              const std::string debug_color,
              const std::string warn_color)
  {
#if HAVE_DUNE_FEM
    try {
      int argc = 0;
      char** argv = new char* [0];
      Dune::Fem::MPIManager::initialize(argc, argv);
    } catch (...) {}
#endif // HAVE_DUNE_FEM
    DSC::TimedLogger().create(info_log_levels,
                              debug_log_levels,
                              enable_warnings,
                              enable_colors,
                              info_color,
                              debug_color,
                              warn_color);
    DSC::TimedLogger().get("cg.thermalblock.example").info() << "creating grid and problem... " << std::endl;
  }
}; // class Initializer


} // namespace internal


template< class GridImp, Dune::GDT::ChooseSpaceBackend space_backend, Dune::Stuff::LA::ChooseBackend la_backend >
class CgExample
  : internal::Initializer
{
  typedef GridImp                                        GridType;
  typedef typename GridType::template Codim< 0 >::Entity EntityType;
  static const size_t                                    dimDomain = GridType::dimension;
  typedef typename GridType::ctype                       DomainFieldType;
  typedef double                                         RangeFieldType;
  static const size_t                                    dimRange = 1;
  typedef Dune::HDD::LinearElliptic::Problems::Thermalblock
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, 1 > ProblemType;
  typedef Dune::Stuff::Grid::Providers::Cube< GridType >            GridProviderType;
public:
  typedef Dune::HDD::LinearElliptic::Discretizations::CG< GridType, Dune::Stuff::Grid::ChooseLayer::leaf,
                                                          RangeFieldType, dimRange, 1,
                                                          space_backend, la_backend > DiscretizationType;

  CgExample(const std::string& num_blocks        = "[2 2 2]",
            const std::string& num_grid_elements = "[32 32 32]",
            const DUNE_STUFF_SSIZE_T info_log_levels  = 0,
            const DUNE_STUFF_SSIZE_T debug_log_levels = -1,
            const bool enable_warnings = true,
            const bool enable_colors   = true,
            const std::string info_color  = DSC::TimedLogging::default_info_color(),
            const std::string debug_color = DSC::TimedLogging::default_debug_color(),
            const std::string warn_color  = DSC::TimedLogging::default_warning_color())
    : Initializer(info_log_levels,
                  debug_log_levels,
                  enable_warnings,
                  enable_colors,
                  info_color,
                  debug_color,
                  warn_color)
    , problem_(ProblemType::create(ProblemType::default_config().add(DSC::Configuration("diffusion_factor.num_elements",
                                                                                        num_blocks),
                                                                     "",
                                                                     true)))
    , grid_provider_(GridProviderType::create(GridProviderType::default_config().add(DSC::Configuration("num_elements",
                                                                                                        num_grid_elements),
                                                                                     "",
                                                                                     true)))
    , discretization_(*grid_provider_, Dune::Stuff::Grid::BoundaryInfoConfigs::AllDirichlet::default_config(), *problem_)
  {
    auto logger = DSC::TimedLogger().get("cg.thermalblock.example");
    logger.info() << "initializing discretization... " << std::flush;
    discretization_.init();
    logger.info() << "done (grid has " << discretization_.grid_view().indexSet().size(0)
                  << " elements, discretization has " << discretization_.ansatz_space().mapper().size() << " DoFs)"
                  << std::endl;
  } // ... CgExample(...)

  DiscretizationType& discretization()
  {
    return discretization_;
  }

private:
  std::unique_ptr< ProblemType > problem_;
  std::unique_ptr< GridProviderType > grid_provider_;
  DiscretizationType discretization_;
}; // class CgExample


#endif // DUNE_HDD_EXAMPLES_LINEARELLIPTIC_THERMALBLOCK_HH
