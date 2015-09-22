// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_EXAMPLES_LINEARPARABOLIC_THERMALBLOCK_HH
#define DUNE_HDD_EXAMPLES_LINEARPARABOLIC_THERMALBLOCK_HH

#include <dune/common/parallel/mpihelper.hh>

#if HAVE_ALUGRID
# include <dune/grid/alugrid.hh>
#endif

#if HAVE_DUNE_FEM
# include <dune/fem/misc/mpimanager.hh>
#endif

#include <dune/stuff/common/timedlogging.hh>
#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/la/container.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/spaces/interface.hh>

#include <dune/hdd/linearelliptic/discretizations/cg.hh>
#include <dune/hdd/linearelliptic/problems/thermalblock.hh>

namespace internal {


template< class GridType >
class Initializer
{
  typedef Dune::Stuff::Grid::Providers::Cube< GridType > GridProviderType;

public:
  Initializer(const std::string& num_grid_elements,
              const DUNE_STUFF_SSIZE_T info_log_levels,
              const DUNE_STUFF_SSIZE_T debug_log_levels,
              const bool enable_warnings,
              const bool enable_colors,
              const std::string info_color,
              const std::string debug_color,
              const std::string warn_color)
  {
    try {
      int argc = 0;
      char** argv = new char* [0];
#if HAVE_DUNE_FEM
      Dune::Fem::MPIManager::initialize(argc, argv);
#else
      Dune::MPIHelper::instance(argc, argv);
#endif
    } catch (...) {}
    try {
      DSC::TimedLogger().create(info_log_levels,
                                debug_log_levels,
                                enable_warnings,
                                enable_colors,
                                info_color,
                                debug_color,
                                warn_color);
    } catch (Dune::Stuff::Exceptions::you_are_using_this_wrong&) {}
    auto logger = DSC::TimedLogger().get("cg.thermalblock.example");
    logger.info() << "creating grid... " << std::flush;
    grid_provider_ = GridProviderType::create(GridProviderType::default_config().add(DSC::Configuration("num_elements",
                                                                                                        num_grid_elements),
                                                                                     "",
                                                                                     true));
#if HAVE_ALUGRID
    if (std::is_same< GridType, Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming > >::value)
      grid_provider_->grid().globalRefine(1);
#endif // HAVE_ALUGRID
    logger.info() << "done (has " << grid_provider_->grid().size(0) << " elements)" << std::endl;
  } // Initializer(...)

protected:
  std::unique_ptr< GridProviderType > grid_provider_;
}; // class Initializer


} // namespace internal


template< class GridImp, Dune::GDT::ChooseSpaceBackend space_backend, Dune::Stuff::LA::ChooseBackend la_backend >
class CgThermalblockExample
  : internal::Initializer< GridImp >
{
  typedef internal::Initializer< GridImp >               BaseType;
  typedef GridImp                                        GridType;
  typedef typename GridType::template Codim< 0 >::Entity EntityType;
  static const size_t                                    dimDomain = GridType::dimension;
  typedef typename GridType::ctype                       DomainFieldType;
  typedef double                                         RangeFieldType;
  static const size_t                                    dimRange = 1;
  typedef Dune::HDD::LinearElliptic::Problems::Thermalblock
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, 1 > ProblemType;
public:
  typedef Dune::HDD::LinearElliptic::Discretizations::CG< GridType, Dune::Stuff::Grid::ChooseLayer::leaf,
                                                          RangeFieldType, dimRange,
                                                          1, space_backend, la_backend > DiscretizationType;
  typedef typename DiscretizationType::VectorType                                        VectorType;

public:
  CgThermalblockExample(const std::string& num_blocks = "[1 1 1]",
                        const std::string& num_grid_elements = "[8 8 8]",
                        const std::vector< std::string >& only_these_products = {},
                        const DUNE_STUFF_SSIZE_T info_log_levels  = 0,
                        const DUNE_STUFF_SSIZE_T debug_log_levels = -1,
                        const bool enable_warnings = true,
                        const bool enable_colors   = true,
                        const std::string info_color  = DSC::TimedLogging::default_info_color(),
                        const std::string debug_color = DSC::TimedLogging::default_debug_color(),
                        const std::string warn_color  = DSC::TimedLogging::default_warning_color())
    : BaseType(num_grid_elements,
               info_log_levels,
               debug_log_levels,
               enable_warnings,
               enable_colors,
               info_color,
               debug_color,
               warn_color)
    , problem_(ProblemType::create(ProblemType::default_config().add(
                                     DSC::Configuration({"diffusion_factor.num_elements", "force.value"},
                                                        {num_blocks,                      "0"}),
                                     "",
                                     true)))
    , boundary_info_cfg_(Dune::Stuff::Grid::BoundaryInfoConfigs::AllDirichlet::default_config())
    , discretization_(*grid_provider_, boundary_info_cfg_, *problem_, -1, only_these_products)
  {
    auto logger = DSC::TimedLogger().get("cg.parabolic.example");
    logger.info() << "initializing discretization... " << std::flush;
    discretization_.init();
    logger.info() << "done (has " << discretization_.ansatz_space().mapper().size() << " DoFs)" << std::endl;
  } // ... CgThermalblockExample(...)

  DiscretizationType& discretization()
  {
    return discretization_;
  }

  void visualize(const std::string& filename_prefix) const
  {
    grid_provider_->visualize(filename_prefix + ".grid", boundary_info_cfg_);
    problem_->visualize(grid_provider_->leaf_view(), filename_prefix + ".problem", /*subsampling=*/ false);
  }

  VectorType project(const std::string& expression) const
  {
    using namespace Dune;
    if (expression.empty())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Given 'expression' is empty!");
    auto logger = DSC::TimedLogger().get("cg.parabolic.project");
    logger.info() << "projecting '" << expression << "'... " << std::flush;
    Stuff::Functions::Expression< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >
        function("x", expression);
    auto discrete_function = GDT::make_discrete_function< VectorType >(discretization_.ansatz_space());
    GDT::project(function, discrete_function);
    logger.info() << "done" << std::endl;
    return discrete_function.vector();
  } // ... project(...)

private:
  using BaseType::grid_provider_;
  std::unique_ptr< ProblemType > problem_;
  DSC::Configuration boundary_info_cfg_;
  DiscretizationType discretization_;
}; // class CgThermalblockExample


#endif DUNE_HDD_EXAMPLES_LINEARPARABOLIC_THERMALBLOCK_HH
