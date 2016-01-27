#ifndef DUNE_HDD_EXAMPLES_LINEARELLIPTIC_GENERIC_HH
#define DUNE_HDD_EXAMPLES_LINEARELLIPTIC_GENERIC_HH

#include <dune/common/parallel/mpihelper.hh>

#if HAVE_DUNE_FEM
# include <dune/fem/misc/mpimanager.hh>
#endif

#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/timedlogging.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/grid/provider.hh>
#include <dune/stuff/la/container/container-interface.hh>
#include <dune/stuff/la/solver.hh>

#include <dune/gdt/operators/projections.hh>
#if HAVE_DUNE_FEM
#include <dune/gdt/spaces/cg.hh>
#include <dune/gdt/operators/darcy.hh>
#endif

#include <dune/hdd/linearelliptic/discretizations/cg.hh>
#include <dune/hdd/linearelliptic/problems.hh>


template< class GridType, Dune::GDT::ChooseSpaceBackend space_backend, Dune::Stuff::LA::ChooseBackend la_backend >
class GenericLinearellipticExample
{
  typedef typename GridType::template Codim< 0 >::Entity E;
  typedef typename GridType::ctype D;
  static const size_t d = GridType::dimension;
  typedef double R;
  static const size_t r = 1;
public:
  typedef Dune::HDD::LinearElliptic::Discretizations::CG< GridType, Dune::Stuff::Grid::ChooseLayer::leaf, R, 1, 1,
                                                          space_backend, la_backend > DiscretizationType;
  typedef typename DiscretizationType::VectorType                                     VectorType;
private:
  typedef typename DiscretizationType::MatrixType                                        MatrixType;
  typedef Dune::Stuff::GridProviders< GridType >                                         GridProvider;
  typedef Dune::Stuff::Grid::BoundaryInfoProvider< typename GridType::LeafIntersection > BoundaryProvider;
  typedef Dune::HDD::LinearElliptic::ProblemsProvider< E, D, d, R, r >                   ProblemProvider;
  typedef Dune::Stuff::LA::Solver< MatrixType >                                          SolverProvider;
#if HAVE_DUNE_FEM
  typedef Dune::GDT::Spaces::CGProvider< GridType, Dune::Stuff::Grid::ChooseLayer::leaf,
                                         Dune::GDT::ChooseSpaceBackend::fem, 1, R, d > VelocitySpaceProvider;
  typedef typename VelocitySpaceProvider::Type VelocitySpaceType;
#endif

public:
  static Dune::Stuff::Common::Configuration logger_options()
  {
    Dune::Stuff::Common::Configuration ret;
    ret["info"]            = "0";
    ret["debug"]           = "-1";
    ret["enable_warnings"] = "true";
    ret["enable_colors"]   = "true";
    ret["info_color"]      = Dune::Stuff::Common::TimedLogging::default_info_color();
    ret["debug_color"]     = Dune::Stuff::Common::TimedLogging::default_debug_color();
    ret["warn_color"]      = Dune::Stuff::Common::TimedLogging::default_warning_color();
    return ret;
  }

  static std::vector< std::string > grid_options()
  {
    return GridProvider::available();
  }

  static Dune::Stuff::Common::Configuration grid_options(const std::string& type)
  {
    return GridProvider::default_config(type);
  }

  static std::vector< std::string > boundary_options()
  {
    return BoundaryProvider::available();
  }

  static Dune::Stuff::Common::Configuration boundary_options(const std::string& type)
  {
    return BoundaryProvider::default_config(type);
  }

  static std::vector< std::string > problem_options()
  {
    return ProblemProvider::available();
  }

  static Dune::Stuff::Common::Configuration problem_options(const std::string& type)
  {
    return ProblemProvider::default_config(type);
  }

  static std::vector< std::string > solver_options()
  {
    return SolverProvider::types();
  }

  static Dune::Stuff::Common::Configuration solver_options(const std::string& type)
  {
    return SolverProvider::options(type);
  }

  GenericLinearellipticExample(const Dune::Stuff::Common::Configuration& logger_cfg = Dune::Stuff::Common::Configuration(),
                               const Dune::Stuff::Common::Configuration& grid_cfg = Dune::Stuff::Common::Configuration(),
                               const Dune::Stuff::Common::Configuration& boundary_cfg = Dune::Stuff::Common::Configuration(),
                               const Dune::Stuff::Common::Configuration& problem_cfg = Dune::Stuff::Common::Configuration())
    : boundary_cfg_(boundary_cfg)
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
      DSC::TimedLogger().create(logger_cfg.get("info",            logger_options().template get< ssize_t     >("info")),
                                logger_cfg.get("debug",           logger_options().template get< ssize_t     >("debug")),
                                logger_cfg.get("enable_warnings", logger_options().template get< bool        >("enable_warnings")),
                                logger_cfg.get("enable_colors",   logger_options().template get< bool        >("enable_colors")),
                                logger_cfg.get("info_color",      logger_options().template get< std::string >("info_color")),
                                logger_cfg.get("debug_color",     logger_options().template get< std::string >("debug_color")),
                                logger_cfg.get("warn_color",      logger_options().template get< std::string >("warn_color")));
    } catch (Dune::Stuff::Exceptions::you_are_using_this_wrong&) {}
    auto logger = DSC::TimedLogger().get("example.linearelliptic.generic");
    logger.info() << "creating grid (" << grid_cfg.get< std::string >("type") << ")... " << std::flush;
    grid_ = GridProvider::create(grid_cfg.get< std::string >("type"), grid_cfg);
    logger.info() << "done (has " << grid_->grid().size(0) << " elements)" << std::endl;

    logger.info() << "creating problem (" << problem_cfg.get< std::string >("type") << ")... " << std::flush;
    problem_= ProblemProvider::create(problem_cfg.get< std::string >("type"), problem_cfg);
    logger.info() << "done" << std::endl;

    logger.info() << "creating discretization... " << std::flush;
    discretization_ = DSC::make_unique< DiscretizationType >(*grid_,
                                                             boundary_cfg_,
                                                             *problem_,
                                                             -1);
    discretization_->init();
#if HAVE_DUNE_FEM
    velocity_space_ = DSC::make_unique< VelocitySpaceType >(VelocitySpaceProvider::create(*grid_));
#endif
    logger.info() << "done (has " << discretization_->ansatz_space().mapper().size() << " DoFs)" << std::endl;
  } // GenericLinearellipticExample(...)

  ~GenericLinearellipticExample()
  {
    DSC::TimedLogger().get("example.linearelliptic.generic").info() << "finished" << std::endl;
  }

  DiscretizationType& discretization()
  {
    return *discretization_;
  }

  void visualize_grid(const std::string& filename) const
  {
    auto logger = DSC::TimedLogger().get("example.linearelliptic.generic");
    logger.info() << "visualizing grid... " << std::flush;
    if (filename.empty())
      DUNE_THROW(Dune::Stuff::Exceptions::wrong_input_given, "Given filename prefix must not be empty!");
    grid_->visualize(filename, boundary_cfg_);
    logger.info() << "done" << std::endl;
  }

  void visualize_problem(const std::string& filename,
                         const Dune::Pymor::Parameter& mu = Dune::Pymor::Parameter()) const
  {
    auto logger = DSC::TimedLogger().get("example.linearelliptic.generic");
    logger.info() << "visualizing problem (mu = " << mu << ")... " << std::flush;
    if (filename.empty())
      DUNE_THROW(Dune::Stuff::Exceptions::wrong_input_given, "Given filename prefix must not be empty!");
    if (mu.empty()) {
      problem_->visualize(grid_->leaf_view(), filename, /*subsampling=*/ false);
    } else {
      problem_->with_mu(mu)->visualize(grid_->leaf_view(), filename, /*subsampling=*/ false);
    }
    logger.info() << "done" << std::endl;
  }

#if HAVE_DUNE_FEM
  void visualize_darcy_velocity(const VectorType& vector,
                                const std::string& filename,
                                const std::string& name,
                                const Dune::Pymor::Parameter& mu = Dune::Pymor::Parameter()) const
  {
    using namespace Dune;
    using namespace Dune::GDT;
    auto logger = DSC::TimedLogger().get("example.linearelliptic.generic");
    logger.info() << "reconstructing darcy velocity... " << std::flush;
    auto problem_mu = problem_->with_mu(mu);
    auto diffusion = *problem_mu->diffusion_factor()->affine_part() * *problem_mu->diffusion_tensor()->affine_part();
    auto pressure = make_const_discrete_function(discretization_->ansatz_space(), vector);
    auto velocity = make_discrete_function< VectorType >(*velocity_space_, name);
    const auto& grid_view = discretization_->grid_view();
    auto darcy_operator = Operators::make_darcy(grid_view, diffusion);
    darcy_operator->apply(pressure, velocity);
    logger.info() << "done" << std::endl;
    logger.info() << "visualizing darcy velocity... " << std::flush;
    velocity.visualize(filename);
    logger.info() << "done" << std::endl;
  } // ... visualize_darcy_velocity(...)
#endif // HAVE_DUNE_FEM

  VectorType project(const std::string& expression) const
  {
    using namespace Dune;
    if (expression.empty())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Given expression must not be empty!");
    auto logger = DSC::TimedLogger().get("example.linearelliptic.generic.project");
    logger.info() << "projecting '" << expression << "'... " << std::flush;
    auto discrete_function = GDT::make_discrete_function< VectorType >(discretization_->ansatz_space());
    GDT::Operators::apply_projection(Stuff::Functions::Expression< E, D, d, R, r >("x", expression), discrete_function);
    logger.info() << "done" << std::endl;
    return discrete_function.vector();
  } // ... project(...)

private:
  const Dune::Stuff::Common::Configuration boundary_cfg_;
  std::unique_ptr< Dune::Stuff::Grid::ProviderInterface< GridType > > grid_;
  std::unique_ptr< Dune::HDD::LinearElliptic::ProblemInterface< E, D, d, R, r > > problem_;
  std::unique_ptr< DiscretizationType > discretization_;
#if HAVE_DUNE_FEM
  std::unique_ptr< VelocitySpaceType > velocity_space_;
#endif
}; // class GenericLinearellipticExample


#endif // DUNE_HDD_EXAMPLES_LINEARELLIPTIC_GENERIC_HH
