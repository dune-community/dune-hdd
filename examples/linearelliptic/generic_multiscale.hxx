#ifndef DUNE_HDD_EXAMPLES_LINEARELLIPTIC_GENERIC_MULTISCALE_HXX
#define DUNE_HDD_EXAMPLES_LINEARELLIPTIC_GENERIC_MULTISCALE_HXX

#include <dune/gdt/operators/prolongations.hh>

#include "generic_multiscale.hh"


template< class G, Dune::GDT::ChooseSpaceBackend sp, Dune::Stuff::LA::ChooseBackend la >
    Dune::Stuff::Common::Configuration GenericLinearellipticMultiscaleExample< G, sp, la >::
logger_options()
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

template< class G, Dune::GDT::ChooseSpaceBackend sp, Dune::Stuff::LA::ChooseBackend la >
    std::vector< std::string > GenericLinearellipticMultiscaleExample< G, sp, la >::
grid_options()
{
  return GridProvider::available();
}

template< class G, Dune::GDT::ChooseSpaceBackend sp, Dune::Stuff::LA::ChooseBackend la >
    Dune::Stuff::Common::Configuration GenericLinearellipticMultiscaleExample< G, sp, la >::
grid_options(const std::string& type)
{
  return GridProvider::default_config(type);
}

template< class G, Dune::GDT::ChooseSpaceBackend sp, Dune::Stuff::LA::ChooseBackend la >
    std::vector< std::string > GenericLinearellipticMultiscaleExample< G, sp, la >::
boundary_options()
{
  return BoundaryProvider::available();
}

template< class G, Dune::GDT::ChooseSpaceBackend sp, Dune::Stuff::LA::ChooseBackend la >
    Dune::Stuff::Common::Configuration GenericLinearellipticMultiscaleExample< G, sp, la >::
boundary_options(const std::string& type)
{
  return BoundaryProvider::default_config(type);
}

template< class G, Dune::GDT::ChooseSpaceBackend sp, Dune::Stuff::LA::ChooseBackend la >
    std::vector< std::string > GenericLinearellipticMultiscaleExample< G, sp, la >::
problem_options()
{
  return ProblemProvider::available();
}

template< class G, Dune::GDT::ChooseSpaceBackend sp, Dune::Stuff::LA::ChooseBackend la >
    Dune::Stuff::Common::Configuration GenericLinearellipticMultiscaleExample< G, sp, la >::
problem_options(const std::string& type)
{
  return ProblemProvider::default_config(type);
}

template< class G, Dune::GDT::ChooseSpaceBackend sp, Dune::Stuff::LA::ChooseBackend la >
    std::vector< std::string > GenericLinearellipticMultiscaleExample< G, sp, la >::
solver_options()
{
  return SolverProvider::types();
}

template< class G, Dune::GDT::ChooseSpaceBackend sp, Dune::Stuff::LA::ChooseBackend la >
    Dune::Stuff::Common::Configuration GenericLinearellipticMultiscaleExample< G, sp, la >::
solver_options(const std::string& type)
{
  return SolverProvider::options(type);
}

template< class G, Dune::GDT::ChooseSpaceBackend sp, Dune::Stuff::LA::ChooseBackend la >
    GenericLinearellipticMultiscaleExample< G, sp, la >::
GenericLinearellipticMultiscaleExample(const Dune::Stuff::Common::Configuration& logger_cfg,
                                       const Dune::Stuff::Common::Configuration& grid_cfg,
                                       const Dune::Stuff::Common::Configuration& boundary_cfg,
                                       const Dune::Stuff::Common::Configuration& problem_cfg)
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
  auto logger = DSC::TimedLogger().get("example.linearelliptic.genericmultiscale");
  logger.info() << "creating grid (" << grid_cfg.get< std::string >("type") << "):" << std::endl;
  grid_ = GridProvider::create(grid_cfg.get< std::string >("type"), grid_cfg);
  logger.info() << "  done (has " << grid_->grid().size(0) << " elements)" << std::endl;

  logger.info() << "creating problem (" << problem_cfg.get< std::string >("type") << ")... " << std::endl;
  problem_= ProblemProvider::create(problem_cfg.get< std::string >("type"), problem_cfg);

  logger.info() << "creating discretization:" << std::endl;
  discretization_ = DSC::make_unique< DiscretizationType >(*grid_,
                                                           boundary_cfg_,
                                                           *problem_,
                                                           /*only_these_products=*/std::vector<std::string>({"l2", "h1", "elliptic"}));
  discretization_->init(/*prune=*/false);
  logger.info() << "  done (has " << discretization_->ansatz_space().mapper().size() << " DoFs)" << std::endl;
} // GenericLinearellipticMultiscaleExample(...)

template< class G, Dune::GDT::ChooseSpaceBackend sp, Dune::Stuff::LA::ChooseBackend la >
    typename GenericLinearellipticMultiscaleExample< G, sp, la >::DiscretizationType& GenericLinearellipticMultiscaleExample< G, sp, la >::
discretization()
{
  return *discretization_;
}

template< class G, Dune::GDT::ChooseSpaceBackend sp, Dune::Stuff::LA::ChooseBackend la >
    void GenericLinearellipticMultiscaleExample< G, sp, la >::
visualize(const std::string& filename_prefix) const
{
  if (filename_prefix.empty())
    DUNE_THROW(Dune::Stuff::Exceptions::wrong_input_given, "Given filename prefix must not be empty!");
  grid_->visualize(filename_prefix + ".grid", boundary_cfg_);
  problem_->visualize(grid_->leaf_view(), filename_prefix + ".problem", /*subsampling=*/ false);
}

template< class G, Dune::GDT::ChooseSpaceBackend sp, Dune::Stuff::LA::ChooseBackend la >
    typename GenericLinearellipticMultiscaleExample< G, sp, la >::VectorType GenericLinearellipticMultiscaleExample< G, sp, la >::
project(const std::string& expression) const
{
  using namespace Dune;
  if (expression.empty())
    DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Given expression must not be empty!");
  auto logger = DSC::TimedLogger().get("example.linearelliptic.genericmultiscale.project");
  logger.info() << "projecting '" << expression << "'... " << std::endl;
  auto discrete_function = GDT::make_discrete_function< VectorType >(discretization_->ansatz_space());
  GDT::Operators::apply_projection(Stuff::Functions::Expression< E, D, d, R, r >("x", expression), discrete_function);
  return discrete_function.vector();
} // ... project(...)

template< class G, Dune::GDT::ChooseSpaceBackend sp, Dune::Stuff::LA::ChooseBackend la >
    typename GenericLinearellipticMultiscaleExample< G, sp, la >::VectorType GenericLinearellipticMultiscaleExample< G, sp, la >::
prolong(const typename GenericLinearellipticMultiscaleExample< G, sp, la >::DiscretizationType& source_disc,
        const typename GenericLinearellipticMultiscaleExample< G, sp, la >::VectorType& source_vec)
{
  using namespace Dune;
  auto logger = DSC::TimedLogger().get("example.linearelliptic.genericmultiscale.project");
  auto source_func = GDT::make_const_discrete_function(source_disc.ansatz_space(), source_vec);
  auto range_func = GDT::make_discrete_function< VectorType >(discretization_->ansatz_space());
  if (source_vec.size() >= range_func.vector().size())
    logger.warn() << "prolonging from space of size " << source_vec.size() << " onto space of size "
                  << range_func.vector().size() << ", this might not be what you want!" << std::endl;
  else
    logger.info() << "prolonging from space of size " << source_vec.size() << " onto space of size "
                  << range_func.vector().size() << "... " << std::endl;
  GDT::Operators::prolong(source_func, range_func);
  return range_func.vector();
} // ... project(...)


#endif // DUNE_HDD_EXAMPLES_LINEARELLIPTIC_GENERIC_MULTISCALE_HXX
