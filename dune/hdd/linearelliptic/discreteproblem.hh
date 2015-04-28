// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_DISCRETEPROBLEM_HH
#define DUNE_HDD_LINEARELLIPTIC_DISCRETEPROBLEM_HH

#include <memory>
#include <limits>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/stuff/common/disable_warnings.hh>
# include <dune/common/parallel/mpihelper.hh>
# include <dune/common/timer.hh>

# if HAVE_DUNE_FEM
#   include <dune/fem/misc/mpimanager.hh>
# endif
#include <dune/stuff/common/reenable_warnings.hh>

#if HAVE_DUNE_GRID_MULTISCALE
# include <dune/grid/multiscale/provider.hh>
#endif

#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/string.hh>
#include <dune/stuff/grid/provider.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/grid/layers.hh>
#include <dune/stuff/common/memory.hh>

#include "problems.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {


template< class GridImp >
class DiscreteProblem
{
public:
  typedef GridImp GridType;
  typedef Stuff::Grid::ProviderInterface< GridType >  GridProviderType;
private:
  typedef Stuff::GridProviders< GridType > GridProviders;
  typedef typename GridType::LeafIntersection IntersectionType;
  typedef Stuff::Grid::BoundaryInfoProvider< IntersectionType > BoundaryInfoProvider;
  typedef typename GridType::template Codim< 0 >::Entity EntityType;
  typedef typename GridType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridType::dimension;
public:
  typedef double RangeFieldType;
  static const unsigned int dimRange = 1;
  typedef LinearElliptic::ProblemInterface< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > ProblemType;
  typedef LinearElliptic::ProblemsProvider< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > ProblemsProvider;

  static void write_config(const std::string filename, const std::string id)
  {
    std::ofstream file;
    file.open(filename);
    file << "[" << id << "]" << std::endl;
    write_keys_to_file("gridprovider", GridProviders::available(), file);
    write_keys_to_file("boundaryinfo", BoundaryInfoProvider::available(), file);
    write_keys_to_file("problem", ProblemsProvider::available(), file);
    file << "[logging]" << std::endl;
    file << "info  = true" << std::endl;
    file << "debug = true" << std::endl;
    file << "file  = false" << std::endl;
    file << "visualize = true" << std::endl;
    file << "[parameter]" << std::endl;
    file << "0.diffusion_factor = [0.1 0.1 1.0 1.0]" << std::endl;
    file << "1.diffusion_factor = [1.0 1.0 0.1 0.1]" << std::endl;
    write_config_to_file< GridProviders >(file);
    write_config_to_file< BoundaryInfoProvider >(file);
    write_config_to_file< ProblemsProvider >(file);
    file.close();
  } // ... write_config(...)

  DiscreteProblem(const std::string id, const std::vector< std::string >& arguments)
  {
    // mpi
    assert(arguments.size() < std::numeric_limits< int >::max());
    int argc = boost::numeric_cast< int >(arguments.size());
    char** argv = Stuff::Common::vectorToMainArgs(arguments);
#if HAVE_DUNE_FEM
    Fem::MPIManager::initialize(argc, argv);
#else
    MPIHelper::instance(argc, argv);
#endif

    // configuration
    config_ = Stuff::Common::Configuration(argc, argv, id + ".cfg");
    if (!config_.has_sub(id))
      DUNE_THROW(Stuff::Exceptions::configuration_error,
                 "Missing sub '" << id << "' in the following Configuration:\n\n" << config_);
    filename_ = config_.get(id + ".filename", id);

    // logger
    const Stuff::Common::Configuration logger_config = config_.sub("logging");
    int log_flags = Stuff::Common::LOG_CONSOLE;
    debug_logging_ = logger_config.get< bool >("debug", false);
    if (logger_config.get< bool >("info"))
      log_flags = log_flags | Stuff::Common::LOG_INFO;
    if (debug_logging_)
      log_flags = log_flags | Stuff::Common::LOG_DEBUG;
    if (logger_config.get< bool >("file", false))
      log_flags = log_flags | Stuff::Common::LOG_FILE;
    Stuff::Common::Logger().create(log_flags, id, "", "");
    auto& info  = Stuff::Common::Logger().info();

    Timer timer;
    const std::string griprovider_type = config_.get< std::string >(id + ".gridprovider");
    info << "creating grid with '" << griprovider_type << "'... " << std::flush;
    grid_provider_ = GridProviders::create(griprovider_type, config_);
    const auto grid_view = grid_provider_->leaf_view();
    info << " done (took " << timer.elapsed()
         << "s, has " << grid_view.indexSet().size(0) << " element";
    if (grid_view.indexSet().size(0) > 1)
      info << "s";
    info << ")" << std::endl;

    const std::string boundary_info_type = config_.get< std::string >(id + ".boundaryinfo");
    if (config_.has_sub(boundary_info_type))
      boundary_info_ = config_.sub(boundary_info_type);
    else
      boundary_info_ = Stuff::Common::Configuration("type", boundary_info_type);

    info << "setting up ";
    timer.reset();
    const std::string problem_type = config_.get< std::string >(id + ".problem");
    if (!debug_logging_)
      info << "'" << problem_type << "'... " << std::flush;
    problem_ = ProblemsProvider::create(problem_type, config_);
    if (debug_logging_)
      info << *problem_ << std::endl;
    else
      info << "done (took " << timer.elapsed() << "s)" << std::endl;

    if (logger_config.get("visualize", true)) {
      info << "visualizing grid and problem... " << std::flush;
      timer.reset();
      grid_provider_->visualize(boundary_info_, filename_ + ".grid");
      problem_->visualize(grid_view, filename_ + ".problem");
      info << "done (took " << timer.elapsed() << "s)" << std::endl;
    } // if (visualize)
  } // DiscreteProblem

  std::string filename() const
  {
    return filename_;
  }

  const Stuff::Common::Configuration& config() const
  {
    return config_;
  }

  bool debug_logging() const
  {
    return debug_logging_;
  }

  GridProviderType& grid_provider()
  {
    return *grid_provider_;
  }

  const GridProviderType& grid_provider() const
  {
    return *grid_provider_;
  }

  const Stuff::Common::Configuration& boundary_info() const
  {
    return boundary_info_;
  }

  const ProblemType& problem() const
  {
    return *problem_;
  }

private:
  static void write_keys_to_file(const std::string name, const std::vector< std::string > keys, std::ofstream& file)
  {
    std::string whitespace = Stuff::Common::whitespaceify(name + " = ");
    file << name << " = " << keys[0] << std::endl;
    for (size_t ii = 1; ii < keys.size(); ++ii)
      file << whitespace << keys[ii] << std::endl;
  } // ... write_keys_to_file(...)

  template< class ConfigProvider >
  static void write_config_to_file(std::ofstream& file)
  {
    for (const auto& type : ConfigProvider::available()) {
      auto config = ConfigProvider::default_config(type, type);
      if (!config.empty())
      file << config;
    }
  } // ... write_config_to_file(...)

  std::string filename_;
  Stuff::Common::Configuration config_;
  bool debug_logging_;
  std::unique_ptr< GridProviderType > grid_provider_;
  Stuff::Common::Configuration boundary_info_;
  std::unique_ptr< const ProblemType > problem_;
}; // class DiscreteProblem


#if HAVE_DUNE_GRID_MULTISCALE


template< class GridImp >
class DiscreteBlockProblem
{
public:
  typedef GridImp GridType;
  typedef grid::Multiscale::ProviderInterface< GridType > GridProviderType;
  typedef grid::Multiscale::MsGridProviders< GridType > GridProviders;
  typedef typename GridType::LeafIntersection IntersectionType;
  typedef Stuff::Grid::BoundaryInfoProvider< IntersectionType > BoundaryInfoProvider;
  typedef typename GridType::template Codim< 0 >::Entity EntityType;
  typedef typename GridType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridType::dimension;
  typedef double RangeFieldType;
  static const unsigned int dimRange = 1;
  typedef LinearElliptic::ProblemInterface
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > ProblemType;
  typedef LinearElliptic::ProblemsProvider
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > ProblemsProvider;

  static void write_config(const std::string filename,
                           const std::string id,
                           const Stuff::Common::Configuration& problem_cfg = Stuff::Common::Configuration(),
                           const Stuff::Common::Configuration& grid_cfg = Stuff::Common::Configuration(),
                           const Stuff::Common::Configuration& additional_cfg = Stuff::Common::Configuration())
  {
    std::ofstream file;
    file.open(filename);
    file << "[" << id << "]" << std::endl;
    if (grid_cfg.has_key("type"))
      file << "gridprovider = " << grid_cfg.get< std::string >("type") << std::endl;
    else
      write_keys_to_file("gridprovider", GridProviders::available(), file);
    if (problem_cfg.has_key("type"))
      file << "problem = " << problem_cfg.get< std::string >("type") << std::endl;
    else
      write_keys_to_file("problem", ProblemsProvider::available(), file);
    file << "[logging]" << std::endl;
    file << "info  = true" << std::endl;
    file << "debug = false" << std::endl;
    file << "file  = false" << std::endl;
    file << "visualize = true" << std::endl;
    if (grid_cfg.has_key("type"))
      file << Stuff::Common::Configuration(grid_cfg, grid_cfg.get< std::string >("type"));
    else
      write_config_to_file< GridProviders >(file);
    if (problem_cfg.has_key("type"))
      file << Stuff::Common::Configuration(problem_cfg, problem_cfg.get< std::string >("type"));
    else
      write_config_to_file< ProblemsProvider >(file);
    file << additional_cfg;
    file.close();
  } // ... write_config(...)

  DiscreteBlockProblem(const std::string id, const std::vector< std::string >& arguments)
  {
    // mpi
    int argc = boost::numeric_cast< int >(arguments.size());
    char** argv = Stuff::Common::vectorToMainArgs(arguments);
#if HAVE_DUNE_FEM
    Fem::MPIManager::initialize(argc, argv);
#else
    MPIHelper::instance(argc, argv);
#endif

    // configuration
    config_ = Stuff::Common::Configuration(argc, argv, id + ".cfg");
    if (!config_.has_sub(id))
      DUNE_THROW(Stuff::Exceptions::configuration_error,
                 "Missing sub '" << id << "' in the following Configuration:\n\n" << config_);
    filename_ = config_.get(id + ".filename", id);

    // logger
    const Stuff::Common::Configuration logger_config = config_.sub("logging");
    int log_flags = Stuff::Common::LOG_CONSOLE;
    debug_logging_ = logger_config.get< bool >("debug", false);
    if (logger_config.get< bool >("info"))
      log_flags = log_flags | Stuff::Common::LOG_INFO;
    if (debug_logging_)
      log_flags = log_flags | Stuff::Common::LOG_DEBUG;
    if (logger_config.get< bool >("file", false))
      log_flags = log_flags | Stuff::Common::LOG_FILE;
    Stuff::Common::Logger().create(log_flags, id, "", "");
    auto& info  = Stuff::Common::Logger().info();

    Timer timer;
    const std::string griprovider_type = config_.get< std::string >(id + ".gridprovider");
    info << "creating grid with '" << griprovider_type << "'..." << std::flush;
    grid_provider_ = GridProviders::create(griprovider_type, config_);
    const auto grid_view = grid_provider_->leaf_view();
    info << " done (took " << timer.elapsed()
         << "s, has " << grid_view.indexSet().size(0) << " element";
    if (grid_view.indexSet().size(0) > 1)
      info << "s";
    info << ")" << std::endl;

    boundary_info_ = Stuff::Grid::BoundaryInfoConfigs::AllDirichlet::default_config();

//    info << "setting up ";
//    timer.reset();
    const std::string problem_type = config_.get< std::string >(id + ".problem");
//    if (!debug_logging_)
//      info << "'" << problem_type << "'... " << std::flush;
    problem_ = ProblemsProvider::create(problem_type, config_);
//    if (debug_logging_)
//      info << *problem_ << std::endl;
//    else
//      info << "done (took " << timer.elapsed() << "s)" << std::endl;

    if (logger_config.get("visualize", true)) {
      info << "visualizing grid and problem... " << std::flush;
      timer.reset();
      grid_provider_->visualize(filename_ + ".ms_grid", false);
      grid_provider_->visualize(boundary_info_, filename_ + ".grid");
      problem_->visualize(grid_view, filename_ + ".problem");
      info << "done (took " << timer.elapsed() << "s)" << std::endl;
    } // if (visualize)
  } // DiscreteBlockProblem

  std::string filename() const
  {
    return filename_;
  }

  const Stuff::Common::Configuration& config() const
  {
    return config_;
  }

  bool debug_logging() const
  {
    return debug_logging_;
  }

  GridProviderType& grid_provider()
  {
    return *grid_provider_;
  }

  const GridProviderType& grid_provider() const
  {
    return *grid_provider_;
  }

  const Stuff::Common::Configuration& boundary_info() const
  {
    return boundary_info_;
  }

  const ProblemType& problem() const
  {
    return *problem_;
  }

private:
  static void write_keys_to_file(const std::string name, const std::vector< std::string > keys, std::ofstream& file)
  {
    std::string whitespace = Stuff::Common::whitespaceify(name + " = ");
    file << name << " = " << keys[0] << std::endl;
    for (size_t ii = 1; ii < keys.size(); ++ii)
      file << whitespace << keys[ii] << std::endl;
  } // ... write_keys_to_file(...)

  template< class ConfigProvider >
  static void write_config_to_file(std::ofstream& file)
  {
    for (const auto& type : ConfigProvider::available()) {
      auto config = ConfigProvider::default_config(type, type);
      if (!config.empty())
      file << config;
    }
  } // ... write_config_to_file(...)

  std::string filename_;
  Stuff::Common::Configuration config_;
  bool debug_logging_;
  std::unique_ptr< GridProviderType > grid_provider_;
  Stuff::Common::Configuration boundary_info_;
  std::unique_ptr< const ProblemType > problem_;
}; // class DiscreteBlockProblem


#endif // HAVE_DUNE_GRID_MULTISCALE


} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_DISCRETEPROBLEM_HH
