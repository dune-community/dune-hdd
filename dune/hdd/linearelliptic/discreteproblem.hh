// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_DISCRETEPROBLEM_HH
#define DUNE_HDD_LINEARELLIPTIC_DISCRETEPROBLEM_HH

#include <memory>
#include <limits>

#include <dune/common/mpihelper.hh>
#include <dune/common/timer.hh>

#if HAVE_ALUGRID
# include <dune/grid/alugrid.hh>
#endif

#if HAVE_DUNE_FEM
# include <dune/stuff/common/disable_warnings.hh>
#   include <dune/fem/misc/mpimanager.hh>
# include <dune/stuff/common/reenable_warnings.hh>
#endif

#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/string.hh>
#include <dune/stuff/grid/provider.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/common/configtree.hh>
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
  typedef LinearElliptic::ProblemProvider< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > ProblemProvider;

  static void write_config(const std::string filename, const std::string id)
  {
    std::ofstream file;
    file.open(filename);
    file << "[" << id << "]" << std::endl;
    write_keys_to_file("gridprovider", GridProviders::available(), file);
    write_keys_to_file("boundaryinfo", BoundaryInfoProvider::available(), file);
    write_keys_to_file("problem", ProblemProvider::available(), file);
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
    write_config_to_file< ProblemProvider >(file);
    file.close();
  } // ... write_config(...)

  DiscreteProblem(const std::string id, const std::vector< std::string >& arguments)
  {
    // mpi
    assert(arguments.size() < std::numeric_limits< int >::max());
    int argc = int(arguments.size());
    char** argv = Stuff::Common::String::vectorToMainArgs(arguments);
#if HAVE_DUNE_FEM
    Fem::MPIManager::initialize(argc, argv);
#else
    MPIHelper::instance(argc, argv);
#endif

    // configuration
    config_ = Stuff::Common::ConfigTree(argc, argv, id + ".cfg");
    if (!config_.has_sub(id))
      DUNE_THROW_COLORFULLY(Stuff::Exceptions::configuration_error,
                            "Missing sub '" << id << "' in the following ConfigTree:\n\n" << config_);
    filename_ = config_.get(id + ".filename", id);

    // logger
    const Stuff::Common::ConfigTree logger_config = config_.sub("logging");
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
         << "s, has " << grid_view->indexSet().size(0) << " element";
    if (grid_view->indexSet().size(0) > 1)
      info << "s";
    info << ")" << std::endl;

    const std::string boundary_info_type = config_.get< std::string >(id + ".boundaryinfo");
    if (config_.has_sub(boundary_info_type))
      boundary_info_ = config_.sub(boundary_info_type);
    else
      boundary_info_ = Stuff::Common::ConfigTree("type", boundary_info_type);

    info << "setting up ";
    timer.reset();
    const std::string problem_type = config_.get< std::string >(id + ".problem");
    if (!debug_logging_)
      info << "'" << problem_type << "'... " << std::flush;
    problem_ = ProblemProvider::create(problem_type, config_);
    if (debug_logging_)
      info << *problem_ << std::endl;
    else
      info << "done (took " << timer.elapsed() << "s)" << std::endl;

    if (logger_config.get("visualize", true)) {
      info << "visualizing grid and problem... " << std::flush;
      timer.reset();
      grid_provider_->visualize(boundary_info_, filename_ + ".grid");
      problem_->visualize(*grid_view, filename_ + ".problem");
      info << "done (took " << timer.elapsed() << "s)" << std::endl;
    } // if (visualize)
  } // DiscreteProblem

  std::string filename() const
  {
    return filename_;
  }

  const Stuff::Common::ConfigTree& config() const
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

  const Stuff::Common::ConfigTree& boundary_info() const
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
  Stuff::Common::ConfigTree config_;
  bool debug_logging_;
  std::unique_ptr< GridProviderType > grid_provider_;
  Stuff::Common::ConfigTree boundary_info_;
  std::unique_ptr< const ProblemType > problem_;
}; // class DiscreteProblem

#if HAVE_ALUGRID


extern template class DiscreteProblem< ALUConformGrid< 2, 2 > >;


#endif // HAVE_ALUGRID

} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_DISCRETEPROBLEM_HH
