// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_PROBLEMS_THERMALBLOCK_HH
#define DUNE_HDD_LINEARELLIPTIC_PROBLEMS_THERMALBLOCK_HH

#include <memory>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/static_assert.hh>
#include <dune/common/timer.hh>

#include <dune/grid/multiscale/provider/cube.hh>

#include <dune/stuff/common/string.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/expression.hh>
#include <dune/stuff/playground/functions/indicator.hh>
#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/grid/boundaryinfo.hh>

#include <dune/pymor/functions/default.hh>
#include <dune/pymor/functions/checkerboard.hh>

#include "default.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Problems {


template< class E, class D, int d, class R, int r = 1 >
class Thermalblock
  : public ProblemInterface< E, D, d, R, r >
{
  Thermalblock() { static_assert(AlwaysFalse< E >::value, "Not available for these dimensions!"); }
};


template< class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp >
class Thermalblock< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 >
  : public Default< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 >
{
  typedef Default< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 > BaseType;
  typedef Thermalblock< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 > ThisType;

public:
  typedef Pymor::Functions::Checkerboard< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 >
      CheckerboardFunctionType;
  using typename BaseType::DiffusionTensorType;
  using typename BaseType::FunctionType;

  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::BaseType::static_id() + ".thermalblock";
  }

  static Stuff::Common::Configuration default_config(const std::string sub_name = "")
  {
    typedef Stuff::Functions::Constant< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 >
        ConstantFunctionType;
    typedef Stuff::Functions::Constant< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, domainDim, domainDim >
        ConstantMatrixFunctionType;
    Stuff::Common::Configuration config;
    Stuff::Common::Configuration checkerboard_config = CheckerboardFunctionType::default_config();
    checkerboard_config["name"] = "diffusion_factor";
    checkerboard_config["type"] = CheckerboardFunctionType::static_id();
    checkerboard_config["num_elements"] = "[4 4 4]";
    checkerboard_config["parameter_name"] = "diffusion_factor";
    config.add(checkerboard_config, "diffusion_factor");
    Stuff::Common::Configuration diffusion_tensor_config = ConstantMatrixFunctionType::default_config();
    diffusion_tensor_config["name"] = "diffusion_tensor";
    diffusion_tensor_config["type"] = ConstantMatrixFunctionType::static_id();
    config.add(diffusion_tensor_config, "diffusion_tensor");
    Stuff::Common::Configuration constant_config = ConstantFunctionType::default_config();
    constant_config["type"] = ConstantFunctionType::static_id();
    constant_config["name"] = "force";
    constant_config["value"] = "1";
    config.add(constant_config, "force");
    constant_config["name"] = "dirichlet";
    constant_config["value"] = "0";
    config.add(constant_config, "dirichlet");
    constant_config["name"] = "neumann";
    config.add(constant_config, "neumann");
    if (sub_name.empty())
      return config;
    else {
      Stuff::Common::Configuration tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  static std::unique_ptr< ThisType > create(const Stuff::Common::Configuration config = default_config(),
                                            const std::string sub_name = static_id())
  {
    const Stuff::Common::Configuration cfg = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    std::shared_ptr< CheckerboardFunctionType >
        checkerboard_function(CheckerboardFunctionType::create(cfg.sub("diffusion_factor")));
    return Stuff::Common::make_unique< ThisType >(checkerboard_function,
                                                  BaseType::create_matrix_function("diffusion_tensor", cfg),
                                                  BaseType::create_vector_function("force", cfg),
                                                  BaseType::create_vector_function("dirichlet", cfg),
                                                  BaseType::create_vector_function("neumann", cfg));
  } // ... create(...)

  Thermalblock(const std::shared_ptr< const CheckerboardFunctionType >& checkerboard_function,
               const std::shared_ptr< const DiffusionTensorType >& diffusion_tensor,
               const std::shared_ptr< const FunctionType >& force,
               const std::shared_ptr< const FunctionType >& dirichlet,
               const std::shared_ptr< const FunctionType >& neumann)
    : BaseType(checkerboard_function, diffusion_tensor, force, dirichlet, neumann)
  {}

  virtual std::string type() const override
  {
    return BaseType::BaseType::static_id() + ".thermalblock";
  }
}; // class Thermalblock< ..., 1 >


template< class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim = 1 >
class LocalThermalblock
{
  static_assert(AlwaysFalse< EntityImp >::value, "Not available for dimRange > 1!");
};


template< class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp >
class LocalThermalblock< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 >
  : public Default< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 >
{
  typedef Default< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 > BaseType;
  typedef LocalThermalblock< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 > ThisType;

  typedef Pymor::Functions::AffinelyDecomposableDefault< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 >
      AffinelyDecomposableDefaultFunctionType;
  typedef Stuff::Functions::Constant< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, domainDim, domainDim >
      ConstantMatrixFunctionType;
public:
  using typename BaseType::FunctionType;

  static std::string static_id()
  {
    return BaseType::BaseType::static_id() + ".localthermalblock";
  }

  static Stuff::Common::Configuration default_config(const std::string sub_name = "")
  {
    typedef Stuff::Functions::Constant< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 >
        ConstantFunctionType;
    Stuff::Common::Configuration config;
    Stuff::Common::Configuration constant_config = ConstantFunctionType::default_config();
    constant_config["type"] = ConstantFunctionType::static_id();
    constant_config["name"] = "force";
    constant_config["value"] = "1";
    config.add(constant_config, "force");
    constant_config["name"] = "dirichlet";
    constant_config["value"] = "0";
    config.add(constant_config, "dirichlet");
    constant_config["name"] = "neumann";
    config.add(constant_config, "neumann");
    if (sub_name.empty())
      return config;
    else {
      Stuff::Common::Configuration tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  static std::unique_ptr< ThisType > create(const Stuff::Common::Configuration config = default_config(),
                                            const std::string sub_name = static_id())
  {
    const Stuff::Common::Configuration cfg = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    return Stuff::Common::make_unique< ThisType >(BaseType::create_vector_function("force", cfg),
                                                  BaseType::create_vector_function("dirichlet", cfg),
                                                  BaseType::create_vector_function("neumann", cfg));
  } // ... create(...)

  LocalThermalblock(const std::shared_ptr< const FunctionType >& force,
                    const std::shared_ptr< const FunctionType >& dirichlet,
                    const std::shared_ptr< const FunctionType >& neumann)
    : BaseType(create_diffusion_factor(),
               create_diffusion_tensor(),
               force,
               dirichlet,
               neumann)
  {}

  virtual std::string type() const override
  {
    return BaseType::BaseType::static_id() + ".localthermalblock";
  }

private:
  static std::shared_ptr< AffinelyDecomposableDefaultFunctionType > create_diffusion_factor()
  {
    typedef Stuff::Functions::Indicator< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 > IndicatorFunctionType;
    const Pymor::ParameterType mu("diffusion_factor", 3);

    auto ret = std::make_shared< AffinelyDecomposableDefaultFunctionType >("diffusion_factor");
    ret->register_component(new IndicatorFunctionType({{{{0.0, 0.0}, {0.5, 0.16}}, 1.0},
                                                       {{{0.0, 0.16}, {0.16, 0.33}}, 1.0},
                                                       {{{0.33, 0.16}, {0.5, 0.33}}, 1.0},
                                                       {{{0.0, 0.33}, {0.5, 1.0}}, 1.0}},
                                                      "left_block"),
                            new Pymor::ParameterFunctional(mu, "diffusion_factor[0]"));
    ret->register_component(new IndicatorFunctionType({{{{0.5, 0.0}, {1.0, 1.0}}, 1.0}},
                                                      "right_block"),
                            new Pymor::ParameterFunctional(mu, "diffusion_factor[1]"));
    ret->register_component(new IndicatorFunctionType({{{{0.16, 0.16}, {0.33, 0.33}}, 1.0}},
                                                      "small_block"),
                            new Pymor::ParameterFunctional(mu, "diffusion_factor[2]"));
    return ret;
  } // ... create_diffusion_factor()

  typedef Pymor::Functions::NonparametricDefault
      < EntityImp, DomainFieldImp, domainDim, RangeFieldImp, domainDim, domainDim > NonparametricFunctionType;

  static std::shared_ptr< NonparametricFunctionType > create_diffusion_tensor()
  {
    Stuff::Common::Configuration diffusion_tensor_config = ConstantMatrixFunctionType::default_config();
    diffusion_tensor_config["name"] = "diffusion_tensor";
    return std::make_shared< NonparametricFunctionType >(ConstantMatrixFunctionType::create(diffusion_tensor_config));
  } // ... create_diffusion_tensor()
}; // class LocalThermalblock< ..., 1 >


} // namespace Problems
namespace DiscreteProblems {


template< class GridImp >
class Thermalblock
{
public:
  typedef GridImp GridType;
  typedef Stuff::Grid::Providers::Cube< GridType > GridProviderType;
private:
  typedef Stuff::Grid::BoundaryInfos::AllDirichlet< typename GridType::LeafIntersection > BoundaryInfoType;
  typedef typename GridType::template Codim< 0 >::Entity EntityType;
  typedef typename GridType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridType::dimension;
public:
  typedef double RangeFieldType;
  static const unsigned int dimRange = 1;
  typedef Problems::Thermalblock< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > ProblemType;

  static void write_config(const std::string filename, const std::string id)
  {
    std::ofstream file;
    file.open(filename);
    file << "[" << id << "]" << std::endl;
    file << "filename = " << id << std::endl;
    file << "[logging]" << std::endl;
    file << "info  = true" << std::endl;
    file << "debug = true" << std::endl;
    file << "file  = false" << std::endl;
    file << "[parameter]" << std::endl;
    file << "0.diffusion_factor = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]" << std::endl;
    file << GridProviderType::default_config(GridProviderType::static_id());
    file << ProblemType::default_config(ProblemType::static_id());
    file << "[pymor]" << std::endl;
    file << "training_set = random" << std::endl;
    file << "num_training_samples = 100" << std::endl;
    file << "reductor = generic" << std::endl;
    file << "extension_algorithm = gram_schmidt" << std::endl;
    file << "extension_algorithm_product = h1_semi" << std::endl;
    file << "greedy_error_norm = h1_semi" << std::endl;
    file << "use_estimator = False" << std::endl;
    file << "max_rb_size = 100" << std::endl;
    file << "target_error = 0.01" << std::endl;
    file << "final_compression = False" << std::endl;
    file << "compression_product = None" << std::endl;
    file << "test_set = training" << std::endl;
    file << "num_test_samples = 100" << std::endl;
    file << "test_error_norm = h1_semi" << std::endl;
    file.close();
  } // ... write_config(...)

  Thermalblock(const std::string id, const std::vector< std::string >& arguments, const bool visualize = true)
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
    const Stuff::Common::Configuration& logger_config = config_.sub("logging");
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
    info << "creating grid with '" << GridProviderType::static_id() << "'... " << std::flush;
    grid_provider_ = GridProviderType::create(config_);
    const auto grid_view = grid_provider_->leaf_view();
    info << " done (took " << timer.elapsed()
         << "s, has " << grid_view->indexSet().size(0) << " element";
    if (grid_view->indexSet().size(0) > 1)
      info << "s";
    info << ")" << std::endl;

    boundary_info_ = Stuff::Common::Configuration("type", BoundaryInfoType::static_id());

    info << "setting up ";
    info << "'" << ProblemType::static_id() << "'... " << std::flush;
    timer.reset();
    problem_ = ProblemType::create(config_);
    info << "done (took " << timer.elapsed() << "s)" << std::endl;

    if (visualize) {
      info << "visualizing grid and problem... " << std::flush;
      timer.reset();
      grid_provider_->visualize(boundary_info_, filename_ + ".grid");
      problem_->visualize(*grid_view, filename_ + ".problem");
      info << "done (took " << timer.elapsed() << "s)" << std::endl;
    } // if (visualize)
  } // Thermalblock

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
  std::string filename_;
  Stuff::Common::Configuration config_;
  bool debug_logging_;
  std::unique_ptr< GridProviderType > grid_provider_;
  Stuff::Common::Configuration boundary_info_;
  std::unique_ptr< const ProblemType > problem_;
}; // class Thermalblock


template< class GridImp >
class ThermalblockMultiScale
{
public:
  typedef GridImp GridType;
  typedef grid::Multiscale::Providers::Cube< GridType > GridProviderType;
private:
  typedef typename GridType::template Codim< 0 >::Entity EntityType;
  typedef typename GridType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridType::dimension;
public:
  typedef double RangeFieldType;
  static const unsigned int dimRange = 1;
  typedef Problems::Thermalblock< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > ProblemType;

  static void write_config(const std::string filename, const std::string id)
  {
    std::ofstream file;
    file.open(filename);
    file << "[" << id << "]" << std::endl;
    file << "filename = " << id << std::endl;
    file << "[logging]" << std::endl;
    file << "info  = true" << std::endl;
    file << "debug = false" << std::endl;
    file << "file  = false" << std::endl;
    file << "visualize = false" << std::endl;
    file << "[parameter]" << std::endl;
    file << "0.diffusion_factor = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]" << std::endl;
    file << GridProviderType::default_config(GridProviderType::static_id());
    file << ProblemType::default_config(ProblemType::static_id());
    file << "[pymor]" << std::endl;
    file << "training_set = random" << std::endl;
    file << "num_training_samples = 100" << std::endl;
    file << "reductor = generic" << std::endl;
    file << "extension_algorithm = pod" << std::endl;
    file << "extension_algorithm_product = h1_semi" << std::endl;
    file << "greedy_error_norm = h1_semi" << std::endl;
    file << "use_estimator = False" << std::endl;
    file << "max_rb_size = 100" << std::endl;
    file << "target_error = 0.01" << std::endl;
    file << "final_compression = False" << std::endl;
    file << "compression_product = None" << std::endl;
    file << "test_set = training" << std::endl;
    file << "num_test_samples = 100" << std::endl;
    file << "test_error_norm = h1_semi" << std::endl;
    file.close();
  } // ... write_config(...)

  ThermalblockMultiScale(const std::string id, const std::vector< std::string >& arguments, const bool visualize = true)
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
    const Stuff::Common::Configuration& logger_config = config_.sub("logging");
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
    info << "creating grid with '" << GridProviderType::static_id() << "'... " << std::flush;
    grid_provider_ = GridProviderType::create(config_);
    const auto grid_view = grid_provider_->leaf_view();
    info << " done (took " << timer.elapsed()
         << "s, has " << grid_view.indexSet().size(0) << " element";
    if (grid_view.indexSet().size(0) > 1)
      info << "s";
    info << ")" << std::endl;

    boundary_info_ = Stuff::Grid::BoundaryInfoConfigs::NormalBased::default_config();
    boundary_info_["dirichlet.0"] = "[0.0 -1.0]";
    boundary_info_["dirichlet.1"] = "[0.0 1.0]";
    boundary_info_["neumann.0"] = "[1.0 0.0]";
    boundary_info_["neumann.1"] = "[-1.0 0.0]";

    info << "setting up ";
    info << "'" << ProblemType::static_id() << "'... " << std::flush;
    timer.reset();
    problem_ = ProblemType::create(config_);
    info << "done (took " << timer.elapsed() << "s)" << std::endl;

    if (visualize && logger_config.get("visualize", false)) {
      info << "visualizing grid and problem... " << std::flush;
      timer.reset();
      grid_provider_->visualize(boundary_info_, filename_ + ".grid");
      grid_provider_->visualize(filename_ + ".msgrid");
      problem_->visualize(grid_view, filename_ + ".problem");
      info << "done (took " << timer.elapsed() << "s)" << std::endl;
    } // if (visualize)
  } // Thermalblock

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
  std::string filename_;
  Stuff::Common::Configuration config_;
  bool debug_logging_;
  std::unique_ptr< GridProviderType > grid_provider_;
  Stuff::Common::Configuration boundary_info_;
  std::unique_ptr< const ProblemType > problem_;
}; // class ThermalblockMultiScale


template< class GridImp >
class LocalThermalblockMultiScale
{
public:
  typedef GridImp GridType;
  typedef grid::Multiscale::Providers::Cube< GridType > GridProviderType;
private:
  typedef typename GridType::template Codim< 0 >::Entity EntityType;
  typedef typename GridType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridType::dimension;
public:
  typedef double RangeFieldType;
  static const unsigned int dimRange = 1;
  typedef Problems::LocalThermalblock< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > ProblemType;

  static void write_config(const std::string filename, const std::string id)
  {
    std::ofstream file;
    file.open(filename);
    file << "[" << id << "]" << std::endl;
    file << "filename = " << id << std::endl;
    file << "[logging]" << std::endl;
    file << "info  = true" << std::endl;
    file << "debug = false" << std::endl;
    file << "file  = false" << std::endl;
    file << "visualize = false" << std::endl;
    file << "[parameter]" << std::endl;
    file << "0.diffusion_factor = [1 2 3]" << std::endl;
    file << GridProviderType::default_config(GridProviderType::static_id());
    file << ProblemType::default_config(ProblemType::static_id());
    file << "[pymor]" << std::endl;
    file << "training_set = random" << std::endl;
    file << "num_training_samples = 100" << std::endl;
    file << "reductor = generic" << std::endl;
    file << "extension_algorithm = pod" << std::endl;
    file << "extension_algorithm_product = h1_semi" << std::endl;
    file << "greedy_error_norm = h1_semi" << std::endl;
    file << "use_estimator = False" << std::endl;
    file << "max_rb_size = 100" << std::endl;
    file << "target_error = 5e-3" << std::endl;
    file << "final_compression = False" << std::endl;
    file << "compression_product = None" << std::endl;
    file << "test_set = training" << std::endl;
    file << "num_test_samples = 100" << std::endl;
    file << "test_error_norm = h1_semi" << std::endl;
    file.close();
  } // ... write_config(...)

  LocalThermalblockMultiScale(const std::string id,
                              const std::vector< std::string >& arguments,
                              const bool visualize = true)
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
    const Stuff::Common::Configuration& logger_config = config_.sub("logging");
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
    info << "creating grid with '" << GridProviderType::static_id() << "'... " << std::flush;
    grid_provider_ = GridProviderType::create(config_);
    const auto grid_view = grid_provider_->leaf_view();
    info << " done (took " << timer.elapsed()
         << "s, has " << grid_view->indexSet().size(0) << " element";
    if (grid_view->indexSet().size(0) > 1)
      info << "s";
    info << ")" << std::endl;

    boundary_info_ = Stuff::Grid::BoundaryInfoConfigs::NormalBased::default_config();
    boundary_info_["dirichlet.0"] = "[0.0 -1.0]";
    boundary_info_["dirichlet.1"] = "[0.0 1.0]";
    boundary_info_["neumann.0"] = "[1.0 0.0]";
    boundary_info_["neumann.1"] = "[-1.0 0.0]";

    info << "setting up ";
    info << "'" << ProblemType::static_id() << "'... " << std::flush;
    timer.reset();
    problem_ = ProblemType::create(config_);
    info << "done (took " << timer.elapsed() << "s)" << std::endl;

    if (visualize && logger_config.get("visualize", false)) {
      info << "visualizing grid and problem... " << std::flush;
      timer.reset();
      grid_provider_->visualize(boundary_info_, filename_ + ".grid");
      grid_provider_->visualize(filename_ + ".msgrid");
      problem_->visualize(*grid_view, filename_ + ".problem");
      info << "done (took " << timer.elapsed() << "s)" << std::endl;
    } // if (visualize)
  } // Thermalblock

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
  std::string filename_;
  Stuff::Common::Configuration config_;
  bool debug_logging_;
  std::unique_ptr< GridProviderType > grid_provider_;
  Stuff::Common::Configuration boundary_info_;
  std::unique_ptr< const ProblemType > problem_;
}; // class LocalThermalblockMultiScale


} // namespace DiscreteProblems
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_PROBLEMS_THERMALBLOCK_HH
