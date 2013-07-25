// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/misc/gridwidth.hh>

#include <dune/grid/part/leaf.hh>
#include <dune/grid/multiscale/provider.hh>

#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/string.hh>
#include <dune/stuff/grid/provider.hh>
#include <dune/stuff/grid/boundaryinfo.hh>

#include <dune/hdd/linearelliptic/problems.hh>
#include <dune/hdd/linearelliptic/discretizations.hh>


#if defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
template< class GridImp = Dune::GridSelector::GridType >
#else
template< class GridImp = Dune::SGrid< 2, 2 > >
#endif
class DiscreteProblem
{
public:
  typedef GridImp                                         GridType;
  typedef Dune::Stuff::GridProviderInterface< GridType >  GridProviderType;
  typedef Dune::Stuff::GridProviders< GridType >          GridProviders;
  typedef Dune::grid::Part::Leaf::Const< GridType >       GridPartType;
  typedef Dune::Stuff::GridboundaryInterface< typename GridPartType::GridViewType > GridboundariesType;
  typedef Dune::Stuff::Gridboundaries< typename GridPartType::GridViewType >        Gridboundaries;
  typedef typename GridPartType::ctype   DomainFieldType;
  static const int DUNE_UNUSED( dimDomain) = GridProviderType::dim;
  typedef double                RangeFieldType;
  static const int DUNE_UNUSED( dimRange) = 1;
  typedef Dune::HDD::LinearElliptic::ProblemInterface< DomainFieldType, dimDomain, RangeFieldType, dimRange >  ProblemType;
  typedef Dune::HDD::LinearElliptic::Problems< DomainFieldType, dimDomain, RangeFieldType, dimRange > Problems;
  typedef Dune::Stuff::Common::ExtendedParameterTree SettingsType;

  static void writeSettingsFile(const std::string filename, const std::string _id)
  {
    std::ofstream file;
    file.open(filename);
    file << "[" << _id << "]" << std::endl;
    writeKeysToFile("gridprovider", GridProviders::available(), file);
    writeKeysToFile("boundaryinfo", Gridboundaries::available(), file);
    writeKeysToFile("problem", Problems::available(), file);
    file << "[logging]" << std::endl;
    file << "info  = true" << std::endl;
    file << "debug = true" << std::endl;
    file << "file  = false" << std::endl;
    file << "[parameter]" << std::endl;
    file << "0.diffusion = [0.1; 0.1; 1.0; 1.0]" << std::endl;
    file << "1.diffusion = [1.0; 1.0; 0.1; 0.1]" << std::endl;
//    file << "[linearsolver]" << std::endl;
//    file << "type      = bicgstab.ilut" << std::endl;
//    file << "            bicgstab.diagonal" << std::endl;
//    file << "precision = 1e-12"  << std::endl;
//    file << "maxIter   = 5000"  << std::endl;
    writeSettingsToFile< GridProviders >(file);
    writeSettingsToFile< Gridboundaries >(file);
    writeSettingsToFile< Problems >(file);
    file.close();
  } // ... writewriteSettingsFileFile(...)

  DiscreteProblem(const std::string id, const std::vector< std::string >& arguments)
  {
    // mpi
    int argc = arguments.size();
    char** argv = Dune::Stuff::Common::String::vectorToMainArgs(arguments);
    Dune::Fem::MPIManager::initialize(argc, argv);

    // parameter
    settings_ = SettingsType(argc, argv, id + ".settings");
    settings_.assertSub(id);
    filename_ = settings_.get(id + ".filename", id);

    // logger
    const SettingsType& logDescription = settings_.sub("logging");
    int logFlags = Dune::Stuff::Common::LOG_CONSOLE;
    debugLogging_ = logDescription.get< bool >("debug", false);
    if (logDescription.get< bool >("info"))
      logFlags = logFlags | Dune::Stuff::Common::LOG_INFO;
    if (debugLogging_)
      logFlags = logFlags | Dune::Stuff::Common::LOG_DEBUG;
    if (logDescription.get< bool >("file", false))
      logFlags = logFlags | Dune::Stuff::Common::LOG_FILE;
    Dune::Stuff::Common::Logger().create(logFlags, id, "", "");
    Dune::Stuff::Common::LogStream& info  = Dune::Stuff::Common::Logger().info();

    Dune::Timer timer;

    const std::string griproviderType = settings_.get< std::string >(id + ".gridprovider");
    info << "setting up grid with '" << griproviderType << "': " << std::endl;
    const std::shared_ptr< const GridProviderType > gridProvider(GridProviders::create(griproviderType, settings_));
    grid_ = gridProvider->grid();
    gridPart_ = std::make_shared< const GridPartType >(*grid_);
    info << "  done (took " << timer.elapsed()
         << " sec, has " << grid_->size(0) << " element";
    if (grid_->size(0) > 1)
      info << "s";
    info << " and a width of "
         << Dune::Fem::GridWidth::calcGridWidth(*gridPart_) << ")" << std::endl;

    const std::string gridbboundaryType = settings_.get< std::string >(id + ".boundaryinfo");
    info << "setting up gridbboundary '" << gridbboundaryType << "'... " << std::flush;
    timer.reset();
    boundaryInfo_ = std::shared_ptr< const GridboundariesType >(Gridboundaries::create(gridbboundaryType,
                                                                                       settings_));
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    info << "visualizing grid... " << std::flush;
    timer.reset();
    gridProvider->visualize(gridbboundaryType, settings_, filename_ + ".grid");
    info << "done (took " << timer.elapsed() << " sek)" << std::endl;

    info << "setting up problem";
    const std::string problemType = settings_.get< std::string >(id + ".problem");
    if (!debugLogging_)
      info << "... ";
    else {
      info << ":" << std::endl;
      info << "  '" << problemType << "'... ";
    }
    info << std::flush;
    timer.reset();
    problem_ = std::shared_ptr< const ProblemType >(Problems::create(problemType, settings_));
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;
    info << "visualizing problem... " << std::flush;
    timer.reset();
    problem_->visualize(gridPart_->gridView(), filename_ + ".problem");
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // Problem

  std::string filename() const
  {
    return filename_;
  }

  const SettingsType& settings() const
  {
    return settings_;
  }

  bool debugLogging() const
  {
    return debugLogging_;
  }

  std::shared_ptr< const GridType > grid() const
  {
    return grid_;
  }

  std::shared_ptr< const GridPartType > gridPart() const
  {
    return gridPart_;
  }

  std::shared_ptr< const GridboundariesType > boundaryInfo() const
  {
    return boundaryInfo_;
  }

  std::shared_ptr< const ProblemType > problem() const
  {
    return problem_;
  }

private:
  static void writeKeysToFile(const std::string name,
             const std::vector< std::string > keys,
             std::ofstream& file)
  {
    std::string whitespace = Dune::Stuff::Common::whitespaceify(name + " = ");
    file << name << " = " << keys[0] << std::endl;
    for (size_t ii = 1; ii < keys.size(); ++ii)
      file << whitespace << keys[ii] << std::endl;
  }

  template< class DescriptionProvider >
  static void writeSettingsToFile(std::ofstream& file)
  {
    for (auto type : DescriptionProvider::available()) {
      const Dune::Stuff::Common::ExtendedParameterTree& settings = DescriptionProvider::defaultSettings(type, type);
      settings.reportNicely(file);
    }
  }

  std::string filename_;
  SettingsType settings_;
  bool debugLogging_;
  std::shared_ptr< const GridType > grid_;
  std::shared_ptr< const GridPartType > gridPart_;
  std::shared_ptr< const GridboundariesType > boundaryInfo_;
  std::shared_ptr< const ProblemType > problem_;
}; // class DiscreteProblem


//#if defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
//template< class GridImp = Dune::GridSelector::GridType >
//#else
//template< class GridImp = Dune::SGrid< 2, 2 > >
//#endif
//class MultiscaleProblem
//{
//public:
//  typedef GridImp                                               GridType;
//  typedef Dune::grid::Multiscale::ProviderInterface< GridType > MsGridProviderType;
//  typedef Dune::grid::Multiscale::Providers< GridType >         MsGridProviders;
//  typedef typename MsGridProviderType::MsGridType               MsGridType;
//  typedef typename MsGridType::GlobalGridViewType               GlobalGridViewType;
//  typedef Dune::Stuff::GridboundaryInterface< GlobalGridViewType > GridboundariesType;
//  typedef Dune::Stuff::Gridboundaries< GlobalGridViewType >        Gridboundaries;
//  typedef typename GridType::ctype   DomainFieldType;
//  static const int DUNE_UNUSED( dimDomain) = MsGridProviderType::dim;
//  typedef double                RangeFieldType;
//  static const int DUNE_UNUSED( dimRange) = 1;
//  typedef Dune::HDD::LinearElliptic::ModelInterface< DomainFieldType, dimDomain, RangeFieldType, dimRange > ModelType;
//  typedef Dune::HDD::LinearElliptic::Models< DomainFieldType, dimDomain, RangeFieldType, dimRange >         Models;
//  typedef Dune::Stuff::Common::ExtendedParameterTree SettingsType;

//  static void writeDescriptionFile(const std::string filename, const std::string _id = id())
//  {
//    std::ofstream file;
//    file.open(filename);
//    file << "[" << _id << "]" << std::endl;
//    writeKeysToFile("gridprovider", MsGridProviders::available(), file);
//    writeKeysToFile("boundaryinfo", Gridboundaries::available(), file);
//    writeKeysToFile("model", Models::available(), file);
//    file << "[logging]" << std::endl;
//    file << "info  = true" << std::endl;
//    file << "debug = true" << std::endl;
//    file << "file  = false" << std::endl;
//    file << "[parameter]" << std::endl;
//    file << "test.size = 2" << std::endl;
//    file << "test.0    = [0.1; 0.1]" << std::endl;
//    file << "test.1    = [1.0; 1.0]" << std::endl;
//    file << "[linearsolver]" << std::endl;
//    file << "type      = bicgstab.ilut" << std::endl;
//    file << "            bicgstab.diagonal" << std::endl;
//    file << "precision = 1e-12"  << std::endl;
//    file << "maxIter   = 5000"  << std::endl;
//    writeDescriptionsToFile< MsGridProviders >(file);
//    writeDescriptionsToFile< Gridboundaries >(file);
//    writeDescriptionsToFile< Models >(file);
//    file.close();
//  } // ... writeDescriptionFile(...)

//  MultiscaleProblem(int argc, char** argv)
//  {
//    // mpi
//    Dune::MPIManager::initialize(argc, argv);

//    // parameter
//    const std::string descriptionFilename = id() + ".description";
//    settings_ = SettingsType(argc, argv, descriptionFilename);
//    settings_.assertSub(id());
//    filename_ = settings_.get(id() + ".filename", id());

//    // logger
//    const SettingsType& logDescription = settings_.sub("logging");
//    int logFlags = Dune::Stuff::Common::LOG_CONSOLE;
//    debugLogging_ = logDescription.get< bool >("debug", false);
//    if (logDescription.get< bool >("info"))
//      logFlags = logFlags | Dune::Stuff::Common::LOG_INFO;
//    if (debugLogging_)
//      logFlags = logFlags | Dune::Stuff::Common::LOG_DEBUG;
//    if (logDescription.get< bool >("file", false))
//      logFlags = logFlags | Dune::Stuff::Common::LOG_FILE;
//    Dune::Stuff::Common::Logger().create(logFlags, id(), "", "");
//    Dune::Stuff::Common::LogStream& info  = Dune::Stuff::Common::Logger().info();

//    Dune::Timer timer;

//    const std::string griproviderType = settings_.get< std::string >(id() + ".gridprovider");
//    info << "setting up multiscale grid with '" << griproviderType << "': " << std::endl;
//    const std::shared_ptr< const MsGridProviderType > msGridProvider(MsGridProviders::create(griproviderType, settings_));
//    grid_ = msGridProvider->grid();
//    msGrid_ = msGridProvider->msGrid();
//    info << "  done (took " << timer.elapsed()
//         << " sec, has " << grid_->size(0) << " element";
//    if (grid_->size(0) > 1)
//      info << "s";
//    info << ", " << msGrid_->size() << " subdomain";
//    if (msGrid_->size() > 1)
//      info << "s";
//    info << " and a width of "
//         << Dune::GridWidth::calcGridWidth(*(msGrid_->globalGridPart())) << ")" << std::endl;

//    const std::string gridbboundaryType = settings_.get< std::string >(id() + ".boundaryinfo");
//    info << "setting up gridbboundary '" << gridbboundaryType << "'... " << std::flush;
//    timer.reset();
//    boundaryInfo_ = std::shared_ptr< const GridboundariesType >(Gridboundaries::create(gridbboundaryType,
//                                                                                       settings_));
//    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

//    info << "visualizing grid... " << std::flush;
//    timer.reset();
//    msGridProvider->visualize(filename_ + ".grid");
//    info << "done (took " << timer.elapsed() << " sek)" << std::endl;

//    const std::string modelType = settings_.get< std::string >(id() + ".model");
//    info << "setting up model '" << modelType << "'... " << std::flush;
//    timer.reset();
//    model_ = std::shared_ptr< const ModelType >(Models::create(modelType, settings_));
//    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

//    info << "visualizing model... " << std::flush;
//    timer.reset();
//    model_->visualize(*(msGrid_->globalGridView()), filename_ + ".model");
//    info << "done (took " << timer.elapsed() << " sec)" << std::endl;
//  } // Problem

//  std::string filename() const
//  {
//    return filename_;
//  }

//  const SettingsType& description() const
//  {
//    return settings_;
//  }

//  bool debugLogging() const
//  {
//    return debugLogging_;
//  }

//  std::shared_ptr< const MsGridType > msGrid() const
//  {
//    return msGrid_;
//  }

//  std::shared_ptr< const GridType > grid() const
//  {
//    return grid_;
//  }

//  std::shared_ptr< const GridboundariesType > boundaryInfo() const
//  {
//    return boundaryInfo_;
//  }

//  std::shared_ptr< const ModelType > model() const
//  {
//    return model_;
//  }

//private:
//  static void writeKeysToFile(const std::string name,
//             const std::vector< std::string > keys,
//             std::ofstream& file)
//  {
//    std::string whitespace = Dune::Stuff::Common::whitespaceify(name + " = ");
//    file << name << " = " << keys[0] << std::endl;
//    for (size_t ii = 1; ii < keys.size(); ++ii)
//      file << whitespace << keys[ii] << std::endl;
//  }

//  template< class DescriptionProvider >
//  static void writeDescriptionsToFile(std::ofstream& file)
//  {
//    for (auto type : DescriptionProvider::available()) {
//      const Dune::Stuff::Common::ExtendedParameterTree& description = DescriptionProvider::createSampleDescription(type, type);
//      description.reportNicely(file);
//    }
//  }

//  std::string filename_;
//  SettingsType settings_;
//  bool debugLogging_;
//  std::shared_ptr< const GridType > grid_;
//  std::shared_ptr< const MsGridType > msGrid_;
//  std::shared_ptr< const GridboundariesType > boundaryInfo_;
//  std::shared_ptr< const ModelType > model_;
//}; // class MultiscaleProblem
