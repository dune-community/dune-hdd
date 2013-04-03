#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include "config.h"
#endif

#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/misc/gridwidth.hh>

#include <dune/grid/part/leaf.hh>

#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/grid/provider.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/common/vector.hh>

#include <dune/detailed-solvers/linearelliptic/model.hh>


namespace Stuff = Dune::Stuff;
namespace Elliptic = Dune::DetailedSolvers::LinearElliptic;


static const int polOrder = 1;
static const int maxParamDim = 50;
typedef Stuff::GridProviderInterface<>            GridProviderType;
typedef Stuff::GridProviders<>                    GridProviders;
typedef GridProviderType::GridType                GridType;
typedef Dune::grid::Part::Leaf::Const< GridType > GridPartType;
typedef Stuff::GridboundaryInterface< typename GridPartType::GridViewType > GridboundariesType;
typedef Stuff::Gridboundaries< typename GridPartType::GridViewType >        Gridboundaries;
typedef GridPartType::ctype                                                               DomainFieldType;
static const int DUNE_UNUSED(                                                             dimDomain) = GridProviderType::dim;
typedef double                                                                            RangeFieldType;
static const int DUNE_UNUSED(                                                             dimRange) = 1;
typedef Elliptic::ModelInterface< DomainFieldType, dimDomain, RangeFieldType, dimRange >  ModelType;
typedef Elliptic::Models< DomainFieldType, dimDomain, RangeFieldType, dimRange >          Models;
typedef Stuff::Common::ExtendedParameterTree DescriptionType;


class Problem
{
public:
  static void writeDescriptionFile(const std::string filename, const std::string _id = id())
  {
    std::ofstream file;
    file.open(filename);
    file << "[" << _id << "]" << std::endl;
    writeKeysToFile("gridprovider", GridProviders::available(), file);
    writeKeysToFile("boundaryinfo", Gridboundaries::available(), file);
    writeKeysToFile("model", Models::available(), file);
    file << "[logging]" << std::endl;
    file << "info  = true" << std::endl;
    file << "debug = true" << std::endl;
    file << "file  = false" << std::endl;
    file << "[parameter]" << std::endl;
    file << "test.size = 2" << std::endl;
    file << "test.0    = [0.1; 0.1]" << std::endl;
    file << "test.1    = [1.0; 1.0]" << std::endl;
    file << "[linearsolver]" << std::endl;
    file << "type      = bicgstab.ilut" << std::endl;
    file << "            bicgstab.diagonal" << std::endl;
    file << "precision = 1e-12"  << std::endl;
    file << "maxIter   = 5000"  << std::endl;
    writeDescriptionsToFile< GridProviders >(file);
    writeDescriptionsToFile< Gridboundaries >(file);
    writeDescriptionsToFile< Models >(file);
    file.close();
  } // void writeParamFile(const std::string filename)

  Problem(int argc, char** argv)
  {
    // mpi
    Dune::MPIManager::initialize(argc, argv);

    // parameter
    const std::string descriptionFilename = id() + ".description";
    description_ = DescriptionType(argc, argv, descriptionFilename);
    description_.assertSub(id());
    filename_ = description_.get(id() + ".filename", id());

    // logger
    const DescriptionType& logDescription = description_.sub("logging");
    int logFlags = Stuff::Common::LOG_CONSOLE;
    debugLogging_ = logDescription.get< bool >("debug", false);
    if (logDescription.get< bool >("info"))
      logFlags = logFlags | Stuff::Common::LOG_INFO;
    if (debugLogging_)
      logFlags = logFlags | Stuff::Common::LOG_DEBUG;
    if (logDescription.get< bool >("file", false))
      logFlags = logFlags | Stuff::Common::LOG_FILE;
    Stuff::Common::Logger().create(logFlags, id(), "", "");
    Stuff::Common::LogStream& info  = Stuff::Common::Logger().info();

    Dune::Timer timer;

    const std::string griproviderType = description_.get< std::string >(id() + ".gridprovider");
    info << "setting up grid with '" << griproviderType << "': " << std::endl;
    const std::shared_ptr< const GridProviderType > gridProvider(GridProviders::create(griproviderType, description_));
    grid_ = gridProvider->grid();
    const std::shared_ptr< const GridPartType > gridPart(new GridPartType(*grid_));
    info << "  done (took " << timer.elapsed()
         << " sec, has " << grid_->size(0) << " element";
    if (grid_->size(0) > 1)
      info << "s";
    info << " and a width of "
         << Dune::GridWidth::calcGridWidth(*gridPart) << ")" << std::endl;

    const std::string gridbboundaryType = description_.get< std::string >(id() + ".boundaryinfo");
    info << "setting up gridbboundary '" << gridbboundaryType << "'... " << std::flush;
    timer.reset();
    boundaryInfo_ = std::shared_ptr< const GridboundariesType >(Gridboundaries::create(gridbboundaryType,
                                                                                       description_));
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    info << "visualizing grid... " << std::flush;
    timer.reset();
    gridProvider->visualize(gridbboundaryType, description_, filename_ + ".grid");
    info << " done (took " << timer.elapsed() << " sek)" << std::endl;

    info << "setting up model";
    const std::string modelType = description_.get< std::string >(id() + ".model");
    if (!debugLogging_)
      info << "... ";
    else {
      info << ":" << std::endl;
      info << "  '" << modelType << "'... ";
    }
    info << std::flush;
    timer.reset();
    model_ = std::shared_ptr< const ModelType >(Models::create(modelType, description_));
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;
    info << "visualizing model... " << std::flush;
    timer.reset();
    model_->visualize(gridPart->gridView(), filename_ + ".model");
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // Problem

  std::string filename() const
  {
    return filename_;
  }

  const DescriptionType& description() const
  {
    return description_;
  }

  bool debugLogging() const
  {
    return debugLogging_;
  }

  std::shared_ptr< const GridType > grid() const
  {
    return grid_;
  }

  std::shared_ptr< const GridboundariesType > boundaryInfo() const
  {
    return boundaryInfo_;
  }

  std::shared_ptr< const ModelType > model() const
  {
    return model_;
  }

private:
  static void writeKeysToFile(const std::string name,
             const std::vector< std::string > keys,
             std::ofstream& file)
  {
    std::string whitespace = Stuff::Common::whitespaceify(name + " = ");
    file << name << " = " << keys[0] << std::endl;
    for (size_t ii = 1; ii < keys.size(); ++ii)
      file << whitespace << keys[ii] << std::endl;
  }

  template< class DescriptionProvider >
  static void writeDescriptionsToFile(std::ofstream& file)
  {
    for (auto type : DescriptionProvider::available()) {
      const Stuff::Common::ExtendedParameterTree& description = DescriptionProvider::createSampleDescription(type, type);
      description.reportNicely(file);
    }
  }

  std::string filename_;
  DescriptionType description_;
  bool debugLogging_;
  std::shared_ptr< const GridType > grid_;
  std::shared_ptr< const GridboundariesType > boundaryInfo_;
  std::shared_ptr< const ModelType > model_;
}; // class Problem
