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
#include <dune/detailed-solvers/linearelliptic/solver/cg/detailed-discretizations.hh>

namespace Stuff = Dune::Stuff;
namespace Elliptic = Dune::DetailedSolvers::LinearElliptic;

// fix some types here, so we can use them in writeDescriptionFile()
static const int polOrder = 1;
static const int maxParamDim = 50;
typedef Stuff::GridProviderInterface<>            GridProviderType;
typedef Stuff::GridProviders<>                    GridProviders;
typedef GridProviderType::GridType                GridType;
typedef Dune::grid::Part::Leaf::Const< GridType > GridPartType;
typedef Stuff::GridboundaryInterface< typename GridPartType::GridViewType > GridboundariesType;
typedef Stuff::Gridboundaries< typename GridPartType::GridViewType >        Gridboundaries;
typedef GridPartType::ctype                                                               DomainFieldType;
const int DUNE_UNUSED(                                                                    dimDomain) = GridProviderType::dim;
typedef double                                                                            RangeFieldType;
const int DUNE_UNUSED(                                                                    dimRange) = 1;
typedef Elliptic::ModelInterface< DomainFieldType, dimDomain, RangeFieldType, dimRange >  ModelType;
typedef Elliptic::Models< DomainFieldType, dimDomain, RangeFieldType, dimRange >          Models;


std::string id(){
  return "examples.linearelliptic.cg.dd";
}


void writeKeysToFile(const std::string name,
           const std::vector< std::string > keys,
           std::ofstream& file)
{
  std::string whitespace = Stuff::Common::whitespaceify(name + " = ");
  file << name << " = " << keys[0] << std::endl;
  for (size_t ii = 1; ii < keys.size(); ++ii)
    file << whitespace << keys[ii] << std::endl;
}


template< class DescriptionProvider >
void writeDescriptionsToFile(std::ofstream& file)
{
  for (auto type : DescriptionProvider::available()) {
    const Stuff::Common::ExtendedParameterTree& description = DescriptionProvider::createSampleDescription(type, type);
    description.reportNicely(file);
  }
}


void writeDescriptionFile(const std::string filename, const std::string _id = id())
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


int run(int argc, char** argv)
{
  try {
    // mpi
    Dune::MPIManager::initialize(argc, argv);

    // parameter
    const std::string descriptionFilename = id() + ".description";
    typedef Stuff::Common::ExtendedParameterTree DescriptionType;
    const DescriptionType description(argc, argv, descriptionFilename);
    description.assertSub(id());
    const std::string filename = description.get(id() + ".filename", id());

    // logger
    const DescriptionType& logDescription = description.sub("logging");
    int logFlags = Stuff::Common::LOG_CONSOLE;
    const bool debugLogging = logDescription.get< bool >("debug", false);
    if (logDescription.get< bool >("info"))
      logFlags = logFlags | Stuff::Common::LOG_INFO;
    if (debugLogging)
      logFlags = logFlags | Stuff::Common::LOG_DEBUG;
    if (logDescription.get< bool >("file", false))
      logFlags = logFlags | Stuff::Common::LOG_FILE;
    Stuff::Common::Logger().create(logFlags, id(), "", "");
    Stuff::Common::LogStream& info  = Stuff::Common::Logger().info();
    Stuff::Common::LogStream& debug = Stuff::Common::Logger().debug();

    Dune::Timer timer;

    const std::string griproviderType = description.get< std::string >(id() + ".gridprovider");
    info << "setting up grid with '" << griproviderType << "': " << std::endl;
    const std::shared_ptr< const GridProviderType > gridProvider(GridProviders::create(griproviderType, description));
    const std::shared_ptr< const GridType > grid = gridProvider->grid();
    const std::shared_ptr< const GridPartType > gridPart(new GridPartType(*grid));
    info << "  done (took " << timer.elapsed()
         << " sec, has " << grid->size(0) << " element";
    if (grid->size(0) > 1)
      info << "s";
    info << " and a width of "
         << Dune::GridWidth::calcGridWidth(*gridPart) << ")" << std::endl;

    const std::string gridbboundaryType = description.get< std::string >(id() + ".boundaryinfo");
    info << "setting up gridbboundary '" << gridbboundaryType << "'... " << std::flush;
    timer.reset();
    const std::shared_ptr< const GridboundariesType > boundaryInfo(Gridboundaries::create(gridbboundaryType,
                                                                                          description));
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    info << "visualizing grid... " << std::flush;
    timer.reset();
    gridProvider->visualize(gridbboundaryType, description, filename + ".grid");
    info << " done (took " << timer.elapsed() << " sek)" << std::endl;

    info << "setting up model";
    const std::string modelType = description.get< std::string >(id() + ".model");
    if (!debugLogging)
      info << "... ";
    else {
      info << ":" << std::endl;
      info << "  '" << modelType << "'... ";
    }
    info << std::flush;
    timer.reset();
    const std::shared_ptr< const ModelType > model(Models::create(modelType, description));
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;
    if (model->parametric() && !model->affineparametric())
      DUNE_THROW(Dune::NotImplemented,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " only implemented for nonparametric or affineparametric models!");
    info << "visualizing model... " << std::flush;
    timer.reset();
    model->visualize(gridPart->gridView(), filename + ".model");
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    info << "initializing solver";
    if (!debugLogging)
      info << "... " << std::flush;
    else
      info << ":" << std::endl;
    timer.reset();
    typedef Elliptic::SolverContinuousGalerkinDD< GridPartType, RangeFieldType, dimRange, 1 > SolverType;
    SolverType solver(gridPart, boundaryInfo, model);
    solver.init(debug, "  ");
    if (!debugLogging)
      info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    const DescriptionType& linearSolverDescription = description.sub("linearsolver");
    typedef typename SolverType::VectorType VectorType;
    const std::string linearSolverType =  linearSolverDescription.get< std::string >("type",      "bicgstab.ilut");
    const double linearSolverPrecision = linearSolverDescription.get< double >(      "precision", 1e-12);
    const size_t linearSolverMaxIter = linearSolverDescription.get< size_t >(        "maxIter",   5000);
    if (!model->parametric()) {
      info << "solving";
      if (!debugLogging)
        info << "... " << std::flush;
      else
        info << ":" << std::endl;
      timer.reset();
      std::shared_ptr< VectorType > solutionVector = solver.createVector();
      solver.solve(solutionVector,
                   linearSolverType,
                   linearSolverPrecision,
                   linearSolverMaxIter,
                   debug,
                   "  ");
      if (!debugLogging)
        info << "done (took " << timer.elapsed() << " sec)" << std::endl;

      info << "writing solution to disc";
      if (!debugLogging)
        info << "... " << std::flush;
      else
        info << ":" << std::endl;
      timer.reset();
      solver.visualize(solutionVector,
                       filename + ".solution",
                       id() + ".solution",
                       debug,
                       "  ");
      if (!debugLogging)
        info << "done (took " << timer.elapsed() << " sec)" << std::endl;
    } else { // if (!model->parametric())
      typedef typename ModelType::ParamFieldType  ParamFieldType;
      typedef typename ModelType::ParamType       ParamType;
      const size_t paramSize = model->paramSize();
      const DescriptionType& parameterDescription = description.sub("parameter");
      const size_t numTestParams = parameterDescription.get< size_t >("test.size");
      // loop over all test parameters
      for (size_t ii = 0; ii < numTestParams; ++ii) {
        const std::string iiString = Dune::Stuff::Common::toString(ii);
        const ParamType testParameter
            = parameterDescription.getDynVector< ParamFieldType >("test." + iiString, paramSize);
        // after this, testParameter is at least as long as paramSize, but it might be too long
        const ParamType mu = Dune::Stuff::Common::resize(testParameter, paramSize);
        info << "solving for parameter [" << mu << "]";
        if (!debugLogging)
          info << "... " << std::flush;
        else
          info << ":" << std::endl;
        timer.reset();
        std::shared_ptr< VectorType > solutionVector = solver.createVector();
        solver.solve(solutionVector,
                     mu,
                     linearSolverType,
                     linearSolverPrecision,
                     linearSolverMaxIter,
                     debug,
                     "  ");
        if (!debugLogging)
          info << "done (took " << timer.elapsed() << " sec)" << std::endl;

        info << "writing solution for parameter [" << mu << "] to disc";
        if (!debugLogging)
          info << "... " << std::flush;
        else
          info << ":" << std::endl;
        timer.reset();
        std::stringstream name;
        name << id() << ".solution." << iiString << " (parameter [" << mu << "])";
        solver.visualize(solutionVector,
                         filename + ".solution." + iiString,
                         name.str(),
                         debug,
                         "  ");
        if (!debugLogging)
          info << "done (took " << timer.elapsed() << " sec)" << std::endl;
      } // loop over all test parameters
    } // if (!model->parametric())

  } catch(Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e.what() << std::endl;
  } catch(std::exception& e) {
    std::cerr << e.what() << std::endl;
  } catch( ... ) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  } // try

  // if we came that far we can as well be happy about it
  return 0;
} // run
