#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#include <vector>

#include <boost/filesystem.hpp>

#include <Eigen/Dense>

#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dune/grid/multiscale/provider/cube.hh>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/misc/gridwidth.hh>

#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/la/container/eigen.hh>
#include <dune/stuff/la/algorithm/normalize.hh>
#include <dune/stuff/la/container/separable.hh>
#include <dune/stuff/la/solver.hh>
#include <dune/stuff/discretefunction/projection/lagrange.hh>
#include <dune/stuff/discretefunction/projection/heterogenous.hh>
#include <dune/stuff/la/algorithm/gramschmidt.hh>

#include <dune/detailed/solvers/stationary/linear/elliptic/model.hh>
#include <dune/detailed/solvers/stationary/linear/elliptic/ms/semicg/detailed-discretizations.hh>

#ifdef POLORDER
  const int polOrder = POLORDER;
#else
  const int polOrder = 1;
#endif

namespace Stationary = Dune::Detailed::Solvers::Stationary;

// some types
typedef Dune::Stuff::Common::ExtendedParameterTree                                      DescriptionType;
typedef Dune::grid::Multiscale::Provider::Cube<>                                        GridProviderType;
typedef GridProviderType::MsGridType                                                    MsGridType;
typedef typename MsGridType::GlobalGridViewType                                         GlobalGridViewType;
typedef typename MsGridType::LocalGridPartType                                          LocalGridPartType;
typedef Dune::Stuff::Grid::BoundaryInfo::Interface< GlobalGridViewType >                BoundaryInfoType;
typedef Dune::Stuff::Grid::BoundaryInfo::Interface< typename MsGridType::LocalGridViewType > LocalBoundaryInfoType;
const int DUNE_UNUSED(                                                                  dimDomain) = GridProviderType::dim;
const int DUNE_UNUSED(                                                                  dimRange) = 1;
typedef GridProviderType::CoordinateType::value_type                                    DomainFieldType;
typedef double                                                                          RangeFieldType;
typedef double                                                                          ParamFieldType;
typedef Stationary::Linear::Elliptic::Model::Interface< DomainFieldType, dimDomain,
                                                        RangeFieldType, dimRange >      ModelType;
typedef Stationary::Linear::Elliptic::MS::SemiCG::DetailedDiscretizations<  ModelType,
                                                                            MsGridType,
                                                                            BoundaryInfoType,
                                                                            polOrder >  SolverType;
typedef SolverType::LocalSolverType                                                     LocalSolverType;
typedef SolverType::VectorType                                                          VectorType;
typedef typename ModelType::ParamType                                                   ParamType;
typedef SolverType::MatrixType                                                          ScalarProductType;
typedef Dune::Stuff::LA::Container::EigenDenseMatrix< RangeFieldType >                  ReducedBasisType;
typedef typename SolverType::AnsatzMapperType                                           MapperType;
typedef SolverType::SeparableMatrixType                                                 DetailedOperatorType;
typedef VectorType                                                                      DetailedFunctionalType;
typedef Dune::Stuff::LA::Container::Separable< ReducedBasisType >                       ReducedOperatorType;
typedef VectorType                                                                      ReducedFunctionalType;
typedef Dune::Stuff::LA::Container::EigenDenseVector< size_t >                          MarkerType;
typedef typename SolverType::DiscreteFunctionType                                       GlobalDiscreteFunctionType;
typedef typename LocalSolverType::DiscreteAnsatzFunctionType                            LocalDiscreteFunctionType;
typedef Dune::Stuff::LA::Solver::Interface< ScalarProductType, VectorType >             LinearSolverType;


std::string id()
{
  return "stationary.linear.elliptic.ms.semidg.detailed-discretizations";
}


void writeDescriptionFile(std::string filename)
{
  // only write param file if there is none
  std::ofstream file;
  file.open(filename);
  file << "[" << id() << "]" << std::endl;
  file << "model = model.stationary.linear.elliptic.parametric.separable.thermalblock" << std::endl;
  file << "visualize = false" << std::endl;
  file << "[grid.multiscale.provider.cube]" << std::endl;
  file << "lowerLeft = [0.0; 0.0]" << std::endl;
  file << "upperRight = [1.0; 1.0]" << std::endl;
  file << "numElements = [50; 50]" << std::endl;
  file << "boundaryId = 7 # a cube from the factory gets the boundary ids 1 to 4 ind 2d and 1 to 6 in 3d (hopefully)" << std::endl;
  file << "partitions = [2; 2]" << std::endl;
  file << "oversamplingLayers = 10" << std::endl;
  file << "[detailed.solvers.stationary.linear.elliptic.ms.semicg.detailed-discretizations]" << std::endl;
  file << "penaltyFactor = 100.0" << std::endl;
  file << "linearsolver.type = bicgstab.ilut" << std::endl;
  file << "linearsolver.maxIter   = 5000"  << std::endl;
  file << "linearsolver.precision = 1e-12"  << std::endl;
  file << "[model.stationary.linear.elliptic.parametric.separable.thermalblock]" << std::endl;
  file << "diffusion.lowerLeft   = [0.0; 0.0] # this should be a bounding box of the above selected grid!" << std::endl;
  file << "diffusion.upperRight  = [1.0; 1.0] # this should be a bounding box of the above selected grid!" << std::endl;
  file << "diffusion.numElements = [6; 6]"  << std::endl;
  file << "diffusion.paramMin    = [0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1]" << std::endl;
  file << "diffusion.paramMax    = [1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0]" << std::endl;
  file << "diffusion.order       = 0" << std::endl;
  file << "diffusion.name        = checkerboard diffusion" << std::endl;
  file << "force.variable   = x" << std::endl;
  file << "force.expression = [1.0; 1.0; 1.0]"  << std::endl;
  file << "force.order      = 0" << std::endl;
  file << "force.name       = constant force" << std::endl;
  file << "dirichlet.variable   = x" << std::endl;
  file << "dirichlet.expression = [0.0; 0.0; 0.0]"  << std::endl;
  file << "dirichlet.order      = 0"  << std::endl;
  file << "dirichlet.name       = trivial dirichlet"  << std::endl;
  file << "neumann.variable   = x" << std::endl;
  file << "neumann.expression = [0.0; 0.0; 0.0]"  << std::endl;
  file << "neumann.order      = 0"  << std::endl;
  file << "neumann.name       = trivial neumann"  << std::endl;
  file << "[logging]" << std::endl;
  file << "info  = true" << std::endl;
  file << "debug = false" << std::endl;
  file << "file  = false" << std::endl;
  file << "[parameter]" << std::endl;
  file << "fixed = [1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0]" << std::endl;
  file << "train.0 = [1.0; 1.0; 1.0; 0.1; 0.1; 0.1; 1.0; 1.0; 1.0; 0.1; 0.1; 0.1; 1.0; 1.0; 1.0; 0.1; 0.1; 0.1; 1.0; 1.0; 1.0; 0.1; 0.1; 0.1; 1.0; 1.0; 1.0; 0.1; 0.1; 0.1; 1.0; 1.0; 1.0; 0.1; 0.1; 0.1]" << std::endl;
  file << "test.0  = [1.0; 1.0; 1.0; 0.1; 0.1; 0.1; 1.0; 1.0; 1.0; 0.1; 0.1; 0.1; 1.0; 1.0; 1.0; 0.1; 0.1; 0.1; 1.0; 1.0; 1.0; 0.1; 0.1; 0.1; 1.0; 1.0; 1.0; 0.1; 0.1; 0.1; 1.0; 1.0; 1.0; 0.1; 0.1; 0.1]" << std::endl;
  file << "test.1  = [1.0; 1.0; 1.0; 0.1; 0.1; 0.1; 1.0; 0.01; 1.0; 0.1; 0.1; 0.1; 1.0; 1.0; 1.0; 0.1; 0.1; 0.1; 1.0; 1.0; 1.0; 0.1; 0.1; 0.1; 1.0; 0.01; 1.0; 0.1; 0.1; 0.1; 1.0; 1.0; 1.0; 0.1; 0.1; 0.1]" << std::endl;
  file << "[reducedBasis]" << std::endl;
  file << "gramSchmidt = true" << std::endl;
  file << "[enrichment]" << std::endl;
  file << "boundaryInfo = stuff.grid.boundaryinfo.alldirichlet" << std::endl;
  file << "maxLocalError = 5e-4" << std::endl;
  file << "maxIterations = 50" << std::endl;
  file << "gramSchmidt = true" << std::endl;
  file << "gramSchmidtTolerance = 1e-10" << std::endl;
  file << "[markSubdomains]" << std::endl;
  file << "strategy = factorTimesMean" << std::endl;
  file << "markingFactor = 1.0" << std::endl;
  file << "age = 2" << std::endl;
  file.close();
} // void ensureParamFile()


DescriptionType initDescription(int argc, char** argv)
{
  // mpi
  Dune::MPIManager::initialize(argc, argv);

  // description
  const std::string descriptionFilename = id() + ".description";
  if (!boost::filesystem::exists(descriptionFilename))
    writeDescriptionFile(descriptionFilename);
  DescriptionType description(argc, argv, descriptionFilename);
  description["filename"] = id();

  // logger
  const DescriptionType& logDescription = description.sub("logging");
  int logFlags = Dune::Stuff::Common::LOG_CONSOLE;
  const bool debugLogging = logDescription.get< bool >("debug", false);
  if (logDescription.get< bool >("info"))
    logFlags = logFlags | Dune::Stuff::Common::LOG_INFO;
  if (debugLogging)
    logFlags = logFlags | Dune::Stuff::Common::LOG_DEBUG;
  if (logDescription.get< bool >("file", false))
    logFlags = logFlags | Dune::Stuff::Common::LOG_FILE;
  Dune::Stuff::Common::Logger().create(logFlags, id(), "", "");

  return description;
} // ... initDescription(...)


Dune::shared_ptr< const MsGridType > initGrid(const DescriptionType& description)
{
  // logger
  const bool visualize = description.get< bool >(id() + ".visualize", true);
  Dune::Stuff::Common::LogStream& info    = Dune::Stuff::Common::Logger().info();
  Dune::Stuff::Common::LogStream& debug   = Dune::Stuff::Common::Logger().debug();
  // timer
  Dune::Timer timer;

  info << "setting up grid: " << std::endl;
  debug.suspend();
  const GridProviderType gridProvider = GridProviderType::createFromDescription(description);
  Dune::shared_ptr< const MsGridType > msGrid = gridProvider.msGrid();
  info << "  took " << timer.elapsed()
       << " sec (has " << gridProvider.grid()->size(0) << " elements, "
       << msGrid->size() << " subdomain";
  if (msGrid->size() > 1)
    info << "s";
  info << " and a width of "
       << std::scientific << Dune::GridWidth::calcGridWidth(*(msGrid->globalGridPart()))
       << std::fixed << ")" << std::endl;
  debug.resume();
  if (visualize) {
    info << "visualizing grid... " << std::flush;
    timer.reset();
    debug.suspend();
    gridProvider.visualize(/*filename +*/ "grid");
    info << "done (took " << timer.elapsed() << " sek)" << std::endl;
    debug.resume();
  } // if (visualize)

  return msGrid;
} // ... initGrid(...)


Dune::shared_ptr< const BoundaryInfoType > initBoundaryInfo(const DescriptionType& description)
{
  // logger
  const bool debugLogging = description.get< bool >("logging.debug", false);
  Dune::Stuff::Common::LogStream& info    = Dune::Stuff::Common::Logger().info();

  // timer
  Dune::Timer timer;

  info << "setting up boundaryinfo";
  const std::string boundaryInfoType = "stuff.grid.boundaryinfo.alldirichlet";
  if (!debugLogging)
    info << "... ";
  else {
    info << ":" << std::endl;
    info << "  '" << boundaryInfoType << "'... ";
  }
  info << std::flush;
  timer.reset();
  Dune::shared_ptr< const BoundaryInfoType >
      boundaryInfo(Dune::Stuff::Grid::BoundaryInfo::create< GlobalGridViewType >(boundaryInfoType, description));
  info << "done (took " << timer.elapsed() << " sec)" << std::endl;

  return boundaryInfo;
} // ... initBoundaryInfo(...)


Dune::shared_ptr< const ModelType > initModel(const DescriptionType& description,
                                              const MsGridType& msGrid)
{
  // logger
  const bool visualize = description.get< bool >(id() + ".visualize", true);
  const bool debugLogging = description.get< bool >("logging.debug", false);
  Dune::Stuff::Common::LogStream& info    = Dune::Stuff::Common::Logger().info();

  // timer
  Dune::Timer timer;

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
  Dune::shared_ptr< const ModelType >
      model(Stationary::Linear::Elliptic::Model::create<  DomainFieldType, dimDomain,
                                                          RangeFieldType, dimRange >(modelType, description));
  info << "done (took " << timer.elapsed() << " sec)" << std::endl;
  if (visualize) {
    info << "visualizing model... " << std::flush;
    timer.reset();
    model->visualize(msGrid.globalGridPart()->gridView(), /*filename +*/ "model");
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;
  }

  return model;
} // ... initModel(...)


std::pair< Dune::shared_ptr< DetailedOperatorType >,
           Dune::shared_ptr< DetailedFunctionalType > >
initDetailedOperatorAndFunctional(const DescriptionType& description,
                                  const Dune::shared_ptr< const ModelType >& model,
                                  const Dune::shared_ptr< const MsGridType >& msGrid,
                                  const Dune::shared_ptr< const BoundaryInfoType >& boundaryInfo)
{
  // logger
//  const bool debugLogging = description.get< bool >("logging.debug", false);
  Dune::Stuff::Common::LogStream& info    = Dune::Stuff::Common::Logger().info();
  Dune::Stuff::Common::LogStream& devnull = Dune::Stuff::Common::Logger().devnull();

  // timer
  Dune::Timer timer;

  info << "assembling detailed operator and functional... "
       << std::flush;
  timer.reset();
  const DescriptionType& discretizationDescription = description.sub(SolverType::id());
  const RangeFieldType penaltyFactor = discretizationDescription.get< RangeFieldType >("penaltyFactor");
  // therefore, prepare the storage for the operator components
  typedef typename DetailedOperatorType::ComponentType DetailedOperatorComponentType;
  std::vector< Dune::shared_ptr< DetailedOperatorComponentType > > detailedOperatorComponents;
  // then assemble the first operator component and the functional
  // * therefore, create the first model
  typedef Stationary::Linear::Elliptic::Model::Default< DomainFieldType, dimDomain,
                                                          RangeFieldType, dimRange > DefaultModelType;
  Dune::shared_ptr< DefaultModelType >
      fixedModel = Dune::make_shared< DefaultModelType >(model->diffusion()->components()[0],
                                                         model->force(),
                                                         model->dirichlet(),
                                                         model->neumann());
  // * create the first solver
  Dune::shared_ptr< SolverType >
      detailedSolver = Dune::make_shared< SolverType >(fixedModel, msGrid, boundaryInfo, penaltyFactor);
  detailedSolver->init("", devnull);
  // * and set the functional
  Dune::shared_ptr< DetailedFunctionalType > detailedFunctional = detailedSolver->vector("force");
  // * and the first operator component
  assert(detailedSolver->matrix("diffusion")->numComponents() == 1);
  detailedOperatorComponents.push_back(detailedSolver->matrix("diffusion")->components()[0]);
  // then we loop over the rest of the components
  for (size_t qq = 1; qq < model->diffusion()->numComponents(); ++qq) {
    // so we create a fixed model
    fixedModel = Dune::make_shared< DefaultModelType >(model->diffusion()->components()[qq],
                                                       model->force(),
                                                       model->dirichlet(),
                                                       model->neumann());
    // and a fixed solver
    detailedSolver = Dune::make_shared< SolverType >(fixedModel, msGrid, boundaryInfo, penaltyFactor);
    detailedSolver->init("", devnull);
    // and extract the operator component
    assert(detailedSolver->matrix("diffusion")->numComponents() == 1);
    detailedOperatorComponents.push_back(detailedSolver->matrix("diffusion")->components()[0]);
  } // then we loop over the rest of the components

  // now we can create the detailed operator
  Dune::shared_ptr< DetailedOperatorType >
      detailedOperator = Dune::make_shared< DetailedOperatorType >(model->diffusion()->paramSize(),
                                                                   detailedOperatorComponents,
                                                                   model->diffusion()->coefficients());
  info << "done (took " << timer.elapsed() << " sec)" << std::endl;

  return std::make_pair(detailedOperator, detailedFunctional);
} // ... initSolver(...)


std::vector< ParamType > getParameters(const DescriptionType& description,
                                       const ModelType& model,
                                       const std::string type)
{
  // get the paramters
  const DescriptionType& parameterDescription = description.sub("parameter");
  std::vector< ParamType > parameters;
  if (parameterDescription.hasSub(type)) {
    size_t index = 0;
    bool search = true;
    while (search) {
      if (parameterDescription.hasKey(type + "." + Dune::Stuff::Common::toString(index))) {
        parameters.push_back(
              parameterDescription.getDynVector< ParamFieldType >(type + "." + Dune::Stuff::Common::toString(index),
                                                                  model.paramSize()));
        ++index;
      } else
        search = false;
    }
  } else if (parameterDescription.hasKey(type)) {
    parameters.push_back(
          parameterDescription.getDynVector< ParamFieldType >(type,
                                                              model.paramSize()));
  } else
    DUNE_THROW(Dune::IOError,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
               << "no '" << type << "' parameter given!");
  for (auto p : parameters)
    assert(p.size() == model.paramSize());
  return parameters;
} // ... getParameters(...)


void nonparametericRun(const DescriptionType& description,
                       const Dune::shared_ptr< const ModelType >& model,
                       const Dune::shared_ptr< const MsGridType >& msGrid,
                       const Dune::shared_ptr< const BoundaryInfoType >& boundaryInfo)
{
  // logger
  const std::string filename = description.get< std::string >("filename");
  const bool debugLogging = description.get< bool >("logging.debug", false);
  Dune::Stuff::Common::LogStream& info    = Dune::Stuff::Common::Logger().info();
  Dune::Stuff::Common::LogStream& debug   = Dune::Stuff::Common::Logger().debug();
//  Dune::Stuff::Common::LogStream& devnull = Dune::Stuff::Common::Logger().devnull();
  Dune::Timer timer;
  // linear solver
  const DescriptionType& discretizationDescription = description.sub(SolverType::id());
  const std::string    linearSolverType      = discretizationDescription.get< std::string >(   "linearsolver.type",
                                                                                               "bicgstab.ilut");
  const size_t         linearSolverMaxIter   = discretizationDescription.get< size_t >(        "linearsolver.maxIter",
                                                                                               5000u);
  const RangeFieldType linearSolverPrecision = discretizationDescription.get< RangeFieldType >("linearsolver.precision",
                                                                                               1e-12);

  info << "initializing solver";
  if (!debugLogging)
    info << "... " << std::flush;
  else
    info << ":" << std::endl;
  SolverType solver(model, msGrid, boundaryInfo, discretizationDescription.get< RangeFieldType >("penaltyFactor"));
  solver.init("  ", debug);
  if (!debugLogging)
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

  info << "solving";
  if (!debugLogging)
    info << "... " << std::flush;
  else
    info << ":" << std::endl;
  timer.reset();
  std::vector< Dune::shared_ptr< VectorType > > solution = solver.createVector();
  solver.solve(solution,
               linearSolverType,
               linearSolverMaxIter,
               linearSolverPrecision,
               "  ",
               debug);
  if (!debugLogging)
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;
  info << "writing solution to disc";
  if (!debugLogging)
    info << "... " << std::flush;
  else
    info << ":" << std::endl;
  timer.reset();
  solver.visualize(solution,
                   filename + ".solution",
                   id() + ".solution",
                   "  ",
                   debug);
  if (!debugLogging)
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;
} // ... nonparametericRun(...)


Dune::shared_ptr< SolverType > computeLocalScalarProducts(const DescriptionType& description,
                                                          const Dune::shared_ptr< const ModelType >& model,
                                                          const Dune::shared_ptr< const MsGridType >& msGrid,
                                                          const Dune::shared_ptr< const BoundaryInfoType >& boundaryInfo,
                                                          const ParamType& mu,
                                                          std::vector< Dune::shared_ptr< ScalarProductType > >& localScalarProducts)
{
  // logger
  Dune::Stuff::Common::LogStream& devnull = Dune::Stuff::Common::Logger().devnull();
  const RangeFieldType penaltyFactor = description.get< RangeFieldType >(SolverType::id() + ".penaltyFactor");
  assert(mu.size() == model->paramSize());

  Dune::shared_ptr< SolverType > fixedGlobalSolver = Dune::make_shared< SolverType >(model->fix(mu),
                                                                                     msGrid,
                                                                                     boundaryInfo,
                                                                                     penaltyFactor);
  fixedGlobalSolver->init("", devnull);

  localScalarProducts = fixedGlobalSolver->systemMatrices(ParamType());

  return fixedGlobalSolver;
} // ... computeLocalScalarProducts(...)


void computeGlobalSnapshot(const DescriptionType& description,
                           const SolverType& fixedGlobalSolver,
                           const DetailedOperatorType& detailedOperator,
                           const DetailedFunctionalType& detailedFunctional,
                           const ParamType& mu,
                           std::vector< Dune::shared_ptr< VectorType > >& solution,
                           Dune::shared_ptr< VectorType >& tmpGlobalVector)
{
  // logger
  Dune::Stuff::Common::LogStream& devnull = Dune::Stuff::Common::Logger().devnull();
  const bool visualize = description.get< bool >(id() + ".visualize", true);
  // get linear solver
  const DescriptionType& discretizationDescription = description.sub(SolverType::id());
  const std::string    linearSolverType      = discretizationDescription.get< std::string >(   "linearsolver.type",
                                                                                               "bicgstab.ilut");
  const size_t         linearSolverMaxIter   = discretizationDescription.get< size_t >(        "linearsolver.maxIter",
                                                                                               5000u);
  const RangeFieldType linearSolverPrecision = discretizationDescription.get< RangeFieldType >("linearsolver.precision",
                                                                                               1e-12);
  const Dune::shared_ptr< const LinearSolverType >
      linearSolver(Dune::Stuff::LA::Solver::create< ScalarProductType, VectorType >(linearSolverType));
  // create system matrix
  const auto systemMatrix = detailedOperator.fix(mu);
  // solve
  linearSolver->apply(*systemMatrix,
                      detailedFunctional,
                      *tmpGlobalVector,
                      linearSolverMaxIter,
                      linearSolverPrecision);
  // copy global vector to local ones
  fixedGlobalSolver.localizeVector(tmpGlobalVector, solution);
  if (visualize) {
    fixedGlobalSolver.visualize(solution,
                     "detailedSolution.param_" + Dune::Stuff::Common::Parameter::reportWOwhitespace(mu),
                     "detailed solution to parameter " + Dune::Stuff::Common::Parameter::report(mu),
                     "",
                     devnull);
  }
} // ... computeGlobalSnapshot(...)


Dune::shared_ptr< MapperType > extendReducedBases(const DescriptionType& description,
                                                  const MarkerType& markedSubdomains,
                                                  const std::vector< Dune::shared_ptr< ScalarProductType > >& localScalarProducts,
                                                  std::vector< Dune::shared_ptr< VectorType > >& localSnapshots,
                                                  std::vector< Dune::shared_ptr< ReducedBasisType > >& localReducedBases)
{
  // logger
  Dune::Stuff::Common::LogStream& debug   = Dune::Stuff::Common::Logger().debug();
  Dune::Timer timer;
  // gram schmidt
  const bool applyGramSchmidt = description.get< bool >("enrichment.gramSchmidt", true);
  const double gramSchmidtTolerance = description.get< double >("enrichment.gramSchmidtTolerance", 1e-10);
  // loop over all marked subdomains
  for (int ii = 0; ii < markedSubdomains.size(); ++ii) {
    const size_t subdomain = markedSubdomains.get(ii);
    debug << "  on subdomain " << subdomain << "... " << std::flush;
    timer.reset();
    // get the local scalar product
    const ScalarProductType& localScalarProduct = *(localScalarProducts[subdomain]);
    // get the local snapshot
    VectorType& localSnapshot = *(localSnapshots[subdomain]);
    // get the local reduced basis
    ReducedBasisType& localReducedBasis = *(localReducedBases[subdomain]);
    // apply a gram schmidt to the local snapshot (inplace) wrt to the existing reduced basis
    bool extend = false;
    if (applyGramSchmidt)
      extend = Dune::Stuff::LA::Algorithm::gramSchmidt(localScalarProduct,
                                                       localReducedBasis,
                                                       localSnapshot,
                                                       gramSchmidtTolerance);
    if (extend) {
      const size_t oldSize = localReducedBasis.cols();
      // resize the local reduced basis matrix by keeping the old columns
      localReducedBasis.backend().conservativeResize(::Eigen::NoChange, oldSize + 1);
      // and set the new column to the orthogonalized local snapshot
      localReducedBasis.backend().col(oldSize) = localSnapshot.backend();
      debug << "done (from " << oldSize << " to " << localReducedBasis.cols() << ", took " << timer.elapsed() << "sec)" << std::endl;
    } else
      debug << "failed" << std::endl;
  } // loop over all marked subdomains
  // also creeate a new mapper
  Dune::shared_ptr< MapperType > mapper = Dune::make_shared< MapperType >();
  mapper->prepare();
  for (size_t ii = 0; ii < localReducedBases.size(); ++ii)
    mapper->add(ii, localReducedBases[ii]->cols());
  mapper->finalize();
  return mapper;
} // ... extendReducedBases(...)


Dune::shared_ptr< MapperType > initReducedBasis(const DescriptionType& description,
                                                const SolverType& fixedGlobalSolver,
                                                const DetailedOperatorType& detailedOperator,
                                                const DetailedFunctionalType& detailedFunctional,
                                                const std::vector< Dune::shared_ptr< ScalarProductType > >& localScalarProducts,
                                                const std::vector< ParamType >& trainingParameters,
                                                std::vector< Dune::shared_ptr< ReducedBasisType > >& localReducedBases,
                                                std::vector< Dune::shared_ptr< VectorType > >& tmpLocalVectors,
                                                Dune::shared_ptr< VectorType >& tmpGlobalvector)
{
  // logger
  Dune::Stuff::Common::LogStream& info = Dune::Stuff::Common::Logger().debug();
  Dune::Timer timer;
  // prepare
  const bool debugLogging = description.get< bool >("logging.debug", false);
  const bool applyGramSchmidt = description.get< bool >("reducedBasis.gramSchmidt", true);
  // sanity checks
  assert(trainingParameters.size() > 0);
  assert(tmpLocalVectors.size() == fixedGlobalSolver.msGrid()->size());
  assert(tmpLocalVectors.size() == localScalarProducts.size());
  assert(tmpLocalVectors.size() == localReducedBases.size());
  // compute the first snapshot
  info << "initializing reduced basis"
       << " with parameter " << Dune::Stuff::Common::Parameter::report(trainingParameters[0]) << "..." << std::flush;
  timer.reset();
  computeGlobalSnapshot(description,
                        fixedGlobalSolver,
                        detailedOperator,
                        detailedFunctional,
                        trainingParameters[0],
                        tmpLocalVectors,        // <-solution
                        tmpGlobalvector);
  // loop over all subdomains
  Dune::shared_ptr< MapperType > mapper = Dune::make_shared< MapperType >();
  mapper->prepare();
  for (size_t ss = 0; ss < localReducedBases.size(); ++ss) {
    // get local scalar product
    const ScalarProductType& localScalarProduct = *(localScalarProducts[ss]);
    // get local snapshot
    VectorType& localSnapshot = *(tmpLocalVectors[ss]);
    assert(localScalarProduct.rows() == localSnapshot.size());
    assert(localScalarProduct.cols() == localSnapshot.size());
    // normalize local snapshot (inplace), if whished
    if (applyGramSchmidt)
      Dune::Stuff::LA::Algorithm::normalize(localScalarProduct, localSnapshot);
    // set local reduced basis
    ReducedBasisType& localReducedBasis = *(localReducedBases[ss]);
    localReducedBasis.backend().col(0) = localSnapshot.backend();
    mapper->add(ss, localReducedBasis.cols());
  } // loop over all subdomains
  mapper->finalize();
  info << " done (took " << timer.elapsed() << "sec)" << std::endl;
  // now extend the basis for all other parameters
  MarkerType markedSubdomains(fixedGlobalSolver.msGrid()->size());
  for (int ii = 0; ii < markedSubdomains.size(); ++ii)
    markedSubdomains.set(ii, ii);
  for (size_t pp = 1; pp < trainingParameters.size(); ++pp) {
    // compute the next snapshot
    info << "extending reduced basis with parameter " << Dune::Stuff::Common::Parameter::report(trainingParameters[pp]);
    if (!debugLogging)
      info << "... " << std::flush;
    else
      info << ":" << std::endl;
    timer.reset();
    computeGlobalSnapshot(description,
                          fixedGlobalSolver,
                          detailedOperator,
                          detailedFunctional,
                          trainingParameters[pp],
                          tmpLocalVectors,        // <- solution
                          tmpGlobalvector);
    // extend basis and get the new mapper
    mapper = extendReducedBases(description,
                                markedSubdomains,
                                localScalarProducts,
                                tmpLocalVectors,
                                localReducedBases);
    if (!debugLogging)
      info << "done (took " << timer.elapsed() << "sec)" << std::endl;
  } // now extend the basis for all other parameters
  return mapper;
} // ... initReducedBasis(...)


std::pair<  Dune::shared_ptr< ReducedOperatorType >,
            Dune::shared_ptr< ReducedFunctionalType > >
createReducedOperatorFunctional(const SolverType& fixedGlobalSolver,
                                const DetailedOperatorType& detailedOperator,
                                const DetailedFunctionalType& detailedFunctional,
                                const std::vector< Dune::shared_ptr< ReducedBasisType > >& localReducedBases,
                                const MapperType& mapper,
                                Dune::shared_ptr< VectorType >& tmpLocalVector,
                                Dune::shared_ptr< VectorType >& tmpGlobalVectorOne,
                                Dune::shared_ptr< VectorType >& tmpGlobalVectorTwo)
{
  // logger
  Dune::Stuff::Common::LogStream& debug = Dune::Stuff::Common::Logger().debug();
  Dune::Timer timer;
  debug << "  from ";
  debug << detailedOperator.components()[0]->rows() << "x" << detailedOperator.components()[0]->cols()
        << " to " << mapper.size() << "x" << mapper.size() << "... " << std::flush;
  timer.reset();
  // init reduced operator
  std::vector< Dune::shared_ptr< ReducedBasisType > > reducedOperatorComponents(detailedOperator.numComponents());
  for (auto& component : reducedOperatorComponents)
    component = Dune::make_shared< ReducedBasisType >(mapper.size(), mapper.size());
  Dune::shared_ptr< ReducedOperatorType >
      reducedOperator = Dune::make_shared< ReducedOperatorType >(detailedOperator.paramSize(),
                                                                 reducedOperatorComponents,
                                                                 detailedOperator.coefficients());
  // init reduced functional
  Dune::shared_ptr< ReducedFunctionalType >
      reducedFunctional = Dune::make_shared< ReducedFunctionalType >(mapper.size());
  // get the multiscale grid
  const MsGridType& msGrid = *(fixedGlobalSolver.msGrid());
  // loop over all subdomains
  for (size_t subdomain = 0; subdomain < msGrid.size(); ++subdomain) {
    // get local reduced basis
    const ReducedBasisType& reducedBasisSubdomain = *(localReducedBases[subdomain]);
    // loop over all local basis functions
    for (int ii = 0; ii < reducedBasisSubdomain.cols(); ++ii) {
      // convert the iith basis function to a global one
      tmpLocalVector->backend() = reducedBasisSubdomain.backend().col(ii);
      fixedGlobalSolver.globalizeVector(tmpLocalVector, subdomain, tmpGlobalVectorOne);
      // compute global index
      const size_t iiGlobal = mapper.toGlobal(subdomain, ii);
      // project the detailed functional wrt iith local basis function
      reducedFunctional->set(iiGlobal, tmpGlobalVectorOne->backend().transpose() * detailedFunctional.backend());
      // loop over all local basis functions
      for (int jj = 0; jj < reducedBasisSubdomain.cols(); ++jj) {
        // convert the jjth basis function to a global one
        tmpLocalVector->backend() = reducedBasisSubdomain.backend().col(jj);
        fixedGlobalSolver.globalizeVector(tmpLocalVector, subdomain, tmpGlobalVectorTwo);
        // compute global index
        const size_t jjGlobal = mapper.toGlobal(subdomain, jj);
        // loop over all components of the operator
        for (size_t qq = 0; qq < detailedOperator.numComponents(); ++qq) {
          // get the components
          const auto& detailedOperatorComponent = *(detailedOperator.components()[qq]);
          ReducedBasisType& reducedOperatorComponent = *(reducedOperator->components()[qq]);
          // project the detailed component wrt the iith and jjth basis function
          reducedOperatorComponent.set(iiGlobal,
                                       jjGlobal,
                                       tmpGlobalVectorOne->backend().transpose()
                                       * detailedOperatorComponent.backend()
                                       * tmpGlobalVectorTwo->backend());
        } // loop over all components of the operator
      } // loop over all local basis functions
    } // loop over all local basis functions
    // loop over all neighbors
    for (size_t neighbor : msGrid.neighborsOf(subdomain)) {
      // get local reduced basis
      const ReducedBasisType& reducedBasisNeighbor = *(localReducedBases[neighbor]);
      // loop over all local basis functions of the subdomain
      for (int ii = 0; ii < reducedBasisSubdomain.cols(); ++ii) {
        // convert the iith basis function to a global one
        tmpLocalVector->backend() = reducedBasisSubdomain.backend().col(ii);
        fixedGlobalSolver.globalizeVector(tmpLocalVector, subdomain, tmpGlobalVectorOne);
        // compute global index
        const size_t iiGlobal = mapper.toGlobal(subdomain, ii);
        // loop over all the local basis functions of the neighbor
        for (int jj = 0; jj < reducedBasisNeighbor.cols(); ++jj) {
          // convert the jjth basis function to a global one
          tmpLocalVector->backend() = reducedBasisNeighbor.backend().col(jj);
          fixedGlobalSolver.globalizeVector(tmpLocalVector, neighbor, tmpGlobalVectorTwo);
          // compute global index
          const size_t jjGlobal = mapper.toGlobal(neighbor, jj);
          // loop over all components of the operator
          for (size_t qq = 0; qq < detailedOperator.numComponents(); ++qq) {
            // get the components
            const auto& detailedOperatorComponent = *(detailedOperator.components()[qq]);
            ReducedBasisType& reducedOperatorComponent = *(reducedOperator->components()[qq]);
            // project the detailed component wrt the iith and jjth basis function
            reducedOperatorComponent.set(iiGlobal,
                                         jjGlobal,
                                         tmpGlobalVectorOne->backend().transpose()
                                         * detailedOperatorComponent.backend()
                                         * tmpGlobalVectorTwo->backend());
          } // loop over all components of the operator
        } // loop over all the local basis functions of the neighbor
      } // loop over all local basis functions of the subdomain
    } // loop over all neighbors
  } // loop over all subdomains

  // finished
  debug << "done (took " << timer.elapsed() << "sec)" << std::endl;
  return std::make_pair(reducedOperator, reducedFunctional);
} // ... reducedOperator(...)


VectorType computeReducedSolution(const ReducedOperatorType& reducedOperator,
                                  const ReducedFunctionalType& reducedFunctional,
                                  const ParamType& mu)
{
  VectorType solution;
  const auto reducedSystemMatrix = reducedOperator.fix(mu);
  solution.backend() = reducedSystemMatrix->backend().colPivHouseholderQr().solve(reducedFunctional.backend());
  return solution;
} // ... computeReducedSolution(...)


void reconstructReducedSolution(const std::vector< Dune::shared_ptr< ReducedBasisType > >& localReducedBases,
                                const MapperType& mapper,
                                const VectorType& reducedSolution,
                                std::vector< Dune::shared_ptr< VectorType > >& reconstruction)
{
  // loop over all subdomains
  for (size_t ss = 0; ss < mapper.numSubdomains(); ++ss) {
    // get the local reduced basis
    const ReducedBasisType& localReducedBasis = *(localReducedBases[ss]);
    // get the local vector
    VectorType& localVector = *(reconstruction[ss]);
    // clear it
    localVector.backend().setZero();
    // loop over all local basis functions
    for (int ii = 0; ii < localReducedBasis.cols(); ++ii) {
      // and reconstruct
      localVector.backend() += localReducedBasis.backend().col(ii) * reducedSolution.get(mapper.toGlobal(ss, ii));
    } // loop over all local basis functions
  } // loop over all subdomains
} // ... reconstructReducedSolution(...)


VectorType computeLocalErrors(const SolverType& fixedGlobalSolver,
                              const DescriptionType& description,
                              const DetailedOperatorType& detailedOperator,
                              const DetailedFunctionalType& detailedFunctional,
                              const std::vector< Dune::shared_ptr< ScalarProductType > >& localScalarProducts,
                              const std::vector< Dune::shared_ptr< ReducedBasisType > >& localReducedBases,
                              const MapperType& mapper,
                              const VectorType& reducedSolution,
                              const ParamType& mu,
                              std::vector< Dune::shared_ptr< VectorType > >& tmpLocalVectorsA,
                              std::vector< Dune::shared_ptr< VectorType > >& tmpLocalVectorsB,
                              Dune::shared_ptr< VectorType >& tmpGlobalVector)
{
  // logger
  Dune::Stuff::Common::LogStream& devnull = Dune::Stuff::Common::Logger().devnull();
  const bool visualize = description.get< bool >(id() + ".visualize", false);

  // compute detailed solution
  computeGlobalSnapshot(description,
                        fixedGlobalSolver,
                        detailedOperator,
                        detailedFunctional,
                        mu,
                        tmpLocalVectorsA,   // <- solution
                        tmpGlobalVector);

  // reconstruct reudced solution
  reconstructReducedSolution(localReducedBases,
                             mapper,
                             reducedSolution,
                             tmpLocalVectorsB); // <- reconstruction

  // compute difference
  for (size_t ss = 0; ss < mapper.numSubdomains(); ++ss)
    tmpLocalVectorsB[ss]->backend() -= tmpLocalVectorsA[ss]->backend();

  // visualize error (if wished)
  if (visualize) {
    fixedGlobalSolver.visualize(tmpLocalVectorsA,
                                "difference.param_" + Dune::Stuff::Common::Parameter::reportWOwhitespace(mu),
                                "difference between detailed and reduced solution to parameter " + Dune::Stuff::Common::Parameter::report(mu),
                                "",
                                devnull);
  } // visualize error (if wished)

  // compute local errors
  VectorType localErrors(mapper.numSubdomains());
  for (size_t ss = 0; ss < mapper.numSubdomains(); ++ss) {
    // get local scalar product
    const ScalarProductType& localScalarProduct = *(localScalarProducts[ss]);
    // get local reference solution
    const VectorType& localReference = *(tmpLocalVectorsA[ss]);
    // get local difference
    const VectorType& localDifference = *(tmpLocalVectorsB[ss]);
    // compute local relative error
    localErrors.set(ss,
                    std::sqrt(localDifference.backend().transpose()
                              * localScalarProduct.backend()
                              * localDifference.backend())
                    / std::sqrt(localReference.backend().transpose()
                                * localScalarProduct.backend()
                                * localReference.backend()));
  } // compute local errors

  return localErrors;
} // ... computeLocalErrors(...)


MarkerType markSubdomains(const DescriptionType& description,
                          const VectorType& localErrors,
                          MarkerType& ageOfSubdomains)
{
  const std::string strategy = description.get< std::string >("markSubdomains.strategy", "factorTimesMean");
  MarkerType markedSubdomains;
  if (strategy == "factorTimesMean") {
    const RangeFieldType factor = description.get< RangeFieldType >("markSubdomains.markingFactor", 1.0);
    const size_t maxAge = description.get< RangeFieldType >("markSubdomains.age", 2);
    const RangeFieldType mean = localErrors.backend().mean();
    const auto which = (localErrors.backend().array() > factor*mean);
    markedSubdomains = MarkerType(0);
    for (int ss = 0; ss < localErrors.size(); ++ss) {
      if (which[ss]) {
        const int oldSize = markedSubdomains.size();
        markedSubdomains.backend().conservativeResize(oldSize + 1);
        markedSubdomains.set(oldSize, ss);
        // since it has been marked, reset the age
        ageOfSubdomains.set(ss, 0);
      } else {
        // if we do not mark it because of the error, lets look at the age
        const size_t oldAge = ageOfSubdomains.get(ss);
        // mark it, if it is too old
        if (oldAge >= maxAge) {
          const int oldSize = markedSubdomains.size();
          markedSubdomains.backend().conservativeResize(oldSize + 1);
          markedSubdomains.set(oldSize, ss);
          // since it has been marked, reset the age
          ageOfSubdomains.set(ss, 0);
        } else {
          // inrease age
          ageOfSubdomains.set(ss, oldAge + 1);
        }
      }
    }
  } else if (strategy == "all") {
    markedSubdomains = MarkerType(localErrors.size());
    for (int ss = 0; ss < localErrors.size(); ++ss)
      markedSubdomains.set(ss, ss);
  } else
    DUNE_THROW(Dune::InvalidStateException,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
               << " wrong 'markSubdomains.strategy' given (has to be one of 'factorTimesMean', 'all', is'"
               << strategy << "')!");
  return markedSubdomains;
} // ... markSubdomains(...)


void computeLocalCorrections(const SolverType& globalSolver,
                             const ModelType& model,
                             const DetailedOperatorType& /*detailedOperator*/,
                             const DetailedFunctionalType& /*detailedFunctional*/,
                             const DescriptionType& description,
                             const std::vector< Dune::shared_ptr< ReducedBasisType > >& localReducedBases,
                             const MapperType& mapper,
                             const VectorType& reducedSolution,
                             const ParamType& mu,
                             const MarkerType& markedSubdomains,
                             std::vector< Dune::shared_ptr< VectorType > >& tmpLocalVectorsA,
                             std::vector< Dune::shared_ptr< VectorType > >& tmpLocalVectorsB)
{
  // logger
  Dune::Stuff::Common::LogStream& debug   = Dune::Stuff::Common::Logger().debug();
  Dune::Stuff::Common::LogStream& devnull = Dune::Stuff::Common::Logger().devnull();
  const bool visualize = description.get< bool >(id() + ".visualize", true);
  Dune::Timer timer;
  // linear solver
  const DescriptionType& discretizationDescription = description.sub(SolverType::id());
  const std::string    linearSolverType      = discretizationDescription.get< std::string >(   "linearsolver.type",
                                                                                               "bicgstab.ilut");
  const size_t         linearSolverMaxIter   = discretizationDescription.get< size_t >(        "linearsolver.maxIter",
                                                                                               5000u);
  const RangeFieldType linearSolverPrecision = discretizationDescription.get< RangeFieldType >("linearsolver.precision",
                                                                                               1e-12);
  const std::string localBoundaryInfoType = description.get< std::string >("enrichment.boundaryInfo",
                                                                           "stuff.grid.boundaryinfo.alldirichlet");
  debug << "  reconstructing reduced solution... " << std::flush;
  reconstructReducedSolution(localReducedBases,
                             mapper,
                             reducedSolution,
                             tmpLocalVectorsA); // <- reconstruction
  const Dune::shared_ptr< const GlobalDiscreteFunctionType >
      reconstructedSolutionFunction = globalSolver.createDiscreteFunction(tmpLocalVectorsA);
  debug << "done (took " << timer.elapsed() << "sec)" << std::endl;
  debug << "  solving for local corrections... " << std::flush;
  timer.reset();
  // create a fixed model for the current parameter
  const Dune::shared_ptr< const ModelType > fixedModel = model.fix(mu);

//  // create a global fixed solver
//  SolverType fixedSolver(fixedModel, globalSolver.msGrid(), globalSolver.boundaryInfo(), globalSolver.penaltyFactor());
//  fixedSolver.init("", devnull);
//  // correct rhs
//  auto tmpGlobalVector = fixedSolver.globalizeVectors(tmpLocalVectorsA);
//  fixedSolver.vector("force")->backend() -= tmpGlobalVector->backend().transpose()
//                                            * fixedSolver.systemMatrix("diffusion")->components()[0]->backend();
//  // solve
//  fixedSolver.solve(tmpLocalVectorsA,       // <- solution
//                    linearSolverType,
//                    linearSolverMaxIter,
//                    linearSolverPrecision,
//                    "",
//                    devnull);
//  const Dune::shared_ptr< const GlobalDiscreteFunctionType >
//      globalCorrectionFunction = globalSolver.createDiscreteFunction(tmpLocalVectorsA);
//  fixedSolver.visualize(tmpLocalVectorsA,
//                        "globalCorrection_parameter_" + Dune::Stuff::Common::Parameter::reportWOwhitespace(mu),
//                        "global correction to parameter " + Dune::Stuff::Common::Parameter::reportWOwhitespace(mu),
//                        "", devnull);

  // get the multiscale grid
  const MsGridType& msGrid = *(globalSolver.msGrid());
  // loop over all marked subdomains
  for (int ii = 0; ii < markedSubdomains.size(); ++ii) {
    const size_t subdomain = markedSubdomains.get(ii);
    // get the local grid part
    const Dune::shared_ptr< const LocalGridPartType > localOversampledGridPart = msGrid.localGridPart(subdomain, true);
    // initialize a local solver on the local oversampled gripart
    DescriptionType localBoundaryInfoDescription;
    localBoundaryInfoDescription["dirichlet"] = "[1; 2; 3; 4; 5; 6]";
    localBoundaryInfoDescription["neumann"] = "[7]";
    const Dune::shared_ptr< const LocalBoundaryInfoType > localBoundaryInfo(
          Dune::Stuff::Grid::BoundaryInfo::create< typename MsGridType::LocalGridViewType >(localBoundaryInfoType,
                                                                                            localBoundaryInfoDescription));
    LocalSolverType localOversampledSolver(localOversampledGridPart, localBoundaryInfo, fixedModel);
    localOversampledSolver.init("", devnull);
    // get the local solver on the non oversampled gridpart
    const LocalSolverType& localSolver = *(globalSolver.localSolver(subdomain));
    // create a local function to hold the reconstructed solution
    Dune::shared_ptr< LocalDiscreteFunctionType >
        localReconstructedSolutionFunction = localOversampledSolver.createAnsatzFunction();
    // and project the current reconstructed solution onto the local oversampled gripart
    Dune::Stuff::HeterogenousProjection<>::project(*reconstructedSolutionFunction,
                                                 *localReconstructedSolutionFunction);

    // compute the local right hand side f := f - lrs' * A
    // * therefore get the local system matrix and right hand side
    const auto localSystemMatrixPtr = localOversampledSolver.matrix("diffusion");
    assert(localSystemMatrixPtr->numComponents() == 1);
    const auto& localSystemMatrix = *(localSystemMatrixPtr->components()[0]);
    auto forceVectorPtr = localOversampledSolver.vector("force");
    assert(forceVectorPtr->numComponents() == 1);
    auto& localRightHandSide = *(forceVectorPtr->components()[0]);
    // * get the local reconstructed solution
    auto& localReconstructedSolution = *(localReconstructedSolutionFunction->vector());
    // * and correct the right hand side
    localRightHandSide.backend() -=localReconstructedSolution.backend().transpose() * localSystemMatrix.backend();
    // then solve locally on the oversampled domain
    Dune::shared_ptr< VectorType > localOversampledSolution = localOversampledSolver.createAnsatzVector();
    localOversampledSolver.solve(localOversampledSolution,
                                 linearSolverType,
                                 linearSolverMaxIter,
                                 linearSolverPrecision,
                                 "",
                                 devnull);
    // and restrict the oversampled solution to the local (not oversampled) subdomain
    Dune::shared_ptr< VectorType > localSolution = localSolver.createAnsatzVector();
//    Dune::Stuff::DiscreteFunction::Projection::Lagrange::project(*localOversampledSolutionFunction,
//                                                                 *localSolutionFunction);
    localSolution->backend() = localOversampledSolution->backend().head(localSolution->size());
    // visualize if wished
    if (visualize) {
      localOversampledSolver.visualizeAnsatzVector(localOversampledSolution,
                                        "localOversampledCorrection_subdomain_" + Dune::Stuff::Common::toString(subdomain) + ".param_" + Dune::Stuff::Common::Parameter::reportWOwhitespace(mu),
                                        "local oversampled correction on subdomain " + Dune::Stuff::Common::toString(subdomain) + " to parameter " + Dune::Stuff::Common::Parameter::report(mu),
                                        "",
                                        devnull);
      localSolver.visualizeAnsatzVector(localSolution,
                                        "localCorrection_subdomain_" + Dune::Stuff::Common::toString(subdomain) + ".param_" + Dune::Stuff::Common::Parameter::reportWOwhitespace(mu),
                                        "local correction on subdomain " + Dune::Stuff::Common::toString(subdomain) + " to parameter " + Dune::Stuff::Common::Parameter::report(mu),
                                        "",
                                        devnull);
//      localOversampledSolutionFunction->vector()->backend() -= tmpLocalVectorsA[0]->backend();
//      localOversampledSolver.visualizeAnsatzVector(localOversampledSolutionFunction->vector(),
//                            "differenceLocalOversampledCorrections_param_" + Dune::Stuff::Common::Parameter::reportWOwhitespace(mu),
//                            "difference local oversampled corrections to parameter " + Dune::Stuff::Common::Parameter::report(mu),
//                            "", devnull);
//      localSolutionFunction->vector()->backend() -= tmpLocalVectorsA[0]->backend();
//      localSolver.visualizeAnsatzVector(localSolutionFunction->vector(),
//                                        "differenceLocalCorrection_subdomain_" + Dune::Stuff::Common::toString(subdomain) + ".param_" + Dune::Stuff::Common::Parameter::reportWOwhitespace(mu),
//                                        "difference local correction on subdomain " + Dune::Stuff::Common::toString(subdomain) + " to parameter " + Dune::Stuff::Common::Parameter::report(mu),
//                                        "",
//                                        devnull);

//      // for the cheating
//      Dune::shared_ptr< LocalDiscreteFunctionType > localizedGlobalCorrectionFunction
//          = localSolver.createAnsatzFunction();
//      Dune::Stuff::DiscreteFunction::Projection::Lagrange::project(*globalCorrectionFunction,
//                                                                   *localizedGlobalCorrectionFunction);
//      localSolver.visualizeAnsatzVector(localizedGlobalCorrectionFunction->vector(),
//                                        "localizedGlobalCorrection_subdomain_" + Dune::Stuff::Common::toString(subdomain) + ".param_" + Dune::Stuff::Common::Parameter::reportWOwhitespace(mu),
//                                        "localized global correction on subdomain " + Dune::Stuff::Common::toString(subdomain) + " to parameter " + Dune::Stuff::Common::Parameter::report(mu),
//                                        "",
//                                        devnull);
//      // also compute the difference
//      localizedGlobalCorrectionFunction->vector()->backend() -= localSolutionFunction->vector()->backend();
//      localSolver.visualizeAnsatzVector(localizedGlobalCorrectionFunction->vector(),
//                                        "localCorrectionDifference_subdomain_" + Dune::Stuff::Common::toString(subdomain) + ".param_" + Dune::Stuff::Common::Parameter::reportWOwhitespace(mu),
//                                        "local crrection difference on subdomain " + Dune::Stuff::Common::toString(subdomain) + " to parameter " + Dune::Stuff::Common::Parameter::report(mu),
//                                        "",
//                                        devnull);
    } // visualize if wished
    // now just copy the local solution into the return argument
    tmpLocalVectorsB[subdomain]->backend() = localSolution->backend();
  } // loop over all marked subdomains
  debug << "done (took " << timer.elapsed() << "sec)" << std::endl;

//  fixedSolver.visualize(tmpLocalVectorsB,
//                        "localCorrections_parameter_" + Dune::Stuff::Common::Parameter::reportWOwhitespace(mu),
//                        "local corrections to parameter " + Dune::Stuff::Common::Parameter::reportWOwhitespace(mu),
//                        "", devnull);
//  // now A nd B should be the same
//  for (size_t ii = 0; ii < tmpLocalVectorsA.size(); ++ii)
//    tmpLocalVectorsA[ii]->backend() -= tmpLocalVectorsB[ii]->backend();
//  fixedSolver.visualize(tmpLocalVectorsA,
//                        "differenceCorrections_parameter_" + Dune::Stuff::Common::Parameter::reportWOwhitespace(mu),
//                        "difference corrections to parameter " + Dune::Stuff::Common::Parameter::reportWOwhitespace(mu),
//                        "", devnull);
} // ... computeLocalCorrections(...)


int main(int argc, char** argv)
{
  try {
    // description
    const DescriptionType description = initDescription(argc, argv);
    // logger
    Dune::Stuff::Common::LogStream& info    = Dune::Stuff::Common::Logger().info();
    Dune::Stuff::Common::LogStream& debug   = Dune::Stuff::Common::Logger().debug();
    Dune::Stuff::Common::LogStream& devnull = Dune::Stuff::Common::Logger().devnull();
    Dune::Timer offlineTimer;
    Dune::Timer timer;
    const bool debugLogging = description.get< bool >("logging.debug", false);
    const bool visualize = description.get< bool >(id() + ".visualize", false);
    // get grid, model and boundaryinfo
    const Dune::shared_ptr< const MsGridType > msGrid = initGrid(description);
    const Dune::shared_ptr< const ModelType > model = initModel(description, *msGrid);
    const Dune::shared_ptr< const BoundaryInfoType > boundaryInfo = initBoundaryInfo(description);
    // do something
    if (model->parametric()) {
      // report
      std::string msg = "  starting offline phase  ";
      info << Dune::Stuff::Common::whitespaceify(msg, '=') << std::endl;
      info << msg << std::endl;
      info << Dune::Stuff::Common::whitespaceify(msg, '=') << std::endl;

      // get parameters
      const std::vector< ParamType > fixedParameters = getParameters(description, *model, "fixed");
      const std::vector< ParamType > trainingParameters = getParameters(description, *model, "train");
      const std::vector< ParamType > testParameters = getParameters(description, *model, "test");

      info << "computing local scalar products";
      if (!debugLogging)
        info << "... " << std::flush;
      else
        info << " for parameter " << Dune::Stuff::Common::Parameter::report(fixedParameters[0]) << "... " << std::flush;
      timer.reset();
      std::vector< Dune::shared_ptr< ScalarProductType > > localScalarProducts;
      const Dune::shared_ptr< const SolverType >
          fixedGlobalSolver = computeLocalScalarProducts(description,
                                                         model,
                                                         msGrid,
                                                         boundaryInfo,
                                                         fixedParameters[0],
                                                         localScalarProducts);
      info << " done (took " << timer.elapsed() << "sec)" << std::endl;

      info << "initializing tmp storage... " << std::flush;
      timer.reset();
      Dune::shared_ptr< VectorType > tmpLocalVectorA = fixedGlobalSolver->localSolver(0)->createAnsatzVector();
      std::vector< Dune::shared_ptr< VectorType > > tmpLocalVectorsA = fixedGlobalSolver->createVector();
      std::vector< Dune::shared_ptr< VectorType > > tmpLocalVectorsB = fixedGlobalSolver->createVector();
      std::vector< Dune::shared_ptr< ReducedBasisType > > localReducedBases(fixedGlobalSolver->msGrid()->size());
      for (size_t ss = 0; ss < fixedGlobalSolver->msGrid()->size(); ++ss)
        localReducedBases[ss]
            = Dune::make_shared< ReducedBasisType >(fixedGlobalSolver->localSolver(ss)->ansatzSpace()->map().size(), 1);
      std::vector< Dune::shared_ptr< ReducedBasisType > > localReducedBasesOld(fixedGlobalSolver->msGrid()->size());
      for (size_t ss = 0; ss < fixedGlobalSolver->msGrid()->size(); ++ss)
        localReducedBasesOld[ss]
            = Dune::make_shared< ReducedBasisType >(fixedGlobalSolver->localSolver(ss)->ansatzSpace()->map().size(), 1);
      Dune::shared_ptr< VectorType > tmpGlobalVectorA
          = Dune::make_shared< VectorType >(fixedGlobalSolver->ansatzMapper().size());
      Dune::shared_ptr< VectorType > tmpGlobalVectorB
          = Dune::make_shared< VectorType >(fixedGlobalSolver->ansatzMapper().size());
      info << " done (took " << timer.elapsed() << "sec)" << std::endl;

      //  detailed operator and functional
      const auto detailedOperatorAndFunctionalPair = initDetailedOperatorAndFunctional(description,
                                                                                       model,
                                                                                       msGrid,
                                                                                       boundaryInfo);
      const DetailedOperatorType& detailedOperator = *(detailedOperatorAndFunctionalPair.first);
      const DetailedFunctionalType& detailedFunctional = *(detailedOperatorAndFunctionalPair.second);

      if (!debugLogging)
        info << "initializing local reduced bases... " << std::flush;
      timer.reset();
      Dune::shared_ptr< MapperType >
          mapper = initReducedBasis(description,
                                    *fixedGlobalSolver,
                                    detailedOperator,
                                    detailedFunctional,
                                    localScalarProducts,
                                    trainingParameters,
                                    localReducedBases,
                                    tmpLocalVectorsA,
                                    tmpGlobalVectorA);
      if (!debugLogging)
        info << " done (took " << timer.elapsed() << "sec)" << std::endl;

      info << "reducing operator and functional";
      if (!debugLogging)
        info << "... " << std::flush;
      else
        info << ":" << std::endl;
      auto reducedPair = createReducedOperatorFunctional(*fixedGlobalSolver,
                                                         detailedOperator,
                                                         detailedFunctional,
                                                         localReducedBases,
                                                         *mapper,
                                                         tmpLocalVectorA,
                                                         tmpGlobalVectorA,
                                                         tmpGlobalVectorB);
      Dune::shared_ptr< const ReducedOperatorType > reducedOperator = reducedPair.first;
      Dune::shared_ptr< const ReducedFunctionalType > reducedFunctional = reducedPair.second;
      if (!debugLogging)
        info << " done (took " << timer.elapsed() << "sec)" << std::endl;

      // report
      msg = "  offline phase finished, took " + Dune::Stuff::Common::toString(offlineTimer.elapsed()) + "sec  ";
      info << Dune::Stuff::Common::whitespaceify(msg, '=') << std::endl;
      info << msg << std::endl;
      info << "  starting online phase" << std::endl;
      info << Dune::Stuff::Common::whitespaceify(msg, '=') << std::endl;

      // now we start with the online phase
      // * therefore get the enriching params
      const DescriptionType& enrichmentDescription = description.sub("enrichment");
      const double enrichmentMaxLocalError = enrichmentDescription.get< double >("maxLocalError");
      const size_t enrichmentMaxIterations = enrichmentDescription.get< size_t >("maxIterations");
      // * then loop over all test parameters
      for (size_t pp = 0; pp < testParameters.size(); ++pp) {
        const ParamType& mu = testParameters[pp];

        // compute a reduced solution
        info << "solving for parameter " << Dune::Stuff::Common::Parameter::report(mu) << "... " << std::flush;
        timer.reset();
        VectorType reducedSolution = computeReducedSolution(*reducedOperator, *reducedFunctional, mu);
        info << " done (took " << std::scientific << timer.elapsed() << std::fixed << "sec)" << std::endl;
        // reconstruct and visualize if wished
        if (visualize) {
          info << "  reconstructing solution... " << std::flush;
          timer.reset();
          reconstructReducedSolution(localReducedBases,
                                     *mapper,
                                     reducedSolution,
                                     tmpLocalVectorsA); // <- reconstruction
          info << " done (took " << std::scientific << timer.elapsed() << std::fixed << "sec)" << std::endl;
          info << "  visualizing reconstructed solution..." << std::flush;
          timer.reset();
          fixedGlobalSolver->visualize(tmpLocalVectorsA,
                                       "reducedSolution.param_" + Dune::Stuff::Common::Parameter::reportWOwhitespace(mu),
                                       "reduced solution to param " + Dune::Stuff::Common::Parameter::report(mu),
                                       "",
                                       devnull);
          info << " done (took " << timer.elapsed() << "sec)" << std::endl;
        } // if (visualize)

        info << "  computing local errors... " << std::flush;
        timer.reset();
        VectorType localErrors = computeLocalErrors(*fixedGlobalSolver,
                                                    description,
                                                    detailedOperator,
                                                    detailedFunctional,
                                                    localScalarProducts,
                                                    localReducedBases,
                                                    *mapper,
                                                    reducedSolution,
                                                    mu,
                                                    tmpLocalVectorsA,
                                                    tmpLocalVectorsB,
                                                    tmpGlobalVectorA);
        VectorType localErrorsOld = localErrors;
        info << " done (took " << timer.elapsed() << "sec, max local error is "
             << std::scientific << localErrors.backend().maxCoeff() << std::fixed << ")" << std::endl;
        // enrich, if max local error is too large
        if (localErrors.backend().maxCoeff() > enrichmentMaxLocalError) {
          msg = "  starting local offline phase  ";
          info << Dune::Stuff::Common::whitespaceify(msg, '=') << std::endl;
          info << msg << std::endl;
          info << Dune::Stuff::Common::whitespaceify(msg, '=') << std::endl;

          // iterate until error is sufficient or we reach max iterations
          size_t enrichmentIteration = 0;
          size_t oldBasisSize = 0;
//          double oldMinLocalError = localErrors.backend().minCoeff();
          std::vector< VectorType > localErrorsVectors;
          std::vector< size_t > totalBasisSizes;
          localErrorsVectors.push_back(localErrors);
          totalBasisSizes.push_back(mapper->size());
          MarkerType ageOfSubdomains(msGrid->size());
          while (enrichmentIteration < enrichmentMaxIterations
                 && localErrors.backend().maxCoeff() > enrichmentMaxLocalError
//                 && !(oldMinLocalError < localErrors.backend().minCoeff())
                 && mapper->size() > oldBasisSize) {
            if (debugLogging) {
              msg = "enrichment iteration " + Dune::Stuff::Common::toString(enrichmentIteration);
              info << msg << std::endl;
              info << Dune::Stuff::Common::whitespaceify(msg, '=') << std::endl;
            }
            debug << "local errors are: "
                 << std::scientific << localErrors.backend().transpose() << std::fixed << std::endl;
//            // enrich on all subdomains
//            MarkerType markedSubdomains(fixedGlobalSolver->msGrid()->size());
//            for (size_t ii = 0; ii < fixedGlobalSolver->msGrid()->size(); ++ii)
//              markedSubdomains.set(ii, ii);
            MarkerType markedSubdomains = markSubdomains(description, localErrors, ageOfSubdomains);
            info << "enriching";
            if (!debugLogging)
              info << "... " << std::flush;
            else
              info << " on subdomains " << markedSubdomains.backend().transpose() << ":" << std::endl;
            timer.reset();
            // therefore, compute the local corrections
            computeLocalCorrections(*fixedGlobalSolver,
                                    *model,
                                    detailedOperator,
                                    detailedFunctional,
                                    description,
                                    localReducedBases,
                                    *mapper,
                                    reducedSolution,
                                    mu,
                                    markedSubdomains,
                                    tmpLocalVectorsA,
                                    tmpLocalVectorsB); // <- the local corrections will be in here
            // keep a copy of the old basis
            for (size_t ss = 0; ss < localReducedBases.size(); ++ss)
              localReducedBasesOld[ss]->backend() = localReducedBases[ss]->backend();
            oldBasisSize = mapper->size();
            // extend the local reduced bases (and obtain a new mapper)
            mapper = extendReducedBases(description,
                                        markedSubdomains,
                                        localScalarProducts,
                                        tmpLocalVectorsB,     // <- those will be changed, if gram schmidt
                                        localReducedBases);   // <- those will be extended
            if (mapper->size() > oldBasisSize) {
              // and project the operator and functional
              if (debugLogging)
                info << "  reducing operator:" << std::endl;
              reducedPair = createReducedOperatorFunctional(*fixedGlobalSolver,
                                                            detailedOperator,
                                                            detailedFunctional,
                                                            localReducedBases,
                                                            *mapper,
                                                            tmpLocalVectorA,
                                                            tmpGlobalVectorA,
                                                            tmpGlobalVectorB);
              reducedOperator = reducedPair.first;
              reducedFunctional = reducedPair.second;
              if (!debugLogging)
                info << " done (took " << timer.elapsed() << "sec)" << std::endl;
              // keep the old local errors
              localErrorsOld = localErrors;
              // compute the new local errors
              info << "  computing local errors... " << std::flush;
              timer.reset();
              reducedSolution = computeReducedSolution(*reducedOperator, *reducedFunctional, mu);
              if (visualize) {
                reconstructReducedSolution(localReducedBases,
                                           *mapper,
                                           reducedSolution,
                                           tmpLocalVectorsA); // <- reconstruction
                fixedGlobalSolver->visualize(tmpLocalVectorsA,
                                             "reducedSolution_iteration_" + Dune::Stuff::Common::toString(enrichmentIteration) + ".param_" + Dune::Stuff::Common::Parameter::reportWOwhitespace(mu),
                                             "reduced solution in iteration " + Dune::Stuff::Common::toString(enrichmentIteration) + " to param " + Dune::Stuff::Common::Parameter::report(mu),
                                             "",
                                             devnull);
              } // if (visualize)
              localErrors = computeLocalErrors(*fixedGlobalSolver,
                                               description,
                                               detailedOperator,
                                               detailedFunctional,
                                               localScalarProducts,
                                               localReducedBases,
                                               *mapper,
                                               reducedSolution,
                                               mu,
                                               tmpLocalVectorsA,
                                               tmpLocalVectorsB,
                                               tmpGlobalVectorA);
//              // now that we have the new local errors, lets throw away those local bases extensions that are crap!
//              debug << "old local erros were " << std::scientific << localErrorsOld.backend().transpose() << std::endl;
//              debug << "new local erros are  " << std::scientific << localErrors.backend().transpose() << std::endl;
//              size_t thrown_away = 0;
//              debug << "  rejecting local bassis extension on subdomains: " << std::flush;
//              for (size_t subdomain = 0; subdomain < msGrid->size(); ++subdomain) {
//                // compare local erros
//                if (localErrorsOld.get(subdomain) < localErrors.get(subdomain)) {
//                  debug << subdomain << " " << std::flush;
//                  // throw away new local basis
//                  localReducedBases[subdomain]->backend() = localReducedBasesOld[subdomain]->backend();
//                  ++thrown_away;
//                }
//              }
//              if (thrown_away) {
//                debug << "  updating reduced quantities:" << std::endl;
//                mapper = Dune::shared_ptr< MapperType >();
//                mapper->prepare();
//                for (size_t subdomain = 0; subdomain < msGrid->size(); ++subdomain)
//                  mapper->add(subdomain, localReducedBases[subdomain]->cols());
//                mapper->finalize();
//                reducedPair = createReducedOperatorFunctional(*fixedGlobalSolver,
//                                                              detailedOperator,
//                                                              detailedFunctional,
//                                                              localReducedBases,
//                                                              *mapper,
//                                                              tmpLocalVectorA,
//                                                              tmpGlobalVectorA,
//                                                              tmpGlobalVectorB);
//                reducedOperator = reducedPair.first;
//                reducedFunctional = reducedPair.second;
//                debug << "  updated local errors are: " << std::flush;
//                reducedSolution = computeReducedSolution(*reducedOperator, *reducedFunctional, mu);
//                localErrors = computeLocalErrors(*fixedGlobalSolver,
//                                                 description,
//                                                 detailedOperator,
//                                                 detailedFunctional,
//                                                 localScalarProducts,
//                                                 localReducedBases,
//                                                 *mapper,
//                                                 reducedSolution,
//                                                 mu,
//                                                 tmpLocalVectorsA,
//                                                 tmpLocalVectorsB,
//                                                 tmpGlobalVectorA);
//                debug << localErrors.backend().transpose() << std::endl;
//              } else
//                debug << "nowhere" << std::endl;

              localErrorsVectors.push_back(localErrors);
              totalBasisSizes.push_back(mapper->size());
              info << " done (took " << timer.elapsed() << "sec, max local error is "
                   << std::scientific << localErrors.backend().maxCoeff() << std::fixed << ")" << std::endl;
            } else {
              if (!debugLogging)
                info << " failed" << std::endl;
            } // if (mapper->size() > oldBasisSize)
            // increase iteration count
            ++enrichmentIteration;
          } // iteate until error is sufficient or we reach max iterations
          info << "final local errors are: "
               << std::scientific << localErrors.backend().transpose() << std::fixed << std::endl;

          const size_t oversamplingLayers = description.get< size_t >("grid.multiscale.provider.cube.oversamplingLayers");
          const size_t numSubdomainElements = msGrid->localGridPart(0)->indexSet().size(0);
          const size_t oversamplingRate = std::floor(100 * (oversamplingLayers / std::floor(std::sqrt(numSubdomainElements))));
          std::ofstream localErrorFile("localErrors_oversampling_"
                                       + Dune::Stuff::Common::toString(oversamplingLayers)
                                       + "_param_" + Dune::Stuff::Common::Parameter::reportWOwhitespace(mu));
          if (localErrorFile.is_open()) {
            for (size_t ii = 0; ii < localErrorsVectors.size(); ++ii) {
              localErrorFile << "% basisSize: " << totalBasisSizes[ii] << "\t || ";
              localErrorFile << "min: " << std::scientific << localErrorsVectors[ii].backend().minCoeff() << " | "
                             << "max: " << std::scientific << localErrorsVectors[ii].backend().maxCoeff() << " | "
                             << "mean: " << std::scientific << localErrorsVectors[ii].backend().mean() << " || ";
              for (int jj = 0; jj < localErrorsVectors[ii].size(); ++jj) {
                localErrorFile << std::scientific << localErrorsVectors[ii].get(jj) << " | ";
              }
              localErrorFile << std::endl;
            }
            localErrorFile << std::endl;
            localErrorFile << "% final local basis sizes: ";
            for (size_t ss = 0; ss < localReducedBases.size(); ++ss)
              localErrorFile << localReducedBases[ss]->cols() << " ";
            localErrorFile << std::endl;
            localErrorFile << "\\addplot[color=black,mark=square] coordinates {" << std::endl;
            for (size_t ii = 0; ii < localErrorsVectors.size(); ++ii)
              localErrorFile << "(" << totalBasisSizes[ii] << ", " << std::scientific << localErrorsVectors[ii].backend().maxCoeff() << ")" << std::endl;
            localErrorFile << "};\n\\addlegendentry{oversampling: " << std::fixed << oversamplingLayers << /*"\\%" <<*/ "}" << std::endl;
//            localErrorFile << "\\addplot[color=gray] coordinates {" << std::endl;
//            for (size_t ii = 0; ii < localErrorsVectors.size(); ++ii)
//              localErrorFile << "(" << totalBasisSizes[ii] << ", " << std::scientific << localErrorsVectors[ii].backend().minCoeff() << ")" << std::endl;
//            localErrorFile << "};" << std::endl;
//            localErrorFile << "\\addplot[color=gray] coordinates {" << std::endl;
//            for (size_t ii = 0; ii < localErrorsVectors.size(); ++ii)
//              localErrorFile << "(" << totalBasisSizes[ii] << ", " << std::scientific << localErrorsVectors[ii].backend().maxCoeff() << ")" << std::endl;
//            localErrorFile << "};" << std::endl;
          } // if (localErrorFile.is_open())

          msg = "  local offline phase finished  ";
          info << Dune::Stuff::Common::whitespaceify(msg, '=') << std::endl;
          info << msg << std::endl;
          info << Dune::Stuff::Common::whitespaceify(msg, '=') << std::endl;

        } // enrich, if max local error is too large
      } // loop over all parameters

    } else {
      nonparametericRun(description, model, msGrid, boundaryInfo);
    } // if (model->parametric())

  } catch(Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e.what() << std::endl;
  } catch(std::exception& e) {
    std::cerr << e.what() << std::endl;
  } catch( ... ) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  } // try

  // if we came that far we can as well be happy about it
  return 0;
} // main
