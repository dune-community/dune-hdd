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
//#include <dune/stuff/discretefunction/norm.hh>
//#include <dune/stuff/function/expression.hh>
#include <dune/stuff/la/container/eigen.hh>
#include <dune/stuff/la/algorithm/normalize.hh>
#include <dune/stuff/la/container/separable.hh>
#include <dune/stuff/la/solver.hh>
#include <dune/stuff/discretefunction/projection/lagrange.hh>
#include <dune/stuff/la/algorithm/gramschmidt.hh>

#include <dune/detailed/solvers/stationary/linear/elliptic/model.hh>
#include <dune/detailed/solvers/stationary/linear/elliptic/ms/semicg/detailed-discretizations.hh>
//#include <dune/detailed/solvers/stationary/linear/elliptic/cg/detailed-discretizations.hh>

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
  file << "model = model.stationary.linear.elliptic.default" << std::endl;
  file << "        model.stationary.linear.elliptic.thermalblock" << std::endl;
  file << "        model.stationary.linear.elliptic.parametric.separable.thermalblock" << std::endl;
//  file << "exact_solution.order = 2" << std::endl;
//  file << "exact_solution.variable = x" << std::endl;
//  file << "exact_solution.expression.0 = -0.5*x[0]*x[0] + 0.5*x[0]" << std::endl;
  file << "visualize = true" << std::endl;
  file << "[grid.multiscale.provider.cube]" << std::endl;
  file << "lowerLeft = [0.0; 0.0; 0.0]" << std::endl;
  file << "upperRight = [1.0; 1.0; 1.0]" << std::endl;
  file << "numElements = [32; 32; 32]" << std::endl;
  file << "boundaryId = 7 # a cube from the factory gets the boundary ids 1 to 4 ind 2d and 1 to 6 in 3d (hopefully)" << std::endl;
  file << "partitions = [2; 2; 2]" << std::endl;
  file << "[detailed.solvers.stationary.linear.elliptic.ms.semicg.detailed-discretizations]" << std::endl;
  file << "penaltyFactor = 100.0" << std::endl;
  file << "linearsolver.type = bicgstab.ilut" << std::endl;
  file << "                    bicgstab.diagonal" << std::endl;
  file << "                    bicgstab" << std::endl;
  file << "                    cg.diagonal" << std::endl;
  file << "                    cg" << std::endl;
  file << "linearsolver.maxIter   = 5000"  << std::endl;
  file << "linearsolver.precision = 1e-12"  << std::endl;
  file << "[model.stationary.linear.elliptic.default]" << std::endl;
  file << "diffusion.order      = 0"  << std::endl;
  file << "diffusion.variable   = x" << std::endl;
  file << "diffusion.expression = [1.0; 1.0; 1.0]"  << std::endl;
  file << "diffusion.name       = constant diffusion"  << std::endl;
  file << "force.order      = 0"  << std::endl;
  file << "force.variable   = x" << std::endl;
  file << "force.expression = [1.0; 1.0; 1.0]"  << std::endl;
  file << "force.name       = constant force"  << std::endl;
  file << "dirichlet.order      = 1"  << std::endl;
  file << "dirichlet.variable   = x" << std::endl;
  file << "dirichlet.expression = [0.0; 0.0; 0.0]"  << std::endl;
  file << "dirichlet.name       = trivial dirichlet"  << std::endl;
  file << "neumann.order      = 0"  << std::endl;
  file << "neumann.variable   = x" << std::endl;
  file << "neumann.expression = [0.0; 0.0; 0.0]"  << std::endl;
  file << "neumann.name       = trivial neumann"  << std::endl;
  file << "[model.stationary.linear.elliptic.thermalblock]" << std::endl;
  file << "diffusion.order = 0"  << std::endl;
  file << "diffusion.lowerLeft   = [0.0; 0.0; 0.0] # this should be a bounding box of the above selected grid!" << std::endl;
  file << "diffusion.upperRight  = [1.0; 1.0; 1.0] # this should be a bounding box of the above selected grid!" << std::endl;
  file << "diffusion.numElements = [2; 2; 2]"  << std::endl;
  file << "diffusion.values      = [1.0; 10.0; 3.0; 2.1]"  << std::endl;
  file << "diffusion.name        = checkerboard diffusion"  << std::endl;
  file << "force.name       = constant force"  << std::endl;
  file << "force.order      = 0"  << std::endl;
  file << "force.variable   = x" << std::endl;
  file << "force.expression = [1.0; 1.0; 1.0]"  << std::endl;
  file << "dirichlet.name       = trivial dirichlet"  << std::endl;
  file << "dirichlet.order      = 0"  << std::endl;
  file << "dirichlet.variable   = x" << std::endl;
  file << "dirichlet.expression = [0.0; 0.0; 0.0]"  << std::endl;
  file << "neumann.name       = trivial neumann"  << std::endl;
  file << "neumann.order      = 0"  << std::endl;
  file << "neumann.variable   = x" << std::endl;
  file << "neumann.expression = [0.0; 0.0; 0.0]"  << std::endl;
//  file << "[model.stationary.linear.elliptic.parametric.separable.default]" << std::endl;
//  file << "diffusion.type = function.parametric.separable.default" << std::endl;
//  file << "diffusion.name = diffusion" << std::endl;
//  file << "diffusion.order = 0" << std::endl;
//  file << "diffusion.component.0 = x[0]" << std::endl;
//  file << "diffusion.component.1 = 2*x[0]" << std::endl;
//  file << "diffusion.coefficient.0 = mu[0]" << std::endl;
//  file << "diffusion.coefficient.1 = mu[1]" << std::endl;
//  file << "diffusion.paramSize = 2" << std::endl;
//  file << "diffusion.paramExplanation = [first diffusion param; second diffusion param]" << std::endl;
//  file << "diffusion.paramMin = [0.0; 0.0]" << std::endl;
//  file << "diffusion.paramMax = [1.0; 1.0]" << std::endl;
//  file << "force.type = function.parametric.separable.default" << std::endl;
//  file << "force.name = force" << std::endl;
//  file << "force.order = 0" << std::endl;
//  file << "force.component.0 = x[0]" << std::endl;
//  file << "force.component.1 = 2*x[0]" << std::endl;
//  file << "force.coefficient.0 = mu[0]" << std::endl;
//  file << "force.coefficient.1 = mu[1]" << std::endl;
//  file << "force.paramSize = 2" << std::endl;
//  file << "force.paramExplanation = [first force param; second force param]" << std::endl;
//  file << "force.paramMin = [0.0; 0.0]" << std::endl;
//  file << "force.paramMax = [1.0; 1.0]" << std::endl;
//  file << "dirichlet.type = function.parametric.separable.default" << std::endl;
//  file << "dirichlet.name = neumann" << std::endl;
//  file << "dirichlet.order = 0" << std::endl;
//  file << "dirichlet.component.0 = x[0]" << std::endl;
//  file << "dirichlet.component.1 = 2*x[0]" << std::endl;
//  file << "dirichlet.coefficient.0 = mu[0]" << std::endl;
//  file << "dirichlet.coefficient.1 = mu[1]" << std::endl;
//  file << "dirichlet.paramSize = 2" << std::endl;
//  file << "dirichlet.paramExplanation = [first dirichlet param; second dirichlet param]" << std::endl;
//  file << "dirichlet.paramMin = [0.0; 0.0]" << std::endl;
//  file << "dirichlet.paramMax = [1.0; 1.0]" << std::endl;
//  file << "neumann.type = function.parametric.separable.default" << std::endl;
//  file << "neumann.name = neumann" << std::endl;
//  file << "neumann.order = 0" << std::endl;
//  file << "neumann.component.0 = x[0]" << std::endl;
//  file << "neumann.component.1 = 2*x[0]" << std::endl;
//  file << "neumann.coefficient.0 = mu[0]" << std::endl;
//  file << "neumann.coefficient.1 = mu[1]" << std::endl;
//  file << "neumann.paramSize = 2" << std::endl;
//  file << "neumann.paramExplanation = [first neumann param; second neumann param]" << std::endl;
//  file << "neumann.paramMin = [0.0; 0.0]" << std::endl;
//  file << "neumann.paramMax = [1.0; 1.0]" << std::endl;
  file << "[model.stationary.linear.elliptic.parametric.separable.thermalblock]" << std::endl;
  file << "diffusion.lowerLeft   = [0.0; 0.0; 0.0] # this should be a bounding box of the above selected grid!" << std::endl;
  file << "diffusion.upperRight  = [1.0; 1.0; 1.0] # this should be a bounding box of the above selected grid!" << std::endl;
  file << "diffusion.numElements = [2; 2; 2]"  << std::endl;
  file << "diffusion.paramMin    = [0.1; 0.1; 0.1; 0.1]" << std::endl;
  file << "diffusion.paramMax    = [10.0; 10.0; 10.0; 10.0]" << std::endl;
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
  file << "debug = true" << std::endl;
  file << "file  = false" << std::endl;
  file << "[parameter]" << std::endl;
  file << "fixed = [1.0; 1.0; 1.0; 1.0]" << std::endl;
  file << "train = [0.1; 0.1; 1.0; 1.0]" << std::endl;
  file << "test  = [1.0; 1.0; 0.1; 0.1]" << std::endl;
  file << "[enrichment]" << std::endl;
  file << "maxLocalError = 1e-5" << std::endl;
  file << "maxIterations = 5" << std::endl;
  file << "gramSchmidtTolerance = 1e-10" << std::endl;
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


Dune::shared_ptr< SolverType > initSolver(const DescriptionType& description)
{
  // logger
  const std::string filename = description.get< std::string >("filename");
  const bool visualize = description.get< bool >(id() + ".visualize", true);
  const bool debugLogging = description.get< bool >("logging.debug", false);
  Dune::Stuff::Common::LogStream& info    = Dune::Stuff::Common::Logger().info();
  Dune::Stuff::Common::LogStream& debug   = Dune::Stuff::Common::Logger().debug();
//  Dune::Stuff::Common::LogStream& devnull = Dune::Stuff::Common::Logger().devnull();

  // timer
  Dune::Timer timer;

  info << "setting up grid: " << std::endl;
  debug.suspend();
  const GridProviderType gridProvider = GridProviderType::createFromDescription(description);
  const Dune::shared_ptr< const MsGridType > msGrid = gridProvider.msGrid();
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
  }

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
  const Dune::shared_ptr< const BoundaryInfoType >
      boundaryInfo(Dune::Stuff::Grid::BoundaryInfo::create< GlobalGridViewType >(boundaryInfoType, description));
  info << "done (took " << timer.elapsed() << " sec)" << std::endl;

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
  const Dune::shared_ptr< const ModelType >
      model(Stationary::Linear::Elliptic::Model::create<  DomainFieldType, dimDomain,
                                                          RangeFieldType, dimRange >(modelType, description));
  info << "done (took " << timer.elapsed() << " sec)" << std::endl;
  if (visualize) {
    info << "visualizing model... " << std::flush;
    timer.reset();
    model->visualize(msGrid->globalGridPart()->gridView(), /*filename +*/ "model");
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;
  }

  info << "setting up solver";
  if (!debugLogging)
    info << "... " << std::flush;
  else
    info << " '" << SolverType::id() << "':" << std::endl;
  timer.reset();
  const DescriptionType& discretizationDescription = description.sub(SolverType::id());
  Dune::shared_ptr< SolverType > solver = Dune::make_shared< SolverType >(model,
                                                                          msGrid,
                                                                          boundaryInfo,
                                                                          discretizationDescription.get< RangeFieldType >("penaltyFactor"));
  solver->init("  ", debug);
  if (!debugLogging)
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

  return solver;
} // ... initSolver(...)


void nonparametericRun(const SolverType& solver, const DescriptionType& description)
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


std::vector< Dune::shared_ptr< ScalarProductType > > computeLocalScalarProducts(const SolverType& solver,
                                                                                const DescriptionType& description)
{
  // logger
//  Dune::Stuff::Common::LogStream& info    = Dune::Stuff::Common::Logger().info();
  Dune::Stuff::Common::LogStream& debug   = Dune::Stuff::Common::Logger().debug();
//  Dune::Stuff::Common::LogStream& devnull = Dune::Stuff::Common::Logger().devnull();
  Dune::Timer timer;
  // paramters
  const DescriptionType& parameterDescription = description.sub("parameter");
  const ParamType fixedParam = parameterDescription.getDynVector< ParamFieldType >(   "fixed", solver.model()->paramSize());
//  const ParamType trainingParam = parameterDescription.getDynVector< ParamFieldType >("train", model->paramSize());
//  const ParamType testParam = parameterDescription.getDynVector< ParamFieldType >(    "test",  model->paramSize());
  assert(fixedParam.size()    == solver.model()->paramSize());
//  assert(trainingParam.size() == model->paramSize());
//  assert(testParam.size()     == model->paramSize());
  debug << "  for parameter " << Dune::Stuff::Common::Parameter::report(fixedParam) << "... " << std::flush;
  const std::vector< Dune::shared_ptr< ScalarProductType > > localScalarProducts = solver.systemMatrices(fixedParam);
  debug << "done (took " << timer.elapsed() << " sec)" << std::endl;
  return localScalarProducts;
} // ... computeLocalScalarProducts(...)


void computeGlobalSnapshot(const SolverType& solver, const DescriptionType& description,
                           const ParamType& parameter, std::vector< Dune::shared_ptr< VectorType > >& targetVectors)
{
  // logger
  Dune::Stuff::Common::LogStream& debug = Dune::Stuff::Common::Logger().debug();
  const bool visualize = description.get< bool >(id() + ".visualize", true);
  // linear solver
  const DescriptionType& discretizationDescription = description.sub(SolverType::id());
  const std::string    linearSolverType      = discretizationDescription.get< std::string >(   "linearsolver.type",
                                                                                               "bicgstab.ilut");
  const size_t         linearSolverMaxIter   = discretizationDescription.get< size_t >(        "linearsolver.maxIter",
                                                                                               5000u);
  const RangeFieldType linearSolverPrecision = discretizationDescription.get< RangeFieldType >("linearsolver.precision",
                                                                                               1e-12);
  // solve
  solver.solve(targetVectors,
               parameter,
               linearSolverType,
               linearSolverMaxIter,
               linearSolverPrecision,
               "  ",
               debug);
  if (visualize) {
    solver.visualize(targetVectors,
                     /*filename +*/ "snapshot.param_" + Dune::Stuff::Common::Parameter::reportWOwhitespace(parameter),
                     /*id() +*/ "snapshot for parameter " + Dune::Stuff::Common::Parameter::report(parameter),
                     "  ",
                     debug);
  }
} // ... computeGlobalSnapshot(...)


void initReducedBasis(const std::vector< Dune::shared_ptr< ScalarProductType > >& localScalarProducts,
                      std::vector< Dune::shared_ptr< VectorType > >& localSnapshots,
                      std::vector< Dune::shared_ptr< ReducedBasisType > >& localReducedBases,
                      MapperType& mapper)
{
  // sanity checks
  assert(localSnapshots.size() == localScalarProducts.size());
  assert(localSnapshots.size() == localReducedBases.size());
  // prepare
  mapper.prepare();
  // loop over all subdomains
  for (size_t ss = 0; ss < localReducedBases.size(); ++ss) {
    // get local scalar product
    const ScalarProductType& localScalarProduct = *(localScalarProducts[ss]);
    // get local snapshot
    VectorType& localSnapshot = *(localSnapshots[ss]);
    assert(localScalarProduct.rows() == localSnapshot.size());
    assert(localScalarProduct.cols() == localSnapshot.size());
    // normalize local snapshot (inplace)
    Dune::Stuff::LA::Algorithm::normalize(localScalarProduct, localSnapshot);
    // set local reduced basis
    ReducedBasisType& localReducedBasis = *(localReducedBases[ss]);
    localReducedBasis.backend().col(0) = localSnapshot.backend();
    mapper.add(ss, 1);
  } // loop over all subdomains
  // finalize
  mapper.finalize();
} // ... initReducedBasis(...)


std::pair<  Dune::shared_ptr< ReducedOperatorType >,
            Dune::shared_ptr< ReducedFunctionalType > >
createReducedOperatorFunctional(const SolverType& solver,
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
  // get detailed operator
  const DetailedOperatorType& detailedOperator = *(solver.systemMatrix("diffusion"));
  debug << detailedOperator.components()[0]->rows() << "x" << detailedOperator.components()[0]->cols()
        << " to " << mapper.size() << "x" << mapper.size() << "... " << std::flush;
  timer.reset();
  // get detailed functional
  const DetailedFunctionalType& detailedFunctional = *(solver.systemVector("force"));
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
  const MsGridType& msGrid = *(solver.msGrid());
  // loop over all subdomains
  for (size_t subdomain = 0; subdomain < msGrid.size(); ++subdomain) {
    // get local reduced basis
    const ReducedBasisType& reducedBasisSubdomain = *(localReducedBases[subdomain]);
    // loop over all local basis functions
    for (int ii = 0; ii < reducedBasisSubdomain.cols(); ++ii) {
      // convert the iith basis function to a global one
      tmpLocalVector->backend() = reducedBasisSubdomain.backend().col(ii);
      solver.globalizeVector(tmpLocalVector, subdomain, tmpGlobalVectorOne);
      // compute global index
      const size_t iiGlobal = mapper.toGlobal(subdomain, ii);
      // project the detailed functional wrt iith local basis function
      reducedFunctional->set(iiGlobal, tmpGlobalVectorOne->backend().transpose() * detailedFunctional.backend());
      // loop over all local basis functions
      for (int jj = 0; jj < reducedBasisSubdomain.cols(); ++jj) {
        // convert the jjth basis function to a global one
        tmpLocalVector->backend() = reducedBasisSubdomain.backend().col(jj);
        solver.globalizeVector(tmpLocalVector, subdomain, tmpGlobalVectorTwo);
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
        solver.globalizeVector(tmpLocalVector, subdomain, tmpGlobalVectorOne);
        // compute global index
        const size_t iiGlobal = mapper.toGlobal(subdomain, ii);
        // loop over all the local basis functions of the neighbor
        for (int jj = 0; jj < reducedBasisNeighbor.cols(); ++jj) {
          // convert the jjth basis function to a global one
          tmpLocalVector->backend() = reducedBasisNeighbor.backend().col(jj);
          solver.globalizeVector(tmpLocalVector, neighbor, tmpGlobalVectorTwo);
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


void computeReducedSolution(const SolverType& solver,
                            const DescriptionType& description,
                            const ReducedOperatorType& reducedOperator,
                            const ReducedFunctionalType& reducedFunctional,
                            const ParamType& mu,
                            VectorType& reducedSolution)
{
  // logger
  Dune::Stuff::Common::LogStream& debug   = Dune::Stuff::Common::Logger().debug();
  const bool visualize = description.get< bool >(id() + ".visualize", false);
  Dune::Timer timer;
  if (visualize)
    debug << "  computing reduced solution... " << std::flush;
  timer.reset();
  const auto reducedSystemMatrix = reducedOperator.fix(solver.model()->mapParam(mu, "diffusion"));
  reducedSolution.backend() = reducedSystemMatrix->backend().colPivHouseholderQr().solve(reducedFunctional.backend());
  if (visualize) {
    debug << " done (took " << std::scientific << timer.elapsed() << std::fixed << "sec)" << std::endl;
  }
} // ... computeReducedSolution(...)


void reconstructReducedSolution(const SolverType& solver,
                                const std::vector< Dune::shared_ptr< ReducedBasisType > >& localReducedBases,
                                const MapperType& mapper,
                                const VectorType& reducedSolution,
                                std::vector< Dune::shared_ptr< VectorType > >& localVectors)
{
  // loop over all subdomains
  for (size_t ss = 0; ss < solver.msGrid()->size(); ++ss) {
    // get the local reduced basis
    const ReducedBasisType& localReducedBasis = *(localReducedBases[ss]);
    // get the local vector
    VectorType& localVector = *(localVectors[ss]);
    // clear it
    localVector.backend().setZero();
    // loop over all local basis functions
    for (int ii = 0; ii < localReducedBasis.cols(); ++ii) {
      // and reconstruct
      localVector.backend() += localReducedBasis.backend().col(ii) * reducedSolution.get(mapper.toGlobal(ss, ii));
    } // loop over all local basis functions
  } // loop over all subdomains
} // ... reconstructReducedSolution(...)


Dune::shared_ptr< VectorType > computeLocalErrors(const SolverType& solver,
                                                  const DescriptionType& description,
                                                  const std::vector< Dune::shared_ptr< ScalarProductType > >& localScalarProducts,
                                                  const std::vector< Dune::shared_ptr< ReducedBasisType > >& localReducedBases,
                                                  const MapperType& mapper,
                                                  const VectorType& reducedSolution,
                                                  const ParamType& mu,
                                                  std::vector< Dune::shared_ptr< VectorType > >& tmpLocalVectorsOne,
                                                  std::vector< Dune::shared_ptr< VectorType > >& tmpLocalVectorsTwo)
{
  // logger
  Dune::Stuff::Common::LogStream& debug   = Dune::Stuff::Common::Logger().debug();
  Dune::Stuff::Common::LogStream& devnull = Dune::Stuff::Common::Logger().devnull();
  const bool visualize = description.get< bool >(id() + ".visualize", false);
  Dune::Timer timer;
  // linear solver
  const DescriptionType& discretizationDescription = description.sub(SolverType::id());
  const std::string    linearSolverType      = discretizationDescription.get< std::string >(   "linearsolver.type",
                                                                                               "bicgstab.ilut");
  const size_t         linearSolverMaxIter   = discretizationDescription.get< size_t >(        "linearsolver.maxIter",
                                                                                               5000u);
  const RangeFieldType linearSolverPrecision = discretizationDescription.get< RangeFieldType >("linearsolver.precision",
                                                                                               1e-12);
  debug << "  computing detailed solution... " << std::flush;
  timer.reset();
  solver.solve(tmpLocalVectorsOne,
               mu,
               linearSolverType,
               linearSolverMaxIter,
               linearSolverPrecision,
               "", devnull);
  debug << "done (took " << timer.elapsed() << " sec)" << std::endl;

  debug << "  reconstructing reduced solution... " << std::flush;
  timer.reset();
  reconstructReducedSolution(solver, localReducedBases, mapper, reducedSolution, tmpLocalVectorsTwo);
  debug << "done (took " << timer.elapsed() << " sec)" << std::endl;

  debug << "  computing difference... " << std::flush;
  timer.reset();
  for (size_t ss = 0; ss < solver.msGrid()->size(); ++ss)
    tmpLocalVectorsTwo[ss]->backend() -= tmpLocalVectorsOne[ss]->backend();
  debug << "done (took " << timer.elapsed() << " sec)" << std::endl;

  // visualize error (if wished)
  if (visualize) {
    solver.visualize(tmpLocalVectorsTwo,
                     "difference.param_" + Dune::Stuff::Common::Parameter::reportWOwhitespace(mu),
                     "difference for parameter " + Dune::Stuff::Common::Parameter::report(mu),
                     "  ",
                     debug);
  } // visualize error (if wished)

  debug << "  computing local errors... " << std::flush;
  timer.reset();
  Dune::shared_ptr< VectorType > localErrors = Dune::make_shared< VectorType >(solver.msGrid()->size());
  // loop over all subdomains
  for (size_t ss = 0; ss < solver.msGrid()->size(); ++ss) {
    // get local scalar product
    const ScalarProductType& localScalarProduct = *(localScalarProducts[ss]);
    // get local reference solution
    const VectorType& localReference = *(tmpLocalVectorsOne[ss]);
    // get local difference
    const VectorType& localDifference = *(tmpLocalVectorsTwo[ss]);
    // compute local relative error
    localErrors->set(ss,
                     std::sqrt(localDifference.backend().transpose()
                               * localScalarProduct.backend()
                               * localDifference.backend())
                     / std::sqrt(localReference.backend().transpose()
                                 * localScalarProduct.backend()
                                 * localReference.backend()));
  } // loop over all subdomains
  debug << "done (took " << timer.elapsed() << " sec)" << std::endl;

  return localErrors;
} // ... computeLocalErrors(...)


void computeLocalCorrections(const SolverType& globalSolver,
                             const DescriptionType& description,
                             const std::vector< Dune::shared_ptr< ReducedBasisType > >& localReducedBases,
                             const MapperType& mapper,
                             const VectorType& reducedSolution,
                             const ParamType& mu,
                             const MarkerType& markedSubdomains,
                             std::vector< Dune::shared_ptr< VectorType > >& localVectorsOne,
                             std::vector< Dune::shared_ptr< VectorType > >& localVectorsTwo)
{
  // linear solver
  const DescriptionType& discretizationDescription = description.sub(SolverType::id());
  const std::string    linearSolverType      = discretizationDescription.get< std::string >(   "linearsolver.type",
                                                                                               "bicgstab.ilut");
  const size_t         linearSolverMaxIter   = discretizationDescription.get< size_t >(        "linearsolver.maxIter",
                                                                                               5000u);
  const RangeFieldType linearSolverPrecision = discretizationDescription.get< RangeFieldType >("linearsolver.precision",
                                                                                               1e-12);
  // logger
  Dune::Stuff::Common::LogStream& debug   = Dune::Stuff::Common::Logger().debug();
  Dune::Stuff::Common::LogStream& devnull = Dune::Stuff::Common::Logger().devnull();
  const bool visualize = description.get< bool >(id() + ".visualize", true);
  Dune::Timer timer;
  debug << "  reconstructing reduced solution... " << std::flush;
  reconstructReducedSolution(globalSolver, localReducedBases, mapper, reducedSolution, localVectorsTwo);
  const Dune::shared_ptr< const GlobalDiscreteFunctionType >
      reconstructedSolutionFunction = globalSolver.createDiscreteFunction(localVectorsTwo);
  debug << "done (took " << timer.elapsed() << "sec)" << std::endl;
  debug << "  solving for local corrections... " << std::flush;
  timer.reset();
  // get the multiscale grid
  const MsGridType& msGrid = *(globalSolver.msGrid());
  // create a fixed model for the current parameter
  const Dune::shared_ptr< const ModelType > fixedModel = globalSolver.model()->fix(mu);
  // loop over all marked subdomains
  for (int ii = 0; ii < markedSubdomains.size(); ++ii) {
    const size_t subdomain = markedSubdomains.get(ii);
    // get the local grid part
    const Dune::shared_ptr< const LocalGridPartType > localOversampledGridPart = msGrid.localGridPart(subdomain, true);
    // initialize a local solver on the local oversampled gripart
    const Dune::shared_ptr< const LocalBoundaryInfoType > localBoundaryInfo(
          Dune::Stuff::Grid::BoundaryInfo::create< typename MsGridType::LocalGridViewType >("stuff.grid.boundaryinfo.alldirichlet"));
    LocalSolverType localOversampledSolver(localOversampledGridPart, localBoundaryInfo, fixedModel);
    localOversampledSolver.init("", devnull);
    // get the local solver on the non oversampled gridpart
    const LocalSolverType& localSolver = *(globalSolver.localSolver(subdomain));
    // create a local function to hold the reconstructed solution
    Dune::shared_ptr< LocalDiscreteFunctionType >
        localReconstructedSolutionFunction = localOversampledSolver.createAnsatzFunction();
    // and project the current reconstructed solution onto the local oversampled gripart
    Dune::Stuff::DiscreteFunction::Projection::Lagrange::project(*reconstructedSolutionFunction,
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
    Dune::shared_ptr< VectorType > localOversampledSolutionVector = localOversampledSolver.createAnsatzVector();
    localOversampledSolver.solve(localOversampledSolutionVector,
                      linearSolverType,
                      linearSolverMaxIter,
                      linearSolverPrecision,
                      "",
                      devnull);
    // and restrict the oversampled solution to the local (not oversampled) subdomain
    Dune::shared_ptr< LocalDiscreteFunctionType >
        localOversampledSolutionFunction = localOversampledSolver.createAnsatzFunction(localOversampledSolutionVector);
    Dune::shared_ptr< LocalDiscreteFunctionType > localSolutionFunction
        = localSolver.createAnsatzFunction();
    Dune::Stuff::DiscreteFunction::Projection::Lagrange::project(*localOversampledSolutionFunction,
                                                                 *localSolutionFunction);
    // visualize if wished
    if (visualize) {
      localSolver.visualizeAnsatzVector(localSolutionFunction->vector(),
                                        "local_correction.param_" + Dune::Stuff::Common::Parameter::reportWOwhitespace(mu),
                                        "local correction for parameter " + Dune::Stuff::Common::Parameter::report(mu),
                                        "",
                                        devnull);
    } // visualize if wished
    // now just copy the local solution into the return argument
    localVectorsOne[subdomain]->backend() = localSolutionFunction->vector()->backend();
  } // loop over all marked subdomains
  debug << "done (took " << timer.elapsed() << "sec)" << std::endl;
} // ... computeLocalCorrections(...)


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
  const double gramSchmidtTolerance = description.get< double >("enrichment.gramSchmidtTolerance", 1e-10);
  debug << "  extending reduced bases:" << std::endl;
  // loop over all marked subdomains
  for (int ii = 0; ii < markedSubdomains.size(); ++ii) {
    const size_t subdomain = markedSubdomains.get(ii);
    debug << "    on subdomain " << subdomain << "... " << std::flush;
    timer.reset();
    // get the local scalar product
    const ScalarProductType& localScalarProduct = *(localScalarProducts[subdomain]);
    // get the local snapshot
    VectorType& localSnapshot = *(localSnapshots[subdomain]);
    // get the local reduced basis
    ReducedBasisType& localReducedBasis = *(localReducedBases[subdomain]);
    // apply a gram schmidt to the local snapshot (inplace) wrt to the existing reduced basis
    const bool success = Dune::Stuff::LA::Algorithm::gramSchmidt(localScalarProduct,
                                                                 localReducedBasis,
                                                                 localSnapshot,
                                                                 gramSchmidtTolerance);
    if (success) {
      const size_t oldSize = localReducedBasis.cols();
      // resize the local reduced basis matrix by keeping the old columns
      localReducedBasis.backend().conservativeResize(::Eigen::NoChange, oldSize + 1);
      // and set the new column to the orthogonalized local snapshot
      localReducedBasis.backend().col(oldSize) = localSnapshot.backend();
      debug << "done (from " << oldSize << " to " << localReducedBasis.cols() << ")" << std::endl;
    } else
      debug << "failed" << std::endl;
  } // loop over all marked subdomains
  // also creeate a new mapper
  Dune::shared_ptr< MapperType > mapper = Dune::make_shared< MapperType >();
  mapper->prepare();
  for (size_t ss = 0; ss < localReducedBases.size(); ++ss)
    mapper->add(ss, localReducedBases[ss]->cols());
  mapper->finalize();
  return mapper;
} // ... extendReducedBases(...)


int main(int argc, char** argv)
{
  try {
    // description
    const DescriptionType description = initDescription(argc, argv);
    // logger
    Dune::Stuff::Common::LogStream& info    = Dune::Stuff::Common::Logger().info();
    Dune::Stuff::Common::LogStream& debug   = Dune::Stuff::Common::Logger().debug();
    Dune::Stuff::Common::LogStream& devnull = Dune::Stuff::Common::Logger().devnull();
    const bool debugLogging = description.get< bool >("logging.debug", false);
    const bool visualize = description.get< bool >(id() + ".visualize", false);
    // report
    std::string msg = "  starting offline phase  ";
    info << Dune::Stuff::Common::whitespaceify(msg, '=') << std::endl;
    info << msg << std::endl;
    info << Dune::Stuff::Common::whitespaceify(msg, '=') << std::endl;
    Dune::Timer offlineTimer;
    Dune::Timer timer;
    // detailed solver
    const Dune::shared_ptr< const SolverType > detailedSolver = initSolver(description);
    // do something
    if (detailedSolver->parametric()) {
      // generate common storage
      info << "initializing common storage... " << std::flush;
      timer.reset();
      Dune::shared_ptr< VectorType > localVector = detailedSolver->localSolver(0)->createAnsatzVector();
      std::vector< Dune::shared_ptr< VectorType > > localVectors = detailedSolver->createVector();
      std::vector< Dune::shared_ptr< VectorType > > localVectorsTwo = detailedSolver->createVector();
      std::vector< Dune::shared_ptr< ReducedBasisType > > localReducedBases(detailedSolver->msGrid()->size());
      VectorType reducedSolution(detailedSolver->msGrid()->size());
      for (size_t ss = 0; ss < detailedSolver->msGrid()->size(); ++ss)
        localReducedBases[ss] = Dune::make_shared< ReducedBasisType >(detailedSolver->localSolver(ss)->ansatzSpace()->map().size(), 1);
      Dune::shared_ptr< VectorType > globalVectorOne = Dune::make_shared< VectorType >(detailedSolver->ansatzMapper().size());
      Dune::shared_ptr< VectorType > globalVectorTwo = Dune::make_shared< VectorType >(detailedSolver->ansatzMapper().size());
      info << " done (took " << timer.elapsed() << "sec)" << std::endl;

      // compute the local scalar products
      info << "computing local scalar products";
      if (!debugLogging)
        info << "... " << std::flush;
      else
        info << ":" << std::endl;
      timer.reset();
      std::vector< Dune::shared_ptr< ScalarProductType > >
          localScalarProducts = computeLocalScalarProducts(*detailedSolver,
                                                           description);
      if (!debugLogging)
        info << " done (took " << timer.elapsed() << "sec)" << std::endl;

      // get the paramters
      const DescriptionType& parameterDescription = description.sub("parameter");
      const ParamType trainingParam
          = parameterDescription.getDynVector< ParamFieldType >("train", detailedSolver->model()->paramSize());
      const ParamType testParam
          = parameterDescription.getDynVector< ParamFieldType >("test",  detailedSolver->model()->paramSize());
      const std::vector< ParamType > parameters = {trainingParam, testParam};

      // compute the global training snapshot
      info << "computing training snapshot";
      if (!debugLogging)
        info << "... " << std::flush;
      else
        info << " for parameter " << Dune::Stuff::Common::Parameter::report(trainingParam) << ":" << std::endl;
      timer.reset();
      computeGlobalSnapshot(*detailedSolver, description, trainingParam, localVectors);
      if (!debugLogging)
        info << " done (took " << timer.elapsed() << "sec)" << std::endl;

      // initiaize the reduced basis
      info << "initializing local reduced bases... " << std::flush;
      timer.reset();
      Dune::shared_ptr< MapperType > mapper = Dune::make_shared< MapperType >();
      initReducedBasis(localScalarProducts, localVectors, localReducedBases, *mapper);
      info << " done (took " << timer.elapsed() << "sec)" << std::endl;

      // project the detailed operator and functional onto the reduced space
      info << "reducing operator";
      if (!debugLogging)
        info << "... " << std::flush;
      else
        info << ":" << std::endl;
      auto reducedPair = createReducedOperatorFunctional(*detailedSolver,
                                                         localReducedBases,
                                                         *mapper,
                                                         localVector,
                                                         globalVectorOne,
                                                         globalVectorTwo);
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
      // * then loop over all parameters
      for (size_t pp = 0; pp < parameters.size(); ++pp) {
        const ParamType& mu = parameters[pp];
        // compute a reduced solution
        info << "solving for parameter " << Dune::Stuff::Common::Parameter::report(mu);
        if (!visualize)
          info << "... " << std::flush;
        else
          info << ":" << std::endl;
        timer.reset();
        computeReducedSolution(*detailedSolver, description, *reducedOperator, *reducedFunctional, mu, reducedSolution);
        if (!visualize)
          info << " done (took " << std::scientific << timer.elapsed() << std::fixed << "sec)" << std::endl;
        // reconstruct and visualize if wished
        if (visualize) {
          info << "  reconstructing solution... " << std::flush;
          timer.reset();
          reconstructReducedSolution(*detailedSolver, localReducedBases, *mapper, reducedSolution, localVectors);
          info << " done (took " << timer.elapsed() << "sec)" << std::endl;
          if (debugLogging) {
            detailedSolver->visualize(localVectors,
                                     "reconstructedSolution.param_" + Dune::Stuff::Common::Parameter::reportWOwhitespace(mu),
                                     "reconstructedSolution to param " + Dune::Stuff::Common::Parameter::report(mu),
                                     "  ",
                                     debug);
          } else {
            info << "visualizing reconstructed solution..." << std::flush;
            timer.reset();
            detailedSolver->visualize(localVectors,
                                     "reconstructedSolution.param_" + Dune::Stuff::Common::Parameter::reportWOwhitespace(mu),
                                     "reconstructedSolution to param " + Dune::Stuff::Common::Parameter::report(mu),
                                     "",
                                     devnull);
            info << " done (took " << timer.elapsed() << "sec)" << std::endl;
          } // if (debuglogging)
        } // if (visualize)
        // estimate error
        info << "estimating error";
        if (!debugLogging)
          info << "... " << std::flush;
        else
          info << ":" << std::endl;
        timer.reset();
        Dune::shared_ptr< VectorType > localErrors = computeLocalErrors(*detailedSolver,
                                                                        description,
                                                                        localScalarProducts,
                                                                        localReducedBases,
                                                                        *mapper,
                                                                        reducedSolution,
                                                                        mu,
                                                                        localVectors,
                                                                        localVectorsTwo);
        if (!debugLogging)
          info << " done (took " << timer.elapsed() << "sec)" << std::endl;
        // enrich, if max local error is too large
        if (localErrors->backend().maxCoeff() > enrichmentMaxLocalError) {
          msg = "  starting local offline phase  ";
          info << Dune::Stuff::Common::whitespaceify(msg, '=') << std::endl;
          info << msg << std::endl;
          info << Dune::Stuff::Common::whitespaceify(msg, '=') << std::endl;

          // iteate until error is sufficient or we reach max iterations
          size_t enrichmentIteration = 0;
          while (enrichmentIteration < enrichmentMaxIterations
                 && localErrors->backend().maxCoeff() > enrichmentMaxLocalError) {
            if (debugLogging) {
              msg = "enrichment iteration " + Dune::Stuff::Common::toString(enrichmentIteration);
              info << msg << std::endl;
              info << Dune::Stuff::Common::whitespaceify(msg, '=') << std::endl;
            }
            info << "local errors are: "
                 << std::scientific << localErrors->backend().transpose() << std::fixed << std::endl;
            // enrich on all subdomains
            MarkerType markedSubdomains(detailedSolver->msGrid()->size());
            for (size_t ii = 0; ii < detailedSolver->msGrid()->size(); ++ii)
              markedSubdomains.set(ii, ii);
            info << "enriching";
            if (!debugLogging)
              info << "... " << std::flush;
            else
              info << " on subdomains " << markedSubdomains.backend().transpose() << ":" << std::endl;
            timer.reset();
            // therefore, compute the local corrections
            computeLocalCorrections(*detailedSolver,
                                    description,
                                    localReducedBases,
                                    *mapper,
                                    reducedSolution,
                                    mu,
                                    markedSubdomains,
                                    localVectors, // <- the local corrections will be in here
                                    localVectorsTwo);
            // extend the local reduced bases (and obtain a new mapper)
            mapper = extendReducedBases(description,
                                        markedSubdomains,
                                        localScalarProducts,
                                        localVectors, // <- this will change the local vector bc of gram schmidt
                                        localReducedBases);
            // and project the operator and functional
            if (debugLogging)
              info << "  reducing operator:" << std::endl;
            reducedPair = createReducedOperatorFunctional(*detailedSolver,
                                                          localReducedBases,
                                                          *mapper,
                                                          localVector,
                                                          globalVectorOne,
                                                          globalVectorTwo);
            reducedOperator = reducedPair.first;
            reducedFunctional = reducedPair.second;
            if (!debugLogging)
              info << " done (took " << timer.elapsed() << "sec)" << std::endl;
            // compute the new local errors
            info << "computing local errors... " << std::flush;
            timer.reset();
            computeReducedSolution(*detailedSolver, description, *reducedOperator, *reducedFunctional, mu, reducedSolution);
            reconstructReducedSolution(*detailedSolver, localReducedBases, *mapper, reducedSolution, localVectors);
            localErrors = computeLocalErrors(*detailedSolver,
                                             description,
                                             localScalarProducts,
                                             localReducedBases,
                                             *mapper,
                                             reducedSolution,
                                             mu,
                                             localVectors,
                                             localVectorsTwo);
            info << " done (took " << timer.elapsed() << "sec)" << std::endl;
            // increase iteration count
            ++enrichmentIteration;
          } // iteate until error is sufficient or we reach max iterations

        } // enrich, if max local error is too large
      } // loop over all parameters

    } else {
      nonparametericRun(*detailedSolver, description);
    } // if (detailedSolver->parametric())

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
