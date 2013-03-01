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
typedef Dune::Stuff::Grid::BoundaryInfo::Interface< GlobalGridViewType >                BoundaryInfoType;
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
typedef SolverType::VectorType                                                          VectorType;
typedef typename ModelType::ParamType                                                   ParamType;
typedef SolverType::MatrixType                                                          ScalarProductType;
typedef Dune::Stuff::LA::Container::EigenDenseMatrix< RangeFieldType >                  ReducedBasisType;
typedef typename SolverType::AnsatzMapperType                                           MapperType;
typedef SolverType::SeparableMatrixType                                                 DetailedOperatorType;
typedef VectorType                                                                      DetailedFunctionalType;
typedef Dune::Stuff::LA::Container::Separable< ReducedBasisType >                       ReducedOperatorType;
typedef VectorType                                                                      ReducedFunctionalType;


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
//  const std::string filename = description.get< std::string >("filename");
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


int main(int argc, char** argv)
{
  try {
    // description
    const DescriptionType description = initDescription(argc, argv);
    // logger
    Dune::Stuff::Common::LogStream& info = Dune::Stuff::Common::Logger().info();
    const bool debugLogging = description.get< bool >("logging.debug", false);
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
      std::vector< Dune::shared_ptr< ReducedBasisType > > localReducedBases(detailedSolver->msGrid()->size());
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

      // compute the global training snapshot
      const DescriptionType& parameterDescription = description.sub("parameter");
      const ParamType trainingParam
          = parameterDescription.getDynVector< ParamFieldType >("train", detailedSolver->model()->paramSize());
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
      MapperType mapper;
      initReducedBasis(localScalarProducts, localVectors, localReducedBases, mapper);
      info << " done (took " << timer.elapsed() << "sec)" << std::endl;

      // project the detailed operator and functional onto the reduced space
      info << "reducing operator";
      if (!debugLogging)
        info << "... " << std::flush;
      else
        info << ":" << std::endl;
      auto reducedPair = createReducedOperatorFunctional(*detailedSolver,
                                                         localReducedBases,
                                                         mapper,
                                                         localVector,
                                                         globalVectorOne,
                                                         globalVectorTwo);
      Dune::shared_ptr< ReducedOperatorType > reducedOperator = reducedPair.first;
      Dune::shared_ptr< ReducedFunctionalType > reducedFunctional = reducedPair.second;
      if (!debugLogging)
        info << " done (took " << timer.elapsed() << "sec)" << std::endl;




      // report
      msg = "  offline phase finished, took " + Dune::Stuff::Common::toString(offlineTimer.elapsed()) + "sec  ";
      info << Dune::Stuff::Common::whitespaceify(msg, '=') << std::endl;
      info << msg << std::endl;
      info << Dune::Stuff::Common::whitespaceify(msg, '=') << std::endl;

      for (auto& component : reducedOperator->components()) {
        info << "=====================================" << std::endl;
        info << component->backend() << std::endl;
      }
      info << "=====================================" << std::endl;

//      info << "reducing operator and functional... " << std::flush;
//      // therefore get the detailed operator and functional
//      Dune::shared_ptr< const DetailedOperatorType > detailedOperator = solver.systemMatrix("diffusion");
//      Dune::shared_ptr< const VectorType > detailedFunctional = solver.systemVector("force");
//      std::vector< Dune::shared_ptr< BasisMatrixType > > reducedOperatorComponents(detailedOperator->numComponents());
//      for (auto& component : reducedOperatorComponents)
//        component = Dune::make_shared< BasisMatrixType >(msGrid->size(), msGrid->size());
//      ReducedOperatorType reducedOperator(detailedOperator->paramSize(),
//                                          reducedOperatorComponents,
//                                          detailedOperator->coefficients());
//      VectorType reducedFunctional(msGrid->size());
//      for (size_t ii = 0; ii < msGrid->size(); ++ii) {
//        for (size_t jj = 0; jj < msGrid->size(); ++jj) {
//        // project all components of the detailed operator onto the reduced space
//          for (size_t qq = 0; qq < detailedOperator->numComponents(); ++qq) {
//            BasisMatrixType& reducedOperatorComponent = *(reducedOperator.components()[qq]);
//            reducedOperatorComponent.set(ii, jj, reducedBasis.backend().col(ii).transpose()
//                                                 * detailedOperator->components()[qq]->backend()
//                                                 * reducedBasis.backend().col(jj));
//          } // project all components of the detailed operator onto the reduced space
//        }
//        // project the detailed functional onto the reduced space
//        reducedFunctional.set(ii, detailedFunctional->backend().transpose() * reducedBasis.backend().col(ii));
//      }
//      info << "done (took " << timer.elapsed() << " sec)" << std::endl;

//      info << "=============================================" << std::endl;
//      info << "offline-phase finished, starting online-phase" << std::endl;
//      info << "=============================================" << std::endl;

//      // this is for the training parameter and should show now error
//      info << "solving for " << Dune::Stuff::Common::Parameter::report(trainingParam) << "... " << std::flush;
//      timer.reset();
//      Dune::shared_ptr< const BasisMatrixType > reducedSystemMatrix
//          = reducedOperator.fix(model->mapParam(trainingParam, "diffusion"));
//      VectorType reducedSolution;
//      reducedSolution.backend() = reducedSystemMatrix->backend().colPivHouseholderQr().solve(reducedFunctional.backend());
//      info << "done (took " << std::scientific << timer.elapsed() << std::fixed << " sec)" << std::endl;
//      info << "estimating error for " << Dune::Stuff::Common::Parameter::report(trainingParam);
//      if (!debugLogging)
//        info << "... " << std::flush;
//      else
//        info << ":" << std::endl;
//      timer.reset();
//      Dune::Timer debugTimer;
//      debug << "  reconstructing solution... " << std::flush;
//      std::vector< Dune::shared_ptr< VectorType > > reconstructedSolution = solver.createVector();
//      tmpGlobalVector->backend().setZero();
//      for (size_t ss = 0; ss < msGrid->size(); ++ss) {
//        tmpGlobalVector->backend() += reducedBasis.backend().col(ss) * reducedSolution.get(ss);
//      }
//      solver.localizeVector(tmpGlobalVector, reconstructedSolution);
//      solver.visualize(reconstructedSolution,
//                       filename + ".reconstructed.training.solution",
//                       id() + ".reconstructed.training.solution",
//                       "  ",
//                       devnull);
//      debug << "done (took " << debugTimer.elapsed() << " sec)" << std::endl;
//      debug << "  computing detailed solution for "
//            << Dune::Stuff::Common::Parameter::report(trainingParam) << "... " << std::flush;
//      std::vector< Dune::shared_ptr< VectorType > > detailedSolution = solver.createVector();
//      solver.solve(detailedSolution, trainingParam,
//                   linearSolverType,
//                   linearSolverMaxIter,
//                   linearSolverPrecision,
//                   "", devnull);
//      debug << "done (took " << debugTimer.elapsed() << " sec)" << std::endl;
//      debug << "  computing local estimates:  ";
//      Dune::DynamicVector< RangeFieldType > localEstimates(msGrid->size());
//      for (size_t ss = 0; ss < msGrid->size(); ++ss) {
//        localEstimates[ss] = std::sqrt((detailedSolution[ss]->backend().transpose() - reconstructedSolution[ss]->backend().transpose())
//                             * localScalarProducts[ss]->backend()
//                             * (detailedSolution[ss]->backend() - reconstructedSolution[ss]->backend()));
//        debug << std::scientific << localEstimates[ss] << std::fixed << " " << std::flush;
//      }
//      debug << std::endl;
//      if (!debugLogging)
//        info << "done (took " << timer.elapsed() << " sec)" << std::endl;

//      // this is for the test parameter and should invoke the postprocessing
//      info << "solving for " << Dune::Stuff::Common::Parameter::report(testParam) << "... " << std::flush;
//      timer.reset();
//      reducedSystemMatrix = reducedOperator.fix(model->mapParam(testParam, "diffusion"));
//      reducedSolution.backend() = reducedSystemMatrix->backend().colPivHouseholderQr().solve(reducedFunctional.backend());
//      info << "done (took " << std::scientific << timer.elapsed() << std::fixed << " sec)" << std::endl;
//      info << "estimating error for " << Dune::Stuff::Common::Parameter::report(testParam);
//      if (!debugLogging)
//        info << "... " << std::flush;
//      else
//        info << ":" << std::endl;
//      timer.reset();
//      debug << "  reconstructing solution... " << std::flush;
//      tmpGlobalVector->backend().setZero();
//      for (size_t ss = 0; ss < msGrid->size(); ++ss) {
//        tmpGlobalVector->backend() += reducedBasis.backend().col(ss) * reducedSolution.get(ss);
//      }
//      solver.localizeVector(tmpGlobalVector, reconstructedSolution);
//      solver.visualize(reconstructedSolution,
//                       filename + ".reconstructed.test.solution",
//                       id() + ".reconstructed.test.solution",
//                       "  ",
//                       devnull);
//      debug << "done (took " << debugTimer.elapsed() << " sec)" << std::endl;
//      debug << "  computing detailed solution for "
//            << Dune::Stuff::Common::Parameter::report(testParam) << "... " << std::flush;
//      solver.solve(detailedSolution, testParam,
//                   linearSolverType,
//                   linearSolverMaxIter,
//                   linearSolverPrecision,
//                   "", devnull);
//      solver.visualize(detailedSolution,
//                       filename + ".test.snapshot",
//                       id() + ".test.snapshot",
//                       "  ",
//                       devnull);
//      debug << "done (took " << debugTimer.elapsed() << " sec)" << std::endl;
//      debug << "  computing local estimates:  ";
//      for (size_t ss = 0; ss < msGrid->size(); ++ss) {
//        localEstimates[ss] = std::sqrt((detailedSolution[ss]->backend().transpose() - reconstructedSolution[ss]->backend().transpose())
//                             * localScalarProducts[ss]->backend()
//                             * (detailedSolution[ss]->backend() - reconstructedSolution[ss]->backend()));
//        debug << std::scientific << localEstimates[ss] << std::fixed << " " << std::flush;
//      }
//      debug << std::endl;
//      if (!debugLogging)
//        info << "done (took " << timer.elapsed() << " sec)" << std::endl;
//      // carry out online enrichment
//      info << "enriching local bases:" << std::endl;
//      BasisMatrixType extendedReducedBasis(solver.ansatzMapper().size(), 2*msGrid->size());
//      // * therefore, fix a model for the online parameter,
//      const Dune::shared_ptr< const ModelType > fixedModel = model->fix(testParam);
//      // * loop over all subdomains
//      for (size_t ss = 0; ss < msGrid->size(); ++ss) {
//        info << "  solving on subdomain " << ss;
////        if (!debugLogging)
//          info << "... " << std::flush;
////        else
////          info << ":" << std::endl;
//        timer.reset();
//        // * get the local grid parts
//        typedef MsGridType::LocalGridPartType LocalGridPartType;
////        const Dune::shared_ptr< const LocalGridPartType > localGridPart = msGrid->localGridPart(ss);
//        const Dune::shared_ptr< const LocalGridPartType > oversampledLocalGridPart = msGrid->localGridPart(ss, true);
//        // * initialize a solver on the local oversampled gripart
//        const Dune::shared_ptr< const Dune::Stuff::Grid::BoundaryInfo::Interface< typename MsGridType::LocalGridViewType > >
//            localBoundaryInfo(Dune::Stuff::Grid::BoundaryInfo::create< typename MsGridType::LocalGridViewType >("stuff.grid.boundaryinfo.alldirichlet"));
//        typedef SolverType::LocalSolverType LocalSolverType;
//        LocalSolverType localSolver(oversampledLocalGridPart, localBoundaryInfo, fixedModel);
//        localSolver.init("    ", devnull);
//        // * project the current reconstructed solution onto the local oversampled gripart
//        typedef typename SolverType::DiscreteFunctionType GlobalDiscreteFunctionType;
//        Dune::shared_ptr< GlobalDiscreteFunctionType >
//            reconstructedSolutionFunction = solver.createDiscreteFunction(reconstructedSolution);
//        typedef typename LocalSolverType::DiscreteAnsatzFunctionType LocalDiscreteFunctionType;
//        Dune::shared_ptr< LocalDiscreteFunctionType >
//            localReconstructedSolutionFunction = localSolver.createAnsatzFunction();
//        Dune::Stuff::DiscreteFunction::Projection::Lagrange::project(*reconstructedSolutionFunction,
//                                                                     *localReconstructedSolutionFunction);
//        // * compute the local right hand side  f := f - lrs' * A
//        auto localSystemMatrixPtr = localSolver.matrix("diffusion");
//        auto forceVectorPtr = localSolver.vector("force");
//        assert(localSystemMatrixPtr->numComponents() == 1);
//        assert(forceVectorPtr->numComponents() == 1);
//        forceVectorPtr->components()[0]->backend()
//            -= localReconstructedSolutionFunction->vector()->backend().transpose()
//               * localSystemMatrixPtr->components()[0]->backend();
//        // * solve locally on the oversampled domain
//        Dune::shared_ptr< VectorType > localOversampledSolution = localSolver.createAnsatzVector();
//        localSolver.solve(localOversampledSolution,
//                          linearSolverType,
//                          linearSolverMaxIter,
//                          linearSolverPrecision,
//                          "    ",
//                          devnull);
//        // * restrict to the local (not oversampled) subdomain
//        Dune::shared_ptr< LocalDiscreteFunctionType >
//            localOversampledSolutionFunction = localSolver.createAnsatzFunction(localOversampledSolution);
//        Dune::shared_ptr< LocalDiscreteFunctionType > localSolutionFunction = solver.localSolver(ss)->createAnsatzFunction();
//        Dune::Stuff::DiscreteFunction::Projection::Lagrange::project(*localOversampledSolutionFunction,
//                                                                     *localSolutionFunction);
////        if (!debugLogging)
//          info << "done (took " << timer.elapsed() << " sec)" << std::endl;
//        // now create the new reduced basis
//        info << "  extending basis on subdomain " << ss << "... " << std::flush;
//        timer.reset();
//        // * therefore copy the old basis function over,
//        extendedReducedBasis.backend().col(2*ss) = reducedBasis.backend().col(ss);
//        // * apply a gram-schmidt to the new local solution
//        const Dune::shared_ptr< ScalarProductType >& localScalarProduct = localScalarProducts[ss];
//        tmpGlobalVector->backend() = reducedBasis.backend().col(ss);
//        solver.localizeVector(tmpGlobalVector, tmpLocalVectors);
//        bool success = Dune::Stuff::LA::Algorithm::gramSchmidt(*localScalarProduct,
//                                                               *(tmpLocalVectors[ss]),
//                                                               *(localSolutionFunction->vector()),
//                                                               linearSolverPrecision);
//        assert(success);
//        solver.globalizeVector(localSolutionFunction->vector(), ss, tmpGlobalVector);
//        extendedReducedBasis.backend().col(2*ss + 1) = tmpGlobalVector->backend();
////        if (!debugLogging)
//          info << "done (took " << timer.elapsed() << " sec)" << std::endl;
//      } // carry out online enrichment

//      info << "reducing operator and functional... " << std::flush;
//      std::vector< Dune::shared_ptr< BasisMatrixType > > extendedReducedOperatorComponents(detailedOperator->numComponents());
//      for (auto& component : extendedReducedOperatorComponents)
//        component = Dune::make_shared< BasisMatrixType >(2*msGrid->size(), 2*msGrid->size());
//      ReducedOperatorType extendedReducedOperator(detailedOperator->paramSize(),
//                                                  extendedReducedOperatorComponents,
//                                                  detailedOperator->coefficients());
//      VectorType extendedReducedFunctional(2*msGrid->size());
//      for (size_t ii = 0; ii < msGrid->size(); ++ii) {
//        for (size_t jj = 0; jj < msGrid->size(); ++jj) {
//        // project all components of the detailed operator onto the reduced space
//          for (size_t qq = 0; qq < detailedOperator->numComponents(); ++qq) {
//            BasisMatrixType& reducedOperatorComponent = *(extendedReducedOperator.components()[qq]);
//            reducedOperatorComponent.set(ii, jj, extendedReducedBasis.backend().col(ii).transpose()
//                                                 * detailedOperator->components()[qq]->backend()
//                                                 * extendedReducedBasis.backend().col(jj));
//          } // project all components of the detailed operator onto the reduced space
//        }
//        // project the detailed functional onto the reduced space
//        extendedReducedFunctional.set(ii, detailedFunctional->backend().transpose() * extendedReducedBasis.backend().col(ii));
//      }
//      info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    } else {
      nonparametericRun(*detailedSolver, description);
    } // if (detailedSolver->parametric())

    // if we came that far we can as well be happy about it
    return 0;
  } catch(Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e.what() << std::endl;
  } catch(std::exception& e) {
    std::cerr << e.what() << std::endl;
  } catch( ... ) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  } // try
} // main
