#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

// system
#include <sstream>

// boost
#include <boost/filesystem.hpp>

// dune-common
#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>
#include <dune/common/timer.hh>

// dune-fem
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/misc/gridwidth.hh>

// dune-grid-multiscale
#include <dune/grid/part/leaf.hh>

// dune-stuff
#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/grid/provider/cube.hh>
#if HAVE_ALUGRID
  #include <dune/stuff/grid/provider/gmsh.hh>
#endif // HAVE_ALUGRID
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/function/expression.hh>
#include <dune/stuff/discretefunction/norm.hh>

// dune-detailed-solvers
#include <dune/detailed/solvers/stationary/linear/elliptic/model/default.hh>
#include <dune/detailed/solvers/stationary/linear/elliptic/model/thermalblock.hh>
#include <dune/detailed/solvers/stationary/linear/elliptic/continuousgalerkin/dune-detailed-discretizations.hh>

#ifdef POLORDER
const int polOrder = POLORDER;
#else
const int polOrder = 1;
#endif

const std::string id = "stationary.linear.elliptic.cg.ddd";

/**
  \brief      Creates a parameter file if it does not exist.

              Nothing is done if the file already exists. If not, a parameter file will be created with all neccessary
              keys and values.
  \param[in]  filename
              (Relative) path to the file.
  **/
void ensureParamFile(std::string filename)
{
  // only write param file if there is none
  if (!boost::filesystem::exists(filename)) {
    std::ofstream file;
    file.open(filename);
    file << "[" << id << "]" << std::endl;
    file << "grid = stuff.grid.provider.cube" << std::endl;
    file << "model = detailed.solvers.stationary.linear.elliptic.model.default" << std::endl;
    file << "exact_solution.order = 2" << std::endl;
    file << "exact_solution.variable = x" << std::endl;
    file << "exact_solution.expression = [-0.5*x[0]*x[0] + 0.5*x[0]]" << std::endl;
    file << "[stuff.grid.provider.cube]" << std::endl;
    file << "lowerLeft = [0.0; 0.0; 0.0]" << std::endl;
    file << "upperRight = [1.0; 1.0; 1.0]" << std::endl;
    file << "numElements = 4" << std::endl;
    file << "filename = " << id << ".grid" << std::endl;
    file << "[detailed.solvers.stationary.linear.elliptic.model.default]" << std::endl;
    file << "diffusion.order = 0"  << std::endl;
    file << "diffusion.variable = x" << std::endl;
    file << "diffusion.expression = [1.0; 1.0; 1.0]"  << std::endl;
    file << "force.order = 0"  << std::endl;
    file << "force.variable = x" << std::endl;
    file << "force.expression = [1.0; 1.0; 1.0]"  << std::endl;
    file << "dirichlet.order = 0"  << std::endl;
    file << "dirichlet.variable = x" << std::endl;
    file << "dirichlet.expression = [1.0; 1.0; 1.0]"  << std::endl;
    file << "neumann.order = 0"  << std::endl;
    file << "neumann.variable = x" << std::endl;
    file << "neumann.expression = [1.0; 1.0; 1.0]"  << std::endl;
    file << "[detailed.solvers.stationary.linear.elliptic.model.thermalblock]" << std::endl;
    file << "diffusion.order = 0"  << std::endl;
    file << "diffusion.lowerLeft = [0.0; 0.0; 0.0]" << std::endl; // should coincide with the grid
    file << "diffusion.upperRight = [1.0; 1.0; 1.0]" << std::endl; // should coincide with the grid
    file << "diffusion.numElements = [2; 2; 2]"  << std::endl;
    file << "diffusion.components = [1.0; 10.0; 3.0; 2.1]"  << std::endl;
    file << "force.order = 0"  << std::endl;
    file << "force.variable = x" << std::endl;
    file << "force.expression = [1.0; 1.0; 1.0]"  << std::endl;
    file << "dirichlet.order = 0"  << std::endl;
    file << "dirichlet.variable = x" << std::endl;
    file << "dirichlet.expression = [1.0; 1.0; 1.0]"  << std::endl;
    file << "neumann.order = 0"  << std::endl;
    file << "neumann.variable = x" << std::endl;
    file << "neumann.expression = [1.0; 1.0; 1.0]"  << std::endl;
    file << "[detailed.solvers.stationary.linear.elliptic.continuousgalerkin]" << std::endl;
    file << "solve.type = eigen.bicgstab.incompletelut" << std::endl;
    file << "solve.maxIter = 5000"  << std::endl;
    file << "solve.precision = 1e-12"  << std::endl;
    file << "visualize.filename = " << id << ".solution"  << std::endl;
    file << "visualize.name = solution"  << std::endl;
    file.close();
  } // only write param file if there is none
} // void ensureParamFile()

template< class DiscreteFunctionType,
          class OutStreamType >
void compute_errors(const Dune::ParameterTree& paramTree,
                    const DiscreteFunctionType& discreteFunction,
                    OutStreamType& out,
                    std::string prefix = "")
{
  // preparations
  typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;
  typedef Dune::Stuff::Function::Expression<
      typename DiscreteFunctionType::DomainFieldType,
      DiscreteFunctionType::dimDomain,
      RangeFieldType,
      DiscreteFunctionType::dimRange >
    FunctionType;
  const FunctionType function = FunctionType::createFromParamTree(paramTree);
  const unsigned int functionOrder = paramTree.get("order", 100);
  // compute norms
  out << prefix << " norm | exact solution | discrete solution | error (abs) | error (rel)" << std::endl;
  out << prefix << "------+----------------+-------------------+-------------+-------------" << std::endl;
  const RangeFieldType L2_reference_norm = Dune::Stuff::DiscreteFunction::Norm::L2(discreteFunction.space().gridPart(),
                                                                                   function,
                                                                                   functionOrder);
  out.precision(2);
  out << prefix << "   L2 | " << std::setw(14) << std::scientific << L2_reference_norm
                                      << " | " << std::flush;
  const RangeFieldType L2_discrete_norm = Dune::Stuff::DiscreteFunction::Norm::L2(discreteFunction);
  out << std::setw(17) << std::scientific << L2_discrete_norm << " | " << std::flush;
  const RangeFieldType L2_difference = Dune::Stuff::DiscreteFunction::Norm::L2_difference(function, functionOrder, discreteFunction);
  out << std::setw(11) << std::scientific << L2_difference << " | " << std::flush;
  out << std::setw(11) << std::scientific << L2_difference/L2_reference_norm << std::endl;
  out << prefix << "------+----------------+-------------------+-------------+-------------" << std::endl;
} // void compute_norms(...)

Dune::Stuff::Grid::Provider::Interface<>* createProvider(const Dune::Stuff::Common::ExtendedParameterTree& paramTree)
{
  // prepare
  if (!paramTree.hasKey(id + ".grid"))
    DUNE_THROW(Dune::RangeError,
               "\nMissing key '" << id << ".grid' in the following Dune::ParameterTree:\n" << paramTree.reportString("  "));
  const std::string providerId = paramTree.get(id + ".grid", "meaningless_default_value");
  // choose provider
  if (providerId == "stuff.grid.provider.cube") {
    typedef Dune::Stuff::Grid::Provider::Cube<> CubeProviderType;
    CubeProviderType* cubeProvider = new CubeProviderType(CubeProviderType::createFromParamTree(paramTree));
    return cubeProvider;
#if HAVE_ALUGRID
  } else if (providerId == "stuff.grid.provider.gmsh") {
    typedef Dune::Stuff::Grid::Provider::Gmsh<> GmshProviderType;
    GmshProviderType* gmshProvider = new GmshProviderType(GmshProviderType::createFromParamTree(paramTree));
    return gmshProvider;
#endif // HAVE_ALUGRID
  } else
    DUNE_THROW(Dune::RangeError,
               "\nError: unknown provider ('" << providerId << "') given in the following Dune::Parametertree:\n" << paramTree.reportString("  "));
} // ... createProvider(...)

template< class DomainFieldType, int dimDomain, class RangeFieldType, int dimRange >
Dune::Detailed::Solvers::Stationary::Linear::Elliptic::Model::Interface< DomainFieldType, dimDomain, RangeFieldType, dimRange >*
createModel(const Dune::Stuff::Common::ExtendedParameterTree& paramTree)
{
  // prepare
  if (!paramTree.hasKey(id + ".model"))
    DUNE_THROW(Dune::RangeError,
               "\nMissing key '" << id << ".model' in the following Dune::ParameterTree:\n" << paramTree.reportString("  "));
  const std::string modelId = paramTree.sub(id).get("model", "meaningless_default_value");
  // choose model
  if (modelId == "detailed.solvers.stationary.linear.elliptic.model.default") {
    typedef Dune::Detailed::Solvers
        ::Stationary
        ::Linear
        ::Elliptic
        ::Model
        ::Default< DomainFieldType, dimDomain, RangeFieldType, dimRange >
      DefaultModelType;
    return new DefaultModelType(DefaultModelType::createFromParamTree(paramTree));
  } else if (modelId == "detailed.solvers.stationary.linear.elliptic.model.thermalblock") {
    typedef Dune::Detailed::Solvers
        ::Stationary
        ::Linear
        ::Elliptic
        ::Model
        ::Thermalblock< DomainFieldType, dimDomain, RangeFieldType, dimRange >
      ThermalblockModelType;
    return new ThermalblockModelType(ThermalblockModelType::createFromParamTree(paramTree));
  } else {
    std::stringstream msg;
    msg << std::endl << "Error in " << id << ": unknown model ('" << modelId << "') given in the following Dune::Parametertree" << std::endl;
    paramTree.report(msg);
    DUNE_THROW(Dune::InvalidStateException, msg.str());
  } // choose model
} // ... createModel(...)

int main(int argc, char** argv)
{
  try {
    // mpi
    Dune::MPIManager::initialize(argc, argv);

    // parameter
    const std::string filename = id + ".param";
    ensureParamFile(filename);
    Dune::Stuff::Common::ExtendedParameterTree paramTree(argc, argv, filename);
    if (!paramTree.hasSub(id))
      DUNE_THROW(Dune::RangeError,
                 "\nError: missing sub '" << id << "' in the following Dune::Parametertree:\n" << paramTree.reportString("  "));

    // logger
    Dune::Stuff::Common::Logger().create(Dune::Stuff::Common::LOG_INFO |
                                         Dune::Stuff::Common::LOG_CONSOLE |
                                         Dune::Stuff::Common::LOG_DEBUG);
    Dune::Stuff::Common::LogStream& info = Dune::Stuff::Common::Logger().info();
    Dune::Stuff::Common::LogStream& debug = Dune::Stuff::Common::Logger().debug();

    // timer
    Dune::Timer timer;

    // grid
    info << "setting up grid: " << std::endl;
    debug.suspend();
    typedef Dune::Stuff::Grid::Provider::Interface<> GridProviderType;
    const GridProviderType* gridProvider = createProvider(paramTree);
    typedef GridProviderType::GridType GridType;
    const Dune::shared_ptr< const GridType > grid = gridProvider->grid();
    typedef Dune::grid::Part::Leaf::Const< GridType > GridPartType;
    const Dune::shared_ptr< const GridPartType > gridPart(new GridPartType(*grid));
    info << "  done (took " << timer.elapsed()
         << " sec, has " << grid->size(0) << " element";
    if (grid->size(0) > 1)
      info << "s";
    info << " and a width of "
         << Dune::GridWidth::calcGridWidth(*gridPart) << ")" << std::endl;
    info << "visualizing grid... " << std::flush;
    timer.reset();
    gridProvider->visualize(id + ".grid");
    info << " done (took " << timer.elapsed() << " sek)" << std::endl;
    debug.resume();

    // model
    info << "setting up model... " << std::flush;
    timer.reset();
    const unsigned int DUNE_UNUSED(dimDomain) = GridProviderType::dim;
    const unsigned int DUNE_UNUSED(dimRange) = 1;
    typedef GridProviderType::CoordinateType::value_type DomainFieldType;
    typedef DomainFieldType RangeFieldType;
    typedef Dune::Detailed::Solvers
        ::Stationary
        ::Linear
        ::Elliptic
        ::Model
        ::Interface< DomainFieldType, dimDomain, RangeFieldType, dimRange >
      ModelType;
    const Dune::shared_ptr< const ModelType > model(createModel< DomainFieldType, dimDomain, RangeFieldType, dimRange >(paramTree));
    typedef Dune::Stuff::Grid::BoundaryInfo::AllDirichlet BoundaryInfoType;
    const Dune::shared_ptr< const BoundaryInfoType > boundaryInfo(new BoundaryInfoType());
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    // solver
    info << "initializing solver:" << std::endl;
//    timer.reset();
    typedef Dune::Detailed::Solvers
        ::Stationary
        ::Linear
        ::Elliptic
        ::ContinuousGalerkin::DuneDetailedDiscretizations< ModelType, GridPartType, BoundaryInfoType, polOrder >
        SolverType;
    SolverType solver(model, gridPart, boundaryInfo);
    solver.init("  ", info);
//    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    info << "solving:" << std::endl;
//    timer.reset();
    if (!paramTree.hasSub(solver.id()))
      DUNE_THROW(Dune::RangeError,
                 "\nError: missing sub '" << solver.id() << "' in the following Dune::Parametertree:\n" << paramTree.reportString("  "));
    if (!paramTree.sub(solver.id()).hasSub("solve"))
      DUNE_THROW(Dune::RangeError,
                 "\nError: missing sub '" << solver.id() << ".solve' in the following Dune::Parametertree:\n" << paramTree.reportString("  "));
    std::string type = "eigen.bicgstab.incompletelut";
    unsigned int maxIter = 5000;
    double precision = 1e-12;
    if (!paramTree.sub(solver.id()).sub("solve").hasKey("type"))
      std::cout << "Warning in " << id << ": no key '" << solver.id() << ".solve.type' given, defaulting to " << type << "!" << std::endl;
    else
      type = paramTree.sub(solver.id()).sub("solve").get("type", type);
    if (!paramTree.sub(solver.id()).sub("solve").hasKey("maxIter"))
      std::cout << "Warning in " << id << ": no key '" << solver.id() << ".solve.maxIter' given, defaulting to " << maxIter << "!" << std::endl;
    else
      maxIter = paramTree.sub(solver.id()).sub("solve").get("maxIter", maxIter);
    if (!paramTree.sub(solver.id()).sub("solve").hasKey("precision"))
      std::cout << "Warning in " << id << ": no key '" << solver.id() << ".solve.precision' given, defaulting to " << precision << "!" << std::endl;
    else
      precision = paramTree.sub(solver.id()).sub("solve").get("precision", precision);
    typedef SolverType::VectorBackendType DofVectorType;
    Dune::shared_ptr< DofVectorType > solution = solver.createVector();
    solver.solve(*solution, type, maxIter, precision, "  ", info);
//    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    info << "postprocessing:" << std::endl;
//    timer.reset();
    if (!paramTree.sub(solver.id()).hasSub("visualize"))
      DUNE_THROW(Dune::RangeError,
                 "\nError: missing sub '" << solver.id() << ".visualize' in the following Dune::Parametertree:\n" << paramTree.reportString("  "));
    std::string visualize_filename = id + ".solution";
    std::string name = "solution";
    if (!paramTree.sub(solver.id()).sub("visualize").hasKey("filename"))
      std::cout << "Warning in " << id << ": no key '" << solver.id() << ".visualize.filename' given, defaulting to " << visualize_filename << "!" << std::endl;
    else
      visualize_filename = paramTree.sub(solver.id()).sub("visualize").get("filename", visualize_filename);
    if (!paramTree.sub(solver.id()).sub("visualize").hasKey("name"))
      std::cout << "Warning in " << id << ": no key '" << solver.id() << ".visualize.name' given, defaulting to " << name << "!" << std::endl;
    else
      name = paramTree.sub(solver.id()).sub("visualize").get("name", name);
    solver.visualize(*solution, visualize_filename, name, "  ", info);
//    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    if (dimDomain == 1) {
      info << "computing errors:" << std::endl;
      const Dune::shared_ptr< const SolverType::DiscreteFunctionType > discreteFunction = solver.createDiscreteFunction(*solution, "discrete_solution");
      compute_errors(paramTree.sub(id).sub("exact_solution"), *discreteFunction, info, "  ");
    }

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
