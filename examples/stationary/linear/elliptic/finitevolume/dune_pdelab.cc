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

// dune-stuff
#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/function/expression.hh>

// dune-detailed-solvers
#include <dune/detailed/solvers/stationary/linear/elliptic/model/default.hh>
#include <dune/detailed/solvers/stationary/linear/elliptic/model/thermalblock.hh>
#include <dune/detailed/solvers/stationary/linear/elliptic/finitevolume/dune-pdelab.hh>

#ifdef POLORDER
const int polOrder = POLORDER;
#else
const int polOrder = 1;
#endif

const std::string id = "stationary.linear.elliptic.fv.dp";

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
    file << "model = detailed.solvers.stationary.linear.elliptic.model.default" << std::endl;
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
    file << "[detailed.solvers.stationary.linear.elliptic.finitevolume.dune-pdelab]" << std::endl;
    file << "solve.type = eigen.bicgstab.incompletelut" << std::endl;
    file << "solve.maxIter = 5000"  << std::endl;
    file << "solve.precision = 1e-12"  << std::endl;
    file << "visualize.filename = " << id << ".solution"  << std::endl;
    file << "visualize.name = solution"  << std::endl;
    file.close();
  } // only write param file if there is none
} // void ensureParamFile()

template< class DomainFieldType, int dimDomain, class RangeFieldType, int dimRange >
Dune::shared_ptr< Dune::Detailed::Solvers::Stationary::Linear::Elliptic::Model::Interface< DomainFieldType, dimDomain, RangeFieldType, dimRange > >
  createModel(const Dune::Stuff::Common::ExtendedParameterTree& paramTree)
{
  // prepare
  paramTree.assertKey(id + ".model");
  const std::string modelId = paramTree.sub(id).get("model", "default_value");
  typedef Dune::Detailed::Solvers::Stationary::Linear::Elliptic::Model::Interface< DomainFieldType, dimDomain, RangeFieldType, dimRange > ModelInterfaceType;
  // choose model
  if (modelId == "detailed.solvers.stationary.linear.elliptic.model.default") {
    typedef Dune::Detailed::Solvers::Stationary::Linear::Elliptic::Model::Default< DomainFieldType, dimDomain, RangeFieldType, dimRange > ModelType;
    paramTree.assertSub(ModelType::id(), id);
    Dune::shared_ptr< ModelInterfaceType > model(new ModelType(paramTree.sub(ModelType::id())));
    return model;
  } else if (modelId == "detailed.solvers.stationary.linear.elliptic.model.thermalblock") {
    typedef Dune::Detailed::Solvers::Stationary::Linear::Elliptic::Model::Thermalblock< DomainFieldType, dimDomain, RangeFieldType, dimRange > ModelType;
    paramTree.assertSub(ModelType::id(), id);
    Dune::shared_ptr< ModelInterfaceType > model(new ModelType(paramTree.sub(ModelType::id())));
    return model;
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
    Dune::MPIHelper::instance(argc, argv);

    // parameter
    const std::string filename = id + ".param";
    ensureParamFile(filename);
    Dune::Stuff::Common::ExtendedParameterTree paramTree(argc, argv, filename);

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
    typedef Dune::Stuff::Grid::Provider::Cube<> GridProviderType;
    paramTree.assertSub(GridProviderType::id(), id);
    const GridProviderType gridProvider(paramTree.sub(GridProviderType::id()));
    typedef GridProviderType::GridType GridType;
    const GridType& grid = gridProvider.grid();
    typedef GridType::LeafGridView GridViewType;
    const Dune::shared_ptr< const GridViewType > gridView(new GridViewType(grid.leafView()));
    info << "  done (took " << timer.elapsed()
         << " sec, has " << grid.size(0) << " element";
    if (grid.size(0) > 1)
      info << "s";
    info << ")" << std::endl;
    info << "visualizing grid... " << std::flush;
    timer.reset();
    gridProvider.visualize(paramTree.sub(GridProviderType::id()).get("filename", id + "_grid"));
    info << " done (took " << timer.elapsed() << " sek)" << std::endl;
    debug.resume();

    // model
    info << "setting up model... " << std::flush;
    timer.reset();
    const unsigned int DUNE_UNUSED(dimDomain) = GridProviderType::dim;
    const unsigned int DUNE_UNUSED(dimRange) = 1;
    typedef GridProviderType::CoordinateType::value_type DomainFieldType;
    typedef DomainFieldType RangeFieldType;
    paramTree.assertSub(id, id);
    typedef Dune::Detailed::Solvers
        ::Stationary
        ::Linear
        ::Elliptic
        ::Model::Interface<  DomainFieldType, dimDomain, RangeFieldType, dimRange > ModelType;
    const Dune::shared_ptr< const ModelType >
        model = createModel< DomainFieldType, dimDomain, RangeFieldType, dimRange >(paramTree);
    typedef Dune::Stuff::Grid::BoundaryInfo::AllDirichlet BoundaryInfoType;
    const Dune::shared_ptr< const BoundaryInfoType > boundaryInfo(new BoundaryInfoType());
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    // solver
    info << "initializing solver:" << std::endl;
    typedef Dune::Detailed::Solvers
        ::Stationary
        ::Linear
        ::Elliptic
        ::FiniteVolume::DunePdelab< ModelType, GridViewType, BoundaryInfoType, polOrder > SolverType;
    paramTree.assertSub(SolverType::id(), id);
    SolverType solver(model, gridView, boundaryInfo);
    solver.init("  ", debug);

    info << "solving:" << std::endl;
    typedef SolverType::TrialVectorType SolutionType;
    typedef Dune::shared_ptr<SolutionType> solution_ptr;
    solution_ptr solution = solution_ptr(new SolutionType(*(solver.gridFunctionSpace()),0.0));
//    evtl. dieses Anlegen des Vektors in eine Methode createVector() auslagern (die in DunePdelab liegt)?!
//    solution_ptr solution = solver.createVector();
    solver.solve(*solution, paramTree.sub(SolverType::id()).sub("solve"), "  ", debug);

    info << "postprocessing:" << std::endl;
    solver.visualize(*solution,
                     paramTree.sub(SolverType::id()).sub("visualize").get("filename", id + "_solution"),
                     paramTree.sub(SolverType::id()).sub("visualize").get("name", "solution"),
                     "  ",
                     debug);

    // if we came that far we can as well be happy about it
    return 0;
  } catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e.what() << std::endl;
  } catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
  } catch ( ... ) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  } // try
} // main
