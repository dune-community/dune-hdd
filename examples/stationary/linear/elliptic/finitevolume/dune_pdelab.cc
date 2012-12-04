#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#include <boost/filesystem.hpp>

#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/grid/provider.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/function/expression.hh>

#include <dune/detailed/solvers/stationary/linear/elliptic/model.hh>
#include <dune/detailed/solvers/stationary/linear/elliptic/finitevolume/dune-pdelab.hh>

#ifdef POLORDER
  const int polOrder = POLORDER;
#else
  const int polOrder = 1;
#endif

const std::string id = "stationary.linear.elliptic.fv.dp";

using namespace Dune;
using namespace Dune::Detailed::Solvers;

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
    file << "filename = " << id << std::endl;
    file << "[stuff.grid.provider.cube]" << std::endl;
    file << "lowerLeft = [0.0; 0.0; 0.0]" << std::endl;
    file << "upperRight = [1.0; 1.0; 1.0]" << std::endl;
    file << "numElements = 4" << std::endl;
    file << "[stuff.boundaryinfo.idbased]" << std::endl;
    file << "dirichlet = [1; 2; 3]" << std::endl;
    file << "neumann = [4]" << std::endl;
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
    file << "visualize.name = solution"  << std::endl;
    file.close();
  } // only write param file if there is none
} // void ensureParamFile()


int main(int argc, char** argv)
{
  try {
    // mpi
    MPIHelper::instance(argc, argv);

    // parameter
    const std::string filename = id + ".param";
    ensureParamFile(filename);
    Stuff::Common::ParameterTreeX paramTree(argc, argv, filename);
    paramTree.assertSub(id);

    // logger
    Stuff::Common::Logger().create(Stuff::Common::LOG_INFO |
                                   Stuff::Common::LOG_CONSOLE |
                                   Stuff::Common::LOG_DEBUG);
    Stuff::Common::LogStream& info = Stuff::Common::Logger().info();
    Stuff::Common::LogStream& debug = Stuff::Common::Logger().debug();

    // timer
    Timer timer;

    // grid
    info << "setting up grid: " << std::endl;
    debug.suspend();
    typedef Stuff::Grid::Provider::Interface<> GridProviderType;
    const std::string gridType = paramTree.get< std::string >(id + ".grid");
    const GridProviderType* gridProvider = Stuff::Grid::Provider::create(gridType,
                                                                         paramTree);
    typedef GridProviderType::GridType GridType;
    const shared_ptr< const GridType > grid = gridProvider->grid();
    typedef GridType::LeafGridView GridViewType;
    const shared_ptr< const GridViewType > gridView(new GridViewType(grid->leafView()));
    info << "  done (took " << timer.elapsed()
         << " sec, has " << grid->size(0) << " element";
    if (grid->size(0) > 1)
      info << "s";
    info << ")" << std::endl;
    info << "visualizing grid... " << std::flush;
    timer.reset();
    gridProvider->visualize(paramTree.get(id + ".filename", id) + ".grid");
    info << " done (took " << timer.elapsed() << " sek)" << std::endl;
    debug.resume();

    // model
    info << "setting up model... " << std::flush;
    timer.reset();
    const unsigned int DUNE_UNUSED(dimDomain) = GridProviderType::dim;
    const unsigned int DUNE_UNUSED(dimRange) = 1;
    typedef GridProviderType::CoordinateType::value_type DomainFieldType;
    typedef DomainFieldType RangeFieldType;
    typedef Stationary::Linear::Elliptic::Model::Interface< DomainFieldType, dimDomain,
                                                            RangeFieldType, dimRange >
        ModelType;
    const std::string modelType = paramTree.get< std::string >(id + ".model");
    const shared_ptr< const ModelType > model(
          Stationary::Linear::Elliptic::Model::create< DomainFieldType, dimDomain, RangeFieldType, dimRange >(modelType,
                                                                                                              paramTree));
    typedef Stuff::Grid::BoundaryInfo::IdBased BoundaryInfoType;
    const shared_ptr< const BoundaryInfoType > boundaryInfo(new BoundaryInfoType(
          BoundaryInfoType::createFromParamTree(paramTree)));
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    // solver
    info << "initializing solver:" << std::endl;
    typedef Stationary::Linear::Elliptic::FiniteVolume::DunePdelab< ModelType,
                                                                    GridViewType,
                                                                    BoundaryInfoType,
                                                                    polOrder >
        SolverType;
    const Stuff::Common::ParameterTreeX solverTree = paramTree.sub(SolverType::id());
    SolverType solver(model, gridView, boundaryInfo);
    solver.init("  ", debug);

    info << "solving:" << std::endl;
    typedef SolverType::TrialVectorType SolutionType;
    Dune::shared_ptr< SolutionType > solution = solver.createVector();
    solver.solve(*solution, paramTree, "  ", debug);

    info << "postprocessing:" << std::endl;
    solver.visualize(*solution,
                     paramTree.get(id + ".filename", id) + ".solution",
                     solverTree.get("visualize.name", "solution"),
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
