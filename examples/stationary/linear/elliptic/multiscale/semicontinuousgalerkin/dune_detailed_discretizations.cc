
#include "config.h"

// boost
#include <boost/filesystem.hpp>

// dune-common
#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>
#include <dune/common/timer.hh>

// dune-grid-multiscale
#include <dune/grid/multiscale/provider/cube.hh>

// dune-stuff
#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/logging.hh>

// dune-detailed-solvers
#include <dune/detailed/solvers/stationary/linear/elliptic/model.hh>
#include <dune/detailed/solvers/stationary/linear/elliptic/multiscale/semicontinuousgalerkin/dune-detailed-discretizations.hh>

#ifdef POLORDER
const int polOrder = POLORDER;
#else
const int polOrder = 1;
#endif

const std::string id = "semicontinuous_ddd_multiscale_solver";

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
    file << "[grid.multiscale.provider.cube]" << std::endl;
    file << "level = 4" << std::endl;
    file << "boundaryId = 7" << std::endl;                  // a cube from the factory gets the boundary ids 1 to 4 ind 2d and 1 to 6 ? in 3d
    file << "partitions.0 = 2" << std::endl;
    file << "partitions.1 = 2" << std::endl;
    file << "partitions.2 = 2" << std::endl;
    file << "filename = " << id << "_msGrid" << std::endl;
    file << "[detailed.solvers.stationary.linear.elliptic.model.default]" << std::endl;
    file << "diffusion.order = 0" << std::endl;
    file << "diffusion.variable = x" << std::endl;
    file << "diffusion.expression.0 = 1.0"  << std::endl;
    file << "diffusion.expression.1 = 1.0"  << std::endl;
    file << "diffusion.expression.2 = 1.0"  << std::endl;
    file << "force.order = 0" << std::endl;
    file << "force.variable = x" << std::endl;
    file << "force.expression.0 = 1.0"  << std::endl;
    file << "force.expression.1 = 1.0"  << std::endl;
    file << "force.expression.2 = 1.0"  << std::endl;
    file << "[detailed.solvers.stationary.linear.elliptic.multiscale.semicontinuousgalerkin]" << std::endl;
    file << "discretization.penaltyFactor = 10.0" << std::endl;
    file << "solve.type = eigen.bicgstab.incompletelut" << std::endl;
    file << "solve.maxIter = 5000"  << std::endl;
    file << "solve.precision = 1e-12"  << std::endl;
    file << "visualize.filename = " << id << "_solution"  << std::endl;
    file << "visualize.name = solution"  << std::endl;
    file.close();
  } // only write param file if there is none
} // void ensureParamFile()

int main(int argc, char** argv)
{
  try {
    // mpi
    Dune::MPIHelper::instance(argc, argv);

    // parameter
    const std::string filename = id + ".param";
    ensureParamFile(filename);
    Dune::ParameterTree paramTree = Dune::Stuff::Common::Parameter::Tree::init(argc, argv, filename);

    // logger
    Dune::Stuff::Common::Logger().create(Dune::Stuff::Common::LOG_INFO |
                                         Dune::Stuff::Common::LOG_CONSOLE |
                                         Dune::Stuff::Common::LOG_DEBUG);
//    Dune::Stuff::Common::LogStream& info = Dune::Stuff::Common::Logger().info();
    std::ostream& info = std::cout;
    Dune::Stuff::Common::LogStream& debug = Dune::Stuff::Common::Logger().debug();

    // timer
    Dune::Timer timer;

    // grid
    info << "setting up grid: " << std::endl;
    debug.suspend();
    typedef Dune::grid::Multiscale::Provider::Cube<> GridProviderType;
    Dune::Stuff::Common::Parameter::Tree::assertSub(paramTree, GridProviderType::id, id);
    const GridProviderType gridProvider(paramTree.sub(GridProviderType::id));
    typedef GridProviderType::MsGridType MsGridType;
    const MsGridType& msGrid = gridProvider.msGrid();
    info << "  took " << timer.elapsed()
         << " sec (has " << msGrid.globalGridPart()->grid().size(0) << " elements, "
         << msGrid.size() << " subdomains)" << std::endl;
    debug.resume();
    info << "visualizing grid... " << std::flush;
    debug.suspend();
    timer.reset();
    msGrid.visualize(paramTree.sub(GridProviderType::id).get("filename", id + "_msGrid"));
    info << " done (took " << timer.elapsed() << " sek)" << std::endl;
    debug.resume();

    // model
    info << "setting up model... " << std::flush;
    timer.reset();
    const unsigned int DUNE_UNUSED(dimDomain) = GridProviderType::dim;
    const unsigned int DUNE_UNUSED(dimRange) = 1;
    typedef GridProviderType::CoordinateType::value_type DomainFieldType;
    typedef DomainFieldType RangeFieldType;
    typedef Dune::Detailed::Solvers::Stationary::Linear::Elliptic::Model::Default< DomainFieldType, dimDomain, RangeFieldType, dimRange > ModelType;
    Dune::Stuff::Common::Parameter::Tree::assertSub(paramTree, ModelType::id, id);
    const ModelType model(paramTree.sub(ModelType::id));
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    // solver
    info << "initializing solver... " << std::endl;
//    debug.suspend();
    timer.reset();
    typedef Dune::Detailed::Solvers::Stationary::Linear::Elliptic::Multiscale::SemicontinuousGalerkin::DuneDetailedDiscretizations< ModelType, MsGridType, polOrder > SolverType;
    Dune::Stuff::Common::Parameter::Tree::assertSub(paramTree, SolverType::id, id);
    SolverType solver(model, msGrid, paramTree.sub(SolverType::id));
    solver.init("  ", std::cout);
//    debug.resume();
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    info << "solving... " << std::endl;
//    debug.suspend();
    timer.reset();
    typedef SolverType::LocalVectorType LocalDofVectorType;
    std::vector< Dune::shared_ptr< LocalDofVectorType > > solution = solver.createVector();
    solver.solve(solution, paramTree.sub(SolverType::id).sub("solve"), "  ", std::cout);
    //    debug.resume();
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    info << "postprocessing..." << std::endl;
//    debug.suspend();
    timer.reset();
    solver.visualize(solution,
                     paramTree.sub(SolverType::id).sub("visualize").get("filename", id + "_solution"),
                     paramTree.sub(SolverType::id).sub("visualize").get("name", "solution"),
                     "  ",
                     std::cout);
//    debug.resume();
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

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
