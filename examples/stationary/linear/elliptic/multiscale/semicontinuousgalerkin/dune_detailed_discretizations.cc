
#include "config.h"

// boost
#include <boost/filesystem.hpp>

// dune-common
#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>
#include <dune/common/timer.hh>

// dune-grid-multiscale
#include <dune/grid/multiscale/provider/cube.hh>

// dune-fem
#include <dune/fem/misc/mpimanager.hh>

// dune-stuff
#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/grid/boundaryinfo.hh>

// dune-detailed-solvers
#include <dune/detailed/solvers/stationary/linear/elliptic/model/default.hh>
//#include <dune/detailed/solvers/stationary/linear/elliptic/model/spe10.hh>
#include <dune/detailed/solvers/stationary/linear/elliptic/multiscale/semicontinuousgalerkin/dune-detailed-discretizations.hh>

#ifdef POLORDER
const int polOrder = POLORDER;
#else
const int polOrder = 1;
#endif

const std::string id = "stationary.linear.elliptic.ms.ddd";

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
    file << "boundaryId = 7" << std::endl;                  // a cube from the factory gets the boundary ids 1 to 4 ind 2d and 1 to 6 in 3d (hopefully)
    file << "partitions.0 = 2" << std::endl;
    file << "partitions.1 = 2" << std::endl;
    file << "partitions.2 = 2" << std::endl;
    file << "filename = " << id << ".msGrid" << std::endl;
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
    file << "dirichlet.order = 0"  << std::endl;
    file << "dirichlet.variable = x" << std::endl;
    file << "dirichlet.expression.0 = x[0]"  << std::endl;
    file << "dirichlet.expression.1 = 1.0"  << std::endl;
    file << "dirichlet.expression.2 = 1.0"  << std::endl;
    file << "[detailed.solvers.stationary.linear.elliptic.multiscale.semicontinuousgalerkin]" << std::endl;
    file << "discretization.penaltyFactor = 10.0" << std::endl;
    file << "solve.type = eigen.bicgstab.incompletelut" << std::endl;
    file << "solve.maxIter = 5000"  << std::endl;
    file << "solve.precision = 1e-12"  << std::endl;
    file << "visualize.filename = " << id << ".solution"  << std::endl;
    file << "visualize.name = solution"  << std::endl;
    file.close();
  } // only write param file if there is none
} // void ensureParamFile()

int main(int argc, char** argv)
{
  try {
    // mpi
    Dune::MPIManager::initialize(argc, argv);

    // parameter
    const std::string filename = id + ".param";
    ensureParamFile(filename);
    Dune::Stuff::Common::ExtendedParameterTree paramTree(argc, argv, filename);

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
    paramTree.assertSub(GridProviderType::id, id);
    const GridProviderType gridProvider(paramTree.sub(GridProviderType::id));
    typedef GridProviderType::MsGridType MsGridType;
    const Dune::shared_ptr< const MsGridType > msGrid = gridProvider.msGridPtr();
    info << "  took " << timer.elapsed()
         << " sec (has " << gridProvider.grid().size(0) << " elements, "
         << msGrid->size() << " subdomains)" << std::endl;
    debug.resume();
    info << "visualizing grid... " << std::flush;
    timer.reset();
    debug.suspend();
    msGrid->visualize(paramTree.sub(GridProviderType::id).get("filename", id + "_msGrid"));
    info << "done (took " << timer.elapsed() << " sek)" << std::endl;
    debug.resume();

    // model
    info << "setting up model... " << std::flush;
    debug.suspend();
    timer.reset();
    const unsigned int DUNE_UNUSED(dimDomain) = GridProviderType::dim;
    const unsigned int DUNE_UNUSED(dimRange) = 1;
    typedef GridProviderType::CoordinateType::value_type DomainFieldType;
    typedef DomainFieldType RangeFieldType;
    typedef Dune::Detailed::Solvers::Stationary::Linear::Elliptic::Model::Default< DomainFieldType, dimDomain, RangeFieldType, dimRange > ModelType;
    paramTree.assertSub(ModelType::id, id);
    const Dune::shared_ptr< const ModelType > model(new ModelType(paramTree.sub(ModelType::id)));
    typedef Dune::Stuff::Grid::BoundaryInfo::AllDirichlet BoundaryInfoType;
    const Dune::shared_ptr< const BoundaryInfoType > boundaryInfo(new BoundaryInfoType());
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;
    debug.resume();

    // solver
    info << "initializing solver";
//    info << ":" << std::endl;
    info << "... " << std::flush;
    debug.suspend();
    timer.reset();
    typedef Dune::Detailed::Solvers::Stationary::Linear::Elliptic::Multiscale::SemicontinuousGalerkin::DuneDetailedDiscretizations<
        ModelType,
        MsGridType,
        BoundaryInfoType,
        polOrder > SolverType;
    paramTree.assertSub(SolverType::id, id);
    SolverType solver(model, msGrid, boundaryInfo, paramTree.sub(SolverType::id));
    solver.init("  ", debug);
    debug.resume();
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    info << "solving";
//    info << ":" << std::endl;
    info << "... " << std::flush;
    debug.suspend();
    timer.reset();
    typedef SolverType::VectorBackendType VectorType;
    Dune::shared_ptr< std::vector< VectorType > > solution = solver.createVector();
    solver.solve(*solution, paramTree.sub(SolverType::id).sub("solve"), "  ", debug);
    debug.resume();
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    info << "postprocessing";
    info << ":" << std::endl;
//    info << "... " << std::flush;
//    debug.suspend();
//    timer.reset();
    solver.visualize(*solution,
                     paramTree.sub(SolverType::id).sub("visualize").get("filename", id + "_solution"),
                     paramTree.sub(SolverType::id).sub("visualize").get("name", "solution"),
                     "  ",
                     std::cout);
//    debug.resume();
//    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

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
