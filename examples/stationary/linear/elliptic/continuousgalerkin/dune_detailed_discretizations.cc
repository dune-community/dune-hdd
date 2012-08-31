
#include "config.h"

// boost
#include <boost/filesystem.hpp>

// dune-common
#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>
#include <dune/common/timer.hh>

// dune-grid-multiscale
#include <dune/grid/part/leaf.hh>

// dune-stuff
#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/grid/provider/cube.hh>

// dune-detailed-solvers
#include <dune/detailed/solvers/stationary/linear/elliptic/model.hh>
#include <dune/detailed/solvers/stationary/linear/elliptic/continuousgalerkin/dune-detailed-discretizations.hh>

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
    file << "[stuff.grid.provider.cube]" << std::endl;
    file << "level = 4" << std::endl;
    file << "[detailed.solvers.stationary.linear.elliptic.model.default]" << std::endl;
    file << "diffusion.order = 0"  << std::endl;
    file << "diffusion.variable = x" << std::endl;
    file << "diffusion.expression.0 = 1.0"  << std::endl;
    file << "diffusion.expression.1 = 1.0"  << std::endl;
    file << "diffusion.expression.2 = 1.0"  << std::endl;
    file << "force.order = 0"  << std::endl;
    file << "force.variable = x" << std::endl;
    file << "force.expression.0 = 1.0"  << std::endl;
    file << "force.expression.1 = 1.0"  << std::endl;
    file << "force.expression.2 = 1.0"  << std::endl;
    file << "dirichlet.order = 0"  << std::endl;
    file << "dirichlet.variable = x" << std::endl;
    file << "dirichlet.expression.0 = 0.0"  << std::endl;
    file << "dirichlet.expression.1 = 0.0"  << std::endl;
    file << "dirichlet.expression.2 = 0.0"  << std::endl;
    file << "[detailed.solvers.stationary.linear.elliptic.continuousgalerkin]" << std::endl;
    file << "init.verbose = true" << std::endl;
    file << "solve.verbose = true" << std::endl;
    file << "solve.type = eigen.bicgstab.incompletelut" << std::endl;
    file << "solve.maxIter = 5000"  << std::endl;
    file << "solve.precision = 1e-12"  << std::endl;
    file << "visualize.verbose = true"  << std::endl;
    file << "visualize.name = solution"  << std::endl;
    file.close();
  } // only write param file if there is none
} // void ensureParamFile()

#ifdef POLORDER
const int polOrder = POLORDER;
#else
const int polOrder = 1;
#endif

int main(int argc, char** argv)
{
  try {
    // mpi
    Dune::MPIHelper::instance(argc, argv);

    // parameter
    const std::string id = "dune_detailed_discretizations";
    const std::string filename = id + ".param";
    ensureParamFile(filename);
    Dune::ParameterTree paramTree = Dune::Stuff::Common::Parameter::Tree::init(argc, argv, filename);

    // timer
    Dune::Timer timer;

    // grid
    std::cout << "setting up grid:" << std::endl;
    typedef Dune::Stuff::Grid::Provider::UnitCube< Dune::GridSelector::GridType > GridProviderType;
    Dune::Stuff::Common::Parameter::Tree::assertSub(paramTree, GridProviderType::id, filename);
    GridProviderType gridProvider(paramTree.sub(GridProviderType::id));
    typedef GridProviderType::GridType GridType;
    GridType& grid = gridProvider.grid();
    typedef Dune::grid::Part::Leaf::Const< GridType > GridPartType;
    const GridPartType gridPart(grid);
    std::cout << "  took " << timer.elapsed() << " sec, has " << gridPart.grid().size(0) << " elements" << std::endl;

    // model
    std::cout << "setting up model... " << std::flush;
    timer.reset();
    const unsigned int DUNE_UNUSED(dimDomain) = GridProviderType::dim;
    const unsigned int DUNE_UNUSED(dimRange) = 1;
    typedef GridProviderType::CoordinateType::value_type DomainFieldType;
    typedef DomainFieldType RangeFieldType;
    typedef Dune::Detailed::Solvers::Stationary::Linear::Elliptic::Model::Default< DomainFieldType, dimDomain, RangeFieldType, dimRange > ModelType;
    Dune::Stuff::Common::Parameter::Tree::assertSub(paramTree, ModelType::id, id);
    const ModelType model(paramTree.sub(ModelType::id));
    std::cout << "done (took " << timer.elapsed() << " sec)" << std::endl;

    // solver
    std::cout << "initializing solver";
    typedef Dune::Detailed::Solvers::Stationary::Linear::Elliptic::ContinuousGalerkin::DuneDetailedDiscretizations< ModelType, GridPartType, polOrder > SolverType;
    if (paramTree.sub(SolverType::id).get("init.verbose", false))
      std::cout << ":" << std::endl;
    else
      std::cout << "... " << std::flush;
    timer.reset();
    Dune::Stuff::Common::Parameter::Tree::assertSub(paramTree, SolverType::id, id);
    paramTree.sub(SolverType::id)["init.prefix"] = "  ";
    SolverType solver(model, gridPart);
    solver.init(paramTree.sub(SolverType::id).sub("init"));
    if (!paramTree.sub(SolverType::id).get("init.verbose", false))
      std::cout << "done (took " << timer.elapsed() << " sec)" << std::endl;

    std::cout << "solving";
    if (paramTree.sub(SolverType::id).get("solve.verbose", false))
      std::cout << ":" << std::endl;
    else
      std::cout << "... " << std::flush;
    timer.reset();
    typedef SolverType::VectorType DofVectorType;
    Dune::shared_ptr< DofVectorType > solution = solver.createVector();
    paramTree.sub(SolverType::id)["solve.prefix"] = "  ";
    solver.solve(solution, paramTree.sub(SolverType::id).sub("solve"));
    if (!paramTree.sub(SolverType::id).get("solve.verbose", false))
      std::cout << "done (took " << timer.elapsed() << " sec)" << std::endl;

    std::cout << "postprocessing";
    if (paramTree.sub(SolverType::id).get("visualize.verbose", false))
      std::cout << ":" << std::endl;
    else
      std::cout << "... " << std::flush;
    timer.reset();
    paramTree.sub(SolverType::id)["visualize.prefix"] = "  ";
    paramTree.sub(SolverType::id)["visualize.filename"] = id + "_solution";
    solver.visualize(solution, paramTree.sub(SolverType::id).sub("visualize"));
    if (!paramTree.sub(SolverType::id).get("visualize.verbose", false))
      std::cout << "done (took " << timer.elapsed() << " sec)" << std::endl;

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
