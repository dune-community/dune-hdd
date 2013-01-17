#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wextra"
#pragma GCC diagnostic ignored "-Wignored-qualifiers"
#pragma GCC diagnostic ignored "-Wtype-limits"
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
  #ifdef HAVE_CMAKE_CONFIG
    #include "cmake_config.h"
  #elif defined (HAVE_CONFIG_H)
    #include "config.h"
  #endif // ifdef HAVE_CMAKE_CONFIG

  #include <dune/common/exceptions.hh>
  #include <dune/common/mpihelper.hh>
  #include <dune/common/timer.hh>

  #include <dune/fem/misc/mpimanager.hh>
  #include <dune/fem/misc/gridwidth.hh>
#pragma GCC diagnostic push

#include <dune/grid/part/leaf.hh>

#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/grid/provider.hh>
#include <dune/stuff/grid/boundaryinfo.hh>

#include <dune/detailed/solvers/stationary/linear/elliptic/model.hh>
#include <dune/detailed/solvers/stationary/linear/elliptic/cg/detailed-discretizations.hh>

#ifdef POLORDER
  const int polOrder = POLORDER;
#else
  const int polOrder = 1;
#endif

const std::string id = "stationary.linear.elliptic.cg.detailed_discretizations";

void writeParamFile(const std::string filename)
{
  std::ofstream file;
  file.open(filename);
  file << "[" << id << "]" << std::endl;
  file << "grid = stuff.grid.provider.cube" << std::endl;
  file << "       stuff.grid.provider.gmsh" << std::endl;
  file << "boundaryinfo = stuff.grid.boundaryinfo.idbased" << std::endl;
  file << "               stuff.grid.boundaryinfo.alldirichlet" << std::endl;
  file << "               stuff.grid.boundaryinfo.allneumann" << std::endl;
  file << "model = detailed.solvers.stationary.linear.elliptic.model.default" << std::endl;
  file << "        detailed.solvers.stationary.linear.elliptic.model.thermalblock" << std::endl;
  file << "filename = " << id << std::endl;
  file << "[logging]" << std::endl;
  file << "info  = true" << std::endl;
  file << "debug = false" << std::endl;
  file << "file  = false" << std::endl;
  file << "[stuff.grid.provider.cube]" << std::endl;
  file << "lowerLeft   = [0.0; 0.0; 0.0]" << std::endl;
  file << "upperRight  = [1.0; 1.0; 1.0]" << std::endl;
  file << "numElements = 32" << std::endl;
  file << "[stuff.grid.provider.gmsh]" << std::endl;
  file << "filename = path_to_g.msh" << std::endl;
  file << "[stuff.grid.boundaryinfo.idbased]" << std::endl;
  file << "dirichlet = [1; 2; 3]" << std::endl;
  file << "neumann   = [4]" << std::endl;
  file << "[detailed.solvers.stationary.linear.elliptic.model.default]" << std::endl;
  file << "diffusion.order      = 0"  << std::endl;
  file << "diffusion.variable   = x" << std::endl;
  file << "diffusion.expression = [1.0; 1.0; 1.0]"  << std::endl;
  file << "force.order      = 0"  << std::endl;
  file << "force.variable   = x" << std::endl;
  file << "force.expression = [1.0; 1.0; 1.0]"  << std::endl;
  file << "dirichlet.order      = 1"  << std::endl;
  file << "dirichlet.variable   = x" << std::endl;
  file << "dirichlet.expression = [0.1*x[0]; 0.0; 0.0]"  << std::endl;
  file << "neumann.order      = 0"  << std::endl;
  file << "neumann.variable   = x" << std::endl;
  file << "neumann.expression = [0.1; 0.0; 0.0]"  << std::endl;
  file << "[detailed.solvers.stationary.linear.elliptic.model.thermalblock]" << std::endl;
  file << "diffusion.order = 0"  << std::endl;
  file << "diffusion.lowerLeft   = [0.0; 0.0; 0.0]" << std::endl; // should coincide with the grid
  file << "diffusion.upperRight  = [1.0; 1.0; 1.0]" << std::endl; // should coincide with the grid
  file << "diffusion.numElements = [2; 2; 2]"  << std::endl;
  file << "diffusion.components  = [1.0; 10.0; 3.0; 2.1]"  << std::endl;
  file << "force.order      = 0"  << std::endl;
  file << "force.variable   = x" << std::endl;
  file << "force.expression = [1.0; 1.0; 1.0]"  << std::endl;
  file << "dirichlet.order      = 0"  << std::endl;
  file << "dirichlet.variable   = x" << std::endl;
  file << "dirichlet.expression = [1.0; 1.0; 1.0]"  << std::endl;
  file << "neumann.order      = 0"  << std::endl;
  file << "neumann.variable   = x" << std::endl;
  file << "neumann.expression = [0.0; 0.0; 0.0]"  << std::endl;
  file << "[solver.linear]" << std::endl;
  file << "type      = eigen.iterative.bicgstab.diagonal" << std::endl;
  file << "            eigen.iterative.bicgstab.incompletelut" << std::endl;
  file << "            eigen.iterative.cg.diagonal" << std::endl;
  file << "            eigen.iterative.cg.incompletelut" << std::endl;
  file << "maxIter   = 5000"  << std::endl;
  file << "precision = 1e-12"  << std::endl;
  file.close();
} // void writeParamFile(const std::string filename)


int run(int argc, char** argv)
{
  try {
    // mpi
    Dune::MPIManager::initialize(argc, argv);

    // parameter
    const std::string paramFilename = id + ".param";
    typedef Dune::Stuff::Common::ExtendedParameterTree ParamTreeType;
    const ParamTreeType paramTree(argc, argv, paramFilename);
    paramTree.assertSub(id);
    const std::string filename = paramTree.get(id + ".filename", id);

    // logger
    const ParamTreeType& logParams = paramTree.sub("logging");
    int logFlags = Dune::Stuff::Common::LOG_CONSOLE;
    const bool debugLogging = logParams.get< bool >("debug", false);
    if (logParams.get< bool >("info"))
      logFlags = logFlags | Dune::Stuff::Common::LOG_INFO;
    if (debugLogging)
      logFlags = logFlags | Dune::Stuff::Common::LOG_DEBUG;
    if (logParams.get< bool >("file", false))
      logFlags = logFlags | Dune::Stuff::Common::LOG_FILE;
    Dune::Stuff::Common::Logger().create(logFlags, id, "", "");
    Dune::Stuff::Common::LogStream& info = Dune::Stuff::Common::Logger().info();
    Dune::Stuff::Common::LogStream& debug = Dune::Stuff::Common::Logger().debug();

    Dune::Timer timer;

    info << "setting up grid: " << std::endl;
    typedef Dune::Stuff::Grid::Provider::Interface<> GridProviderType;
    const GridProviderType* gridProvider
        = Dune::Stuff::Grid::Provider::create(paramTree.get< std::string >(id + ".grid", "stuff.grid.provider.cube"),
                                              paramTree);
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
    gridProvider->visualize(filename + ".grid");
    info << " done (took " << timer.elapsed() << " sek)" << std::endl;

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
        ::Model::Interface< DomainFieldType, dimDomain, RangeFieldType, dimRange >ModelType;
    const Dune::shared_ptr< const ModelType > model(Dune::Detailed::Solvers
        ::Stationary::Linear::Elliptic::Model::create< DomainFieldType, dimDomain, RangeFieldType, dimRange >(
                                                      paramTree.get< std::string >(id + ".model"),
                                                      paramTree));
    typedef Dune::Stuff::Grid::BoundaryInfo::Interface< typename GridPartType::GridViewType > BoundaryInfoType;
    const Dune::shared_ptr< const BoundaryInfoType > boundaryInfo(
          Dune::Stuff::Grid::BoundaryInfo::create< typename GridPartType::GridViewType >(
            paramTree.get< std::string >(id + ".boundaryinfo"),
            paramTree));
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;
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
    typedef Dune::Detailed::Solvers
        ::Stationary
        ::Linear
        ::Elliptic
        ::CG::DetailedDiscretizations< GridPartType, polOrder, RangeFieldType, 1 >
      SolverType;
    SolverType solver(gridPart, boundaryInfo, model);
    solver.init("  ", debug);
    if (!debugLogging)
      info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    info << "solving";
    if (!debugLogging)
      info << "... " << std::flush;
    else
      info << ":" << std::endl;
    timer.reset();
    const ParamTreeType& linearSolverParams = paramTree.sub("solver.linear");
    typedef typename SolverType::VectorType VectorType;
    Dune::shared_ptr< VectorType > solutionVector = solver.createAnsatzVector();
    solver.solve(*solutionVector,
                 linearSolverParams.get< std::string >( "type",      "eigen.iterative.bicgstab.diagonal"),
                 linearSolverParams.get< unsigned int >("maxIter",   5000),
                 linearSolverParams.get< double >(      "precision", 1e-12),
                 "  ",
                 debug);
    if (!debugLogging)
      info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    info << "postprocessing";
    if (!debugLogging)
      info << "... " << std::flush;
    else
      info << ":" << std::endl;
    timer.reset();
    solver.visualizeAnsatzVector(*solutionVector,
                                 filename + ".solution",
                                 id + ".solution",
                                 "  ",
                                 debug);
    if (!debugLogging)
      info << "done (took " << timer.elapsed() << " sec)" << std::endl;

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
