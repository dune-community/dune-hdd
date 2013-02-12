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

#include <dune/grid/part/leaf.hh>

#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/grid/provider.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/common/vector.hh>

#include <dune/detailed/solvers/stationary/linear/elliptic/model.hh>
#include <dune/detailed/solvers/stationary/linear/elliptic/cg/detailed-discretizations.hh>

#ifdef POLORDER
  const int polOrder = POLORDER;
#else
  const int polOrder = 1;
#endif

#ifdef MAX_PARAMS
  const int maxParamDim = MAX_PARMS;
#else
  const int maxParamDim = 50;
#endif

using namespace Dune::Detailed::Solvers;


std::string id(){
  return "stationary.linear.elliptic.cg.detailed_discretizations";
}


void writeDescriptionFile(const std::string filename, const std::string _id = id())
{
  std::ofstream file;
  file.open(filename);
  file << "[" << _id << "]" << std::endl;
  file << "grid = stuff.grid.provider.cube" << std::endl;
  file << "       stuff.grid.provider.gmsh" << std::endl;
  file << "boundaryinfo = stuff.grid.boundaryinfo.idbased" << std::endl;
  file << "               stuff.grid.boundaryinfo.alldirichlet" << std::endl;
  file << "               stuff.grid.boundaryinfo.allneumann" << std::endl;
  file << "model = model.stationary.linear.elliptic.default" << std::endl;
  file << "        model.stationary.linear.elliptic.thermalblock" << std::endl;
  file << "        model.stationary.linear.elliptic.parametric.separable.default" << std::endl;
  file << "        model.stationary.linear.elliptic.parametric.separable.thermalblock" << std::endl;
  file << "filename = " << _id << std::endl;
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
  file << "dirichlet.expression = [0.1*x[0]; 0.0; 0.0]"  << std::endl;
  file << "dirichlet.name       = linear dirichlet"  << std::endl;
  file << "neumann.order      = 0"  << std::endl;
  file << "neumann.variable   = x" << std::endl;
  file << "neumann.expression = [0.1; 0.0; 0.0]"  << std::endl;
  file << "neumann.name       = constant neumann"  << std::endl;
  file << "[model.stationary.linear.elliptic.thermalblock]" << std::endl;
  file << "diffusion.order = 0"  << std::endl;
  file << "diffusion.lowerLeft   = [0.0; 0.0; 0.0] # this should be a bounding box of the above selected grid!" << std::endl;
  file << "diffusion.upperRight  = [1.0; 1.0; 1.0] # this should be a bounding box of the above selected grid!" << std::endl;
  file << "diffusion.numElements = [2; 2; 2]"  << std::endl;
  file << "diffusion.components  = [1.0; 10.0; 3.0; 2.1]"  << std::endl;
  file << "force.name       = constant force"  << std::endl;
  file << "force.order      = 0"  << std::endl;
  file << "force.variable   = x" << std::endl;
  file << "force.expression = [1.0; 1.0; 1.0]"  << std::endl;
  file << "dirichlet.name       = constant dirichlet"  << std::endl;
  file << "dirichlet.order      = 0"  << std::endl;
  file << "dirichlet.variable   = x" << std::endl;
  file << "dirichlet.expression = [1.0; 1.0; 1.0]"  << std::endl;
  file << "neumann.name       = trivial neumann"  << std::endl;
  file << "neumann.order      = 0"  << std::endl;
  file << "neumann.variable   = x" << std::endl;
  file << "neumann.expression = [0.0; 0.0; 0.0]"  << std::endl;
  file << "[model.stationary.linear.elliptic.parametric.separable.default]" << std::endl;
  file << "diffusion.type = function.parametric.separable.default" << std::endl;
  file << "diffusion.name = diffusion" << std::endl;
  file << "diffusion.order = 0" << std::endl;
  file << "diffusion.component.0 = x[0]" << std::endl;
  file << "diffusion.component.1 = 2*x[0]" << std::endl;
  file << "diffusion.coefficient.0 = mu[0]" << std::endl;
  file << "diffusion.coefficient.1 = mu[1]" << std::endl;
  file << "diffusion.paramSize = 2" << std::endl;
  file << "diffusion.paramExplanation = [first diffusion param; second diffusion param]" << std::endl;
  file << "diffusion.paramMin = [0.0; 0.0]" << std::endl;
  file << "diffusion.paramMax = [1.0; 1.0]" << std::endl;
  file << "force.type = function.parametric.separable.default" << std::endl;
  file << "force.name = force" << std::endl;
  file << "force.order = 0" << std::endl;
  file << "force.component.0 = x[0]" << std::endl;
  file << "force.component.1 = 2*x[0]" << std::endl;
  file << "force.coefficient.0 = mu[0]" << std::endl;
  file << "force.coefficient.1 = mu[1]" << std::endl;
  file << "force.paramSize = 2" << std::endl;
  file << "force.paramExplanation = [first force param; second force param]" << std::endl;
  file << "force.paramMin = [0.0; 0.0]" << std::endl;
  file << "force.paramMax = [1.0; 1.0]" << std::endl;
  file << "dirichlet.type = function.parametric.separable.default" << std::endl;
  file << "dirichlet.name = neumann" << std::endl;
  file << "dirichlet.order = 0" << std::endl;
  file << "dirichlet.component.0 = x[0]" << std::endl;
  file << "dirichlet.component.1 = 2*x[0]" << std::endl;
  file << "dirichlet.coefficient.0 = mu[0]" << std::endl;
  file << "dirichlet.coefficient.1 = mu[1]" << std::endl;
  file << "dirichlet.paramSize = 2" << std::endl;
  file << "dirichlet.paramExplanation = [first dirichlet param; second dirichlet param]" << std::endl;
  file << "dirichlet.paramMin = [0.0; 0.0]" << std::endl;
  file << "dirichlet.paramMax = [1.0; 1.0]" << std::endl;
  file << "neumann.type = function.parametric.separable.default" << std::endl;
  file << "neumann.name = neumann" << std::endl;
  file << "neumann.order = 0" << std::endl;
  file << "neumann.component.0 = x[0]" << std::endl;
  file << "neumann.component.1 = 2*x[0]" << std::endl;
  file << "neumann.coefficient.0 = mu[0]" << std::endl;
  file << "neumann.coefficient.1 = mu[1]" << std::endl;
  file << "neumann.paramSize = 2" << std::endl;
  file << "neumann.paramExplanation = [first neumann param; second neumann param]" << std::endl;
  file << "neumann.paramMin = [0.0; 0.0]" << std::endl;
  file << "neumann.paramMax = [1.0; 1.0]" << std::endl;
  file << "[model.stationary.linear.elliptic.parametric.separable.thermalblock]" << std::endl;
  file << "diffusion.lowerLeft   = [0.0; 0.0; 0.0] # this should be a bounding box of the above selected grid!" << std::endl;
  file << "diffusion.upperRight  = [1.0; 1.0; 1.0] # this should be a bounding box of the above selected grid!" << std::endl;
  file << "diffusion.numElements = [2; 2; 2]"  << std::endl;
  file << "diffusion.paramMin    = [0.1; 0.1; 0.1; 0.1]" << std::endl;
  file << "diffusion.paramMax    = [10.0; 10.0; 10.0; 10.0]" << std::endl;
  file << "diffusion.order       = 0" << std::endl;
  file << "force.variable   = x" << std::endl;
  file << "force.expression = [1.0; 1.0; 1.0]"  << std::endl;
  file << "force.order      = 0" << std::endl;
  file << "force.name       = constant force" << std::endl;
  file << "dirichlet.variable   = x" << std::endl;
  file << "dirichlet.expression = [0.1*x[0]; 0.0; 0.0]"  << std::endl;
  file << "dirichlet.order      = 0"  << std::endl;
  file << "dirichlet.name       = dirichlet"  << std::endl;
  file << "neumann.variable   = x" << std::endl;
  file << "neumann.expression = [0.1; 0.0; 0.0]"  << std::endl;
  file << "neumann.order      = 0"  << std::endl;
  file << "neumann.name       = constant neumann"  << std::endl;
  file << "[parameter]" << std::endl;
  file << "fixed         = [1.0; 1.0; 1.0; 1.0]" << std::endl;
  file << "test.size = 2" << std::endl;
  file << "test.0    = [0.1; 0.1; 1.0; 1.0; 0.1; 0.1; 1.0; 1.0]" << std::endl;
  file << "test.1    = [1.0; 1.0; 0.1; 0.1; 1.0; 1.0; 0.1; 0.1]" << std::endl;
  file << "[solver.linear]" << std::endl;
  file << "type      = bicgstab.diagonal" << std::endl;
  file << "            bicgstab.ilut" << std::endl;
  file << "            bicgstab" << std::endl;
  file << "            cg.diagonal" << std::endl;
  file << "            cg" << std::endl;
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
    const std::string descriptionFilename = id() + ".description";
    typedef Dune::Stuff::Common::ExtendedParameterTree DescriptionType;
    const DescriptionType description(argc, argv, descriptionFilename);
    description.assertSub(id());
    const std::string filename = description.get(id() + ".filename", id());

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
    Dune::Stuff::Common::LogStream& info = Dune::Stuff::Common::Logger().info();
    Dune::Stuff::Common::LogStream& debug = Dune::Stuff::Common::Logger().debug();

    Dune::Timer timer;

    info << "setting up grid: " << std::endl;
    typedef Dune::Stuff::Grid::Provider::Interface<> GridProviderType;
    const GridProviderType* gridProvider
        = Dune::Stuff::Grid::Provider::create(description.get< std::string >(id() + ".grid", "stuff.grid.provider.cube"),
                                              description);
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

    info << "setting up boundaryinfo";
    const std::string boundaryInfoType = description.get< std::string >(id() + ".boundaryinfo");
    if (!debugLogging)
      info << "... ";
    else {
      info << ":" << std::endl;
      info << "  '" << boundaryInfoType << "'... ";
    }
    info << std::flush;
    timer.reset();
    typedef typename GridPartType::GridViewType GridViewType;
    typedef Dune::Stuff::Grid::BoundaryInfo::Interface< GridViewType > BoundaryInfoType;
    const Dune::shared_ptr< const BoundaryInfoType >
        boundaryInfo(Dune::Stuff::Grid::BoundaryInfo::create< GridViewType >(boundaryInfoType, description));
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    info << "visualizing grid... " << std::flush;
    timer.reset();
    gridProvider->visualize(boundaryInfoType, description, filename + ".grid");
    info << " done (took " << timer.elapsed() << " sek)" << std::endl;

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
    const int DUNE_UNUSED(dimDomain) = GridProviderType::dim;
    const int DUNE_UNUSED(dimRange) = 1;
    typedef GridProviderType::CoordinateType::value_type DomainFieldType;
    typedef double RangeFieldType;
    typedef double ParamFieldType;
    typedef Stationary::Linear::Elliptic::Model::Interface< DomainFieldType, dimDomain,
                                                            RangeFieldType, dimRange > ModelType;
    const Dune::shared_ptr< const ModelType >
        model(Stationary::Linear::Elliptic::Model::create<  DomainFieldType, dimDomain,
                                                              RangeFieldType, dimRange >(modelType, description));
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;
    if (model->parametric() && !model->separable())
      DUNE_THROW(Dune::NotImplemented,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " only implemented for nonparametric or separable parametric models!");
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
    typedef Stationary::Linear::Elliptic::CG::DetailedDiscretizations<  GridPartType,
                                                                        polOrder,
                                                                        RangeFieldType, dimRange > SolverType;
    SolverType solver(gridPart, boundaryInfo, model);
    solver.init("  ", debug);
    if (!debugLogging)
      info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    const DescriptionType& linearSolverDescription = description.sub("solver.linear");
    typedef typename SolverType::VectorType VectorType;
    const std::string linearSolverType =  linearSolverDescription.get< std::string >("type",      "bicgstab.diagonal");
    const size_t linearSolverMaxIter = linearSolverDescription.get< size_t >(        "maxIter",   5000);
    const double linearSolverPrecision = linearSolverDescription.get< double >(      "precision", 1e-12);
    if (!solver.parametric()) {
      info << "solving";
      if (!debugLogging)
        info << "... " << std::flush;
      else
        info << ":" << std::endl;
      timer.reset();
      Dune::shared_ptr< VectorType > solutionVector = solver.createAnsatzVector();
      solver.solve(solutionVector,
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
      solver.visualizeAnsatzVector(solutionVector,
                                   filename + ".solution",
                                   id() + ".solution",
                                   "  ",
                                   debug);
      if (!debugLogging)
        info << "done (took " << timer.elapsed() << " sec)" << std::endl;
    } else { // if (!solver.parametric())
      typedef typename ModelType::ParamFieldType  ParamFieldType;
      typedef typename ModelType::ParamType       ParamType;
      const size_t paramSize = model->paramSize();
      const DescriptionType& parameterDescription = description.sub("parameter");
      const size_t numTestParams = parameterDescription.get< size_t >("test.size");
      // loop over all test parameters
      for (size_t ii = 0; ii < numTestParams; ++ii) {
        const std::string iiString = Dune::Stuff::Common::toString(ii);
        const ParamType testParameter
            = parameterDescription.getDynVector< ParamFieldType >("test." + iiString, paramSize);
        // after this, testParameter is at least as long as paramSize, but it might be too long
        const ParamType mu = Dune::Stuff::Common::resize(testParameter, paramSize);
        info << "solving for parameter [" << mu << "]";
        if (!debugLogging)
          info << "... " << std::flush;
        else
          info << ":" << std::endl;
        timer.reset();
        Dune::shared_ptr< VectorType > solutionVector = solver.createAnsatzVector();
        solver.solve(solutionVector,
                     mu,
                     linearSolverType,
                     linearSolverMaxIter,
                     linearSolverPrecision,
                     "  ",
                     debug);
        if (!debugLogging)
          info << "done (took " << timer.elapsed() << " sec)" << std::endl;

        info << "writing solution for parameter [" << mu << "] to disc";
        if (!debugLogging)
          info << "... " << std::flush;
        else
          info << ":" << std::endl;
        timer.reset();
        std::stringstream name;
        name << id() << ".solution." << iiString << " (parameter [" << mu << "])";
        solver.visualizeAnsatzVector(solutionVector,
                                     filename + ".solution." + iiString,
                                     name.str(),
                                     "  ",
                                     debug);
        if (!debugLogging)
          info << "done (took " << timer.elapsed() << " sec)" << std::endl;
      } // loop over all test parameters
    } // if (!solver.parametric())

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
