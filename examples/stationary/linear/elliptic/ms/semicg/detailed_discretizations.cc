#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#include <vector>

#include <boost/filesystem.hpp>

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

#include <dune/detailed/solvers/stationary/linear/elliptic/model.hh>
#include <dune/detailed/solvers/stationary/linear/elliptic/ms/semicg/detailed-discretizations.hh>
//#include <dune/detailed/solvers/stationary/linear/elliptic/cg/detailed-discretizations.hh>

#ifdef POLORDER
  const int polOrder = POLORDER;
#else
  const int polOrder = 1;
#endif

namespace Stationary = Dune::Detailed::Solvers::Stationary;

std::string id()
{
  return "stationary.linear.elliptic.ms.semidg.detailed-discretizations";
}

/**
  \brief      Creates a parameter file if it does not exist.

              Nothing is done if the file already exists. If not, a parameter file will be created with all neccessary
              keys and values.
  \param[in]  filename
              (Relative) path to the file.
  **/
void writeDescriptionFile(std::string filename)
{
  // only write param file if there is none
  std::ofstream file;
  file.open(filename);
  file << "[" << id() << "]" << std::endl;
  file << "model = model.stationary.linear.elliptic.default" << std::endl;
//  file << "        model.stationary.linear.elliptic.thermalblock" << std::endl;
//  file << "        model.stationary.linear.elliptic.parametric.separable.default" << std::endl;
//  file << "        model.stationary.linear.elliptic.parametric.separable.thermalblock" << std::endl;
//  file << "exact_solution.order = 2" << std::endl;
//  file << "exact_solution.variable = x" << std::endl;
//  file << "exact_solution.expression.0 = -0.5*x[0]*x[0] + 0.5*x[0]" << std::endl;
  file << "[grid.multiscale.provider.cube]" << std::endl;
  file << "lowerLeft = [0.0; 0.0; 0.0]" << std::endl;
  file << "upperRight = [1.0; 1.0; 1.0]" << std::endl;
  file << "numElements = [4; 4; 4]" << std::endl;
  file << "boundaryId = 7 # a cube from the factory gets the boundary ids 1 to 4 ind 2d and 1 to 6 in 3d (hopefully)" << std::endl;
  file << "partitions = [2; 2; 2]" << std::endl;
  file << "[detailed.solvers.stationary.linear.elliptic.ms.semicg.detailed_discretizations]" << std::endl;
  file << "penaltyFactor = 10.0" << std::endl;
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
  file << "dirichlet.expression = [0.1*x[0]; 0.0; 0.0]"  << std::endl;
  file << "dirichlet.name       = linear dirichlet"  << std::endl;
  file << "neumann.order      = 0"  << std::endl;
  file << "neumann.variable   = x" << std::endl;
  file << "neumann.expression = [0.1; 0.0; 0.0]"  << std::endl;
  file << "neumann.name       = constant neumann"  << std::endl;
//  file << "[model.stationary.linear.elliptic.thermalblock]" << std::endl;
//  file << "diffusion.order = 0"  << std::endl;
//  file << "diffusion.lowerLeft   = [0.0; 0.0; 0.0] # this should be a bounding box of the above selected grid!" << std::endl;
//  file << "diffusion.upperRight  = [1.0; 1.0; 1.0] # this should be a bounding box of the above selected grid!" << std::endl;
//  file << "diffusion.numElements = [2; 2; 2]"  << std::endl;
//  file << "diffusion.components  = [1.0; 10.0; 3.0; 2.1]"  << std::endl;
//  file << "force.name       = constant force"  << std::endl;
//  file << "force.order      = 0"  << std::endl;
//  file << "force.variable   = x" << std::endl;
//  file << "force.expression = [1.0; 1.0; 1.0]"  << std::endl;
//  file << "dirichlet.name       = constant dirichlet"  << std::endl;
//  file << "dirichlet.order      = 0"  << std::endl;
//  file << "dirichlet.variable   = x" << std::endl;
//  file << "dirichlet.expression = [1.0; 1.0; 1.0]"  << std::endl;
//  file << "neumann.name       = trivial neumann"  << std::endl;
//  file << "neumann.order      = 0"  << std::endl;
//  file << "neumann.variable   = x" << std::endl;
//  file << "neumann.expression = [0.0; 0.0; 0.0]"  << std::endl;
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
//  file << "[model.stationary.linear.elliptic.parametric.separable.thermalblock]" << std::endl;
//  file << "diffusion.lowerLeft   = [0.0; 0.0; 0.0] # this should be a bounding box of the above selected grid!" << std::endl;
//  file << "diffusion.upperRight  = [1.0; 1.0; 1.0] # this should be a bounding box of the above selected grid!" << std::endl;
//  file << "diffusion.numElements = [2; 2; 2]"  << std::endl;
//  file << "diffusion.paramMin    = [0.1; 0.1; 0.1; 0.1]" << std::endl;
//  file << "diffusion.paramMax    = [10.0; 10.0; 10.0; 10.0]" << std::endl;
//  file << "diffusion.order       = 0" << std::endl;
//  file << "force.variable   = x" << std::endl;
//  file << "force.expression = [1.0; 1.0; 1.0]"  << std::endl;
//  file << "force.order      = 0" << std::endl;
//  file << "force.name       = constant force" << std::endl;
//  file << "dirichlet.variable   = x" << std::endl;
//  file << "dirichlet.expression = [0.1*x[0]; 0.0; 0.0]"  << std::endl;
//  file << "dirichlet.order      = 0"  << std::endl;
//  file << "dirichlet.name       = dirichlet"  << std::endl;
//  file << "neumann.variable   = x" << std::endl;
//  file << "neumann.expression = [0.1; 0.0; 0.0]"  << std::endl;
//  file << "neumann.order      = 0"  << std::endl;
//  file << "neumann.name       = constant neumann"  << std::endl;
  file << "[logging]" << std::endl;
  file << "info  = true" << std::endl;
  file << "debug = true" << std::endl;
  file << "file  = false" << std::endl;
  file.close();
} // void ensureParamFile()


//template< class MultiscaleDiscreteFunctionType,
//          class OutStreamType >
//void compute_errors(Dune::ParameterTree& paramTree,
//                    const MultiscaleDiscreteFunctionType& discreteFunction,
//                    OutStreamType& out,
//                    std::string prefix = "")
//{
//  // exact solution
//  typedef typename MultiscaleDiscreteFunctionType::LocalDiscreteFunctionType DiscreteFunctionType;
//  typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;
//  typedef typename DiscreteFunctionType::DomainFieldType DomainFieldType;
//  typedef Dune::Stuff::Function::Expression<
//      DomainFieldType,
//      DiscreteFunctionType::dimDomain,
//      RangeFieldType,
//      DiscreteFunctionType::dimRange >
//    FunctionType;
//  const FunctionType function(paramTree);
//  const unsigned int functionOrder = paramTree.get("order", 100);
//  // compute norms
//  out << prefix << " norm | exact solution | discrete solution | error (abs) | error (rel)" << std::endl;
//  out << prefix << "------+----------------+-------------------+-------------+-------------" << std::endl;
//  const RangeFieldType L2_reference_norm = Dune::Stuff::DiscreteFunction::Norm::L2(*(discreteFunction.msGrid().globalGridPart()),
//                                                                               function,
//                                                                               functionOrder);
//  out.precision(2);
//  out << prefix << "   L2 | " << std::setw(14) << std::scientific << L2_reference_norm
//                                      << " | " << std::flush;
//  const RangeFieldType L2_discrete_norm = Dune::Stuff::DiscreteFunction::Norm::L2(discreteFunction);
//  out << std::setw(17) << std::scientific << L2_discrete_norm
//                                                          << " | " << std::flush;
//  const RangeFieldType L2_difference = Dune::Stuff::DiscreteFunction::Norm::L2_difference(function, functionOrder, discreteFunction);
//  out << std::setw(11) << std::scientific << L2_difference
//                                                                        << " | " << std::flush;
//  out << std::setw(11) << std::scientific << L2_difference/L2_reference_norm << std::endl;
//  out << prefix << "------+----------------+-------------------+-------------+-------------" << std::endl;
////  const RangeFieldType h1_reference_norm = Dune::Stuff::DiscreteFunction::Norm::h1(referenceSolution);
////  out.precision(2);
////  out << prefix << "   h1 | " << std::setw(9) << std::scientific << h1_reference_norm << " | " << std::flush;
////  const RangeFieldType h1_multiscale_norm = Dune::Stuff::DiscreteFunction::Norm::h1(multiscaleSolution);
////  out << std::setw(10) << std::scientific << h1_multiscale_norm << " | " << std::flush;
////  const RangeFieldType h1_difference = Dune::Stuff::DiscreteFunction::Norm::h1_difference(referenceSolution, multiscaleSolution);
////  out << std::setw(11) << std::scientific << h1_difference << " | " << std::flush;
////  out << std::setw(11) << std::scientific << h1_difference/h1_reference_norm << std::endl;
////  out << prefix << "------+-----------+------------+-------------+-------------" << std::endl;
//} // void compute_norms(...)


int main(int argc, char** argv)
{
  try {
    // mpi
    Dune::MPIManager::initialize(argc, argv);

    // parameter
    const std::string filename = id();
    const std::string descriptionFilename = filename + ".description";
    if (!boost::filesystem::exists(descriptionFilename))
      writeDescriptionFile(descriptionFilename);
    typedef Dune::Stuff::Common::ExtendedParameterTree DescriptionType;
    const DescriptionType description(argc, argv, descriptionFilename);

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
    Dune::Stuff::Common::LogStream& info  = Dune::Stuff::Common::Logger().info();
    Dune::Stuff::Common::LogStream& debug = Dune::Stuff::Common::Logger().debug();

    // timer
    Dune::Timer timer;

    info << "setting up grid: " << std::endl;
    debug.suspend();
    typedef Dune::grid::Multiscale::Provider::Cube<> GridProviderType;
    const GridProviderType gridProvider = GridProviderType::createFromDescription(description);
    typedef GridProviderType::MsGridType MsGridType;
    const Dune::shared_ptr< const MsGridType > msGrid = gridProvider.msGrid();
    info << "  took " << timer.elapsed()
         << " sec (has " << gridProvider.grid()->size(0) << " elements, "
         << msGrid->size() << " subdomain";
    if (msGrid->size() > 1)
      info << "s";
    info << " and a width of "
         << Dune::GridWidth::calcGridWidth(*(msGrid->globalGridPart())) << ")" << std::endl;
    debug.resume();
    info << "visualizing grid... " << std::flush;
    timer.reset();
    debug.suspend();
    gridProvider.visualize(filename + ".grid");
    info << "done (took " << timer.elapsed() << " sek)" << std::endl;
    debug.resume();

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
    typedef typename MsGridType::GlobalGridViewType GlobalGridViewType;
    typedef Dune::Stuff::Grid::BoundaryInfo::Interface< GlobalGridViewType > BoundaryInfoType;
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
    if (model->parametric())
      DUNE_THROW(Dune::NotImplemented,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " only implemented for nonparametric models at the moment!");
    info << "visualizing model... " << std::flush;
    timer.reset();
    model->visualize(msGrid->globalGridPart()->gridView(), filename + ".model");
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    info << "setting up solver";
    typedef Stationary::Linear::Elliptic::MS::SemiCG::DetailedDiscretizations<  ModelType,
                                                                                MsGridType,
                                                                                BoundaryInfoType,
                                                                                polOrder > SolverType;
    if (!debugLogging)
      info << "... " << std::flush;
    else
      info << " '" << SolverType::id() << "':" << std::endl;
    timer.reset();
    const DescriptionType& discretizationDescription = description.sub(SolverType::id());
    SolverType solver(model, msGrid, boundaryInfo, discretizationDescription.get< RangeFieldType >("penaltyFactor"));
    solver.init("  ", debug);
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    info << "solving";
    if (!debugLogging)
      info << "... " << std::flush;
    else
      info << ":" << std::endl;
    const std::string linearSolverType =  discretizationDescription.get< std::string >("linearsolver.type",      "bicgstab.ilut");
    const size_t linearSolverMaxIter = discretizationDescription.get< size_t >(        "linearsolver.maxIter",   5000u);
    const double linearSolverPrecision = discretizationDescription.get< double >(      "linearsolver.precision", 1e-12);
    timer.reset();
    typedef SolverType::VectorType VectorType;
    std::vector< Dune::shared_ptr< VectorType > > solution = solver.createVector();
    solver.solve(solution,
                 linearSolverType,
                 linearSolverMaxIter,
                 linearSolverPrecision,
                 "  ",
                 debug);
    if (!debugLogging)
      info << "done (took " << timer.elapsed() << " sec)" << std::endl;

////    info << "computing detailed reference solution... " << std::flush;
////    debug.suspend();
////    timer.reset();
////    typedef Dune::Detailed::Solvers::Stationary::Linear::Elliptic::ContinuousGalerkin::DuneDetailedDiscretizations<
////        ModelType,
////        typename MsGridType::GlobalGridPartType,
////        BoundaryInfoType,
////        polOrder >
////      ReferenceSolverType;
////    ReferenceSolverType referenceSolver(model, msGrid->globalGridPart(), boundaryInfo);
////    referenceSolver.init();
////    Dune::shared_ptr< VectorType > referenceSolution = referenceSolver.createVector();
////    referenceSolver.solve(*referenceSolution, paramTree.sub(SolverType::id).sub("solve"), "  ", debug);
////    info << "done (took " << timer.elapsed() << " sec)" << std::endl;
////    debug.resume();

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

////    if (dimDomain == 1) {
////      info << "computing norms:" << std::endl;
////      debug.suspend();
////      referenceSolver.visualize(*referenceSolution,
////                                paramTree.sub(SolverType::id).sub("visualize").get("filename", id + ".solution") + "_reference",
////                                paramTree.sub(SolverType::id).sub("visualize").get("name", "solution") + "_reference",
////                                "  ",
////                                debug);

////      Dune::shared_ptr< typename SolverType::DiscreteFunctionType > discreteSolution = solver.createDiscreteFunction(*solution, "discrete_solution");
////      debug.resume();
////      compute_errors(paramTree.sub(id).sub("exact_solution"),
////                     *discreteSolution,
////                     info, "  ");
////    }
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
