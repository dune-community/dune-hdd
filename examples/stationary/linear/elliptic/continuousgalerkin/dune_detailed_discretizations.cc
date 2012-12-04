#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#include <boost/filesystem.hpp>

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
#include <dune/stuff/function/expression.hh>
#include <dune/stuff/discretefunction/norm.hh>

#include <dune/detailed/solvers/stationary/linear/elliptic/model.hh>
#include <dune/detailed/solvers/stationary/linear/elliptic/continuousgalerkin/dune-detailed-discretizations.hh>

#ifdef POLORDER
  const int polOrder = POLORDER;
#else
  const int polOrder = 1;
#endif

const std::string id = "stationary.linear.elliptic.cg.ddd";

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
    file << "exact_solution.order = 2" << std::endl;
    file << "exact_solution.variable = x" << std::endl;
    file << "exact_solution.expression = [-0.5*x[0]*x[0] + 0.5*x[0]]" << std::endl;
    file << "[stuff.grid.provider.cube]" << std::endl;
    file << "lowerLeft = [0.0; 0.0; 0.0]" << std::endl;
    file << "upperRight = [1.0; 1.0; 1.0]" << std::endl;
    file << "numElements = 4" << std::endl;
    file << "[detailed.solvers.stationary.linear.elliptic.model.default]" << std::endl;
    file << "diffusion.order = 0"  << std::endl;
    file << "diffusion.variable = x" << std::endl;
    file << "diffusion.expression = [1.0; 1.0; 1.0]"  << std::endl;
    file << "force.order = 0"  << std::endl;
    file << "force.variable = x" << std::endl;
    file << "force.expression = [1.0; 1.0; 1.0]"  << std::endl;
    file << "dirichlet.order = 0"  << std::endl;
    file << "dirichlet.variable = x" << std::endl;
    file << "dirichlet.expression = [0.0; 0.0; 0.0]"  << std::endl;
    file << "neumann.order = 0"  << std::endl;
    file << "neumann.variable = x" << std::endl;
    file << "neumann.expression = [0.0; 0.0; 0.0]"  << std::endl;
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


int main(int argc, char** argv)
{
  try {
    // mpi
    MPIManager::initialize(argc, argv);

    // parameter
    const std::string filename = id + ".param";
    ensureParamFile(filename);
    Stuff::Common::ExtendedParameterTree paramTree(argc, argv, filename);
    paramTree.assertSub(id);

    // logger
    Stuff::Common::Logger().create(Stuff::Common::LOG_INFO |
                                   Stuff::Common::LOG_CONSOLE |
                                   Stuff::Common::LOG_DEBUG);
    Stuff::Common::LogStream& info = Stuff::Common::Logger().info();
    Stuff::Common::LogStream& debug = Stuff::Common::Logger().debug();

    // timer
    Dune::Timer timer;

    // grid
    info << "setting up grid: " << std::endl;
    debug.suspend();
    typedef Dune::Stuff::Grid::Provider::Interface<> GridProviderType;
    const std::string gridType = paramTree.get< std::string >(id + ".grid", "stuff.grid.provider.cube");
    const GridProviderType* gridProvider = Stuff::Grid::Provider::create(gridType,
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
    typedef Dune::Stuff::Grid::BoundaryInfo::AllDirichlet BoundaryInfoType;
    const Dune::shared_ptr< const BoundaryInfoType > boundaryInfo(new BoundaryInfoType());
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;
    info << "visualizing model... " << std::flush;
    timer.reset();
    model->visualize(gridPart->gridView(), paramTree.get(id + ".filename", id) + ".model");
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    // solver
    info << "initializing solver... " << std::flush;
    debug.suspend();
    timer.reset();
    typedef Stationary::Linear::Elliptic::ContinuousGalerkin::DuneDetailedDiscretizations<  ModelType,
                                                                                            GridPartType,
                                                                                            BoundaryInfoType,
                                                                                            polOrder >
        SolverType;
    SolverType solver(model, gridPart, boundaryInfo);
    solver.init("  ", debug);
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;
    debug.resume();

    info << "solving... " << std::flush;
    debug.suspend();
    timer.reset();
    const Dune::Stuff::Common::ExtendedParameterTree solverTree = paramTree.sub(solver.id());
    typedef SolverType::VectorBackendType DofVectorType;
    Dune::shared_ptr< DofVectorType > solution = solver.createVector();
    solver.solve(*solution,
                 solverTree.get("solve.type", "eigen.bicgstab.diagonal"),
                 solverTree.get("solve.maxIter", 1000),
                 solverTree.get("solve.precision", 1e-12),
                 "  ",
                 debug);
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;
    debug.resume();

    info << "postprocessing:" << std::endl;
    solver.visualize(*solution,
                     paramTree.get< std::string >(id + ".filename") + ".solution",
                     solverTree.get< std::string >("visualize.name"),
                     "  ",
                     info);

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
