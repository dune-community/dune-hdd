
#include "config.h"

// boost
#include <boost/filesystem.hpp>

// dune-common
#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>
#include <dune/common/timer.hh>

// dune-fem
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/gridpart/gridpartview.hh>

// dune-stuff
#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/grid/provider/cube.hh>

// dune-detailed-solvers
#include <dune/detailed-solvers/stationary/linear/elliptic/model.hh>
#include <dune/detailed-solvers/stationary/linear/elliptic/continuousgalerkin/dune-detailed-discretizations.hh>

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
    file << "[diffusion]" << std::endl;
    file << "variable = x" << std::endl;
    file << "expression.0 = 1.0"  << std::endl;
    file << "expression.1 = 1.0"  << std::endl;
    file << "expression.2 = 1.0"  << std::endl;
    file << "order = 0"  << std::endl;
    file << "[force]" << std::endl;
    file << "variable = x" << std::endl;
    file << "expression.0 = 1.0"  << std::endl;
    file << "expression.1 = 1.0"  << std::endl;
    file << "expression.2 = 1.0"  << std::endl;
    file << "order = 0"  << std::endl;
    file << "[solver]" << std::endl;
    file << "maxIter = 5000"  << std::endl;
    file << "precision = 1e-12"  << std::endl;
    file.close();
  } // only write param file if there is none
} // void ensureParamFile()

#ifndef POLORDER
const int polOrder = 1;
#else
const int polOrder = POLORDER;
#endif

int main(int argc, char** argv)
{
  try {
    // mpi
    Dune::MPIHelper::instance(argc, argv);

    // parameter
    const std::string id = "dune_fem";
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
    typedef Dune::LeafGridPart< GridType > GridPartType;
    GridPartType gridPart(grid);
    typedef Dune::GridPartView< GridPartType > GridViewType;
    GridViewType gridView(gridPart);
    std::cout << "took " << timer.elapsed() << " sec, has " << gridView.size(0) << " entities" << std::endl;

    // model
    const unsigned int dimDomain = GridProviderType::dim;
    const unsigned int dimRange = 1;
    typedef GridProviderType::CoordinateType::value_type DomainFieldType;
    typedef DomainFieldType RangeFieldType;
    typedef Dune::DetailedSolvers::Stationary::Linear::Elliptic::Model< DomainFieldType, dimDomain, RangeFieldType, dimRange > ModelType;
    const ModelType model(paramTree);

























    // done
    return 0;
  } catch(Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e.what() << std::endl;
  } catch(std::exception& e) {
    std::cerr << e.what() << std::endl;
  } catch( ... ) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  } // try
} // main
