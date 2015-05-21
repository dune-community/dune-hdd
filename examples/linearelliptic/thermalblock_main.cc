// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include <string>
#include <vector>

#include <boost/exception/exception.hpp>

#include <dune/grid/yaspgrid.hh>

#include "thermalblock.hh"

using namespace Dune;


int main(int /*argc*/, char** /*argv*/)
{
  try {
    typedef CgExample< YaspGrid< 3 >, GDT::ChooseSpaceBackend::pdelab, Stuff::LA::ChooseBackend::istl_sparse >
        ExampleType;
    ExampleType example;
    auto& disc = example.discretization();
    auto solution = disc.create_vector();
    disc.solve(solution, Pymor::Parameter("diffusion", {1, 1, 1, 1, 1, 1, 1, 1}));
    disc.visualize(solution, "solution", "solution");

    // if we came that far we can as well be happy about it
    return 0;
  } catch (Dune::Exception& e) {
    std::cerr << "\nDune reported error: " << e << std::endl;
    std::abort();
  } catch (boost::exception& e) {
    std::cerr << "\nboost reported error!" << std::endl;
    std::abort();
  } catch (std::exception& e) {
    std::cerr << "\nstl reported error: " << e.what() << std::endl;
    std::abort();
  } catch (...) {
    std::cerr << "\nUnknown exception thrown!" << std::endl;
    std::abort();
  } // try
} // ... main(...)
