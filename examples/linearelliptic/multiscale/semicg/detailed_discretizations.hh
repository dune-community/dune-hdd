// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <string>

std::string id(){
  return "examples.linearelliptic.ms.semicg.dd";
}


#include "../../problem.hh"
#include <dune/detailed-solvers/linearelliptic/solver/multiscale/semicg/detailed-discretizations.hh>


namespace Stuff = Dune::Stuff;
namespace Elliptic = Dune::DetailedSolvers::LinearElliptic;


// define the 'analytical' types
typedef MultiscaleProblem<>           ProblemType;
typedef ProblemType::DescriptionType  DescriptionType;
typedef ProblemType::GridType         GridType;
typedef ProblemType::RangeFieldType   RangeFieldType;
static const int DUNE_UNUSED(         dimDomain) = ProblemType::dimDomain;
static const int DUNE_UNUSED(         dimRange) = ProblemType::dimRange;
typedef ProblemType::ModelType        ModelType;


int run(int argc, char** argv)
{
  try {
    // init problem
    ProblemType problem(argc, argv);
    const bool debugLogging = problem.debugLogging();
    Stuff::Common::LogStream& info  = Stuff::Common::Logger().info();
    Stuff::Common::LogStream& debug = Stuff::Common::Logger().debug();
    Dune::Timer timer;
    const DescriptionType& description = problem.description();
    const std::string filename = problem.filename();
    // check
    if (problem.model()->parametric() && !problem.model()->affineparametric())
      DUNE_THROW(Dune::NotImplemented,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " only implemented for nonparametric or affineparametric models!");

    info << "initializing solver";
    if (!debugLogging)
      info << "... " << std::flush;
    else
      info << ":" << std::endl;
    timer.reset();
    typedef Elliptic::MultiscaleSolverSemiContinuousGalerkinDD< GridType, RangeFieldType, dimRange, 1 > SolverType;
    SolverType solver(problem.msGrid(),
                      problem.boundaryInfo(),
                      problem.model(),
                      description.get(SolverType::id() + ".penaltyParameter", 100.0));
    solver.init(debug, "  ");
    if (!debugLogging)
      info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    const DescriptionType& linearSolverDescription = description.sub("linearsolver");
    typedef typename SolverType::VectorType VectorType;
    const std::string linearSolverType =  linearSolverDescription.get< std::string >("type",      "bicgstab.ilut");
    const double linearSolverPrecision = linearSolverDescription.get< double >(      "precision", 1e-12);
    const size_t linearSolverMaxIter = linearSolverDescription.get< size_t >(        "maxIter",   5000);
    if (!problem.model()->parametric()) {
      info << "solving";
      if (!debugLogging)
        info << "... " << std::flush;
      else
        info << ":" << std::endl;
      timer.reset();
      std::vector< std::shared_ptr< VectorType > > solution = solver.createVector();
      solver.solve(solution,
                   linearSolverType,
                   linearSolverPrecision,
                   linearSolverMaxIter,
                   debug,
                   "  ");
      if (!debugLogging)
        info << "done (took " << timer.elapsed() << " sec)" << std::endl;

      info << "writing solution to disc";
      if (!debugLogging)
        info << "... " << std::flush;
      else
        info << ":" << std::endl;
      timer.reset();
      solver.visualize(solution,
                       filename + ".solution",
                       id() + ".solution",
                       debug,
                       "  ");
      if (!debugLogging)
        info << "done (took " << timer.elapsed() << " sec)" << std::endl;
    } else { // if (!model->parametric())
      typedef typename ModelType::ParamFieldType  ParamFieldType;
      typedef typename ModelType::ParamType       ParamType;
      const size_t paramSize = problem.model()->paramSize();
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
        std::vector< std::shared_ptr< VectorType > > solution = solver.createVector();
        solver.solve(solution,
                     mu,
                     linearSolverType,
                     linearSolverPrecision,
                     linearSolverMaxIter,
                     debug,
                     "  ");
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
        solver.visualize(solution,
                         filename + ".solution." + iiString,
                         name.str(),
                         debug,
                         "  ");
        if (!debugLogging)
          info << "done (took " << timer.elapsed() << " sec)" << std::endl;
      } // loop over all test parameters
    } // if (!model->parametric())

  } catch(Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e.what() << std::endl;
  } catch(std::exception& e) {
    std::cerr << e.what() << std::endl;
  } catch(...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  } // try

  // if we came that far we can as well be happy about it
  return 0;
} // run
