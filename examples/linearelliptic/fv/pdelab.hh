#include <string>

std::string id(){
  return "examples.linearelliptic.fv.pdelab";
}


#include "../problem.hh"
#include <dune/detailed-solvers/linearelliptic/solver/fv/pdelab.hh>


int run(int argc, char** argv)
{
  try {
    // init problem
    Problem problem(argc, argv);
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
    // grid part
    const std::shared_ptr< const GridPartType > gridPart(new GridPartType(*(problem.grid())));

    info << "initializing solver";
    if (!debugLogging)
      info << "... " << std::flush;
    else
      info << ":" << std::endl;
    timer.reset();
    typedef Dune::DetailedSolvers::LinearElliptic::SolverFiniteVolumePdelab< GridPartType, RangeFieldType, 1> SolverType;
    SolverType solver(gridPart, problem.boundaryInfo(), problem.model());
    solver.init("  ", debug);
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
      std::shared_ptr< VectorType > solutionVector = solver.createVector();
      solver.solve(solutionVector,
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
      solver.visualize(solutionVector,
                       filename + ".solution",
                       id() + ".solution",
                       debug,
                       "  ");
      if (!debugLogging)
        info << "done (took " << timer.elapsed() << " sec)" << std::endl;
    }
  } catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e.what() << std::endl;
  } catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
  } catch ( ... ) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  } // try
  // if we came that far we can as well be happy about it
  return 0;
} // main
