#include "detailed_discretizations.hh"

std::string id()
{
  return "examples.linearelliptic.cg.dd";
}

void writeDescriptionFile(const std::string filename, const std::string _id)
{
  LinearEllipticExampleCG::ProblemType::writeDescriptionFile(filename, _id);
}

MyVector::MyVector(double s)
  : data_(s)
{}

double MyVector::data()
{
  return data_;
}

Operator::Operator()
{}

MyVector* Operator::apply(MyVector* vector)
{
  return new MyVector(2.0 * vector->data());
}

LinearEllipticExampleCG::LinearEllipticExampleCG(/*int argc, char** argv*/std::vector< std::string > arguments)
  : problem_(arguments.size(), Dune::Stuff::Common::String::vectorToMainArgs(arguments))
{
  const bool debugLogging = problem_.debugLogging();
  Dune::Stuff::Common::LogStream& info  = Dune::Stuff::Common::Logger().info();
  Dune::Stuff::Common::LogStream& debug = Dune::Stuff::Common::Logger().debug();
  Dune::Timer timer;
//    const DescriptionType& description = problem.description();
//    const std::string filename = problem.filename();
  // check
  if (problem_.model()->parametric() && !problem_.model()->affineparametric())
    DUNE_THROW(Dune::NotImplemented,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
               << " only implemented for nonparametric or affineparametric models!");
  // grid part
  const std::shared_ptr< const GridPartType > gridPart(new GridPartType(*(problem_.grid())));

  info << "initializing solver";
  if (!debugLogging)
    info << "... " << std::flush;
  else
    info << ":" << std::endl;
  timer.reset();
  solver_ = std::make_shared< SolverType >(gridPart,
                                           problem_.boundaryInfo(),
                                           problem_.model());
  solver_->init(debug, "  ");
  if (!debugLogging)
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;
} // LinearEllipticExampleCG

MyVector* LinearEllipticExampleCG::solve()
{
  return new MyVector(10);
}

//#include <boost/filesystem.hpp>

//int main(int argc, char** argv)
//{
//  const std::string descriptionFilename = id() + ".description";
//  if (!boost::filesystem::exists(descriptionFilename))
//    Problem::writeDescriptionFile(descriptionFilename);

//  return run(argc, argv);
//} // main
