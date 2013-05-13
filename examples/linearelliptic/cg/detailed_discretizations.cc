#include "detailed_discretizations.hh"

std::string id()
{
  return "examples.linearelliptic.cg.dd";
}

void writeDescriptionFile(const std::string filename, const std::string _id)
{
  LinearEllipticExampleCG::ProblemType::writeDescriptionFile(filename, _id);
}

DuneVector::DuneVector()
  : backend_(std::make_shared< BackendType >())
{}

DuneVector::DuneVector(const int ss)
  : backend_(std::make_shared< BackendType >(ss))
{
  assert(ss > 0);
}

DuneVector::DuneVector(const DuneVector* other)
  : backend_(std::make_shared< BackendType >(other->backend()))
{}

DuneVector::DuneVector(std::shared_ptr< BackendType > other)
  : backend_(other)
{}

int DuneVector::len() const
{
  return backend_->size();
}

double DuneVector::dot(const DuneVector* other) const
{
  return backend_->backend().transpose() * other->backend().backend();
}

void DuneVector::scale(const double scalar)
{
  backend_->backend() *= scalar;
}

DuneVector* DuneVector::add(const DuneVector* other) const
{
  std::shared_ptr< BackendType > ret = std::make_shared< BackendType >(*backend_);
  ret->backend() += other->backend().backend();
  return new DuneVector(ret);
}

DuneVector::BackendType& DuneVector::backend()
{
  return *backend_;
}

const DuneVector::BackendType& DuneVector::backend() const
{
  return *backend_;
}


DuneOperator::DuneOperator(const std::shared_ptr< const BackendType >& matrix)
  : matrix_(matrix)
{}

DuneVector* DuneOperator::apply(const DuneVector* vector) const
{
  assert(matrix_->cols() == vector->backend().size());
  DuneVector* ret = new DuneVector(vector->backend().size());
  ret->backend().backend() = matrix_->backend() * vector->backend().backend();
  return ret;
}

double DuneOperator::apply2(const DuneVector* vectorOne, const DuneVector* vectorTwo) const
{
  assert(matrix_->rows() == vectorOne->backend().size());
  assert(matrix_->cols() == vectorTwo->backend().size());
  return vectorOne->backend().backend().transpose() * matrix_->backend() * vectorTwo->backend().backend();
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

bool LinearEllipticExampleCG::parametric() const
{
  return problem_.model()->parametric();
}

int LinearEllipticExampleCG::paramSize() const
{
  return problem_.model()->paramSize();
}

DuneVector* LinearEllipticExampleCG::solve() const
{
  // check state
  if (problem_.model()->parametric())
    DUNE_THROW(Dune::InvalidStateException,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
               << " nonparametric solve() called for a parametric problem (check parametric() first)!");
  // prepare
  const bool debugLogging = problem_.debugLogging();
  Dune::Stuff::Common::LogStream& info  = Dune::Stuff::Common::Logger().info();
  Dune::Stuff::Common::LogStream& debug = Dune::Stuff::Common::Logger().debug();
  Dune::Timer timer;
  const DescriptionType& linearSolverDescription = problem_.description().sub("linearsolver");
  const std::string linearSolverType =  linearSolverDescription.get< std::string >("type",      "bicgstab.ilut");
  const double linearSolverPrecision = linearSolverDescription.get< double >(      "precision", 1e-12);
  const size_t linearSolverMaxIter = linearSolverDescription.get< size_t >(        "maxIter",   5000);
  info << "solving";
  if (!debugLogging)
    info << "... " << std::flush;
  else
    info << ":" << std::endl;
  timer.reset();
  std::shared_ptr< VectorType > solutionVector = solver_->createVector();
  solver_->solve(solutionVector,
                 linearSolverType,
                 linearSolverPrecision,
                 linearSolverMaxIter,
                 debug,
                 "  ");
  if (!debugLogging)
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;
  return new DuneVector(solutionVector);
}

DuneVector* LinearEllipticExampleCG::solve(std::vector< double > mu) const
{
  // check state
  if (!problem_.model()->parametric())
    DUNE_THROW(Dune::InvalidStateException,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
               << " parametric solve() called for a nonparametric problem (check parametric() first)!");
  // prepare
  const bool debugLogging = problem_.debugLogging();
  Dune::Stuff::Common::LogStream& info  = Dune::Stuff::Common::Logger().info();
  Dune::Stuff::Common::LogStream& debug = Dune::Stuff::Common::Logger().debug();
  Dune::Timer timer;
  const DescriptionType& linearSolverDescription = problem_.description().sub("linearsolver");
  const std::string linearSolverType =  linearSolverDescription.get< std::string >("type",      "bicgstab.ilut");
  const double linearSolverPrecision = linearSolverDescription.get< double >(      "precision", 1e-12);
  const size_t linearSolverMaxIter = linearSolverDescription.get< size_t >(        "maxIter",   5000);
  const Dune::Stuff::Common::Parameter::Type param = Dune::Stuff::Common::Parameter::create(mu);
  info << "solving";
  if (!debugLogging)
    info << "... " << std::flush;
  else
    info << ":" << std::endl;
  timer.reset();
  std::shared_ptr< VectorType > solutionVector = solver_->createVector();
  solver_->solve(solutionVector,
                 param,
                 linearSolverType,
                 linearSolverPrecision,
                 linearSolverMaxIter,
                 debug,
                 "  ");
  if (!debugLogging)
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;
  return new DuneVector(solutionVector);
}

std::vector< DuneOperator* > LinearEllipticExampleCG::operators() const
{
  std::vector< DuneOperator* > ret;
  const auto matrix = solver_->parametricMatrix();
  const auto& components = matrix->components();
  for (size_t qq = 0; qq < components.size(); ++qq)
    ret.push_back(new DuneOperator(components[qq]));
  return ret;
}

//int run()
//{
//  try {
//    LinearEllipticExampleCG example;
//    DuneVector* solution = example.solve({1.0, 1.0, 1.0, 1.0});

//  } catch(Dune::Exception& e) {
//    std::cerr << "Dune reported error: " << e.what() << std::endl;
//  } catch(std::exception& e) {
//    std::cerr << e.what() << std::endl;
//  } catch( ... ) {
//    std::cerr << "Unknown exception thrown!" << std::endl;
//  } // try

//  // if we came that far we can as well be happy about it
//  return 0;
//} // run
