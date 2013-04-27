#include "detailed_discretizations.hh"

std::string id()
{
  return "examples.linearelliptic.cg.dd";
}

void writeDescriptionFile(const std::string filename, const std::string _id)
{
  LinearEllipticExampleCG::ProblemType::writeDescriptionFile(filename, _id);
}

ColMajorDenseMatrix::ColMajorDenseMatrix()
  : BaseType()
{}

ColMajorDenseMatrix::ColMajorDenseMatrix(const int _rows, const int _cols)
  : BaseType(_rows, _cols)
{
  assert(_rows > 0);
  assert(_cols > 0);
}

ColMajorDenseMatrix::ColMajorDenseMatrix(const EigenDenseVectorType& vector)
  : BaseType(vector.backend().transpose())
{}

ColMajorDenseMatrix::ColMajorDenseMatrix(const ColMajorDenseMatrix& other)
  : BaseType(other.backend())
{}

double *ColMajorDenseMatrix::data()
{
  return BaseType::backend().data();
}

std::vector< int > ColMajorDenseMatrix::shape() const
{
  return {BaseType::rows(), BaseType::cols()};
}

int ColMajorDenseMatrix::len() const
{
  return BaseType::rows();
}

ColMajorDenseMatrix ColMajorDenseMatrix::operator+(const ColMajorDenseMatrix& other) const
{
  ColMajorDenseMatrix ret(*this);
  ret.backend() += other.backend();
  return ret;
}


Operator::Operator(const std::shared_ptr< const MatrixType >& matrix)
  : matrix_(matrix)
{}

ColMajorDenseMatrix Operator::apply(const ColMajorDenseMatrix& vectors) const
{
  ColMajorDenseMatrix ret(vectors.rows(), vectors.cols());
  for (unsigned int ii = 0; ii < vectors.rows(); ++ii)
    ret.backend().row(ii) = matrix_->backend() * vectors.backend().row(ii).transpose();
  return ret;
}

ColMajorDenseMatrix Operator::apply2(const ColMajorDenseMatrix& vectorsOne, const ColMajorDenseMatrix& vectorsTwo, const bool pairwise) const
{
  if (pairwise) {
    const unsigned int numVectors = vectorsOne.len();
    assert(vectorsTwo.len() == numVectors);
    ColMajorDenseMatrix ret(numVectors, 1);
    for (unsigned int ii = 0; ii < numVectors; ++ii)
      ret.set(ii, 0, vectorsOne.backend().row(ii) * matrix_->backend() * vectorsTwo.backend().row(ii).transpose());
    return ret;
  } else {
    std::cout << "Not yet implemented!" << std::endl;
    return ColMajorDenseMatrix();
  }
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

ColMajorDenseMatrix* LinearEllipticExampleCG::solve() const
{
  // check state
  if (parametric())
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
  return new ColMajorDenseMatrix(*solutionVector);
}

Operator* LinearEllipticExampleCG::getOperator() const
{
  // check state
  if (parametric())
    DUNE_THROW(Dune::InvalidStateException,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
               << " nonparametric getOperator() called for a parametric problem (check parametric() first)!");
  return new Operator(solver_->operatorMatrix());
}
