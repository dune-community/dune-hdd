#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include "config.h"
#endif

#include <string>
#include <memory>
#include <vector>

#include <dune/grid/sgrid.hh>

#include <dune/stuff/common/string.hh>
#include <dune/stuff/la/container/eigen.hh>

#include <dune/detailed-solvers/linearelliptic/solver/cg/detailed-discretizations.hh>


std::string id();

#include "../problem.hh"


void writeDescriptionFile(const std::string filename, const std::string _id = id());


class ColMajorDenseMatrix
{
public:
  typedef Dune::Stuff::LA::Container::EigenDenseMatrix< double > EigenDenseMatrixType;
  typedef Dune::Stuff::LA::Container::EigenDenseVector< double > EigenDenseVectorType;

  ColMajorDenseMatrix(const EigenDenseVectorType& vector);

  std::string data() const;

  std::vector< int > shape() const;

  int len() const;

  ColMajorDenseMatrix* copy(int first, int last) const;

  ColMajorDenseMatrix* operator+(const ColMajorDenseMatrix& other) const;
};



class LinearEllipticExampleCG
{
public:
  typedef Problem< Dune::SGrid< 2, 2 > > ProblemType;
//  typedef Dune::Stuff::GridProviderInterface< GridType >  GridProviderType;
//  typedef Dune::Stuff::GridProviders< GridType >          GridProviders;
  typedef typename ProblemType::GridPartType       GridPartType;
//  typedef Dune::Stuff::GridboundaryInterface< typename GridPartType::GridViewType > GridboundariesType;
//  typedef Dune::Stuff::Gridboundaries< typename GridPartType::GridViewType >        Gridboundaries;
//  typedef GridPartType::ctype   DomainFieldType;
//  static const int DUNE_UNUSED( dimDomain) = GridProviderType::dim;
  typedef typename ProblemType::RangeFieldType                RangeFieldType;
  static const int DUNE_UNUSED(dimRange) = ProblemType::dimRange;
//  typedef Dune::DetailedSolvers::LinearElliptic::ModelInterface< DomainFieldType, dimDomain, RangeFieldType, dimRange >  ModelType;
//  typedef Dune::DetailedSolvers::LinearElliptic::Models< DomainFieldType, dimDomain, RangeFieldType, dimRange >          Models;
//  typedef Dune::Stuff::Common::ExtendedParameterTree DescriptionType;
  typedef Dune::DetailedSolvers::LinearElliptic::SolverContinuousGalerkinDD<  GridPartType, RangeFieldType,
                                                                              dimRange, 1 >                 SolverType;
  typedef typename SolverType::VectorType Vectortype;

  LinearEllipticExampleCG(/*int argc, char** argv*/std::vector< std::string > arguments);

  bool parametric() const;

  ColMajorDenseMatrix* solve() const;

private:
  ProblemType problem_;
  std::shared_ptr< SolverType > solver_;
}; // class LinearEllipticExampleCG


//int run(int argc, char** argv)
//{
//  try {

//    const DescriptionType& linearSolverDescription = description.sub("linearsolver");
//    typedef typename SolverType::VectorType VectorType;
//    const std::string linearSolverType =  linearSolverDescription.get< std::string >("type",      "bicgstab.ilut");
//    const double linearSolverPrecision = linearSolverDescription.get< double >(      "precision", 1e-12);
//    const size_t linearSolverMaxIter = linearSolverDescription.get< size_t >(        "maxIter",   5000);
//    if (!problem.model()->parametric()) {
//      info << "solving";
//      if (!debugLogging)
//        info << "... " << std::flush;
//      else
//        info << ":" << std::endl;
//      timer.reset();
//      std::shared_ptr< VectorType > solutionVector = solver.createVector();
//      solver.solve(solutionVector,
//                   linearSolverType,
//                   linearSolverPrecision,
//                   linearSolverMaxIter,
//                   debug,
//                   "  ");
//      if (!debugLogging)
//        info << "done (took " << timer.elapsed() << " sec)" << std::endl;

//      info << "writing solution to disc";
//      if (!debugLogging)
//        info << "... " << std::flush;
//      else
//        info << ":" << std::endl;
//      timer.reset();
//      solver.visualize(solutionVector,
//                       filename + ".solution",
//                       id() + ".solution",
//                       debug,
//                       "  ");
//      if (!debugLogging)
//        info << "done (took " << timer.elapsed() << " sec)" << std::endl;
//    } else { // if (!model->parametric())
//      typedef typename ModelType::ParamFieldType  ParamFieldType;
//      typedef typename ModelType::ParamType       ParamType;
//      const size_t paramSize = problem.model()->paramSize();
//      const DescriptionType& parameterDescription = description.sub("parameter");
//      const size_t numTestParams = parameterDescription.get< size_t >("test.size");
//      // loop over all test parameters
//      for (size_t ii = 0; ii < numTestParams; ++ii) {
//        const std::string iiString = Dune::Stuff::Common::toString(ii);
//        const ParamType testParameter
//            = parameterDescription.getDynVector< ParamFieldType >("test." + iiString, paramSize);
//        // after this, testParameter is at least as long as paramSize, but it might be too long
//        const ParamType mu = Dune::Stuff::Common::resize(testParameter, paramSize);
//        info << "solving for parameter [" << mu << "]";
//        if (!debugLogging)
//          info << "... " << std::flush;
//        else
//          info << ":" << std::endl;
//        timer.reset();
//        std::shared_ptr< VectorType > solutionVector = solver.createVector();
//        solver.solve(solutionVector,
//                     mu,
//                     linearSolverType,
//                     linearSolverPrecision,
//                     linearSolverMaxIter,
//                     debug,
//                     "  ");
//        if (!debugLogging)
//          info << "done (took " << timer.elapsed() << " sec)" << std::endl;

//        info << "writing solution for parameter [" << mu << "] to disc";
//        if (!debugLogging)
//          info << "... " << std::flush;
//        else
//          info << ":" << std::endl;
//        timer.reset();
//        std::stringstream name;
//        name << id() << ".solution." << iiString << " (parameter [" << mu << "])";
//        solver.visualize(solutionVector,
//                         filename + ".solution." + iiString,
//                         name.str(),
//                         debug,
//                         "  ");
//        if (!debugLogging)
//          info << "done (took " << timer.elapsed() << " sec)" << std::endl;
//      } // loop over all test parameters
//    } // if (!model->parametric())

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
