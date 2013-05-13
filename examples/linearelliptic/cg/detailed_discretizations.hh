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
#include <dune/stuff/common/color.hh>
#include <dune/stuff/la/container/eigen.hh>
#include <dune/stuff/common/parameter.hh>

#include <dune/detailed-solvers/linearelliptic/solver/cg/detailed-discretizations.hh>


std::string id();

#include "../problem.hh"


void writeDescriptionFile(const std::string filename = id(), const std::string _id = id());


class DuneVector
{
public:
  typedef Dune::Stuff::LA::Container::EigenDenseVector< double > BackendType;

  DuneVector();

  DuneVector(const int size);

  DuneVector(const DuneVector* other);

  DuneVector(std::shared_ptr< BackendType > other);

  int len() const;

  double dot(const DuneVector* other) const;

  void scale(const double scalar);

  DuneVector* add(const DuneVector* other) const;

  BackendType& backend();

  const BackendType& backend() const;
private:
  std::shared_ptr< BackendType > backend_;
}; // class DuneVector


class DuneOperator
{
public:
  typedef Dune::Stuff::LA::Container::EigenRowMajorSparseMatrix< double > BackendType;

  DuneOperator(const std::shared_ptr< const BackendType >& matrix);

  DuneVector* apply(const DuneVector* vector) const;

  double apply2(const DuneVector* vectorOne, const DuneVector* vectorTwo) const;

private:
  const std::shared_ptr< const BackendType > matrix_;
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
  typedef ProblemType::DescriptionType DescriptionType;
  typedef Dune::DetailedSolvers::LinearElliptic::SolverContinuousGalerkinDD<  GridPartType, RangeFieldType,
                                                                              dimRange, 1 >                 SolverType;
  typedef typename SolverType::VectorType VectorType;

  LinearEllipticExampleCG(/*int argc, char** argv*/std::vector< std::string > arguments = std::vector< std::string >());

  bool parametric() const;

  int paramSize() const;

  DuneVector* solve() const;

  DuneVector* solve(std::vector< double > mu) const;

  std::vector< DuneOperator* > operators() const;

private:
  ProblemType problem_;
  std::shared_ptr< SolverType > solver_;
}; // class LinearEllipticExampleCG


//int run();
