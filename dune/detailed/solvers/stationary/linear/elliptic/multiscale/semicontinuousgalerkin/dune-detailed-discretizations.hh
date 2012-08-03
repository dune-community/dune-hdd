#ifndef DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MULTISCALE_SEMICONTINUOUSGALERKIN_DUNE_DETAILED_DISCRETIZATIONS_HH
#define DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MULTISCALE_SEMICONTINUOUSGALERKIN_DUNE_DETAILED_DISCRETIZATIONS_HH

// system
#include <vector>
#include <sstream>

// dune-common
#include <dune/common/exceptions.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/timer.hh>

// dune-fem
#include <dune/fem/space/common/functionspace.hh>

// dune-grid-multiscale
#include <dune/grid/multiscale/default.hh>

// dune-detailed-discretizations
#include <dune/detailed/discretizations/la/factory/eigen.hh>
#include <dune/detailed/discretizations/mapper/multiscale.hh>
#include <dune/detailed/discretizations/discretefunction/multiscale.hh>

// dune-detailed-solvers
#include <dune/detailed/solvers/stationary/linear/elliptic/continuousgalerkin/dune-detailed-discretizations.hh>

// dune-stuff
#include <dune/stuff/common/logging.hh>

namespace Dune {

namespace Detailed {

namespace Solvers {

namespace Stationary {

namespace Linear {

namespace Elliptic {

namespace Multiscale {

namespace SemicontinuousGalerkin {

template< class ModelImp, class MsGridImp, int polynomialOrder >
class DuneDetailedDiscretizations
{
public:
  typedef ModelImp ModelType;

  typedef MsGridImp MsGridType;

  typedef typename MsGridType::GlobalGridPartType GlobalGridPartType;

  typedef GlobalGridPartType GridPartType;

  static const int polOrder = polynomialOrder;

  typedef DuneDetailedDiscretizations< ModelType, MsGridType, polOrder > ThisType;

  static const std::string id;

private:
  typedef typename ModelType::DomainFieldType DomainFieldType;

  static const int dimDomain = ModelType::dimDomain;

  typedef typename ModelType::RangeFieldType RangeFieldType;

  static const int dimRange = ModelType::dimRange;

  typedef Dune::Detailed::Discretizations::LA::Factory::Eigen< RangeFieldType > ContainerFactory;

  typedef typename ContainerFactory::SparseMatrixType MatrixBackendType;

  typedef typename ContainerFactory::DenseVectorType VectorBackendType;

  typedef Dune::FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, dimRange > FunctionSpaceType;

  typedef typename MsGridType::LocalGridPartType LocalGridPartType;

  typedef typename Dune::Detailed::Solvers::Stationary::Linear::Elliptic::ContinuousGalerkin::DuneDetailedDiscretizations< ModelType, LocalGridPartType, polOrder > LocalSolverType;

  typedef typename LocalSolverType::SparsityPatternType LocalSparsityPatternType;

  typedef LocalSparsityPatternType GlobalSparsityPatternType;

  // should be the index type of the local problem, not unsinged int
  typedef Dune::Detailed::Discretizations::Mapper::Multiscale< unsigned int > AnsatzMapperType;

  typedef AnsatzMapperType TestMapperType;

  typedef typename MatrixBackendType::StorageType MatrixType;

public:
  typedef typename VectorBackendType::StorageType LocalVectorType;

  typedef typename LocalSolverType::DiscreteFunctionType LocalDiscreteFunctionType;

  typedef typename Dune::Detailed::Discretizations::DiscreteFunction::Multiscale< MsGridType, LocalDiscreteFunctionType > DiscreteFunctionType;

  DuneDetailedDiscretizations(const ModelType& model, const MsGridType& msGrid)
    : model_(model),
      msGrid_(msGrid),
      ansatzMapper_(),
      testMapper_()
  {}

  void init(Dune::ParameterTree paramTree = Dune::ParameterTree())
  {
    // logging
    const std::string prefix = paramTree.get("prefix", "");
//    Dune::Stuff::Common::LogStream& debug = Dune::Stuff::Common::Logger().info();
    std::ostream& debug = std::cout;

    // timer
    Dune::Timer timer;

    // create local solvers for each subdomain
    const unsigned int subdomains = msGrid_.size();
    debug << std::endl << prefix << id << ".init:" << std::endl;
    debug << prefix << "initializing " << subdomains << " local continuous galerkin solvers... " << std::flush;
    ansatzMapper_.prepare();
    testMapper_.prepare();
    for (unsigned int subdomain = 0; subdomain < subdomains; ++subdomain) {
      localGridParts_.push_back(msGrid_.localGridPart(subdomain));
      localSolvers_.push_back(Dune::shared_ptr< LocalSolverType >(new LocalSolverType(model_, *(localGridParts_[subdomain]))));
      localSolvers_[subdomain]->init();
      localSparsityPatterns_.push_back(localSolvers_[subdomain]->sparsityPattern());
      ansatzMapper_.add(subdomain, localSolvers_[subdomain]->ansatzMapper().size());
      testMapper_.add(subdomain, localSolvers_[subdomain]->testMapper().size());
    }
    ansatzMapper_.finalize();
    testMapper_.finalize();
    debug << "done (took " << timer.elapsed() << " sek)" << std::endl;

//    // system matrix and right hand side
//    debug << prefix << "setting up global matrix (of size " << ansatzMapper_.size() << "x" << ansatzMapper_.size() << ") and vector container... " << std::flush;
//    timer.reset();
//    computeGlobalSparsityPatternType();
//    matrix_ = Dune::shared_ptr< MatrixBackendType >(new MatrixBackendType(ContainerFactory::createSparseMatrix(ansatzMapper_.size(), testMapper_.size(), *globalSparsityPattern_)));
//    rhs_ = Dune::shared_ptr< VectorBackendType >(new VectorBackendType(ContainerFactory::createDenseVector(testMapper_.size())));
//    debug << "done (took " << timer.elapsed() << " sec)" << std::endl;

//    // assemble system
//    debug << prefix << "assembling system" << std::flush;
//    timer.reset();
//    for (unsigned int subdomain = 0; subdomain < msGrid_.size(); ++subdomain) {
//      debug << "." << std::flush;
//      copyLocalToGlobalMatrix(subdomain);
//      copyLocalToGlobalVector(subdomain);
//    }
//    debug << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void init(Dune::ParameterTree paramTree = Dune::ParameterTree())

  void solve(std::vector< Dune::shared_ptr< LocalVectorType > >& solution, Dune::ParameterTree paramTree = Dune::ParameterTree()) const
  {
    for (unsigned int subdomain = 0; subdomain < msGrid_.size(); ++subdomain) {
      localSolvers_[subdomain]->solve(solution[subdomain], paramTree);
    }
  } // void solve(std::vector< Dune::shared_ptr< LocalVectorType > >& solution, Dune::ParameterTree paramTree = Dune::ParameterTree()) const

  void visualize(const std::vector< Dune::shared_ptr< LocalVectorType > >& vector, Dune::ParameterTree paramTree = Dune::ParameterTree()) const
  {
    assert(vector.size() == msGrid_.size());
    // logging
    const std::string prefix = paramTree.get("prefix", "");
//    Dune::Stuff::Common::LogStream& debug = Dune::Stuff::Common::Logger().info();
    std::ostream& debug = std::cout;
    const std::string name = paramTree.get("name", id);
    const std::string filename = paramTree.get("filename", "visualization");
    Dune::Timer timer;
    debug << std::endl;
    debug << prefix << "writing '" << name << "' to '" << filename;
    if (dimDomain == 1)
      debug << ".vtp";
    else
      debug << ".vtu";
    debug << "'... " << std::flush;
    // create vector of local discrete functions
    std::vector< VectorBackendType > localVectorBackends;
    std::vector< Dune::shared_ptr< LocalDiscreteFunctionType > > localDiscreteFunctions;
    for (unsigned int subdomain = 0; subdomain < msGrid_.size(); ++subdomain) {
      localVectorBackends.push_back(vector[subdomain]);
      localDiscreteFunctions.push_back(Dune::shared_ptr< LocalDiscreteFunctionType >(new LocalDiscreteFunctionType(localSolvers_[subdomain]->ansatzSpace(), localVectorBackends[subdomain])));
    }
    // create multiscale discrete function
    Dune::shared_ptr< DiscreteFunctionType > discreteFunction(new DiscreteFunctionType(msGrid_, localDiscreteFunctions, name));
    typedef Dune::VTKWriter< typename MsGridType::GlobalGridViewType > VTKWriterType;
    VTKWriterType vtkWriter(*(msGrid_.globalGridView()));
    vtkWriter.addVertexData(discreteFunction);
    vtkWriter.write(filename);
    debug << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void visualize(const std::vector< Dune::shared_ptr< LocalVectorType > >& vector, Dune::ParameterTree paramTree = Dune::ParameterTree()) const

  std::vector< Dune::shared_ptr< LocalVectorType > > createVector() const
  {
    typename std::vector< Dune::shared_ptr< LocalVectorType > > ret;
    for (unsigned int subdomain = 0; subdomain < msGrid_.size(); ++subdomain) {
      ret.push_back(localSolvers_[subdomain]->createVector());
    }
    return ret;
  } // std::vector< Dune::shared_ptr< LocalVectorType > > createVector() const

private:
  void computeGlobalSparsityPatternType()
  {
    globalSparsityPattern_ = Dune::shared_ptr< GlobalSparsityPatternType >(new GlobalSparsityPatternType());
    // loop over all local sparsity patterns
    for (unsigned int subdomain = 0; subdomain < msGrid_.size(); ++subdomain) {
      // get local sparsity pattern
      const LocalSparsityPatternType& localSparsityPattern = *(localSparsityPatterns_[subdomain]);
      // loop over all rows in the sparsity pattern
      for (typename LocalSparsityPatternType::const_iterator rowIterator = localSparsityPattern.begin();
           rowIterator != localSparsityPattern.end();
           ++rowIterator) {
        const unsigned int localRow = rowIterator->first;
        const unsigned int globalRow = ansatzMapper_.toGlobal(subdomain, localRow);
        typename GlobalSparsityPatternType::iterator result = globalSparsityPattern_->find(globalRow);
        if (result == globalSparsityPattern_->end())
          globalSparsityPattern_->insert(typename GlobalSparsityPatternType::value_type(globalRow, typename GlobalSparsityPatternType::mapped_type()));
        result = globalSparsityPattern_->find(globalRow);
        assert(result != globalSparsityPattern_->end());
        const typename LocalSparsityPatternType::mapped_type& localRowSet = rowIterator->second;
        // loop over all cols in the current row
        for (typename LocalSparsityPatternType::mapped_type::const_iterator colIterator = localRowSet.begin();
             colIterator != localRowSet.end();
             ++colIterator) {
          const unsigned int localCol = *colIterator;
          const unsigned int globalCol = testMapper_.toGlobal(subdomain, localCol);
          result->second.insert(globalCol);
        } // loop over all cols in the current row
      } // loop over all rows in the sparsity pattern
    } // loop over all local sparsity patterns
    return;
  } // void computeGlobalSparsityPatternType()

  void copyLocalToGlobalMatrix(const unsigned int subdomain)
  {
    // get local system matrix
    const MatrixBackendType localMatrix(localSolvers_[subdomain]->getSystemMatrix());
    // get local sparsity pattern
    const LocalSparsityPatternType& localSparsityPattern = *(localSparsityPatterns_[subdomain]);
    // loop over all local rows
    for (typename LocalSparsityPatternType::const_iterator rowIterator = localSparsityPattern.begin();
         rowIterator != localSparsityPattern.end();
         ++rowIterator) {
      const unsigned int localRow = rowIterator->first;
      const unsigned int globalRow = ansatzMapper_.toGlobal(subdomain, localRow);
      const typename LocalSparsityPatternType::mapped_type& localRowSet = rowIterator->second;
      // loop over all cols in the current row
      for (typename LocalSparsityPatternType::mapped_type::const_iterator colIterator = localRowSet.begin();
           colIterator != localRowSet.end();
           ++colIterator) {
        const unsigned int localCol = *colIterator;
        const unsigned int globalCol = testMapper_.toGlobal(subdomain, localCol);
        // add entry
        matrix_->add(globalRow, globalCol, localMatrix.get(localRow, localCol));
      } // loop over all cols in the current row
    } // loop over all local rows
  } // void copyLocalToGlobal(const unsigned int subdomain)

  void copyLocalToGlobalVector(const unsigned int subdomain)
  {
    // get local rhs
    const VectorBackendType localVector(localSolvers_[subdomain]->getRightHandSide());
    // copy entries
    for (unsigned int localI = 0; localI < localVector.size(); ++localI) {
      const unsigned int globalI = testMapper_.toGlobal(subdomain, localI);
      rhs_->set(globalI, localVector.get(localI));
    } // copy entries
  } // void copyLocalToGlobalVector(const unsigned int subdomain)

  const ModelType& model_;
  const MsGridType& msGrid_;
  AnsatzMapperType ansatzMapper_;
  TestMapperType testMapper_;
  std::vector< Dune::shared_ptr< const LocalGridPartType > > localGridParts_;
  std::vector< Dune::shared_ptr< LocalSolverType > > localSolvers_;
  std::vector< Dune::shared_ptr< LocalSparsityPatternType > > localSparsityPatterns_;
  Dune::shared_ptr< GlobalSparsityPatternType > globalSparsityPattern_;
  Dune::shared_ptr< MatrixBackendType > matrix_;
  Dune::shared_ptr< VectorBackendType > rhs_;
}; // class DuneDetailedDiscretizations

template< class ModelType, class GridPartType, int polOrder >
const std::string DuneDetailedDiscretizations< ModelType, GridPartType, polOrder >::id = "detailed.solvers.stationary.linear.elliptic.multiscale.semicontinuousgalerkin";

} // namespace SemicontinuousGalerkin

} // namespace Multiscale

} // namespace Elliptic

} // namespace Linear

} // namespace Stationary

} // namespace Solvers

} // namespace Detailed

} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MULTISCALE_SEMICONTINUOUSGALERKIN_DUNE_DETAILED_DISCRETIZATIONS_HH
