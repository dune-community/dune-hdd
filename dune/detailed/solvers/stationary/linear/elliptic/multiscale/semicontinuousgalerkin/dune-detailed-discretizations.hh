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
#include <dune/detailed/discretizations/assembler/multiscale/system.hh>

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

    // system matrix and right hand side
    debug << prefix << "setting up global matrix (of size " << ansatzMapper_.size() << "x" << ansatzMapper_.size() << ") and vector container... " << std::flush;
    timer.reset();
    computeGlobalSparsityPattern();
    computeCouplingSparsityPattern();
    matrix_ = Dune::shared_ptr< MatrixBackendType >(new MatrixBackendType(ContainerFactory::createSparseMatrix(ansatzMapper_.size(), testMapper_.size(), *globalSparsityPattern_)));
    rhs_ = Dune::shared_ptr< VectorBackendType >(new VectorBackendType(ContainerFactory::createDenseVector(testMapper_.size())));
    debug << "done (took " << timer.elapsed() << " sec)" << std::endl;

    // compute coupling matrices
    for (unsigned int subdomain = 0; subdomain < msGrid_.size(); ++subdomain) {
      // get coupling sparsity pattern map
      std::map< unsigned int, LocalSparsityPatternType >& couplingSparsityPatternMap = couplingSparsityPatternMaps_[subdomain];
      // init storage
      couplingMatricesMaps_.push_back(std::map< unsigned int, MatrixBackendType >());
      std::map< unsigned int, MatrixBackendType >& couplingMatricesMap = couplingMatricesMaps_[subdomain];
      // loop over all neighboring subdomain
      typedef typename MsGridType::NeighboringSubdomainsSetType NeighborsType;
      const Dune::shared_ptr< const NeighborsType > neighbors = msGrid_.neighborsOf(subdomain);
      for (typename NeighborsType::const_iterator neighborIt = neighbors->begin();
           neighborIt != neighbors->end();
           ++neighborIt) {
        const unsigned int neighbor = *neighborIt;
        // get coupling sparsity pattern
        assert(couplingSparsityPatternMap.find(neighbor) != couplingSparsityPatternMap.end());
        LocalSparsityPatternType& couplingSparsityPattern = couplingSparsityPatternMap[neighbor];
        // create matrix
        couplingMatricesMap.insert(std::pair< unsigned int, MatrixBackendType >(
          neighbor,
          MatrixBackendType(ContainerFactory::createSparseMatrix(localSolvers_[subdomain]->ansatzMapper().size(),
                                                                 localSolvers_[neighbor]->testMapper().size(),
                                                                 couplingSparsityPattern))));
      } // loop over all neighboring subdomain
    } // compute coupling matrices

    // assemble coupling matrices
    debug << prefix << "assembling coupling matrices... " << std::flush;
    timer.reset();
    typedef Dune::Detailed::Discretizations::Assembler::Multiscale::System::Primal< MsGridType, LocalSolverType > CouplingSystemAssemblerType;
    const CouplingSystemAssemblerType couplingSystemAssembler(msGrid_, localSolvers_);
    couplingSystemAssembler.assembleCoupling(couplingMatricesMaps_);
    debug << "done (took " << timer.elapsed() << " sec)" << std::endl;

    // assemble system
    debug << prefix << "assembling system" << std::flush;
    timer.reset();
    for (unsigned int subdomain = 0; subdomain < msGrid_.size(); ++subdomain) {
      debug << "." << std::flush;
      copyLocalToGlobalMatrix(subdomain);
      copyLocalToGlobalVector(subdomain);
    }
    debug << " done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void init(Dune::ParameterTree paramTree = Dune::ParameterTree())

  void solve(std::vector< Dune::shared_ptr< LocalVectorType > >& solution, Dune::ParameterTree paramTree = Dune::ParameterTree()) const
  {
    assert(solution.size() == msGrid_.size());
    // logging
    const std::string prefix = paramTree.get("prefix", "");
//    Dune::Stuff::Common::LogStream& debug = Dune::Stuff::Common::Logger().info();
    std::ostream& debug = std::cout;

    //    // solve locally
//    for (unsigned int subdomain = 0; subdomain < msGrid_.size(); ++subdomain) {
//      localSolvers_[subdomain]->solve(solution[subdomain], paramTree);
//    }

    // preparations
    const std::string type = paramTree.get("type", "eigen.bicgstab.incompletelut");
    const unsigned int maxIter = paramTree.get("maxIter", 5000);
    const double precision = paramTree.get("precision", 1e-12);
    Dune::Timer timer;

    // create global solution vector
    VectorBackendType tmp = ContainerFactory::createDenseVector(testMapper_.size());

    // solve
    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::BicgstabIlut BicgstabIlutSolver;
    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::BicgstabDiagonal BicgstabDiagonalSolver;
    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::CgDiagonalUpper CgDiagonalUpperSolver;
    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::CgDiagonalLower CgDiagonalLowerSolver;
    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::SimplicialcholeskyUpper SimplicialcholeskyUpperSolver;
    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::SimplicialcholeskyLower SimplicialcholeskyLowerSolver;
    debug << std::endl;
    debug << prefix << "solving linear system of size " << matrix_->rows() << "x" << matrix_->cols() << std::endl
                << prefix << "using " << type << "... " << std::flush;
    if (type == "eigen.bicgstab.incompletelut"){
      BicgstabIlutSolver::apply(*matrix_, tmp, *rhs_, maxIter, precision);
    } else if (type == "eigen.bicgstab.diagonal"){
      BicgstabDiagonalSolver::apply(*matrix_, tmp, *rhs_, maxIter, precision);
    } else if (type == "eigen.cg.diagonal.upper"){
      CgDiagonalUpperSolver::apply(*matrix_, tmp, *rhs_, maxIter, precision);
    } else if (type == "eigen.cg.diagonal.lower"){
      CgDiagonalLowerSolver::apply(*matrix_, tmp, *rhs_, maxIter, precision);
    } else if (type == "eigen.simplicialcholesky.upper"){
      SimplicialcholeskyUpperSolver::apply(*matrix_, tmp, *rhs_, maxIter, precision);
    } else if (type == "eigen.simplicialcholesky.lower"){
      SimplicialcholeskyLowerSolver::apply(*matrix_, tmp, *rhs_, maxIter, precision);
    } else {
      std::stringstream msg;
      msg << "Error";
      if (id != "") {
        msg << " in " << id;
      }
      msg << ": solver type '" << type << "not supported!" << std::endl;
      DUNE_THROW(Dune::InvalidStateException, msg.str());
    } // solve
    debug << "done (took " << timer.elapsed() << " sec)" << std::endl;

    // copy global vector to local vector
    debug << prefix << "copying global vector to local vectors...  " << std::flush;
    timer.reset();
    for (unsigned int subdomain = 0; subdomain < msGrid_.size(); ++subdomain) {
      VectorBackendType localVectorBackend(solution[subdomain]);
      for (unsigned int localI = 0; localI < localVectorBackend.size(); ++localI) {
        const unsigned int globalI = testMapper_.toGlobal(subdomain, localI);
        localVectorBackend.set(localI, tmp.get(globalI));
      }
    } // copy global vector to local vector
    debug << "done (took " << timer.elapsed() << " sec)" << std::endl;
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
  void computeGlobalSparsityPattern()
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
  } // void computeGlobalSparsityPatternType()

  //! there are too many assumptions and too much happening here, the sparsity pattern should be provided by the function spaces, not sure how to do this yet
  //! \attention  computeGlobalSparsityPattern() has to be called first
  void computeCouplingSparsityPattern()
  {
    // loop over all subdomains
    for (unsigned int subdomain = 0; subdomain < msGrid_.size(); ++subdomain) {
      // get the subdomains ansatz space
      const typename LocalSolverType::AnsatzSpaceType& ansatzSpace = localSolvers_[subdomain]->ansatzSpace();
      // init storage
      couplingSparsityPatternMaps_.push_back(std::map< unsigned int, LocalSparsityPatternType >());
      std::map< unsigned int, LocalSparsityPatternType >& couplingSparsityPatternMap = couplingSparsityPatternMaps_[subdomain];
      // loop over all neighboring subdomain
      typedef typename MsGridType::NeighboringSubdomainsSetType NeighborsType;
      const Dune::shared_ptr< const NeighborsType > neighbors = msGrid_.neighborsOf(subdomain);
      for (typename NeighborsType::const_iterator neighborIt = neighbors->begin();
           neighborIt != neighbors->end();
           ++neighborIt) {
        const unsigned int neighborSubdomain = *neighborIt;
        // init storage
        couplingSparsityPatternMap.insert(std::pair< unsigned int, LocalSparsityPatternType >(neighborSubdomain, LocalSparsityPatternType()));
        LocalSparsityPatternType& couplingSparsityPattern = couplingSparsityPatternMap.find(neighborSubdomain)->second;
        // get the neighbor subdomains test space
        const typename LocalSolverType::TestSpaceType& testSpace = localSolvers_[neighborSubdomain]->testSpace();
        // get the coupling grid part
        typedef typename MsGridType::CouplingGridPartType CouplingGridPartType;
        const Dune::shared_ptr< const CouplingGridPartType > couplingGridPart = msGrid_.couplingGridPart(subdomain, neighborSubdomain);
        // loop over the coupling grid part
        for (typename CouplingGridPartType::template Codim< 0 >::IteratorType entityIt = couplingGridPart->template begin< 0 >();
             entityIt != couplingGridPart->template end< 0 >();
             ++entityIt) {
          typedef typename CouplingGridPartType::template Codim< 0 >::EntityType EntityType;
          const EntityType& entity = *entityIt;
          // loop over the neighbors
          for (typename CouplingGridPartType::IntersectionIteratorType intersectionIt = couplingGridPart->ibegin(entity);
               intersectionIt != couplingGridPart->iend(entity);
               ++intersectionIt) {
            typedef typename CouplingGridPartType::IntersectionIteratorType::Intersection IntersectionType;
            const IntersectionType& intersection = *intersectionIt;
            const typename IntersectionType::EntityPointer neighborPtr = intersection.outside();
            const EntityType& neighbor = *neighborPtr;
            // loop over all ansatz base functions of the entity
            for (unsigned int i = 0; i < ansatzSpace.baseFunctionSet().local(entity).size(); ++i) {
              // get row in coupling sparsity pattern (and create, if necessary)
              const unsigned int localI = ansatzSpace.map().toGlobal(entity, i);
              std::map< unsigned int, std::set< unsigned int > >::iterator localRow = couplingSparsityPattern.find(localI);
              if (localRow == couplingSparsityPattern.end()) {
                couplingSparsityPattern.insert(std::pair< unsigned int, std::set< unsigned int > >(localI, std::set< unsigned int >()));
                localRow = couplingSparsityPattern.find(localI);
              }
              assert(localRow != couplingSparsityPattern.end());
              // get row in global sparsity pattern (and create, if necessary)
              const unsigned int globalI = ansatzMapper_.toGlobal(subdomain, localI);
              std::map< unsigned int, std::set< unsigned int > >::iterator globalRow = globalSparsityPattern_->find(globalI);
              if (globalRow == globalSparsityPattern_->end()) {
                globalSparsityPattern_->insert(std::pair< unsigned int, std::set< unsigned int > >(globalI, std::set< unsigned int >()));
                globalRow = globalSparsityPattern_->find(globalI);
              }
              assert(globalRow != globalSparsityPattern_->end());
              // loop over all test base functions of the neighbor
              for(unsigned int j = 0; j < testSpace.baseFunctionSet().local(neighbor).size(); ++j) {
                // insert into coupling sparsity pattern
                const unsigned int localJ = testSpace.map().toGlobal(neighbor, j);
                localRow->second.insert(localJ);
                // insert into global sparsity pattern
                const unsigned int globalJ = testMapper_.toGlobal(neighborSubdomain, localJ);
                globalRow->second.insert(globalJ);
              } // loop over all test base functions
            } // loop over all ansatz base functions
          } // loop over the neighbors
        } // loop over the coupling grid part
      } // loop over all neighboring subdomain
    } // loop over all subdomains
  } // void computeCouplingSparsityPattern()

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
  std::vector< Dune::shared_ptr< const LocalSparsityPatternType > > localSparsityPatterns_;
  std::vector< std::map< unsigned int, LocalSparsityPatternType > > couplingSparsityPatternMaps_;
  Dune::shared_ptr< GlobalSparsityPatternType > globalSparsityPattern_;
  Dune::shared_ptr< MatrixBackendType > matrix_;
  Dune::shared_ptr< VectorBackendType > rhs_;
  std::vector< std::map< unsigned int, MatrixBackendType > > couplingMatricesMaps_;
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
