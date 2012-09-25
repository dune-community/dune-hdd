#ifndef DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MULTISCALE_SEMICONTINUOUSGALERKIN_DUNE_DETAILED_DISCRETIZATIONS_HH
#define DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MULTISCALE_SEMICONTINUOUSGALERKIN_DUNE_DETAILED_DISCRETIZATIONS_HH

// system
#include <vector>
#include <sstream>

// dune-common
#include <dune/common/exceptions.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/timer.hh>

// dune-grid-multiscale
#include <dune/grid/multiscale/default.hh>

// dune-detailed-discretizations
#include <dune/detailed/discretizations/mapper/multiscale.hh>
#include <dune/detailed/discretizations/la/factory/eigen.hh>
#include <dune/detailed/discretizations/discretefunction/multiscale.hh>

// dune-detailed-solvers
#include <dune/detailed/solvers/stationary/linear/elliptic/continuousgalerkin/dune-detailed-discretizations.hh>
#include <dune/detailed/solvers/stationary/linear/elliptic/coupling/primal/dune-detailed-discretizations.hh>
#include <dune/detailed/solvers/stationary/linear/elliptic/boundary/dune-detailed-discretizations.hh>

// dune-stuff
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/discretefunction/projection/dirichlet.hh>

namespace Dune {

namespace Detailed {

namespace Solvers {

namespace Stationary {

namespace Linear {

namespace Elliptic {

namespace Multiscale {

namespace SemicontinuousGalerkin {

template< class ModelImp, class MsGridImp, class BoundaryInfoImp, int polynomialOrder >
class DuneDetailedDiscretizations
{
public:
  typedef ModelImp ModelType;

  typedef MsGridImp MsGridType;

  typedef BoundaryInfoImp BoundaryInfoType;

  static const int polOrder = polynomialOrder;

  typedef DuneDetailedDiscretizations< ModelType, MsGridType, BoundaryInfoType, polOrder > ThisType;

  static const std::string id;

  typedef typename MsGridType::GlobalGridPartType GlobalGridPartType;

  typedef typename MsGridType::LocalGridPartType LocalGridPartType;

private:
  typedef typename ModelType::DomainFieldType DomainFieldType;

  static const int dimDomain = ModelType::dimDomain;

  typedef typename ModelType::RangeFieldType RangeFieldType;

  static const int dimRange = ModelType::dimRange;

  typedef Dune::Detailed::Discretizations::LA::Factory::Eigen< RangeFieldType > ContainerFactory;

  typedef typename ContainerFactory::SparseMatrixType MatrixBackendType;

public:
  typedef typename ContainerFactory::DenseVectorType VectorBackendType;

private:
  typedef typename MsGridType::CouplingGridPartType CouplingGridPartType;

  typedef typename MsGridType::BoundaryGridPartType BoundaryGridPartType;

  typedef Dune::Stuff::Grid::BoundaryInfo::AllNeumann LocalBoundaryInfoType;

  typedef Dune::Detailed::Solvers::Stationary::Linear::Elliptic::ContinuousGalerkin::DuneDetailedDiscretizations<
      ModelType,
      LocalGridPartType,
      LocalBoundaryInfoType,
      polOrder >
    LocalSolverType;

  typedef Dune::Detailed::Solvers::Stationary::Linear::Elliptic::Coupling::Primal::DuneDetailedDiscretizations<
      CouplingGridPartType,
      LocalSolverType,
      LocalSolverType,
      ContainerFactory >
    CouplingSolverType;

  typedef Dune::Detailed::Solvers::Stationary::Linear::Elliptic::Boundary::DuneDetailedDiscretizations<
      BoundaryGridPartType,
      BoundaryInfoType,
      LocalSolverType,
      ContainerFactory >
    BoundarySolverType;

  typedef typename LocalSolverType::AnsatzSpaceType LocalAnsatzSpaceType;

  typedef typename LocalSolverType::TestSpaceType LocalTestSpaceType;

  typedef typename LocalAnsatzSpaceType::PatternType PatternType;

  typedef Dune::Detailed::Discretizations::Mapper::Multiscale< typename LocalAnsatzSpaceType::MapperType::IndexType > AnsatzMapperType;

  typedef Dune::Detailed::Discretizations::Mapper::Multiscale< typename LocalTestSpaceType::MapperType::IndexType > TestMapperType;

public:
  typedef typename LocalSolverType::DiscreteFunctionType LocalDiscreteFunctionType;

  typedef typename Dune::Detailed::Discretizations::DiscreteFunction::Multiscale< MsGridType, LocalDiscreteFunctionType > DiscreteFunctionType;

  DuneDetailedDiscretizations(const Dune::shared_ptr< const ModelType > model,
                              const Dune::shared_ptr< const MsGridType > msGrid,
                              const Dune::shared_ptr< const BoundaryInfoType > boundaryInfo,
                              const Dune::ParameterTree& paramTree)
    : model_(model)
    , msGrid_(msGrid)
    , boundaryInfo_(boundaryInfo)
    , paramTree_(paramTree)
    , initialized_(false)
    , penaltyFactor_(-1)
    , ansatzMapper_()
    , testMapper_()
    , localSolvers_(msGrid_->size())
    , pattern_(Dune::shared_ptr< PatternType >(new PatternType()))
  {
    // get penalty factor
    std::string key = "discretization.penaltyFactor";
    Dune::Stuff::Common::Parameter::Tree::assertKey(paramTree, key, id);
    penaltyFactor_ = paramTree.get(key, RangeFieldType(-1.0));
    if (!penaltyFactor_ > 0) {
      std::stringstream msg;
      msg << "Error in " << id << ":"
          << "wrong '" << key << "' given (should be posititve double, is '" << penaltyFactor_ << "') in the following Dune::ParameterTree:" << std::endl;
      paramTree.report(msg);
      DUNE_THROW(Dune::InvalidStateException, msg.str());
    }
    // test orders
    assert(model_->diffusionOrder() >= 0 && "Please provide a nonnegative order for the diffusion!");
    assert(model_->forceOrder() >= 0 && "Please provide a nonnegative order for the force!");
  } // DuneDetailedDiscretizations()

  void init(const std::string prefix = "", std::ostream& out = Dune::Stuff::Common::Logger().debug())
  {
    if (!initialized_) {
      // prepare
      Dune::Timer timer;
      const unsigned int subdomains = msGrid_->size();
      const bool verbose = (subdomains >= 3) && (subdomains <= std::pow(3, 3));
      std::ostream& devnull = Dune::Stuff::Common::Logger().devnull();
      std::vector< std::map< unsigned int, Dune::shared_ptr< CouplingSolverType > > > couplingSolverMaps(msGrid_->size());
      std::vector< Dune::shared_ptr< BoundarySolverType > > boundarySolvers(msGrid_->size());

      // walk the subdomains for the first time
      //   * to initialize the local solvers and
      //   * to build up the multiscale mappers
      out << prefix << "initializing local solvers (on " << subdomains << " subdomains)"  << std::flush;
      if (!verbose)
        out << "..." << std::flush;
      ansatzMapper_.prepare();
      testMapper_.prepare();
      Dune::shared_ptr< const LocalBoundaryInfoType > localBoundaryInfo(new LocalBoundaryInfoType());
      for (unsigned int subdomain = 0; subdomain < msGrid_->size(); ++subdomain) {
        // initialize local solver
        localSolvers_[subdomain] = Dune::shared_ptr< LocalSolverType >(
              new LocalSolverType(model_, msGrid_->localGridPart(subdomain), localBoundaryInfo));
        localSolvers_[subdomain]->init(prefix + "  ", devnull);
        // initilalize the multiscale mappers
        ansatzMapper_.add(subdomain, localSolvers_[subdomain]->ansatzSpace().map().size());
        testMapper_.add(subdomain, localSolvers_[subdomain]->testSpace().map().size());
        if (verbose)
          out << "." << std::flush;
      } // walk the subdomains for the first time
      ansatzMapper_.finalize();
      testMapper_.finalize();
      out<< " done (took " << timer.elapsed() << " sek)" << std::endl;

      // walk the subdomains for the second time
      //   * to initialize the coupling solvers,
      //   * to initialize the boundary solvers and
      //   * to build up the global sparsity pattern
      out << prefix << "initializing coupling and boundary solvers (on " << subdomains << " subdomains)"  << std::flush;
      if (!verbose)
        out << "..." << std::flush;
      timer.reset();
      for (unsigned int subdomain = 0; subdomain < msGrid_->size(); ++subdomain) {
        // copy the local solvers pattern (which has to be done here, bc the multiscale mappers are not ready above)
        addLocalToGlobalPattern(*(localSolvers_[subdomain]->pattern()),
                                subdomain,
                                subdomain,
                                *pattern_);
        // initialize the boundary solvers
        Dune::shared_ptr< BoundarySolverType > boundarySolver(new BoundarySolverType(msGrid_->boundaryGridPart(subdomain),
                                                                                     boundaryInfo_,
                                                                                     localSolvers_[subdomain],
                                                                                     paramTree_));
        boundarySolver->init(prefix + "  ", devnull);
        boundarySolvers[subdomain] = boundarySolver;
        // and copy its pattern
        addLocalToGlobalPattern(*(boundarySolver->pattern()),
                                subdomain,
                                subdomain,
                                *pattern_);
        // walk the neighbors
        std::map< unsigned int, Dune::shared_ptr< CouplingSolverType > >& couplingSolverMap = couplingSolverMaps[subdomain];
        const typename MsGridType::NeighborSetType& neighbors = msGrid_->neighborsOf(subdomain);
        for (typename MsGridType::NeighborSetType::const_iterator neighborIt = neighbors.begin();
             neighborIt != neighbors.end();
             ++neighborIt) {
          const unsigned int neighboringSubdomain = *neighborIt;
          // visit each coupling only once (assemble primaly)
          if (subdomain < neighboringSubdomain) {
            // initialize the coupling solvers
            Dune::shared_ptr< CouplingSolverType > couplingSolver(new CouplingSolverType(msGrid_->couplingGridPart(subdomain, neighboringSubdomain),
                                                                                         msGrid_->couplingGridPart(neighboringSubdomain, subdomain),
                                                                                         localSolvers_[subdomain],
                                                                                         localSolvers_[neighboringSubdomain],
                                                                                         paramTree_));
            couplingSolver->init(prefix + "  ", devnull);
            couplingSolverMap.insert(std::pair< unsigned int, Dune::shared_ptr< CouplingSolverType > >(
                neighboringSubdomain,
                couplingSolver));
            // and copy its patterns
            addLocalToGlobalPattern(*(couplingSolver->innerInnerPattern()),
                                    subdomain,
                                    subdomain,
                                    *pattern_);
            addLocalToGlobalPattern(*(couplingSolver->innerOuterPattern()),
                                    subdomain,
                                    neighboringSubdomain,
                                    *pattern_);
            addLocalToGlobalPattern(*(couplingSolver->outerInnerPattern()),
                                    neighboringSubdomain,
                                    subdomain,
                                    *pattern_);
            addLocalToGlobalPattern(*(couplingSolver->outerOuterPattern()),
                                    neighboringSubdomain,
                                    neighboringSubdomain,
                                    *pattern_);
          } // visit each coupling only once
        } // walk the neighbors
        if (verbose)
          out << "." << std::flush;
      } // walk the subdomains for the second time
      out<< " done (took " << timer.elapsed() << " sek)" << std::endl;

      // walk the subdomains for the third time
      //   * to assemble the global system matrix and right hand side
      out << prefix << "generating global matrix and vector containers (of size "
          << ansatzMapper_.size()
          << "x"
          << testMapper_.size()
          << ")" << std::flush;
      if (!verbose)
        out << "..." << std::flush;
      timer.reset();
      matrix_ = Dune::shared_ptr< MatrixBackendType >(new MatrixBackendType(
          ContainerFactory::createSparseMatrix(ansatzMapper_.size(), testMapper_.size(), *pattern_)));
      rhs_ = Dune::shared_ptr< VectorBackendType >(new VectorBackendType(
          ContainerFactory::createDenseVector(testMapper_.size())));
      for (unsigned int subdomain = 0; subdomain < msGrid_->size(); ++subdomain) {
        // copy the local containers of the subdomain solver
        copyLocalToGlobalMatrix(*(localSolvers_[subdomain]->systemMatrix()),
                                *(localSolvers_[subdomain]->pattern()),
                                subdomain,
                                subdomain,
                                *matrix_);
        copyLocalToGlobalVector(*(localSolvers_[subdomain]->rightHandSide()), subdomain, *rhs_);
        // copy the local containers of the boundary solver
        copyLocalToGlobalMatrix(*(boundarySolvers[subdomain]->systemMatrix()),
                                *(boundarySolvers[subdomain]->pattern()),
                                subdomain,
                                subdomain,
                                *matrix_);
        copyLocalToGlobalVector(*(boundarySolvers[subdomain]->rightHandSide()), subdomain, *rhs_);
        // walk the neighbors
        std::map< unsigned int, Dune::shared_ptr< CouplingSolverType > >& couplingSolverMap = couplingSolverMaps[subdomain];
        const typename MsGridType::NeighborSetType& neighbors = msGrid_->neighborsOf(subdomain);
        for (typename MsGridType::NeighborSetType::const_iterator neighborIt = neighbors.begin();
             neighborIt != neighbors.end();
             ++neighborIt) {
          const unsigned int neighboringSubdomain = *neighborIt;
          // visit each coupling only once (assemble primaly)
          if (subdomain < neighboringSubdomain) {
            const CouplingSolverType& couplingSolver = *(couplingSolverMap.find(neighboringSubdomain)->second);
            // copy the local containers of the subdomain solver
            copyLocalToGlobalMatrix(*(couplingSolver.innerInnerMatrix()),
                                    *(couplingSolver.innerInnerPattern()),
                                    subdomain,
                                    subdomain,
                                    *matrix_);
            copyLocalToGlobalMatrix(*(couplingSolver.innerOuterMatrix()),
                                    *(couplingSolver.innerOuterPattern()),
                                    subdomain,
                                    neighboringSubdomain,
                                    *matrix_);
            copyLocalToGlobalMatrix(*(couplingSolver.outerInnerMatrix()),
                                    *(couplingSolver.outerInnerPattern()),
                                    neighboringSubdomain,
                                    subdomain,
                                    *matrix_);
            copyLocalToGlobalMatrix(*(couplingSolver.outerOuterMatrix()),
                                    *(couplingSolver.outerOuterPattern()),
                                    neighboringSubdomain,
                                    neighboringSubdomain,
                                    *matrix_);
          } // visit each coupling only once
        } // walk the neighbors
        if (verbose)
          out << "." << std::flush;
      } // walk the subdomains for the third time
      out<< " done (took " << timer.elapsed() << " sek)" << std::endl;
      // done
      initialized_ = true;
    } // if (!initialized_)
  } // void init()

  Dune::shared_ptr< std::vector< VectorBackendType > > createVector() const
  {
    typename Dune::shared_ptr< std::vector< VectorBackendType > > retPtr(new std::vector< VectorBackendType >());
    typename std::vector< VectorBackendType >& ret = *retPtr;
    for (unsigned int subdomain = 0; subdomain < msGrid_->size(); ++subdomain) {
      const LocalSolverType& localSolver = *(localSolvers_[subdomain]);
      typename Dune::shared_ptr< VectorBackendType > localVector = localSolver.createVector();
      ret.push_back(*localVector);
    }
    return retPtr;
  } // Dune::shared_ptr< std::vector< VectorBackendType > > createVector() const

  void solve(std::vector< VectorBackendType >& solution,
             const Dune::ParameterTree paramTree = Dune::ParameterTree(),
             const std::string prefix = "",
             std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
  {
    // preparations
    assert(initialized_ && "The system can only be solved after init() has been called! ");
    assert(solution.size() == msGrid_->size() && "Given vector has wrong size!");
    const std::string type = paramTree.get("type", "eigen.bicgstab.incompletelut");
    const unsigned int maxIter = paramTree.get("maxIter", 5000);
    const double precision = paramTree.get("precision", 1e-12);
    Dune::Timer timer;
    // create global solution vector
    VectorBackendType tmpSolution = ContainerFactory::createDenseVector(testMapper_.size());
    // solve
    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::BicgstabIlut BicgstabIlutSolver;
    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::BicgstabDiagonal BicgstabDiagonalSolver;
    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::CgDiagonalUpper CgDiagonalUpperSolver;
    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::CgDiagonalLower CgDiagonalLowerSolver;
    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::SimplicialcholeskyUpper SimplicialcholeskyUpperSolver;
    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::SimplicialcholeskyLower SimplicialcholeskyLowerSolver;
    out << prefix << "solving linear system of size " << matrix_->rows() << "x" << matrix_->cols() << std::endl
        << prefix << "  using " << type << "... " << std::flush;
    if (type == "eigen.bicgstab.incompletelut"){
      BicgstabIlutSolver::apply(*matrix_, tmpSolution, *rhs_, maxIter, precision);
    } else if (type == "eigen.bicgstab.diagonal"){
      BicgstabDiagonalSolver::apply(*matrix_, tmpSolution, *rhs_, maxIter, precision);
    } else if (type == "eigen.cg.diagonal.upper"){
      CgDiagonalUpperSolver::apply(*matrix_, tmpSolution, *rhs_, maxIter, precision);
    } else if (type == "eigen.cg.diagonal.lower"){
      CgDiagonalLowerSolver::apply(*matrix_, tmpSolution, *rhs_, maxIter, precision);
    } else if (type == "eigen.simplicialcholesky.upper"){
      SimplicialcholeskyUpperSolver::apply(*matrix_, tmpSolution, *rhs_, maxIter, precision);
    } else if (type == "eigen.simplicialcholesky.lower"){
      SimplicialcholeskyLowerSolver::apply(*matrix_, tmpSolution, *rhs_, maxIter, precision);
    } else {
      std::stringstream msg;
      msg << "Error";
      if (id != "") {
        msg << " in " << id;
      }
      msg << ": solver type '" << type << "not supported!" << std::endl;
      DUNE_THROW(Dune::InvalidStateException, msg.str());
    } // solve
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
    // copy global vector to local vectors
    out << prefix << "copying global vector to local...  " << std::flush;
    timer.reset();
    for (unsigned int subdomain = 0; subdomain < msGrid_->size(); ++subdomain) {
      VectorBackendType& localVector = solution[subdomain];
      for (unsigned int localI = 0; localI < localVector.size(); ++localI) {
        const unsigned int globalI = testMapper_.toGlobal(subdomain, localI);
        localVector.set(localI, tmpSolution.get(globalI));
      }
    } // copy global vector to local vector
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void solve(...) const

  void visualize(const std::vector< VectorBackendType >& vector,
                 const std::string filename = "solution",
                 const std::string name = "solution",
                 const std::string prefix = "",
                 std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
  {
    // preparations
    assert(initialized_ && "A vector can only be visualized after init() has been called! ");
    assert(vector.size() == msGrid_->size() && "Given vector has wrong size!");
    Dune::Timer timer;
    out << prefix << "writing '" << name << "' to '" << filename;
    if (dimDomain == 1)
      out << ".vtp";
    else
      out << ".vtu";
    out << "'... " << std::flush;
    // create vector of local discrete functions
    std::vector< Dune::shared_ptr< LocalDiscreteFunctionType > > localDiscreteFunctions;
    for (unsigned int subdomain = 0; subdomain < msGrid_->size(); ++subdomain) {
      localDiscreteFunctions.push_back(Dune::shared_ptr< LocalDiscreteFunctionType >(
          new LocalDiscreteFunctionType(localSolvers_[subdomain]->ansatzSpace(),
                                        vector[subdomain])));
    }
    // create multiscale discrete function
    Dune::shared_ptr< DiscreteFunctionType > discreteFunction(new DiscreteFunctionType(*msGrid_, localDiscreteFunctions, name));
    // visualize
    visualize(discreteFunction, filename, "", Dune::Stuff::Common::Logger().devnull());
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void visualize(...) const

  void visualize(const Dune::shared_ptr< const DiscreteFunctionType > discreteFunction,
                 const std::string filename = "discreteFunction",
                 const std::string prefix = "",
                 std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
  {
    // preparations
    assert(initialized_ && "Please call init() before calling visualize()! ");
    Dune::Timer timer;
    out << prefix << "writing '" << discreteFunction->name() << "' to '" << filename;
    if (dimDomain == 1)
      out << ".vtp";
    else
      out << ".vtu";
    out << "'... " << std::flush;
    typedef Dune::VTKWriter< typename MsGridType::GlobalGridViewType > VTKWriterType;
    VTKWriterType vtkWriter(*(msGrid_->globalGridView()));
    vtkWriter.addVertexData(discreteFunction);
    vtkWriter.write(filename);
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void visualize(...) const

  const Dune::shared_ptr< const LocalSolverType > localSolver(const unsigned int subdomain) const
  {
    assert(initialized_ && "Please call init() before calling localSolver()!");
    assert(subdomain < msGrid_->size());
    return localSolvers_[subdomain];
  }

  Dune::shared_ptr< LocalSolverType > localSolver(const unsigned int subdomain)
  {
    assert(initialized_ && "Please call init() before calling localSolver()!");
    assert(subdomain < msGrid_->size());
    return localSolvers_[subdomain];
  }

private:
  void addLocalToGlobalPattern(const PatternType& local,
                               const unsigned int ansatzSubdomain,
                               const unsigned int testSubdomain,
                               PatternType& global) const
  {
    // loop over all rows of the local pattern
    for (typename PatternType::const_iterator localRowIt = local.begin();
         localRowIt != local.end();
         ++localRowIt) {
      const typename PatternType::key_type localRowIndex = localRowIt->first;
      const typename PatternType::key_type globalRowIndex = ansatzMapper_.toGlobal(ansatzSubdomain, localRowIndex);
      const typename PatternType::mapped_type localRowSet = localRowIt->second;
      // get corresponding row of the global pattern (make use of automatic creation of the map with [])
      typename PatternType::mapped_type& globalRowSet = global[globalRowIndex];
      // loop over all column entries in local row
      for (typename PatternType::mapped_type::const_iterator localColumnIt = localRowSet.begin();
           localColumnIt != localRowSet.end();
           ++localColumnIt) {
        // map local to global
        const typename PatternType::mapped_type::key_type localColumnIndex = *localColumnIt;
        const typename PatternType::mapped_type::key_type globalColumnIndex = testMapper_.toGlobal(testSubdomain, localColumnIndex);
        // add colum
        globalRowSet.insert(globalColumnIndex);
      } // loop over all column entries in args row
    } // loop pver all rows of the input pattern
  } // void addLocalToGlobalPattern() const

  void copyLocalToGlobalMatrix(const MatrixBackendType& localMatrix,
                               const PatternType& localPattern,
                               const unsigned int ansatzSubdomain,
                               const unsigned int testSubdomain,
                               MatrixBackendType& global) const
  {
    // loop over all local rows
    for (typename PatternType::const_iterator rowIterator = localPattern.begin();
         rowIterator != localPattern.end();
         ++rowIterator) {
      const unsigned int localRow = rowIterator->first;
      const unsigned int globalRow = ansatzMapper_.toGlobal(ansatzSubdomain, localRow);
      const typename PatternType::mapped_type& localRowSet = rowIterator->second;
      // loop over all cols in the current row
      for (typename PatternType::mapped_type::const_iterator colIterator = localRowSet.begin();
           colIterator != localRowSet.end();
           ++colIterator) {
        const unsigned int localCol = *colIterator;
        const unsigned int globalCol = testMapper_.toGlobal(testSubdomain, localCol);
        // add entry
        global.add(globalRow, globalCol, localMatrix.get(localRow, localCol));
      } // loop over all cols in the current row
    } // loop over all local rows
  } // void copyLocalToGlobalMatrix() const

  void copyLocalToGlobalVector(const VectorBackendType& local, const unsigned int subdomain, VectorBackendType& global)
  {
    for (unsigned int localI = 0; localI < local.size(); ++localI) {
      const unsigned int globalI = testMapper_.toGlobal(subdomain, localI);
      global.add(globalI, local.get(localI));
    }
  } // copyLocalToGlobalVector()

  const Dune::shared_ptr< const ModelType > model_;
  const Dune::shared_ptr< const MsGridType > msGrid_;
  const Dune::shared_ptr< const BoundaryInfoType > boundaryInfo_;
  const Dune::ParameterTree& paramTree_;
  bool initialized_;
  RangeFieldType penaltyFactor_;
  AnsatzMapperType ansatzMapper_;
  TestMapperType testMapper_;
  std::vector< Dune::shared_ptr< LocalSolverType > > localSolvers_;
  Dune::shared_ptr< PatternType > pattern_;
  Dune::shared_ptr< MatrixBackendType > matrix_;
  Dune::shared_ptr< VectorBackendType > rhs_;
}; // class DuneDetailedDiscretizations

template< class ModelType, class GridPartType, class BoundaryInfoType, int polOrder >
const std::string DuneDetailedDiscretizations< ModelType, GridPartType, BoundaryInfoType, polOrder >::id = "detailed.solvers.stationary.linear.elliptic.multiscale.semicontinuousgalerkin";

} // namespace SemicontinuousGalerkin

} // namespace Multiscale

} // namespace Elliptic

} // namespace Linear

} // namespace Stationary

} // namespace Solvers

} // namespace Detailed

} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MULTISCALE_SEMICONTINUOUSGALERKIN_DUNE_DETAILED_DISCRETIZATIONS_HH
