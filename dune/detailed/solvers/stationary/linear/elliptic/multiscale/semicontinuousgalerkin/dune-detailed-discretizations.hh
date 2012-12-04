#ifndef DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MULTISCALE_SEMICONTINUOUSGALERKIN_DUNE_DETAILED_DISCRETIZATIONS_HH
#define DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MULTISCALE_SEMICONTINUOUSGALERKIN_DUNE_DETAILED_DISCRETIZATIONS_HH

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#include <vector>
#include <sstream>

#include <dune/common/exceptions.hh>
#include <dune/common/timer.hh>

#include <dune/grid/multiscale/default.hh>

#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/grid/boundaryinfo.hh>

#include <dune/detailed/discretizations/mapper/multiscale.hh>
#include <dune/detailed/discretizations/la/factory/eigen.hh>
#include <dune/detailed/discretizations/discretefunction/multiscale.hh>
#include <dune/detailed/discretizations/evaluation/local/binary/ipdgfluxes.hh>
#include <dune/detailed/discretizations/discreteoperator/local/codim1/boundaryintegral.hh>
#include <dune/detailed/discretizations/assembler/local/codim1/matrix.hh>
#include <dune/detailed/discretizations/assembler/multiscale/boundary.hh>
#include <dune/detailed/discretizations/discreteoperator/local/codim1/innerintegral.hh>
#include <dune/detailed/discretizations/evaluation/local/quaternary/ipdgfluxes.hh>
#include <dune/detailed/discretizations/assembler/local/codim1/matrix.hh>
#include <dune/detailed/discretizations/assembler/multiscale/coupling.hh>

#include <dune/detailed/solvers/stationary/linear/elliptic/model/interface.hh>
#include <dune/detailed/solvers/stationary/linear/elliptic/continuousgalerkin/dune-detailed-discretizations.hh>

namespace Dune {
namespace Detailed {
namespace Solvers {
namespace Stationary {
namespace Linear {
namespace Elliptic {
namespace Multiscale {
namespace SemicontinuousGalerkin {

/**
 *  \todo Implement non-zero dirichlet and neumann values in assembleBoundaryContribution()
 */
template< class ModelImp, class MsGridImp, class BoundaryInfoImp, int polynomialOrder >
class DuneDetailedDiscretizations
{
public:
  typedef ModelImp ModelType;

  typedef MsGridImp MsGridType;

  typedef BoundaryInfoImp BoundaryInfoType;

  static const int polOrder = polynomialOrder;

  typedef DuneDetailedDiscretizations< ModelType, MsGridType, BoundaryInfoType, polOrder > ThisType;

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

  typedef typename LocalSolverType::AnsatzSpaceType LocalAnsatzSpaceType;

  typedef typename LocalSolverType::TestSpaceType LocalTestSpaceType;

  typedef typename LocalAnsatzSpaceType::PatternType PatternType;

  typedef Dune::Detailed::Discretizations::Mapper::Multiscale< typename LocalAnsatzSpaceType::MapperType::IndexType > AnsatzMapperType;

  typedef Dune::Detailed::Discretizations::Mapper::Multiscale< typename LocalTestSpaceType::MapperType::IndexType > TestMapperType;

public:
  typedef typename LocalSolverType::DiscreteFunctionType LocalDiscreteFunctionType;

  typedef typename Dune::Detailed::Discretizations::DiscreteFunction::Multiscale< MsGridType, LocalDiscreteFunctionType > DiscreteFunctionType;

  static const std::string id()
  {
    return "detailed.solvers.stationary.linear.elliptic.multiscale.semicontinuousgalerkin";
  }

  DuneDetailedDiscretizations(const Dune::shared_ptr< const ModelType > model,
                              const Dune::shared_ptr< const MsGridType > msGrid,
                              const Dune::shared_ptr< const BoundaryInfoType > boundaryInfo,
                              const RangeFieldType& _penaltyFactor)
    : model_(model)
    , msGrid_(msGrid)
    , boundaryInfo_(boundaryInfo)
    , penaltyFactor_(_penaltyFactor)
    , initialized_(false)
    , ansatzMapper_()
    , testMapper_()
    , localSolvers_(msGrid_->size())
    , pattern_(Dune::shared_ptr< PatternType >(new PatternType()))
  {
    // test penalty factor
    assert(penaltyFactor_ > 0 && "Please provide a positive penalty factor!");
    // test orders
    assert(model_->diffusionOrder() >= 0 && "Please provide a nonnegative integration order for the diffusion!");
    assert(model_->forceOrder() >= 0 && "Please provide a nonnegative integration order for the force!");
  } // DuneDetailedDiscretizations()

  const Dune::shared_ptr< const ModelType > model() const
  {
    return model_;
  }

  const Dune::shared_ptr< const MsGridType > msGrid() const
  {
    return msGrid_;
  }

  const Dune::shared_ptr< const BoundaryInfoType > boundaryInfo() const
  {
    return boundaryInfo_;
  }

  void init(const std::string prefix = "", std::ostream& out = Dune::Stuff::Common::Logger().debug())
  {
    if (!initialized_) {
      // prepare
      Dune::Timer timer;
      const unsigned int subdomains = msGrid_->size();
      const bool verbose = (subdomains >= 3) && (subdomains <= std::pow(3, 3));
      std::ostream& devnull = Dune::Stuff::Common::Logger().devnull();
      std::map< unsigned int, Dune::shared_ptr< const PatternType > > boundaryPatternMap;
      std::vector< std::map< unsigned int, std::map< std::string, Dune::shared_ptr< const PatternType > > > > couplingPatternMapMaps(msGrid_->size());

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
      //   * to initialize the coupling pattern,
      //   * to initialize the boundary pattern and
      //   * to build up the global sparsity pattern
      out << prefix << "initializing coupling patterns (on " << subdomains << " subdomains)"  << std::flush;
      if (!verbose)
        out << "..." << std::flush;
      timer.reset();
      for (unsigned int subdomain = 0; subdomain < msGrid_->size(); ++subdomain) {
        // copy the local solvers pattern (which has to be done here, bc the multiscale mappers are not ready yet above)
        addLocalToGlobalPattern(*(localSolvers_[subdomain]->pattern()),
                                subdomain,
                                subdomain,
                                *pattern_);
        // create the boundary pattern
        const typename LocalSolverType::AnsatzSpaceType& innerAnsatzSpace = localSolvers_[subdomain]->ansatzSpace();
        const typename LocalSolverType::TestSpaceType& innerTestSpace = localSolvers_[subdomain]->testSpace();
        if (msGrid_->boundary(subdomain)) {
          const Dune::shared_ptr< const PatternType > boundaryPattern
              = innerAnsatzSpace.computeLocalPattern(*(msGrid_->boundaryGridPart(subdomain)),
                                                     innerTestSpace);
          boundaryPatternMap.insert(std::pair< unsigned int, Dune::shared_ptr< const PatternType > >(
                                      subdomain,
                                      boundaryPattern));
          // and copy it
          addLocalToGlobalPattern(*boundaryPattern,
                                  subdomain,
                                  subdomain,
                                  *pattern_);
        } // if (msGrid->boundary(subdomain))
        // walk the neighbors
        std::map< unsigned int, std::map< std::string, Dune::shared_ptr< const PatternType > > >& couplingPatternMapMap
            = couplingPatternMapMaps[subdomain];
        const typename MsGridType::NeighborSetType& neighbors = msGrid_->neighborsOf(subdomain);
        for (typename MsGridType::NeighborSetType::const_iterator neighborIt = neighbors.begin();
             neighborIt != neighbors.end();
             ++neighborIt) {
          const unsigned int neighboringSubdomain = *neighborIt;
          // visit each coupling only once (assemble primaly)
          if (subdomain < neighboringSubdomain) {
            // create the coupling patterns
            const typename LocalSolverType::TestSpaceType& outerAnsatzSpace = localSolvers_[neighboringSubdomain]->ansatzSpace();
            const typename LocalSolverType::TestSpaceType& outerTestSpace = localSolvers_[neighboringSubdomain]->testSpace();
            const typename MsGridType::CouplingGridPartType& insideOutsideGridPart = *(msGrid_->couplingGridPart(subdomain, neighboringSubdomain));
            const typename MsGridType::CouplingGridPartType& outsideInsideGridPart = *(msGrid_->couplingGridPart(neighboringSubdomain, subdomain));
            std::map< std::string, Dune::shared_ptr< const PatternType > >& couplingPatternMap
                = couplingPatternMapMap[neighboringSubdomain];
            couplingPatternMap.insert(std::pair< std::string, Dune::shared_ptr< const PatternType > >(
                "inside/inside",
                innerAnsatzSpace.computeLocalPattern(insideOutsideGridPart,
                                                     innerTestSpace)));
            couplingPatternMap.insert(std::pair< std::string, Dune::shared_ptr< const PatternType > >(
                "inside/outside",
                innerAnsatzSpace.computeCouplingPattern(insideOutsideGridPart,
                                                        outerTestSpace)));
            couplingPatternMap.insert(std::pair< std::string, Dune::shared_ptr< const PatternType > >(
                "outside/inside",
                outerAnsatzSpace.computeCouplingPattern(outsideInsideGridPart,
                                                        innerTestSpace)));
            couplingPatternMap.insert(std::pair< std::string, Dune::shared_ptr< const PatternType > >(
                "outside/outside",
                outerAnsatzSpace.computeLocalPattern(outsideInsideGridPart,
                                                     outerTestSpace)));
            // and copy them
            addLocalToGlobalPattern(*(couplingPatternMap["inside/inside"]),
                                    subdomain,
                                    subdomain,
                                    *pattern_);
            addLocalToGlobalPattern(*(couplingPatternMap["inside/outside"]),
                                    subdomain,
                                    neighboringSubdomain,
                                    *pattern_);
            addLocalToGlobalPattern(*(couplingPatternMap["outside/inside"]),
                                    neighboringSubdomain,
                                    subdomain,
                                    *pattern_);
            addLocalToGlobalPattern(*(couplingPatternMap["outside/outside"]),
                                    neighboringSubdomain,
                                    neighboringSubdomain,
                                    *pattern_);
          } // visit each coupling only once (assemble primaly)
        } // walk the neighbors
        if (verbose)
          out << "." << std::flush;
      } // walk the subdomains for the second time
      out<< " done (took " << timer.elapsed() << " sek)" << std::endl;

      // walk the subdomains for the third time
      //   * to assemble the coupling matrices,
      //   * to assemble the boundary matrices and vectors and
      //   * to build up the global matrix and vector
      out << prefix << "generating global matrix and vector containers (of size "
          << ansatzMapper_.size()
          << "x"
          << testMapper_.size()
          << ")" << std::flush;
      if (!verbose)
        out << "..." << std::flush;
      timer.reset();
      // initialize the global matrix and vector container
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
        // for the boundary contribution
        const typename LocalSolverType::AnsatzSpaceType& innerAnsatzSpace = localSolvers_[subdomain]->ansatzSpace();
        const typename LocalSolverType::TestSpaceType& innerTestSpace = localSolvers_[subdomain]->testSpace();
        if (msGrid_->boundary(subdomain)) {
          //   * initialize the boundary matrix and vector,
          typename std::map< unsigned int, Dune::shared_ptr< const PatternType > >::const_iterator result = boundaryPatternMap.find(subdomain);
          assert(result != boundaryPatternMap.end());
          const Dune::shared_ptr< const PatternType > boundaryPattern = result->second;
          Dune::shared_ptr< MatrixBackendType > boundaryMatrix(new MatrixBackendType(
              ContainerFactory::createSparseMatrix(innerAnsatzSpace.map().size(),
                                                   innerTestSpace.map().size(),
                                                   *boundaryPattern)));
          Dune::shared_ptr< VectorBackendType > boundaryRhs(new VectorBackendType(
              ContainerFactory::createDenseVector(innerTestSpace.map().size())));
          //   * assemble them
          assembleBoundaryContribution(subdomain, *boundaryMatrix, *boundaryRhs);
          //   * and copy them into the global matrix and vector
          copyLocalToGlobalMatrix(*boundaryMatrix,
                                  *boundaryPattern,
                                  subdomain,
                                  subdomain,
                                  *matrix_);
          copyLocalToGlobalVector(*boundaryRhs, subdomain, *rhs_);
        } // if (msGrid_->boundary(subdomain))
        // walk the neighbors
        std::map< unsigned int, std::map< std::string, Dune::shared_ptr< const PatternType > > >& couplingPatternMapMap
            = couplingPatternMapMaps[subdomain];
        const typename MsGridType::NeighborSetType& neighbors = msGrid_->neighborsOf(subdomain);
        for (typename MsGridType::NeighborSetType::const_iterator neighborIt = neighbors.begin();
             neighborIt != neighbors.end();
             ++neighborIt) {
          const unsigned int neighboringSubdomain = *neighborIt;
          // visit each coupling only once (assemble primaly)
          if (subdomain < neighboringSubdomain) {
            // for the coupling contribution
            const typename LocalSolverType::TestSpaceType& outerAnsatzSpace = localSolvers_[neighboringSubdomain]->ansatzSpace();
            const typename LocalSolverType::TestSpaceType& outerTestSpace = localSolvers_[neighboringSubdomain]->testSpace();
            std::map< std::string, Dune::shared_ptr< const PatternType > >& couplingPatternMap
                = couplingPatternMapMap[neighboringSubdomain];
            const PatternType& insideInsidePattern = *(couplingPatternMap["inside/inside"]);
            const PatternType& insideOutsidePattern = *(couplingPatternMap["inside/outside"]);
            const PatternType& outsideInsidePattern = *(couplingPatternMap["outside/inside"]);
            const PatternType& outsideOutsidePattern = *(couplingPatternMap["outside/outside"]);
            //   * initialize the coupling matrizes,
            Dune::shared_ptr< MatrixBackendType > insideInsideMatrix(new MatrixBackendType(
                ContainerFactory::createSparseMatrix(innerAnsatzSpace.map().size(),
                                                     innerTestSpace.map().size(),
                                                     insideInsidePattern)));
            Dune::shared_ptr< MatrixBackendType > insideOutsideMatrix(new MatrixBackendType(
                ContainerFactory::createSparseMatrix(innerAnsatzSpace.map().size(),
                                                     outerTestSpace.map().size(),
                                                     insideOutsidePattern)));
            Dune::shared_ptr< MatrixBackendType > outsideInsideMatrix(new MatrixBackendType(
                ContainerFactory::createSparseMatrix(outerAnsatzSpace.map().size(),
                                                     innerTestSpace.map().size(),
                                                     outsideInsidePattern)));
            Dune::shared_ptr< MatrixBackendType > outsideOutsideMatrix(new MatrixBackendType(
                ContainerFactory::createSparseMatrix(outerAnsatzSpace.map().size(),
                                                     outerTestSpace.map().size(),
                                                     outsideOutsidePattern)));
            //   * assemble them
            assembleCouplingContribution(subdomain,
                                         neighboringSubdomain,
                                         *insideInsideMatrix,
                                         *insideOutsideMatrix,
                                         *outsideInsideMatrix,
                                         *outsideOutsideMatrix);
            //   * and copy them into the global matrix
            copyLocalToGlobalMatrix(*insideInsideMatrix,
                                    insideInsidePattern,
                                    subdomain,
                                    subdomain,
                                    *matrix_);
            copyLocalToGlobalMatrix(*insideOutsideMatrix,
                                    insideOutsidePattern,
                                    subdomain,
                                    neighboringSubdomain,
                                    *matrix_);
            copyLocalToGlobalMatrix(*outsideInsideMatrix,
                                    outsideInsidePattern,
                                    neighboringSubdomain,
                                    subdomain,
                                    *matrix_);
            copyLocalToGlobalMatrix(*outsideOutsideMatrix,
                                    outsideOutsidePattern,
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
    assert(initialized_ && "Please call init() beafore calling createDiscreteFunction()!");
    typename Dune::shared_ptr< std::vector< VectorBackendType > > retPtr(new std::vector< VectorBackendType >());
    typename std::vector< VectorBackendType >& ret = *retPtr;
    for (unsigned int subdomain = 0; subdomain < msGrid_->size(); ++subdomain) {
      const LocalSolverType& localSolver = *(localSolvers_[subdomain]);
      typename Dune::shared_ptr< VectorBackendType > localVector = localSolver.createVector();
      ret.push_back(*localVector);
    }
    return retPtr;
  } // Dune::shared_ptr< std::vector< VectorBackendType > > createVector() const

  Dune::shared_ptr< DiscreteFunctionType > createDiscreteFunction(const std::string name = "discrete_function") const
  {
    assert(initialized_ && "Please call init() beafore calling createDiscreteFunction()!");
    // create vector of local discrete functions
    std::vector< Dune::shared_ptr< LocalDiscreteFunctionType > > localDiscreteFunctions;
    for (unsigned int subdomain = 0; subdomain < msGrid_->size(); ++subdomain) {
      localDiscreteFunctions.push_back(Dune::shared_ptr< LocalDiscreteFunctionType >(
          new LocalDiscreteFunctionType(localSolvers_[subdomain]->ansatzSpace())));
    }
    // create multiscale discrete function
    Dune::shared_ptr< DiscreteFunctionType > discreteFunction(new DiscreteFunctionType(*msGrid_,
                                                                                       localDiscreteFunctions,
                                                                                       name));
    return discreteFunction;
  }

  Dune::shared_ptr< DiscreteFunctionType > createDiscreteFunction(std::vector< VectorBackendType >& vectors,
                                                                  const std::string name = "discrete_function") const
  {
    assert(initialized_ && "Please call init() beafore calling createDiscreteFunction()!");
    // create vector of local discrete functions
    std::vector< Dune::shared_ptr< LocalDiscreteFunctionType > > localDiscreteFunctions;
    for (unsigned int subdomain = 0; subdomain < msGrid_->size(); ++subdomain) {
      localDiscreteFunctions.push_back(Dune::shared_ptr< LocalDiscreteFunctionType >(
          new LocalDiscreteFunctionType(localSolvers_[subdomain]->ansatzSpace(),
                                        vectors[subdomain])));
    }
    // create multiscale discrete function
    Dune::shared_ptr< DiscreteFunctionType > discreteFunction(new DiscreteFunctionType(*msGrid_,
                                                                                       localDiscreteFunctions,
                                                                                       name));
    return discreteFunction;
  }

  void solve(std::vector< VectorBackendType >& solution,
             const std::string type = "eigen.bicgstab.diagonal",
             const unsigned int maxIter = 5000,
             const double precision = 1e-12,
             const std::string prefix = "",
             std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
  {
    // preparations
    assert(initialized_ && "The system can only be solved after init() has been called! ");
    assert(solution.size() == msGrid_->size() && "Given vector has wrong size!");
    assert(precision > 0 && "Please provide a positive precision!");
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
      if (id() != "") {
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

  const Dune::shared_ptr< const MatrixBackendType > systemMatrix() const
  {
    assert(initialized_);
    return matrix_;
  }

  Dune::shared_ptr< MatrixBackendType > systemMatrix()
  {
    assert(initialized_);
    return matrix_;
  }

  const Dune::shared_ptr< const VectorBackendType > rightHandSide() const
  {
    assert(initialized_);
    return rhs_;
  }

  Dune::shared_ptr< VectorBackendType > rightHandSide()
  {
    assert(initialized_);
    return rhs_;
  }

  const Dune::shared_ptr< const PatternType > pattern() const
  {
    assert(initialized_);
    return pattern_;
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
  } // void addLocalToGlobalPattern(...) const

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
  } // void copyLocalToGlobalMatrix(...) const

  void copyLocalToGlobalVector(const VectorBackendType& local, const unsigned int subdomain, VectorBackendType& global)
  {
    for (unsigned int localI = 0; localI < local.size(); ++localI) {
      const unsigned int globalI = testMapper_.toGlobal(subdomain, localI);
      global.add(globalI, local.get(localI));
    }
  } // void copyLocalToGlobalVector(...)

  void assembleBoundaryContribution(const unsigned int subdomain,
                                    MatrixBackendType& boundaryMatrix,
                                    VectorBackendType& /*boundaryRhs*/) const
  {
    typedef Dune::FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, dimRange > FunctionSpaceType;
    // operator
    typedef typename ModelType::DiffusionType DiffusionType;
    typedef Dune::Detailed::Discretizations::Evaluation::Local::Binary::IPDGfluxes::Dirichlet< FunctionSpaceType, DiffusionType > IPDGfluxType;
    const IPDGfluxType ipdgFlux(model_->diffusion(), model_->diffusionOrder(), penaltyFactor_);
    typedef Dune::Detailed::Discretizations::DiscreteOperator::Local::Codim1::BoundaryIntegral< IPDGfluxType > DirichletOperatorType;
    const DirichletOperatorType dirichletOperator(ipdgFlux);
    // local matrix assembler
    typedef Dune::Detailed::Discretizations::Assembler::Local::Codim1::Boundary< DirichletOperatorType > DirichletMatrixAssemblerType;
    const DirichletMatrixAssemblerType dirichletMatrixAssembler(dirichletOperator);
    // functional
    // ...
    // local vector assmebler
    // ...
    // boundary assembler
    typedef Dune::Detailed::Discretizations::Assembler::Multiscale::Boundary<
        BoundaryGridPartType,
        BoundaryInfoType,
        typename LocalSolverType::AnsatzSpaceType,
        typename LocalSolverType::TestSpaceType >
      BoundaryAssemblerType;
    const BoundaryAssemblerType boundaryAssembler(*(msGrid_->boundaryGridPart(subdomain)),
                                                  boundaryInfo_,
                                                  localSolvers_[subdomain]->ansatzSpace(),
                                                  localSolvers_[subdomain]->testSpace());
    boundaryAssembler.assemble(dirichletMatrixAssembler,
                               boundaryMatrix/*,
                               dirichletVectorAssembler,
                               boundaryRhs*/);
  } // void assembleBoundaryContribution(...)

  void assembleCouplingContribution(const unsigned int subdomain,
                                    const unsigned int neighboringSubdomain,
                                    MatrixBackendType& insideInsideMatrix,
                                    MatrixBackendType& insideOutsideMatrix,
                                    MatrixBackendType& outsideInsideMatrix,
                                    MatrixBackendType& outsideOutsideMatrix) const
  {
    // operator
    typedef Dune::FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, dimRange > FunctionSpaceType;
    typedef typename ModelType::DiffusionType DiffusionType;
    typedef Dune::Detailed::Discretizations::Evaluation::Local::Quaternary::IPDGfluxes::Inner< FunctionSpaceType, DiffusionType > IPDGfluxType;
    const IPDGfluxType ipdgFlux(model_->diffusion(), model_->diffusionOrder(), penaltyFactor_);
    typedef Dune::Detailed::Discretizations::DiscreteOperator::Local::Codim1::InnerIntegral< IPDGfluxType > IPDGoperatorType;
    const IPDGoperatorType ipdgOperator(ipdgFlux);
    // local matrix assembler
    typedef Dune::Detailed::Discretizations::Assembler::Local::Codim1::Inner< IPDGoperatorType > CouplingMatrixAssemblerType;
    const CouplingMatrixAssemblerType couplingMatrixAssembler(ipdgOperator);
    // coupling assembler
    typedef Dune::Detailed::Discretizations::Assembler::Multiscale::Coupling::Primal<
        CouplingGridPartType,
        typename LocalSolverType::AnsatzSpaceType,
        typename LocalSolverType::TestSpaceType,
        typename LocalSolverType::AnsatzSpaceType,
        typename LocalSolverType::TestSpaceType >
      CouplingAssemblerType;
    const CouplingAssemblerType couplingAssembler(*(msGrid_->couplingGridPart(subdomain, neighboringSubdomain)),
                                                  localSolvers_[subdomain]->ansatzSpace(),
                                                  localSolvers_[subdomain]->testSpace(),
                                                  localSolvers_[neighboringSubdomain]->ansatzSpace(),
                                                  localSolvers_[neighboringSubdomain]->testSpace());
    couplingAssembler.assembleMatrices(couplingMatrixAssembler,
                                       insideInsideMatrix,
                                       insideOutsideMatrix,
                                       outsideInsideMatrix,
                                       outsideOutsideMatrix);
  } // void assembleCouplingContribution(...)

  const Dune::shared_ptr< const ModelType > model_;
  const Dune::shared_ptr< const MsGridType > msGrid_;
  const Dune::shared_ptr< const BoundaryInfoType > boundaryInfo_;
  RangeFieldType penaltyFactor_;
  bool initialized_;
  AnsatzMapperType ansatzMapper_;
  TestMapperType testMapper_;
  std::vector< Dune::shared_ptr< LocalSolverType > > localSolvers_;
  Dune::shared_ptr< PatternType > pattern_;
  Dune::shared_ptr< MatrixBackendType > matrix_;
  Dune::shared_ptr< VectorBackendType > rhs_;
}; // class DuneDetailedDiscretizations

} // namespace SemicontinuousGalerkin
} // namespace Multiscale
} // namespace Elliptic
} // namespace Linear
} // namespace Stationary
} // namespace Solvers
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MULTISCALE_SEMICONTINUOUSGALERKIN_DUNE_DETAILED_DISCRETIZATIONS_HH
