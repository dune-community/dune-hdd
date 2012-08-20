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
#include <dune/detailed/discretizations/discretefunctionspace/continuous/lagrange.hh>
#include <dune/detailed/discretizations/discretefunctionspace/sub/linear.hh>
#include <dune/detailed/discretizations/evaluation/local/binary/elliptic.hh>
#include <dune/detailed/discretizations/discreteoperator/local/codim0/integral.hh>
#include <dune/detailed/discretizations/evaluation/local/unary/scale.hh>
#include <dune/detailed/discretizations/discretefunctional/local/codim0/integral.hh>
#include <dune/detailed/discretizations/la/factory/eigen.hh>
#include <dune/detailed/discretizations/mapper/multiscale.hh>
#include <dune/detailed/discretizations/discretefunction/multiscale.hh>
#include <dune/detailed/discretizations/assembler/multiscale/coupling.hh>
#include <dune/detailed/discretizations/evaluation/local/quaternary/ipdgfluxes.hh>
#include <dune/detailed/discretizations/discreteoperator/local/codim1/innerintegral.hh>
#include <dune/detailed/discretizations/assembler/local/codim1/matrix.hh>

// dune-detailed-solvers
#include <dune/detailed/solvers/stationary/linear/elliptic/continuousgalerkin/dune-detailed-discretizations.hh>

// dune-stuff
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/grid/boundaryinfo.hh>

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

  typedef typename MsGridType::LocalGridPartType LocalGridPartType;

  typedef GlobalGridPartType GridPartType;

  static const int polOrder = polynomialOrder;

  typedef DuneDetailedDiscretizations< ModelType, MsGridType, polOrder > ThisType;

  static const std::string id;

private:
  typedef typename ModelType::DomainFieldType DomainFieldType;

  static const int dimDomain = ModelType::dimDomain;

  typedef typename ModelType::RangeFieldType RangeFieldType;

  static const int dimRange = ModelType::dimRange;

  typedef Dune::FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, dimRange > FunctionSpaceType;

  typedef Dune::Detailed::Discretizations::DiscreteFunctionSpace::Continuous::Lagrange< FunctionSpaceType, LocalGridPartType, polOrder > LocalDiscreteH1Type;

  typedef Dune::Stuff::Grid::BoundaryInfo::IdBased BoundaryInfoType;

  typedef Dune::Detailed::Discretizations::DiscreteFunctionSpace::Sub::Linear::Dirichlet< LocalDiscreteH1Type, BoundaryInfoType > LocalAnsatzSpaceType;

  typedef LocalAnsatzSpaceType LocalTestSpaceType;

  typedef typename LocalAnsatzSpaceType::PatternType SparsityPatternType;

  typedef Dune::Detailed::Discretizations::Mapper::Multiscale< typename LocalAnsatzSpaceType::MapperType::IndexType > AnsatzMapperType;

  typedef AnsatzMapperType TestMapperType;

  typedef typename MsGridType::CouplingGridPartType CouplingGridPartType;

  typedef Dune::Detailed::Discretizations::LA::Factory::Eigen< RangeFieldType > ContainerFactory;

  typedef typename ContainerFactory::SparseMatrixType MatrixBackendType;

  typedef typename ContainerFactory::DenseVectorType VectorBackendType;

public:
  typedef typename VectorBackendType::StorageType LocalVectorType;

  typedef Dune::Detailed::Discretizations::DiscreteFunction::Default< LocalAnsatzSpaceType, VectorBackendType > LocalDiscreteFunctionType;

  typedef typename Dune::Detailed::Discretizations::DiscreteFunction::Multiscale< MsGridType, LocalDiscreteFunctionType > DiscreteFunctionType;

  DuneDetailedDiscretizations(const ModelType& model, const MsGridType& msGrid, const Dune::ParameterTree& paramTree)
    : model_(model)
    , msGrid_(msGrid)
    , initialized_(false)
    , penaltyFactor_(-1)
    , ansatzMapper_()
    , testMapper_()
    , localGridParts_(msGrid_.size())
    , localDiscreteH1s_(msGrid_.size())
    , localAnsatzSpaces_(msGrid_.size())
    , localTestSpaces_(msGrid_.size())
    , localSparsityPatterns_(msGrid_.size())
    , couplingSpartityPatternMaps_(msGrid_.size())
    , localMatrices_(msGrid_.size())
    , localVectors_(msGrid_.size())
    , couplingMatricesMaps_(msGrid_.size())
    , globalSparsityPattern_(Dune::shared_ptr< SparsityPatternType >(new SparsityPatternType()))
  {
    // get penalty factor
    std::string key = "discretization.penaltyFactor";
    Dune::Stuff::Common::Parameter::Tree::assertKey(paramTree, key, id);
    penaltyFactor_ = paramTree.get(key, -1);
    if (!penaltyFactor_ > 0) {
      std::stringstream msg;
      msg << "Error in " << id << ":"
          << "wrong '" << key << "' given (should be posititve double) in the following Dune::ParameterTree:" << std::endl;
      paramTree.report(msg);
      DUNE_THROW(Dune::InvalidStateException, msg.str());
    }
  } // DuneDetailedDiscretizations()

  void init(const std::string prefix = "", std::ostream& out = Dune::Stuff::Common::Logger().debug())
  {
    if (!initialized_) {
      // preparations
      Dune::Timer timer;
      const unsigned int subdomains = msGrid_.size();
      const bool verbose = subdomains < std::pow(3, 3);
      ansatzMapper_.prepare();
      testMapper_.prepare();

      out << prefix << "initializing locally (on " << subdomains << " subdomains)"  << std::flush;
      if (!verbose)
        out << "...";
      // boundary info
      typedef typename BoundaryInfoType::IdSetType BoundaryIdSetType;
      typedef typename BoundaryInfoType::IdSetMapType BoundaryIdSetMapType;
      Dune::shared_ptr< BoundaryIdSetMapType > boundaryIdSetMap = Dune::shared_ptr< BoundaryIdSetMapType >(new BoundaryIdSetMapType());
      boundaryIdSetMap->insert(std::pair< std::string, BoundaryIdSetType >("dirichlet", BoundaryIdSetType()));
      for (int id = 1; id < 7; ++id) boundaryIdSetMap->operator[]("dirichlet").insert(id); // because the inner boundaries will hopefully have id 7
      boundaryInfo_ = Dune::shared_ptr< BoundaryInfoType >(new BoundaryInfoType(boundaryIdSetMap));
      // walk the subdomains to initialize locally
      for (unsigned int subdomain = 0; subdomain < subdomains; ++subdomain) {
        // local grid part
        localGridParts_[subdomain] = Dune::shared_ptr< const LocalGridPartType >(msGrid_.localGridPart(subdomain));
        // function spaces
        localDiscreteH1s_[subdomain]
            = Dune::shared_ptr< const LocalDiscreteH1Type >(new LocalDiscreteH1Type(*(localGridParts_[subdomain])));
        localAnsatzSpaces_[subdomain]
            = Dune::shared_ptr< const LocalAnsatzSpaceType >(new LocalAnsatzSpaceType(*(localDiscreteH1s_[subdomain]), boundaryInfo_));
        localTestSpaces_[subdomain]
            = Dune::shared_ptr< const LocalTestSpaceType >(new LocalTestSpaceType(*(localDiscreteH1s_[subdomain]), boundaryInfo_));
        // multiscale mapper
        ansatzMapper_.add(subdomain, localAnsatzSpaces_[subdomain]->map().size());
        testMapper_.add(subdomain, localTestSpaces_[subdomain]->map().size());
        // sparsity pattern
        localSparsityPatterns_[subdomain] = localAnsatzSpaces_[subdomain]->computePattern(*(localTestSpaces_[subdomain]));
        if (verbose)
          out << "." << std::flush;
      } // walk the subdomains to initialize locally
      ansatzMapper_.finalize();
      testMapper_.finalize();
      out<< " done (took " << timer.elapsed() << " sek)" << std::endl;

      out << prefix << "initializing coupling (on " << subdomains << " subdomains)"  << std::flush;
      if (!verbose)
        out << "...";
      timer.reset();
      // walk the subdomains to initialize coupling
      for (unsigned int subdomain = 0; subdomain < subdomains; ++subdomain) {
        // sparsity pattern map
        std::map< unsigned int, Dune::shared_ptr< const SparsityPatternType > >&
            couplingSparsityPatternMap = couplingSpartityPatternMaps_[subdomain];
        // init inside / inside pattern
        Dune::shared_ptr< SparsityPatternType > insideInsidePatterns(new SparsityPatternType());
        // walk the neighbors
        const Dune::shared_ptr< const typename MsGridType::NeighboringSubdomainsSetType > neighbors = msGrid_.neighborsOf(subdomain);
        for (typename MsGridType::NeighboringSubdomainsSetType::const_iterator neighborIt = neighbors->begin();
             neighborIt != neighbors->end();
             ++neighborIt) {
          const unsigned int neighboringSubdomain = *neighborIt;
          // get coupling grid part
          const Dune::shared_ptr< const typename MsGridType::CouplingGridPartType >
              couplingGridPart = msGrid_.couplingGridPart(subdomain, neighboringSubdomain);
          // get inside / inside pattern for this coupling grid part (the outside / outside pattern is computed when the neighbor is the subdomain)
          const Dune::shared_ptr< const SparsityPatternType >
              insideInsidePattern = localAnsatzSpaces_[subdomain]->computeLocalPattern(*couplingGridPart);
          // add to the overall inside / inside pattern
          addPatternToPattern(*insideInsidePattern, *insideInsidePatterns);
          // compute the inside / outside pattern (the outside / inside pattern is computed when the neighbor is the subdomain)
          const Dune::shared_ptr< const SparsityPatternType >
              insideOutsidePattern = localAnsatzSpaces_[subdomain]->computeCouplingPattern(*couplingGridPart,
                                                                                           *(localTestSpaces_[neighboringSubdomain]));
          // add it to the map
          couplingSparsityPatternMap.insert(std::pair< unsigned int, Dune::shared_ptr< const SparsityPatternType > >(
            neighboringSubdomain,
            insideOutsidePattern));
        } // walk the neighbors
        // finally add the overall inside / inside pattern
        addPatternToPattern(*(localSparsityPatterns_[subdomain]), *insideInsidePatterns);
        localSparsityPatterns_[subdomain] = insideInsidePatterns;
        if (verbose)
          out << "." << std::flush;
      } // walk the subdomains to initialize coupling
      out<< " done (took " << timer.elapsed() << " sek)" << std::endl;

      out << prefix << "initializing local matrices (on " << subdomains << " subdomains)"  << std::flush;
      if (!verbose)
        out << "...";
      timer.reset();
      // walk the subdomains to initialize storage
      for (unsigned int subdomain = 0; subdomain < subdomains; ++subdomain) {
        // create local matrix and vector (inside / inside)
        localMatrices_[subdomain] = Dune::shared_ptr< MatrixBackendType >(new MatrixBackendType(
            ContainerFactory::createSparseMatrix(localAnsatzSpaces_[subdomain]->map().size(),
                                                 localTestSpaces_[subdomain]->map().size(),
                                                 *(localSparsityPatterns_[subdomain]))));
        localVectors_[subdomain] = Dune::shared_ptr< VectorBackendType >(new VectorBackendType(
            ContainerFactory::createDenseVector(localTestSpaces_[subdomain]->map().size())));
        // get coupling matrices map
        typename std::map< unsigned int, Dune::shared_ptr< MatrixBackendType > >& couplingMatricesMap = couplingMatricesMaps_[subdomain];
        // walk the neighbors
        const Dune::shared_ptr< const typename MsGridType::NeighboringSubdomainsSetType > neighbors = msGrid_.neighborsOf(subdomain);
        for (typename MsGridType::NeighboringSubdomainsSetType::const_iterator neighborIt = neighbors->begin();
             neighborIt != neighbors->end();
             ++neighborIt) {
          const unsigned int neighboringSubdomain = *neighborIt;
          // create coupling matrix (inside / outside) (we assume here, that the correct coupling pattern exists, thus no checking of the map!)
          couplingMatricesMap.insert(std::pair< unsigned int, Dune::shared_ptr< MatrixBackendType > >(
              neighboringSubdomain,
              Dune::shared_ptr< MatrixBackendType >(new MatrixBackendType(
                  ContainerFactory::createSparseMatrix(localAnsatzSpaces_[subdomain]->map().size(),
                                                       localTestSpaces_[neighboringSubdomain]->map().size(),
                                                       *(couplingSpartityPatternMaps_[subdomain][neighboringSubdomain]))))));
        } // walk the neighbors
        if (verbose)
          out << "." << std::flush;
      } // walk the subdomains to initialize storage
      out<< " done (took " << timer.elapsed() << " sek)" << std::endl;

      out << prefix << "initializing global matrix (of size " << ansatzMapper_.size() << "x" << testMapper_.size() << ") and vector" << std::flush;
      if (!verbose)
        out << "...";
      // walk the subdomains to create the global sparsity pattern
      for (unsigned int subdomain = 0; subdomain < subdomains; ++subdomain) {
        // add the local (inside / inside) pattern
        addLocalToGlobalPattern(*(localSparsityPatterns_[subdomain]), subdomain, subdomain, *globalSparsityPattern_);
        // walk the neighbors
        const Dune::shared_ptr< const typename MsGridType::NeighboringSubdomainsSetType > neighbors = msGrid_.neighborsOf(subdomain);
        for (typename MsGridType::NeighboringSubdomainsSetType::const_iterator neighborIt = neighbors->begin();
             neighborIt != neighbors->end();
             ++neighborIt) {
          const unsigned int neighboringSubdomain = *neighborIt;
          // add the coupling (inside / outside) pattern
          addLocalToGlobalPattern(*(couplingSpartityPatternMaps_[subdomain][neighboringSubdomain]),
                                  subdomain,
                                  neighboringSubdomain,
                                  *globalSparsityPattern_);
        } // walk the neighbors
        if (verbose)
          out << "." << std::flush;
      } // walk the subdomains to create the global sparsity pattern
      matrix_ = Dune::shared_ptr< MatrixBackendType >(new MatrixBackendType(
          ContainerFactory::createSparseMatrix(ansatzMapper_.size(),
                                               testMapper_.size(),
                                               *globalSparsityPattern_)));
      rhs_ = Dune::shared_ptr< VectorBackendType >(new VectorBackendType(
          ContainerFactory::createDenseVector(testMapper_.size())));
      out<< " done (took " << timer.elapsed() << " sek)" << std::endl;

      out << prefix << "assembling local and coupling matrices and vectors (on " << subdomains << " subdomains)"  << std::flush;
      if (!verbose)
        out << "...";
      timer.reset();
      // walk the subdomains to assemble locally
      for (unsigned int subdomain = 0; subdomain < subdomains; ++subdomain) {
        // left hand side
        //   operator
        typedef typename ModelType::DiffusionType DiffusionType;
        typedef Dune::Detailed::Discretizations::Evaluation::Local::Binary::Elliptic< FunctionSpaceType, DiffusionType > EllipticEvaluationType;
        assert(model_.diffusionOrder() >= 0 && "Please provide a nonnegative integration order!");
        const EllipticEvaluationType ellipticEvaluation(model_.diffusion(), model_.diffusionOrder());
        typedef Dune::Detailed::Discretizations::DiscreteOperator::Local::Codim0::Integral< EllipticEvaluationType > EllipticOperatorType;
        const EllipticOperatorType ellipticOperator(ellipticEvaluation);
        //   matrix assembler
        typedef Dune::Detailed::Discretizations::Assembler::Local::Codim0::Matrix< EllipticOperatorType > LocalMatrixAssemblerType;
        const LocalMatrixAssemblerType localMatrixAssembler(ellipticOperator);
        // right hand side
        //   functional
        typedef typename ModelType::ForceType ForceType;
        typedef Dune::Detailed::Discretizations::Evaluation::Local::Unary::Scale< FunctionSpaceType, ForceType > ProductEvaluationType;
        assert(model_.forceOrder() >= 0 && "Please provide a nonnegative integration order!");
        const ProductEvaluationType productEvaluation(model_.force(), model_.forceOrder());
        typedef Dune::Detailed::Discretizations::DiscreteFunctional::Local::Codim0::Integral< ProductEvaluationType > L2FunctionalType;
        const L2FunctionalType l2Functional(productEvaluation);
        //   vector assembler
        typedef Dune::Detailed::Discretizations::Assembler::Local::Codim0::Vector< L2FunctionalType > LocalVectorAssemblerType;
        const LocalVectorAssemblerType localVectorAssembler(l2Functional);
        // assemble system
        typedef Dune::Detailed::Discretizations::Assembler::System::Constrained< LocalAnsatzSpaceType, LocalTestSpaceType > LocalSystemAssemblerType;
        const LocalSystemAssemblerType localSystemAssembler(*(localAnsatzSpaces_[subdomain]), *(localTestSpaces_[subdomain]));
        localSystemAssembler.assembleSystem(localMatrixAssembler,
                                            *(localMatrices_[subdomain]),
                                            localVectorAssembler,
                                            *(localVectors_[subdomain]));
  //      out << std::endl << "===== subdomain " << subdomain << " matrix (before) =====" << std::endl;
  //      out << *(localMatrices_[subdomain]->storage());
  //      out << std::endl << "=============== " << subdomain << " =====================" << std::endl;
        // walk the neighbors
        const Dune::shared_ptr< const typename MsGridType::NeighboringSubdomainsSetType > neighbors = msGrid_.neighborsOf(subdomain);
        for (typename MsGridType::NeighboringSubdomainsSetType::const_iterator neighborIt = neighbors->begin();
             neighborIt != neighbors->end();
             ++neighborIt) {
          const unsigned int neighboringSubdomain = *neighborIt;
          // assemble primaly
          if (subdomain < neighboringSubdomain) {
            // evaluation
            typedef Dune::Detailed::Discretizations::Evaluation::Local::Quaternary::IPDGfluxes::Inner< FunctionSpaceType, DiffusionType > IPDGfluxType;
            const IPDGfluxType ipdgFlux(model_.diffusion(), model_.diffusionOrder(), penaltyFactor_);
            // operator
            typedef Dune::Detailed::Discretizations::DiscreteOperator::Local::Codim1::InnerIntegral< IPDGfluxType > IPDGoperatorType;
            const IPDGoperatorType ipdgOperator(ipdgFlux);
            // local assembler
            typedef Dune::Detailed::Discretizations::Assembler::Local::Codim1::Inner< IPDGoperatorType > CouplingMatrixAssemblerType;
            const CouplingMatrixAssemblerType couplingMatrixAssembler(ipdgOperator);
            // assemble coupling matrix
            typedef Dune::Detailed::Discretizations::Assembler::Multiscale::Coupling::Primal<
                typename MsGridType::CouplingGridPartType,
                LocalAnsatzSpaceType,
                LocalTestSpaceType,
                LocalAnsatzSpaceType,
                LocalTestSpaceType > CouplingAssemblerType;
            const CouplingAssemblerType couplingAssembler(*(msGrid_.couplingGridPart(subdomain, neighboringSubdomain)),
                                                          *(localAnsatzSpaces_[subdomain]),
                                                          *(localTestSpaces_[subdomain]),
                                                          *(localAnsatzSpaces_[neighboringSubdomain]),
                                                          *(localTestSpaces_[neighboringSubdomain]));
            couplingAssembler.assembleMatrices(couplingMatrixAssembler,
                                               *(localMatrices_[subdomain]),
                                               *(localMatrices_[neighboringSubdomain]),
                                               *(couplingMatricesMaps_[subdomain][neighboringSubdomain]),
                                               *(couplingMatricesMaps_[neighboringSubdomain][subdomain]));
          } // assemble primaly
        } // walk the neighbors
  //      out << std::endl << "===== subdomain " << subdomain << " matrix (after) =====" << std::endl;
  //      out << *(localMatrices_[subdomain]->storage());
  //      out << std::endl << "=============== " << subdomain << " ====================" << std::endl;
        if (verbose)
          out << "." << std::flush;
      } // walk the subdomains to assemble locally
      out<< " done (took " << timer.elapsed() << " sek)" << std::endl;

      out << prefix << "copying local matrices and vectors to global (on " << subdomains << " subdomains)"  << std::flush;
      if (!verbose)
        out << "...";
      timer.reset();
      // walk the subdomains to copy from local to global
      for (unsigned int subdomain = 0; subdomain < subdomains; ++subdomain) {
        // copy local matrix
        copyLocalToGlobalMatrix(*(localMatrices_[subdomain]),
                                *(localSparsityPatterns_[subdomain]),
                                subdomain,
                                subdomain,
                                *matrix_);
        copyLocalToGlobalVector(*(localVectors_[subdomain]), subdomain, *rhs_);
        // walk the neighbors
        const Dune::shared_ptr< const typename MsGridType::NeighboringSubdomainsSetType > neighbors = msGrid_.neighborsOf(subdomain);
        for (typename MsGridType::NeighboringSubdomainsSetType::const_iterator neighborIt = neighbors->begin();
             neighborIt != neighbors->end();
             ++neighborIt) {
          const unsigned int neighboringSubdomain = *neighborIt;
          copyLocalToGlobalMatrix(*(couplingMatricesMaps_[subdomain][neighboringSubdomain]),
                                  *(couplingSpartityPatternMaps_[subdomain][neighboringSubdomain]),
                                  subdomain,
                                  neighboringSubdomain,
                                  *matrix_);
        } // walk the neighbors
        if (verbose)
          out << "." << std::flush;
      } // walk the subdomains to copy from local to global
      out<< " done (took " << timer.elapsed() << " sek)" << std::endl;
      // done
      initialized_ = true;
    } // if (!initialized_)
  } // init()

  std::vector< Dune::shared_ptr< LocalVectorType > > createVector() const
  {
    typename std::vector< Dune::shared_ptr< LocalVectorType > > ret;
    for (unsigned int subdomain = 0; subdomain < msGrid_.size(); ++subdomain) {
      VectorBackendType tmp = ContainerFactory::createDenseVector(*(localTestSpaces_[subdomain]));
      ret.push_back(tmp.storage());
    }
    return ret;
  } // std::vector< Dune::shared_ptr< LocalVectorType > > createVector() const

  void solve(std::vector< Dune::shared_ptr< LocalVectorType > >& solution,
             const Dune::ParameterTree paramTree = Dune::ParameterTree(),
             const std::string prefix = "",
             std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
  {
    // preparations
    assert(initialized_);
    assert(solution.size() == msGrid_.size());
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
    out << prefix << "solving linear system of size " << matrix_->rows() << "x" << matrix_->cols() << std::endl
        << prefix << "  using " << type << "... " << std::flush;
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
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;

    // copy global vector to local vector
    out << prefix << "copying global vector to local...  " << std::flush;
    timer.reset();
    for (unsigned int subdomain = 0; subdomain < msGrid_.size(); ++subdomain) {
      VectorBackendType localVectorBackend(solution[subdomain]);
      for (unsigned int localI = 0; localI < localVectorBackend.size(); ++localI) {
        const unsigned int globalI = testMapper_.toGlobal(subdomain, localI);
        localVectorBackend.set(localI, tmp.get(globalI));
      }
    } // copy global vector to local vector
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // solve() const

  void visualize(const std::vector< Dune::shared_ptr< LocalVectorType > >& vector,
                 const std::string filename = "solution",
                 const std::string name = "solution",
                 const std::string prefix = "",
                 std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
  {
    // preparations
    assert(initialized_);
    assert(vector.size() == msGrid_.size());
    Dune::Timer timer;
    out << prefix << "writing '" << name << "' to '" << filename;
    if (dimDomain == 1)
      out << ".vtp";
    else
      out << ".vtu";
    out << "'... " << std::flush;
    // create vector of local discrete functions
    std::vector< VectorBackendType > localVectorBackends;
    typedef Dune::Detailed::Discretizations::DiscreteFunction::Default< LocalAnsatzSpaceType, VectorBackendType > LocalDiscreteFunctionType;
    std::vector< Dune::shared_ptr< LocalDiscreteFunctionType > > localDiscreteFunctions;
    for (unsigned int subdomain = 0; subdomain < msGrid_.size(); ++subdomain) {
      localVectorBackends.push_back(vector[subdomain]);
      localDiscreteFunctions.push_back(Dune::shared_ptr< LocalDiscreteFunctionType >(new LocalDiscreteFunctionType(
          *(localAnsatzSpaces_[subdomain]),
          localVectorBackends[subdomain])));
    }
    // create multiscale discrete function
    typedef Dune::Detailed::Discretizations::DiscreteFunction::Multiscale< MsGridType, LocalDiscreteFunctionType > DiscreteFunctionType;
    Dune::shared_ptr< DiscreteFunctionType > discreteFunction(new DiscreteFunctionType(msGrid_, localDiscreteFunctions, name));
    typedef Dune::VTKWriter< typename MsGridType::GlobalGridViewType > VTKWriterType;
    VTKWriterType vtkWriter(*(msGrid_.globalGridView()));
    vtkWriter.addVertexData(discreteFunction);
    vtkWriter.write(filename);
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // visualize() const

  std::vector< Dune::shared_ptr< const LocalAnsatzSpaceType > > ansatzSpace() const
  {
    return localAnsatzSpaces_;
  }

  std::vector< Dune::shared_ptr< const LocalTestSpaceType > > testSpace() const
  {
    return localTestSpaces_;
  }

private:
  template< class SparsityPatternType >
  void addPatternToPattern(const SparsityPatternType& arg, SparsityPatternType& ret) const
  {
    // loop over all rows of the input pattern
    for (typename SparsityPatternType::const_iterator argRowIt = arg.begin();
         argRowIt != arg.end();
         ++argRowIt) {
      const typename SparsityPatternType::key_type argRowIndex = argRowIt->first;
      const typename SparsityPatternType::mapped_type argRowSet = argRowIt->second;
      // get corresponding row of the output pattern (make use of automatic creation of the map with [])
      typename SparsityPatternType::mapped_type& retRowSet = ret[argRowIndex];
      // loop over all column entries in args row
      for (typename SparsityPatternType::mapped_type::const_iterator argColumnIt = argRowSet.begin();
           argColumnIt != argRowSet.end();
           ++argColumnIt) {
        // add colum
        retRowSet.insert(*argColumnIt);
      } // loop over all column entries in args row
    } // loop pver all rows of the input pattern
  } // addPatternToPattern()

  template< class SparsityPatternType >
  void addLocalToGlobalPattern(const SparsityPatternType& local,
                               const unsigned int ansatzSubdomain,
                               const unsigned int testSubdomain,
                               SparsityPatternType& global) const
  {
    // loop over all rows of the local pattern
    for (typename SparsityPatternType::const_iterator localRowIt = local.begin();
         localRowIt != local.end();
         ++localRowIt) {
      const typename SparsityPatternType::key_type localRowIndex = localRowIt->first;
      const typename SparsityPatternType::key_type globalRowIndex = ansatzMapper_.toGlobal(ansatzSubdomain, localRowIndex);
      const typename SparsityPatternType::mapped_type localRowSet = localRowIt->second;
      // get corresponding row of the global pattern (make use of automatic creation of the map with [])
      typename SparsityPatternType::mapped_type& globalRowSet = global[globalRowIndex];
      // loop over all column entries in local row
      for (typename SparsityPatternType::mapped_type::const_iterator localColumnIt = localRowSet.begin();
           localColumnIt != localRowSet.end();
           ++localColumnIt) {
        // map local to global
        const typename SparsityPatternType::mapped_type::key_type localColumnIndex = *localColumnIt;
        const typename SparsityPatternType::mapped_type::key_type globalColumnIndex = testMapper_.toGlobal(testSubdomain, localColumnIndex);
        // add colum
        globalRowSet.insert(globalColumnIndex);
      } // loop over all column entries in args row
    } // loop pver all rows of the input pattern
  } // addLocalToGlobalPattern()

  void copyLocalToGlobalMatrix(const MatrixBackendType& localMatrix,
                               const SparsityPatternType& localSparsityPattern,
                               const unsigned int ansatzSubdomain,
                               const unsigned int testSubdomain,
                               MatrixBackendType& global)
  {
    // loop over all local rows
    for (typename SparsityPatternType::const_iterator rowIterator = localSparsityPattern.begin();
         rowIterator != localSparsityPattern.end();
         ++rowIterator) {
      const unsigned int localRow = rowIterator->first;
      const unsigned int globalRow = ansatzMapper_.toGlobal(ansatzSubdomain, localRow);
      const typename SparsityPatternType::mapped_type& localRowSet = rowIterator->second;
      // loop over all cols in the current row
      for (typename SparsityPatternType::mapped_type::const_iterator colIterator = localRowSet.begin();
           colIterator != localRowSet.end();
           ++colIterator) {
        const unsigned int localCol = *colIterator;
        const unsigned int globalCol = testMapper_.toGlobal(testSubdomain, localCol);
        // add entry
        global.add(globalRow, globalCol, localMatrix.get(localRow, localCol));
      } // loop over all cols in the current row
    } // loop over all local rows
  } // copyLocalToGlobal()

  void copyLocalToGlobalVector(const VectorBackendType& local, const unsigned int subdomain, VectorBackendType& global)
  {
    for (unsigned int localI = 0; localI < local.size(); ++localI) {
      const unsigned int globalI = testMapper_.toGlobal(subdomain, localI);
      global.set(globalI, local.get(localI));
    }
  } // copyLocalToGlobalVector()

  // initialized
  const ModelType& model_;
  const MsGridType& msGrid_;
  bool initialized_;
  RangeFieldType penaltyFactor_;
  AnsatzMapperType ansatzMapper_;
  TestMapperType testMapper_;
  std::vector< Dune::shared_ptr< const LocalGridPartType > > localGridParts_;
  std::vector< Dune::shared_ptr< const LocalDiscreteH1Type > > localDiscreteH1s_;
  std::vector< Dune::shared_ptr< const LocalAnsatzSpaceType > > localAnsatzSpaces_;
  std::vector< Dune::shared_ptr< const LocalTestSpaceType > > localTestSpaces_;
  std::vector< Dune::shared_ptr< const SparsityPatternType > > localSparsityPatterns_;
  std::vector< std::map< unsigned int, Dune::shared_ptr< const SparsityPatternType > > > couplingSpartityPatternMaps_;
  std::vector< Dune::shared_ptr< MatrixBackendType > > localMatrices_;
  std::vector< Dune::shared_ptr< VectorBackendType > > localVectors_;
  std::vector< std::map< unsigned int, Dune::shared_ptr< MatrixBackendType > > > couplingMatricesMaps_;
  Dune::shared_ptr< SparsityPatternType > globalSparsityPattern_;
  // uninitialized
  Dune::shared_ptr< BoundaryInfoType > boundaryInfo_;
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
