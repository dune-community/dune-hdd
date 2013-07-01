#ifndef DUNE_HDD_LINEARELLIPTIC_SOLVER_MULTISCALE_SEMICG_DETAILED_DISCRETIZATIONS_HH
#define DUNE_HDD_LINEARELLIPTIC_SOLVER_MULTISCALE_SEMICG_DETAILED_DISCRETIZATIONS_HH

#include <memory>
#include <vector>
#include <sstream>
#include <ostream>

#include <dune/common/exceptions.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/timer.hh>

#include <dune/grid/multiscale/default.hh>

#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/la/container/affineparametric.hh>

#include <dune/detailed/discretizations/mapper/multiscale.hh>
#include <dune/detailed/discretizations/la/container/factory/eigen.hh>
#include <dune/detailed/discretizations/discretefunction/multiscale.hh>
#include <dune/detailed/discretizations/evaluation/local/binary/ipdgfluxes.hh>
#include <dune/detailed/discretizations/discreteoperator/local/codim1/boundaryintegral.hh>
#include <dune/detailed/discretizations/assembler/local/codim1/matrix.hh>
#include <dune/detailed/discretizations/assembler/multiscale/boundary.hh>
#include <dune/detailed/discretizations/discreteoperator/local/codim1/innerintegral.hh>
#include <dune/detailed/discretizations/evaluation/local/quaternary/ipdgfluxes.hh>
#include <dune/detailed/discretizations/assembler/local/codim1/matrix.hh>
#include <dune/detailed/discretizations/assembler/multiscale/coupling.hh>

#include "../../cg/detailed-discretizations.hh"
#include "../../../model/interface.hh"


namespace Dune {
namespace HDD {
namespace LinearElliptic {


// forward of the solver, to be used in the traits and allow for specialization
template< class GridImp, class RangeFieldImp, int rangeDim, int polynomialOrder >
class MultiscaleSolverSemiContinuousGalerkinDD
{
public:
  MultiscaleSolverSemiContinuousGalerkinDD() = delete;
};


/**
 *  \todo Implement non-zero dirichlet and neumann values in assembleBoundaryContribution()
 */
template< class GridImp, class RangeFieldImp, int polynomialOrder >
class MultiscaleSolverSemiContinuousGalerkinDD< GridImp, RangeFieldImp, 1, polynomialOrder >
{
public:
  typedef MultiscaleSolverSemiContinuousGalerkinDD< GridImp, RangeFieldImp, 1, polynomialOrder > ThisType;

  typedef Dune::grid::Multiscale::Default< GridImp >  MsGridType;
  typedef typename MsGridType::GridType               GridType;
  typedef typename MsGridType::GlobalGridPartType     GlobalGridPartType;
  typedef typename MsGridType::GlobalGridViewType     GlobalGridViewType;
  typedef typename MsGridType::LocalGridPartType      LocalGridPartType;
  typedef typename MsGridType::CouplingGridPartType   CouplingGridPartType;
  typedef typename MsGridType::BoundaryGridPartType   BoundaryGridPartType;

  static const int polOrder = polynomialOrder;

  typedef typename GridType::ctype  DomainFieldType;
  static const int                  dimDomain = GridType::dimension;
  typedef RangeFieldImp             RangeFieldType;
  static const int                  dimRange = 1;

  typedef ModelInterface< DomainFieldType, dimDomain, RangeFieldType, dimRange >  ModelType;
  typedef Dune::Stuff::GridboundaryInterface< GlobalGridViewType >                BoundaryInfoType;

  typedef typename ModelType::ParamType       ParamType;

  typedef Dune::Stuff::GridboundaryAllNeumann< typename LocalGridPartType::GridViewType >     LocalBoundaryInfoType;
  typedef SolverContinuousGalerkinDD< LocalGridPartType, RangeFieldType, dimRange, polOrder > LocalSolverType;

  typedef typename LocalSolverType::PatternType PatternType;
  typedef typename LocalSolverType::MatrixType  MatrixType;
  typedef typename LocalSolverType::VectorType  VectorType;

  typedef Dune::Stuff::LA::Container::AffineParametric< MatrixType > AffineParametricMatrixType;
  typedef Dune::Stuff::LA::Container::AffineParametric< VectorType > AffineParametricVectorType;

  typedef Dune::Detailed::Discretizations::Mapper::Multiscale<> AnsatzMapperType;
  typedef Dune::Detailed::Discretizations::Mapper::Multiscale<> TestMapperType;

  typedef typename LocalSolverType::DiscreteAnsatzFunctionConstType LocalDiscreteFunctionType;
  typedef typename Dune::Detailed::Discretizations::DiscreteFunction::Multiscale< MsGridType, LocalDiscreteFunctionType > DiscreteFunctionType;

  static const std::string id()
  {
    return "solver.linearelliptic.ms.semicg.dd";
  }

  MultiscaleSolverSemiContinuousGalerkinDD(const std::shared_ptr< const MsGridType > _msGrid,
                                           const std::shared_ptr< const BoundaryInfoType > _boundaryInfo,
                                           const std::shared_ptr< const ModelType > _model,
                                           const RangeFieldType _penaltyFactor)
    : msGrid_(_msGrid)
    , boundaryInfo_(_boundaryInfo)
    , model_(_model)
    , penaltyFactor_(_penaltyFactor)
    , initialized_(false)
    , localSolvers_(msGrid_->size())
    , ansatzMapper_()
    , testMapper_()
  {
    // sanity checks
    std::stringstream msg;
    unsigned int throw_up = 0;
    // * penalty factor
    if (!(penaltyFactor_ > 0)) {
      msg << "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
          << " nonpositive penalty factor given!";
      ++throw_up;
    }
    // * integration orders
    if (model_->diffusion()->order() < 0) {
      msg << "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
          << " negative integration order given for the diffusion!";
      ++throw_up;
    }
    if (model_->force()->order() < 0) {
      msg << "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
          << " negative integration order given for the force!";
      ++throw_up;
    }
    if (model_->neumann()->order() < 0) {
      msg << "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
          << " negative integration order given for the neumann values!";
      ++throw_up;
    }
    // * parametrization
    if (model_->parametric() && !model_->affineparametric()) {
      msg << "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
          << " only implemented for nonparametric or affineparametric models!";
      ++throw_up;
    }
    if (model_->dirichlet()->parametric()) {
      msg << "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
          << " parametric dirichlet given (not yet implemented)!";
      ++throw_up;
    }
    if (model_->neumann()->parametric()) {
      msg << "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
          << " parametric neumann given (not yet implemented)!";
      ++throw_up;
    }
    if (throw_up)
      DUNE_THROW(Dune::InvalidStateException, msg.str());

    // create the local solvers and build the multiscale the mappers
    ansatzMapper_.prepare();
    testMapper_.prepare();
    // walk all subdomains
    for (unsigned int subdomain = 0; subdomain < msGrid_->size(); ++subdomain) {
      localSolvers_[subdomain] = std::make_shared< LocalSolverType >(msGrid_->localGridPart(subdomain),
                                                                     std::make_shared< LocalBoundaryInfoType >(),
                                                                     model_);
      ansatzMapper_.add(subdomain, localSolvers_[subdomain]->ansatzSpace()->map().size());
      testMapper_.add(subdomain, localSolvers_[subdomain]->testSpace()->map().size());
    } // walk all subdomains
    ansatzMapper_.finalize();
    testMapper_.finalize();
  } // DetailedDiscretizations()

  const std::shared_ptr< const MsGridType > msGrid() const
  {
    return msGrid_;
  }

  const std::shared_ptr< const BoundaryInfoType > boundaryInfo() const
  {
    return boundaryInfo_;
  }

  const std::shared_ptr< const ModelType > model() const
  {
    return model_;
  }

  const RangeFieldType penaltyFactor() const
  {
    return penaltyFactor_;
  }

  std::vector< std::shared_ptr< VectorType > > createVector() const
  {
    std::vector< std::shared_ptr< VectorType > > ret(msGrid_->size());
    for (unsigned int subdomain = 0; subdomain < msGrid_->size(); ++subdomain) {
      ret[subdomain] = localSolvers_[subdomain]->createVector();
    }
    return ret/*std::make_shared< VectorType >(ansatzMapper_.size())*/;
  } // ... createVector() const

  void visualize(const std::vector< std::shared_ptr< VectorType > >& vector,
                 const std::string filename = "solution",
                 const std::string name = "solution",
                 std::ostream& out = Dune::Stuff::Common::Logger().devnull(),
                 const std::string prefix = "") const
  {
    // preparations
    assert(vector.size() == msGrid_->size() && "Given vector has wrong size!");
    Dune::Timer timer;
    out << prefix << "writing '" << name << "'" << std::endl
        << prefix << "     to '" << filename;
    if (dimDomain == 1)
      out << ".vtp";
    else
      out << ".vtu";
    out << "'... " << std::flush;
    // create vector of local discrete functions
    std::vector< std::shared_ptr< LocalDiscreteFunctionType > > localDiscreteFunctions;
    for (unsigned int subdomain = 0; subdomain < msGrid_->size(); ++subdomain) {
      localDiscreteFunctions.push_back(std::shared_ptr< LocalDiscreteFunctionType >(
          new LocalDiscreteFunctionType(*(localSolvers_[subdomain]->ansatzSpace()),
                                        vector[subdomain])));
    }
    // visualize
    visualizeFunction(std::make_shared< DiscreteFunctionType >(*msGrid_, localDiscreteFunctions, name),
                      filename,
                      Dune::Stuff::Common::Logger().devnull(),
                      "");
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void visualize(...) const

  /**
   *  \todo check ansatz vs. test space for pattern!
   *  \attention Only trivial dirichlet atm, no neumann!
   */
  void init(std::ostream& out = Dune::Stuff::Common::Logger().devnull(), const std::string prefix = "")
  {
    if (!initialized_) {
      // prepare
      Dune::Timer timer;
      const size_t subdomains = msGrid_->size();
      const bool verbose = (subdomains >= 3) && (subdomains <= std::pow(3, 3));
      std::map< unsigned int, std::shared_ptr< const PatternType > > boundaryPatternMap;
      std::vector< std::map< unsigned int, std::map< std::string, std::shared_ptr< const PatternType > > > >
          couplingPatternMapMaps(msGrid_->size());
      // initialize the global pattern
      pattern_ = std::make_shared< PatternType >(ansatzMapper_.size());
      out<< " done (took " << timer.elapsed() << " sek)" << std::endl;
      // walk the subdomains for the first time
      //   * to initialize the coupling pattern,
      //   * to initialize the boundary pattern and
      //   * to build up the global sparsity pattern
      out << prefix << "initializing coupling patterns (on " << subdomains << " subdomains)"  << std::flush;
      if (!verbose)
        out << "..." << std::flush;
      timer.reset();
      for (unsigned int subdomain = 0; subdomain < msGrid_->size(); ++subdomain) {
        // init the local solver (assembles matrices and patterns)
        localSolvers_[subdomain]->init();
        // copy the local solvers pattern
        addLocalToGlobalPattern(localSolvers_[subdomain]->pattern("diffusion"),
                                subdomain,
                                subdomain,
                                pattern_);
        // create the boundary pattern
        const typename LocalSolverType::AnsatzSpaceType& innerAnsatzSpace = *(localSolvers_[subdomain]->ansatzSpace());
        const typename LocalSolverType::TestSpaceType& innerTestSpace = *(localSolvers_[subdomain]->testSpace());
        if (msGrid_->boundary(subdomain)) {
          const std::shared_ptr< const PatternType > boundaryPattern
              = innerAnsatzSpace.computeLocalPattern(*(msGrid_->boundaryGridPart(subdomain)),
                                                     innerTestSpace);
          boundaryPatternMap.insert(std::make_pair(subdomain, boundaryPattern));
          // and copy it
          addLocalToGlobalPattern(boundaryPattern,
                                  subdomain,
                                  subdomain,
                                  pattern_);
        } // if (msGrid->boundary(subdomain))
        // walk the neighbors
        std::map< unsigned int, std::map< std::string, std::shared_ptr< const PatternType > > >& couplingPatternMapMap
            = couplingPatternMapMaps[subdomain];
        for (size_t neighboringSubdomain : msGrid_->neighborsOf(subdomain)) {
          // visit each coupling only once (assemble primaly)
          if (subdomain < neighboringSubdomain) {
            // create the coupling patterns
            const typename LocalSolverType::AnsatzSpaceType& outerAnsatzSpace = *(localSolvers_[neighboringSubdomain]->ansatzSpace());
            const typename LocalSolverType::TestSpaceType& outerTestSpace = *(localSolvers_[neighboringSubdomain]->testSpace());
            const typename MsGridType::CouplingGridPartType& insideOutsideGridPart = *(msGrid_->couplingGridPart(subdomain, neighboringSubdomain));
            const typename MsGridType::CouplingGridPartType& outsideInsideGridPart = *(msGrid_->couplingGridPart(neighboringSubdomain, subdomain));
            std::map< std::string, std::shared_ptr< const PatternType > >& couplingPatternMap
                = couplingPatternMapMap[neighboringSubdomain];
            couplingPatternMap.insert(std::pair< std::string, std::shared_ptr< const PatternType > >(
                "inside/inside",
                innerAnsatzSpace.computeLocalPattern(insideOutsideGridPart,
                                                     innerTestSpace)));
            couplingPatternMap.insert(std::pair< std::string, std::shared_ptr< const PatternType > >(
                "inside/outside",
                innerAnsatzSpace.computeCouplingPattern(insideOutsideGridPart,
                                                        outerTestSpace)));
            couplingPatternMap.insert(std::pair< std::string, std::shared_ptr< const PatternType > >(
                "outside/inside",
                outerAnsatzSpace.computeCouplingPattern(outsideInsideGridPart,
                                                        innerTestSpace)));
            couplingPatternMap.insert(std::pair< std::string, std::shared_ptr< const PatternType > >(
                "outside/outside",
                outerAnsatzSpace.computeLocalPattern(outsideInsideGridPart,
                                                     outerTestSpace)));
            // and copy them
            addLocalToGlobalPattern(couplingPatternMap["inside/inside"],
                                    subdomain,
                                    subdomain,
                                    pattern_);
            addLocalToGlobalPattern(couplingPatternMap["inside/outside"],
                                    subdomain,
                                    neighboringSubdomain,
                                    pattern_);
            addLocalToGlobalPattern(couplingPatternMap["outside/inside"],
                                    neighboringSubdomain,
                                    subdomain,
                                    pattern_);
            addLocalToGlobalPattern(couplingPatternMap["outside/outside"],
                                    neighboringSubdomain,
                                    neighboringSubdomain,
                                    pattern_);
          } // visit each coupling only once (assemble primaly)
        } // walk the neighbors
        if (verbose)
          out << "." << std::flush;
      } // walk the subdomains for the first time
      out<< " done (took " << timer.elapsed() << " sek)" << std::endl;
      // walk the subdomains for the second time
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
      if (!model_->diffusion()->parametric()) {
        matrix_ = std::make_shared< AffineParametricMatrixType >(std::make_shared< MatrixType >(ansatzMapper_.size(),
                                                                                                testMapper_.size(),
                                                                                                *pattern_));
      } else {
        // create one matrix for each component
        std::vector< std::shared_ptr< MatrixType > > diffusionMatrices;
        for (size_t qq = 0; qq < model_->diffusion()->numComponents(); ++qq)
          diffusionMatrices.push_back(std::make_shared< MatrixType >(ansatzMapper_.size(),
                                                                     testMapper_.size(),
                                                                     *pattern_));
        // and one for the affine part
        diffusionMatrices.push_back(std::make_shared< MatrixType >(ansatzMapper_.size(),
                                                                   testMapper_.size(),
                                                                   *pattern_));
        matrix_ = std::make_shared< AffineParametricMatrixType >(model_->diffusion()->paramSize(),
                                                                 diffusionMatrices,
                                                                 model_->diffusion()->coefficients());
      } // if (!model_->diffusion()->parametric())
      if (!model_->force()->parametric()) {
        rhs_ = std::make_shared< AffineParametricVectorType >(std::make_shared< VectorType >(testMapper_.size()));
      } else {
        // create one vector for each component
        std::vector< std::shared_ptr< VectorType > > forceVectors;
        for (size_t qq = 0; qq < model_->force()->numComponents(); ++qq)
          forceVectors.push_back(std::make_shared< VectorType >(testMapper_.size()));
        rhs_ = std::make_shared< AffineParametricVectorType >(model_->force()->paramSize(),
                                                              forceVectors,
                                                              model_->force()->coefficients());
      } // if (!model_->force()->parametric())
      for (unsigned int subdomain = 0; subdomain < msGrid_->size(); ++subdomain) {
        // copy the local containers of the subdomain solver
        const std::shared_ptr< const AffineParametricMatrixType >
            localMatrix = localSolvers_[subdomain]->matrix("diffusion");
        assert(matrix_->numComponents() == localMatrix->numComponents());
        for (size_t qq = 0; qq < matrix_->numComponents(); ++qq) {
          copyLocalToGlobalMatrix(localMatrix->components()[qq],
                                  localSolvers_[subdomain]->pattern("diffusion"),
                                  subdomain,
                                  subdomain,
                                  matrix_->components()[qq]);
        }
        const std::shared_ptr< const AffineParametricVectorType >
            localVector = localSolvers_[subdomain]->vector("force");
        assert(rhs_->numComponents() == localVector->numComponents());
        for (size_t qq = 0; qq < rhs_->numComponents(); ++qq)
          copyLocalToGlobalVector(localVector->components()[qq], subdomain, rhs_->components()[qq]);
        // for the boundary contribution
        const auto& innerAnsatzMapper = localSolvers_[subdomain]->ansatzSpace()->map();
        const auto& innerTestMapper = localSolvers_[subdomain]->testSpace()->map();
        if (msGrid_->boundary(subdomain)) {
          //   * initialize the boundary matrix and vector,
          typename std::map< unsigned int, std::shared_ptr< const PatternType > >::const_iterator result = boundaryPatternMap.find(subdomain);
          assert(result != boundaryPatternMap.end());
          const std::shared_ptr< const PatternType > boundaryPattern = result->second;
          std::shared_ptr< AffineParametricMatrixType > boundaryMatrix;
          if (!model_->diffusion()->parametric()) {
            boundaryMatrix
                = std::make_shared< AffineParametricMatrixType >(std::make_shared< MatrixType >(innerAnsatzMapper.size(),
                                                                                                innerTestMapper.size(),
                                                                                                *boundaryPattern));
          } else {
            // create one matrix for each component
            std::vector< std::shared_ptr< MatrixType > > boundaryMatrices;
            for (size_t qq = 0; qq < model_->diffusion()->numComponents(); ++qq)
              boundaryMatrices.push_back(std::make_shared< MatrixType >(innerAnsatzMapper.size(),
                                                                        innerTestMapper.size(),
                                                                        *boundaryPattern));
            boundaryMatrix = std::make_shared< AffineParametricMatrixType >(model_->diffusion()->paramSize(),
                                                                            boundaryMatrices,
                                                                            model_->diffusion()->coefficients());
          } // if (!model_->diffusion()->parametric())
          // we only work with trivial dirichlet at the moment
//          std::shared_ptr< VectorType > boundaryRhs = std::make_shared< VectorType >(innerTestSpace.map().size());
          //   * assemble them
          assembleBoundaryContribution(subdomain, boundaryMatrix/*, boundaryRhs*/);
          //   * and copy them into the global matrix and vector
          for (size_t qq = 0; qq < boundaryMatrix->numComponents(); ++qq) {
            copyLocalToGlobalMatrix(boundaryMatrix->components()[qq],
                                    boundaryPattern,
                                    subdomain,
                                    subdomain,
                                    matrix_->components()[qq]);
          }
//          copyLocalToGlobalVector(boundaryRhs, subdomain, rhs_);
        } // if (msGrid_->boundary(subdomain))
        // walk the neighbors
        std::map< unsigned int, std::map< std::string, std::shared_ptr< const PatternType > > >& couplingPatternMapMap
            = couplingPatternMapMaps[subdomain];
        for (size_t neighboringSubdomain : msGrid_->neighborsOf(subdomain)) {
          // visit each coupling only once (assemble primaly)
          if (subdomain < neighboringSubdomain) {
            // for the coupling contribution
            const auto& outerAnsatzMapper = localSolvers_[neighboringSubdomain]->ansatzSpace()->map();
            const auto& outerTestMapper = localSolvers_[neighboringSubdomain]->testSpace()->map();
            std::map< std::string, std::shared_ptr< const PatternType > >& couplingPatternMap
                = couplingPatternMapMap[neighboringSubdomain];
            const std::shared_ptr< const PatternType >& insideInsidePattern = couplingPatternMap["inside/inside"];
            const std::shared_ptr< const PatternType >& insideOutsidePattern = couplingPatternMap["inside/outside"];
            const std::shared_ptr< const PatternType >& outsideInsidePattern = couplingPatternMap["outside/inside"];
            const std::shared_ptr< const PatternType >& outsideOutsidePattern = couplingPatternMap["outside/outside"];
            // initialize the coupling matrices
            // * inside/inside
            std::shared_ptr< AffineParametricMatrixType > insideInsideMatrix;
            if (!model_->diffusion()->parametric()) {
              insideInsideMatrix
                  = std::make_shared< AffineParametricMatrixType >(std::make_shared< MatrixType >(innerAnsatzMapper.size(),
                                                                                                  innerTestMapper.size(),
                                                                                                  *insideInsidePattern));
            } else {
              // create one matrix for each component
              std::vector< std::shared_ptr< MatrixType > > insideInsideMatrices;
              for (size_t qq = 0; qq < model_->diffusion()->numComponents(); ++qq)
                insideInsideMatrices.push_back(std::make_shared< MatrixType >(innerAnsatzMapper.size(),
                                                                              innerTestMapper.size(),
                                                                              *insideInsidePattern));
              insideInsideMatrix = std::make_shared< AffineParametricMatrixType >(model_->diffusion()->paramSize(),
                                                                                  insideInsideMatrices,
                                                                                  model_->diffusion()->coefficients());
            } // if (!model_->diffusion()->parametric())
            // * inside/outside
            std::shared_ptr< AffineParametricMatrixType > insideOutsideMatrix;
            if (!model_->diffusion()->parametric()) {
              insideOutsideMatrix
                  = std::make_shared< AffineParametricMatrixType >(std::make_shared< MatrixType >(innerAnsatzMapper.size(),
                                                                                                  outerTestMapper.size(),
                                                                                                  *insideOutsidePattern));
            } else {
              // create one matrix for each component
              std::vector< std::shared_ptr< MatrixType > > insideOutsideMatrices;
              for (size_t qq = 0; qq < model_->diffusion()->numComponents(); ++qq)
                insideOutsideMatrices.push_back(std::make_shared< MatrixType >(innerAnsatzMapper.size(),
                                                                               outerTestMapper.size(),
                                                                               *insideOutsidePattern));
              insideOutsideMatrix = std::make_shared< AffineParametricMatrixType >(model_->diffusion()->paramSize(),
                                                                                   insideOutsideMatrices,
                                                                                   model_->diffusion()->coefficients());
            } // if (!model_->diffusion()->parametric())
            // * outside/inside
            std::shared_ptr< AffineParametricMatrixType > outsideInsideMatrix;
            if (!model_->diffusion()->parametric()) {
              outsideInsideMatrix
                  = std::make_shared< AffineParametricMatrixType >(std::make_shared< MatrixType >(outerAnsatzMapper.size(),
                                                                                                  innerTestMapper.size(),
                                                                                                  *outsideInsidePattern));
            } else {
              // create one matrix for each component
              std::vector< std::shared_ptr< MatrixType > > outsideInsideMatrices;
              for (size_t qq = 0; qq < model_->diffusion()->numComponents(); ++qq)
                outsideInsideMatrices.push_back(std::make_shared< MatrixType >(outerAnsatzMapper.size(),
                                                                               innerTestMapper.size(),
                                                                               *outsideInsidePattern));
              outsideInsideMatrix = std::make_shared< AffineParametricMatrixType >(model_->diffusion()->paramSize(),
                                                                                   outsideInsideMatrices,
                                                                                   model_->diffusion()->coefficients());
            } // if (!model_->diffusion()->parametric())
            // * outside/outside
            std::shared_ptr< AffineParametricMatrixType > outsideOutsideMatrix;
            if (!model_->diffusion()->parametric()) {
              outsideOutsideMatrix
                  = std::make_shared< AffineParametricMatrixType >(std::make_shared< MatrixType >(outerAnsatzMapper.size(),
                                                                                                  outerTestMapper.size(),
                                                                                                  *outsideOutsidePattern));
            } else {
              // create one matrix for each component
              std::vector< std::shared_ptr< MatrixType > > outsideOutsideMatrices;
              for (size_t qq = 0; qq < model_->diffusion()->numComponents(); ++qq)
                outsideOutsideMatrices.push_back(std::make_shared< MatrixType >(outerAnsatzMapper.size(),
                                                                                outerTestMapper.size(),
                                                                                *outsideOutsidePattern));
              outsideOutsideMatrix = std::make_shared< AffineParametricMatrixType >(model_->diffusion()->paramSize(),
                                                                                    outsideOutsideMatrices,
                                                                                    model_->diffusion()->coefficients());
            } // if (!model_->diffusion()->parametric())
            //   * assemble them
            assembleCouplingContribution(subdomain,
                                         neighboringSubdomain,
                                         insideInsideMatrix,
                                         insideOutsideMatrix,
                                         outsideInsideMatrix,
                                         outsideOutsideMatrix);
            //   * and copy them into the global matrix
            for (size_t qq = 0; qq < insideInsideMatrix->numComponents(); ++qq)
              copyLocalToGlobalMatrix(insideInsideMatrix->components()[qq],
                                      insideInsidePattern,
                                      subdomain,
                                      subdomain,
                                      matrix_->components()[qq]);
            for (size_t qq = 0; qq < insideInsideMatrix->numComponents(); ++qq)
              copyLocalToGlobalMatrix(insideOutsideMatrix->components()[qq],
                                      insideOutsidePattern,
                                      subdomain,
                                      neighboringSubdomain,
                                      matrix_->components()[qq]);
            for (size_t qq = 0; qq < insideInsideMatrix->numComponents(); ++qq)
              copyLocalToGlobalMatrix(outsideInsideMatrix->components()[qq],
                                      outsideInsidePattern,
                                      neighboringSubdomain,
                                      subdomain,
                                      matrix_->components()[qq]);
            for (size_t qq = 0; qq < insideInsideMatrix->numComponents(); ++qq)
              copyLocalToGlobalMatrix(outsideOutsideMatrix->components()[qq],
                                      outsideOutsidePattern,
                                      neighboringSubdomain,
                                      neighboringSubdomain,
                                      matrix_->components()[qq]);
          } // visit each coupling only once
        } // walk the neighbors
        if (verbose)
          out << "." << std::flush;
      } // walk the subdomains for the second time
      out<< " done (took " << timer.elapsed() << " sek)" << std::endl;
      // done
      initialized_ = true;
    } // if (!initialized_)
  } // void init(...)

  void solve(std::vector< std::shared_ptr< VectorType > >& solution,
             const std::string linearSolverType = "bicgstab.ilut",
             const double linearSolverPrecision = 1e-12,
             const size_t linearSolverMaxIter = 5000,
             std::ostream& out = Dune::Stuff::Common::Logger().debug(),
             const std::string prefix = "") const
  {
    if (!initialized_)
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " call init() before calling solve()!");
    if (model_->parametric())
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " nonparametric solve() called for a parametric model!");
    generic_solve(solution,
                  ParamType(),
                  linearSolverType, linearSolverPrecision, linearSolverMaxIter,
                  out, prefix);
  } // ... solve(...)

  void solve(std::vector< std::shared_ptr< VectorType > >& solutionVector,
             const ParamType& mu,
             const std::string linearSolverType = "bicgstab.diagonal",
             const double linearSolverPrecision = 1e-12,
             const size_t linearSolverMaxIter = 5000,
             std::ostream& out = Dune::Stuff::Common::Logger().devnull(),
             const std::string prefix = "") const
  {
    if (!initialized_)
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " call init() before calling solve()!");
    // check, that we are really in the nonparametric setting!
    if (!model_->parametric())
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " parametric solve() called for a nonparametric model!");
    generic_solve(solutionVector,
                  mu,
                  linearSolverType, linearSolverPrecision, linearSolverMaxIter,
                  out, prefix);
  } // ... solve(..., mu, ...)

private:
  void addLocalToGlobalPattern(const std::shared_ptr< const PatternType >& local,
                               const unsigned int ansatzSubdomain,
                               const unsigned int testSubdomain,
                               std::shared_ptr< PatternType >& global) const
  {
    // loop over all rows of the local pattern
    for (size_t localRowIndex = 0; localRowIndex < local->size(); ++localRowIndex) {
      const size_t globalRowIndex = ansatzMapper_.toGlobal(ansatzSubdomain, localRowIndex);
      const auto& localRowSet = local->set(localRowIndex);
      // get corresponding row of the global pattern (make use of automatic creation of the map with [])
      auto& globalRowSet = global->set(globalRowIndex);
      // loop over all column entries in local row
      for (auto localColIndex : localRowSet) {
        // map local to global
        const size_t globalColIndex = testMapper_.toGlobal(testSubdomain, localColIndex);
        // add colum
        globalRowSet.insert(globalColIndex);
      } // loop over all column entries in args row
    } // loop pver all rows of the input pattern
  } // void addLocalToGlobalPattern(...) const

  void copyLocalToGlobalMatrix(const std::shared_ptr< const MatrixType >& localMatrix,
                               const std::shared_ptr< const PatternType >& localPattern,
                               const unsigned int ansatzSubdomain,
                               const unsigned int testSubdomain,
                               std::shared_ptr< MatrixType >& global) const
  {
    // loop over all local rows
    for (size_t localRow = 0; localRow < localPattern->size(); ++localRow) {
      const unsigned int globalRow = ansatzMapper_.toGlobal(ansatzSubdomain, localRow);
      // loop over all cols in the current row
      for (auto localCol : localPattern->set(localRow)) {
        const size_t globalCol = testMapper_.toGlobal(testSubdomain, localCol);
        // add entry
        global->add(globalRow, globalCol, localMatrix->get(localRow, localCol));
      } // loop over all cols in the current row
    } // loop over all local rows
  } // void copyLocalToGlobalMatrix(...) const

  void copyLocalToGlobalVector(const std::shared_ptr< const VectorType >& local,
                               const unsigned int subdomain,
                               std::shared_ptr< VectorType >& global) const
  {
    for (unsigned int localI = 0; localI < local->size(); ++localI) {
      const unsigned int globalI = testMapper_.toGlobal(subdomain, localI);
      global->add(globalI, local->get(localI));
    }
  } // void copyLocalToGlobalVector(...)

  void assembleBoundaryContribution(const unsigned int subdomain,
                                    std::shared_ptr< AffineParametricMatrixType >& boundaryMatrix/*,
                                    std::shared_ptr< VectorType >& boundaryRhs*/) const
  {
    typedef Dune::FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, dimRange > FunctionSpaceType;
    // operator
    typedef typename ModelType::FunctionType DiffusionType;
    typedef Dune::Detailed::Discretizations::Evaluation::Local::Binary::IPDGfluxes::Dirichlet<  FunctionSpaceType,
                                                                                                DiffusionType >
        IPDGfluxType;
    std::vector< IPDGfluxType* > ipdgFluxes;
    typedef Dune::Detailed::Discretizations::DiscreteOperator::Local::Codim1::BoundaryIntegral< IPDGfluxType >
        DirichletOperatorType;
    std::vector< DirichletOperatorType* > dirichletOperators;
    if (!model_->diffusion()->parametric()) {
      ipdgFluxes.push_back(new IPDGfluxType(model_->diffusion(), model_->diffusion()->order(), penaltyFactor_));
      dirichletOperators.push_back(new DirichletOperatorType(*(ipdgFluxes[0])));
    } else {
      // we are separable (see constructor), loop over all components
      for (size_t qq = 0; qq < model_->diffusion()->numComponents(); ++qq) {
        ipdgFluxes.push_back(new IPDGfluxType(model_->diffusion()->components()[qq],
                                              model_->diffusion()->order(),
                                              penaltyFactor_));
        dirichletOperators.push_back(new DirichletOperatorType(*(ipdgFluxes[qq])));
      } // loop over all components
    } // if (!model_->diffusion()->parametric())
    // local matrix assembler
    typedef Dune::Detailed::Discretizations::Assembler::Local::Codim1::Boundary< DirichletOperatorType >
        LocalDirichletMatrixAssemblerType;
    std::vector< std::shared_ptr< const LocalDirichletMatrixAssemblerType > > localDirichletMatrixAssembler;
    for (size_t qq = 0; qq < dirichletOperators.size(); ++qq)
      localDirichletMatrixAssembler.push_back(
            std::make_shared< LocalDirichletMatrixAssemblerType >(*(dirichletOperators[qq])));
    // functional
    // ...
    // local vector assembler
    // ...
    // boundary assembler
    typedef Dune::Detailed::Discretizations::Assembler::Multiscale::Boundary<
        BoundaryGridPartType,
        BoundaryInfoType,
        typename LocalSolverType::AnsatzSpaceType,
        typename LocalSolverType::TestSpaceType >
      BoundaryAssemblerType;
    BoundaryAssemblerType boundaryAssembler(*(msGrid_->boundaryGridPart(subdomain)),
                                            boundaryInfo_,
                                            *(localSolvers_[subdomain]->ansatzSpace()),
                                            *(localSolvers_[subdomain]->testSpace()));
    assert(boundaryMatrix->numComponents() == localDirichletMatrixAssembler.size());
    for (size_t qq = 0; qq < localDirichletMatrixAssembler.size(); ++qq)
      boundaryAssembler.addLocalMatrixAssembler(localDirichletMatrixAssembler[qq], boundaryMatrix->components()[qq]);
    boundaryAssembler.assemble(/*dirichletMatrixAssembler,
                               *boundaryMatrix*//*,
                               dirichletVectorAssembler,
                               boundaryRhs*/);
    // clean up
    for (size_t ii = 0; ii < dirichletOperators.size(); ++ii)
      delete dirichletOperators[ii];
    for (size_t ii = 0; ii < ipdgFluxes.size(); ++ii)
      delete ipdgFluxes[ii];
  } // void assembleBoundaryContribution(...)

  void assembleCouplingContribution(const unsigned int subdomain,
                                    const unsigned int neighboringSubdomain,
                                    std::shared_ptr< AffineParametricMatrixType >& insideInsideMatrix,
                                    std::shared_ptr< AffineParametricMatrixType >& insideOutsideMatrix,
                                    std::shared_ptr< AffineParametricMatrixType >& outsideInsideMatrix,
                                    std::shared_ptr< AffineParametricMatrixType >& outsideOutsideMatrix) const
  {
    // operator
    typedef Dune::FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, dimRange > FunctionSpaceType;
    typedef typename ModelType::FunctionType DiffusionType;
    typedef Dune::Detailed::Discretizations::Evaluation::Local::Quaternary::IPDGfluxes::Inner<  FunctionSpaceType,
                                                                                                DiffusionType >
        IPDGfluxType;
    std::vector< IPDGfluxType* > ipdgFluxes;
//    const IPDGfluxType ipdgFlux(model_->diffusion(), model_->diffusion()->order(), penaltyFactor_);
    typedef Dune::Detailed::Discretizations::DiscreteOperator::Local::Codim1::InnerIntegral< IPDGfluxType >
        IPDGoperatorType;
    std::vector< IPDGoperatorType* > ipdgOperators;
//    const IPDGoperatorType ipdgOperator(ipdgFlux);
    if (!model_->diffusion()->parametric()) {
      ipdgFluxes.push_back(new IPDGfluxType(model_->diffusion(), model_->diffusion()->order(), penaltyFactor_));
      ipdgOperators.push_back(new IPDGoperatorType(*(ipdgFluxes[0])));
    } else {
      // we are separable (see constructor), loop over all components
      for (size_t qq = 0; qq < model_->diffusion()->numComponents(); ++qq) {
        ipdgFluxes.push_back(new IPDGfluxType(model_->diffusion()->components()[qq],
                                              model_->diffusion()->order(),
                                              penaltyFactor_));
        ipdgOperators.push_back(new IPDGoperatorType(*(ipdgFluxes[qq])));
      } // loop over all components
    } // if (!model_->diffusion()->parametric())
    // local matrix assembler
    typedef Dune::Detailed::Discretizations::Assembler::Local::Codim1::Inner< IPDGoperatorType >
        LocalCouplingMatrixAssemblerType;
    std::vector< std::shared_ptr< const LocalCouplingMatrixAssemblerType > > localCouplingMatrixAssembler;
    for (size_t qq = 0; qq < ipdgOperators.size(); ++qq)
      localCouplingMatrixAssembler.push_back(
            std::make_shared< LocalCouplingMatrixAssemblerType >(*(ipdgOperators[qq])));
    // coupling assembler
    typedef Dune::Detailed::Discretizations::Assembler::Multiscale::Coupling::Primal<
        CouplingGridPartType,
        typename LocalSolverType::AnsatzSpaceType,
        typename LocalSolverType::TestSpaceType,
        typename LocalSolverType::AnsatzSpaceType,
        typename LocalSolverType::TestSpaceType >
      CouplingAssemblerType;
    CouplingAssemblerType couplingAssembler(*(msGrid_->couplingGridPart(subdomain, neighboringSubdomain)),
                                            *(localSolvers_[subdomain]->ansatzSpace()),
                                            *(localSolvers_[subdomain]->testSpace()),
                                            *(localSolvers_[neighboringSubdomain]->ansatzSpace()),
                                            *(localSolvers_[neighboringSubdomain]->testSpace()));
    assert(insideInsideMatrix->numComponents() == localCouplingMatrixAssembler.size());
    assert(insideOutsideMatrix->numComponents() == localCouplingMatrixAssembler.size());
    assert(outsideInsideMatrix->numComponents() == localCouplingMatrixAssembler.size());
    assert(outsideOutsideMatrix->numComponents() == localCouplingMatrixAssembler.size());
    for (size_t qq = 0; qq < localCouplingMatrixAssembler.size(); ++qq)
      couplingAssembler.addLocalMatrixAssembler(localCouplingMatrixAssembler[qq],
                                                insideInsideMatrix->components()[qq],
                                                outsideOutsideMatrix->components()[qq],
                                                insideOutsideMatrix->components()[qq],
                                                outsideInsideMatrix->components()[qq]);
    couplingAssembler.assemble();/*Matrices(couplingMatrixAssembler,
                                       *insideInsideMatrix,
                                       *insideOutsideMatrix,
                                       *outsideInsideMatrix,
                                       *outsideOutsideMatrix);*/
    // clean up
    for (auto ipdgOperator : ipdgOperators)
      delete ipdgOperator;
    for (auto ipdgFlux : ipdgFluxes)
      delete ipdgFlux;
  } // void assembleCouplingContribution(...)

  void generic_solve(std::vector< std::shared_ptr< VectorType > >& solution,
                     const ParamType& mu,
                     const std::string& linearSolverType,
                     const double linearSolverPrecision,
                     const size_t& linearSolverMaxIter,
                     std::ostream& out,
                     const std::string prefix) const
  {
    // first of all, get the corect parameters (the model returns empty ones for nonparametric functions)
    const ParamType muDiffusion = model_->mapParam(mu, "diffusion");
    const ParamType muForce = model_->mapParam(mu, "force");
//    const ParamType muDirichlet = model_->mapParam(mu, "dirichlet");
//    const ParamType muNeumann = model_->mapParam(mu, "neumann");
    Dune::Timer timer;
    out << prefix << "computing system matrix...   " << std::flush;
    std::shared_ptr< const MatrixType > systemMatrix;
    if (model_->diffusion()->parametric())
      systemMatrix = matrix_->fix(muDiffusion);
    else {
      assert(matrix_->numComponents() == 1);
      systemMatrix = matrix_->components()[0];
    }
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
    out << prefix << "computing right hand side...   " << std::flush;
    std::shared_ptr< const VectorType > rightHandSide;
    if (model_->force()->parametric())
      rightHandSide = rhs_->fix(muForce);
    else {
      assert(rhs_->numComponents() == 1);
      rightHandSide = rhs_->components()[0];
    }
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
    out << prefix << "solving linear system (of size " << systemMatrix->rows()
        << "x" << systemMatrix->cols() << ")" << std::endl;
    out << prefix << "  using '" << linearSolverType << "'... " << std::flush;
    timer.reset();
    // create global solution vector
    VectorType tmpSolutionVector(testMapper_.size());
    typedef typename Dune::Stuff::LA::Solver::Interface< MatrixType, VectorType > SolverType;
    const std::shared_ptr< const SolverType >
        solver(Dune::Stuff::LA::Solver::create< MatrixType, VectorType >(linearSolverType));
    const unsigned int failure = solver->apply(*systemMatrix,
                                               *rightHandSide,
                                               tmpSolutionVector,
                                               linearSolverMaxIter,
                                               linearSolverPrecision);
    if (failure)
      DUNE_THROW(Dune::MathError,
                 "\n"
                 << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " linear solver '" << linearSolverType << "' reported error code " << failure << "!\n"
                 << "  1: did not converge\n"
                 << "  2: had numerical issues\n"
                 << "  3: dude, I have no idea");
    if (tmpSolutionVector.size() != int(ansatzMapper_.size()))
      DUNE_THROW(Dune::MathError,
                 "\n"
                 << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " linear solver '" << linearSolverType << "' produced a solution of wrong size (is "
                 << tmpSolutionVector.size() << ", should be " << ansatzMapper_.size() << ")!");
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
    out << prefix << "copying global vector to local...  " << std::flush;
    timer.reset();
    for (unsigned int subdomain = 0; subdomain < msGrid_->size(); ++subdomain) {
      std::shared_ptr< VectorType >& localVector = solution[subdomain];
      assert(localVector->size() == localSolvers_[subdomain]->ansatzSpace()->map().size());
      for (unsigned int localI = 0; localI < localVector->size(); ++localI) {
        const unsigned int globalI = testMapper_.toGlobal(subdomain, localI);
        localVector->set(localI, tmpSolutionVector.get(globalI));
      }
    } // copy global vector to local vector
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void generic_solve(...)

  void visualizeFunction(const std::shared_ptr< const DiscreteFunctionType > discreteFunction,
                 const std::string filename = "discreteFunction",
                 std::ostream& out = Dune::Stuff::Common::Logger().devnull(),
                         const std::string prefix = "") const
  {
    // preparations
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

  const std::shared_ptr< const MsGridType > msGrid_;
  const std::shared_ptr< const BoundaryInfoType > boundaryInfo_;
  const std::shared_ptr< const ModelType > model_;
  const RangeFieldType penaltyFactor_;
  bool initialized_;
  std::vector< std::shared_ptr< LocalSolverType > > localSolvers_;
  AnsatzMapperType ansatzMapper_;
  TestMapperType testMapper_;
  std::shared_ptr< PatternType > pattern_;
  std::shared_ptr< AffineParametricMatrixType > matrix_;
  std::shared_ptr< AffineParametricVectorType > rhs_;
}; // class MultiscaleSolverSemiContinuousGalerkinDD

} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_SOLVER_MULTISCALE_SEMICG_DETAILED_DISCRETIZATIONS_HH
