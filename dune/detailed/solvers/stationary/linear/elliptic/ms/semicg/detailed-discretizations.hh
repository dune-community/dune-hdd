#ifndef DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MS_SEMICG_DETAILED_DISCRETIZATIONS_HH
#define DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MS_SEMICG_DETAILED_DISCRETIZATIONS_HH

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

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
#include <dune/stuff/la/container/separable.hh>

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

#include <dune/detailed/solvers/stationary/linear/elliptic/model/interface.hh>
#include <dune/detailed/solvers/stationary/linear/elliptic/cg/detailed-discretizations.hh>

namespace Dune {
namespace Detailed {
namespace Solvers {
namespace Stationary {
namespace Linear {
namespace Elliptic {
namespace MS {
namespace SemiCG {

/**
 *  \todo Implement non-zero dirichlet and neumann values in assembleBoundaryContribution()
 */
template< class ModelImp, class MsGridImp, class BoundaryInfoImp, int polynomialOrder >
class DetailedDiscretizations
{
public:
  typedef ModelImp        ModelType;
  typedef MsGridImp       MsGridType;
  typedef BoundaryInfoImp BoundaryInfoType;
  static const int        polOrder = polynomialOrder;

  typedef DetailedDiscretizations< ModelType, MsGridType, BoundaryInfoType, polOrder > ThisType;

  typedef typename MsGridType::GlobalGridPartType GlobalGridPartType;
  typedef typename MsGridType::LocalGridPartType  LocalGridPartType;

private:
  typedef typename ModelType::DomainFieldType DomainFieldType;
  static const int                            dimDomain = ModelType::dimDomain;
  typedef typename ModelType::RangeFieldType  RangeFieldType;
  static const int                            dimRange = ModelType::dimRange;
public:
  typedef typename ModelType::ParamType       ParamType;

private:
  typedef typename MsGridType::CouplingGridPartType CouplingGridPartType;
  typedef typename MsGridType::BoundaryGridPartType BoundaryGridPartType;

  typedef Dune::Stuff::Grid::BoundaryInfo::AllNeumann< typename LocalGridPartType::GridViewType > LocalBoundaryInfoType;
public:
  typedef Elliptic::CG::DetailedDiscretizations<  LocalGridPartType,
                                                  polOrder,
                                                  RangeFieldType,
                                                  dimRange >                                      LocalSolverType;
private:

  typedef typename LocalSolverType::AnsatzSpaceType   LocalAnsatzSpaceType;
  typedef typename LocalSolverType::TestSpaceType     LocalTestSpaceType;

public:
  typedef typename LocalSolverType::PatternType PatternType;
  typedef typename LocalSolverType::MatrixType  MatrixType;
  typedef typename LocalSolverType::VectorType  VectorType;

  typedef Dune::Stuff::LA::Container::Separable< MatrixType > SeparableMatrixType;
  typedef Dune::Stuff::LA::Container::Separable< VectorType > SeparableVectorType;

//private:
  typedef Dune::Detailed::Discretizations::Mapper::Multiscale<> AnsatzMapperType;
  typedef Dune::Detailed::Discretizations::Mapper::Multiscale<> TestMapperType;

//public:
  typedef typename LocalSolverType::DiscreteAnsatzFunctionConstType LocalDiscreteFunctionType;
  typedef typename Dune::Detailed::Discretizations::DiscreteFunction::Multiscale< MsGridType, LocalDiscreteFunctionType > DiscreteFunctionType;

  static const std::string id()
  {
    return "detailed.solvers.stationary.linear.elliptic.ms.semicg.detailed-discretizations";
  }

  DetailedDiscretizations(const Dune::shared_ptr< const ModelType > _model,
                          const Dune::shared_ptr< const MsGridType > _msGrid,
                          const Dune::shared_ptr< const BoundaryInfoType > _boundaryInfo,
                          const RangeFieldType _penaltyFactor)
    : model_(_model)
    , msGrid_(_msGrid)
    , boundaryInfo_(_boundaryInfo)
    , penaltyFactor_(_penaltyFactor)
    , initialized_(false)
    , ansatzMapper_()
    , testMapper_()
    , localSolvers_(msGrid_->size())
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
    if (model_->parametric() && !model_->separable()) {
      msg << "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
          << " only implemented for nonparametric or separable-parametric models!";
      ++throw_up;
    }
    if (model_->force()->parametric()) {
      msg << "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
          << " parametric force given!";
      ++throw_up;
    }
    if (model_->dirichlet()->parametric()) {
      msg << "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
          << " parametric dirichlet given!";
      ++throw_up;
    }
    if (model_->neumann()->parametric()) {
      msg << "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
          << " parametric neumann given!";
      ++throw_up;
    }
    if (throw_up)
      DUNE_THROW(Dune::InvalidStateException, msg.str());
  } // DetailedDiscretizations()

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

  const RangeFieldType penaltyFactor() const
  {
    return penaltyFactor_;
  }

  bool parametric() const
  {
    return model_->parametric();
  }

  /**
   *  \todo check ansatz vs. test space for pattern!
   *  \attention Only trivial dirichlet atm, no neumann!
   */
  void init(const std::string prefix = "", std::ostream& out = Dune::Stuff::Common::Logger().debug())
  {
    if (!initialized_) {
      // prepare
      Dune::Timer timer;
      const size_t subdomains = msGrid_->size();
      const bool verbose = (subdomains >= 3) && (subdomains <= std::pow(3, 3));
      std::ostream& devnull = Dune::Stuff::Common::Logger().devnull();
      std::map< unsigned int, Dune::shared_ptr< const PatternType > > boundaryPatternMap;
      std::vector< std::map< unsigned int, std::map< std::string, Dune::shared_ptr< const PatternType > > > >
          couplingPatternMapMaps(msGrid_->size());
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
        localSolvers_[subdomain] = Dune::make_shared< LocalSolverType >(msGrid_->localGridPart(subdomain),
                                                                        localBoundaryInfo,
                                                                        model_);
        localSolvers_[subdomain]->init(prefix + "  ", devnull);
        // initilalize the multiscale mappers
        ansatzMapper_.add(subdomain, localSolvers_[subdomain]->ansatzSpace()->map().size());
        testMapper_.add(subdomain, localSolvers_[subdomain]->testSpace()->map().size());
        if (verbose)
          out << "." << std::flush;
      } // walk the subdomains for the first time
      ansatzMapper_.finalize();
      testMapper_.finalize();
      // initialize the global pattern
      pattern_ = Dune::make_shared< PatternType >(ansatzMapper_.size());
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
        addLocalToGlobalPattern(localSolvers_[subdomain]->systemPattern(),
                                subdomain,
                                subdomain,
                                pattern_);
        // create the boundary pattern
        const typename LocalSolverType::AnsatzSpaceType& innerAnsatzSpace = *(localSolvers_[subdomain]->ansatzSpace());
        const typename LocalSolverType::TestSpaceType& innerTestSpace = *(localSolvers_[subdomain]->testSpace());
        if (msGrid_->boundary(subdomain)) {
          const Dune::shared_ptr< const PatternType > boundaryPattern
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
        std::map< unsigned int, std::map< std::string, Dune::shared_ptr< const PatternType > > >& couplingPatternMapMap
            = couplingPatternMapMaps[subdomain];
        for (size_t neighboringSubdomain : msGrid_->neighborsOf(subdomain)) {
          // visit each coupling only once (assemble primaly)
          if (subdomain < neighboringSubdomain) {
            // create the coupling patterns
            const typename LocalSolverType::TestSpaceType& outerAnsatzSpace = *(localSolvers_[neighboringSubdomain]->ansatzSpace());
            const typename LocalSolverType::TestSpaceType& outerTestSpace = *(localSolvers_[neighboringSubdomain]->testSpace());
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
      if (!model_->diffusion()->parametric()) {
        matrix_ = Dune::make_shared< SeparableMatrixType >(Dune::make_shared< MatrixType >(ansatzMapper_.size(),
                                                                                           testMapper_.size(),
                                                                                           *pattern_));
      } else {
        // create one matrix for each component
        std::vector< Dune::shared_ptr< MatrixType > > diffusionMatrices;
        for (size_t qq = 0; qq < model_->diffusion()->numComponents(); ++qq)
          diffusionMatrices.push_back(Dune::make_shared< MatrixType >(ansatzMapper_.size(),
                                                                      testMapper_.size(),
                                                                      *pattern_));
        matrix_ = Dune::make_shared< SeparableMatrixType >(model_->diffusion()->paramSize(),
                                                           diffusionMatrices,
                                                           model_->diffusion()->coefficients());
      } // if (!model_->diffusion()->parametric())
      rhs_ = Dune::make_shared< VectorType >(testMapper_.size());
      for (unsigned int subdomain = 0; subdomain < msGrid_->size(); ++subdomain) {
        // copy the local containers of the subdomain solver
        const Dune::shared_ptr< const SeparableMatrixType > localMatrix = localSolvers_[subdomain]->systemMatrix("diffusion");
        assert(matrix_->numComponents() == localMatrix->numComponents());
        for (size_t qq = 0; qq < matrix_->numComponents(); ++qq) {
//          writeMatrixToDisc(*(localMatrix->components()[qq]),
//                            *(localSolvers_[subdomain]->systemPattern("diffusion")),
//                            "localMatrix_subdomain_"
//                            + Dune::Stuff::Common::toString(subdomain)
//                            + "_component_"
//                            + Dune::Stuff::Common::toString(qq));
//          std::cout << "== localMatrix_subdomain_"
//                       + Dune::Stuff::Common::toString(subdomain)
//                       + "_component_"
//                       + Dune::Stuff::Common::toString(qq) << " =======================" << std::endl;
//          std::cout << localMatrix->components()[qq]->backend();
          copyLocalToGlobalMatrix(localMatrix->components()[qq],
                                  localSolvers_[subdomain]->systemPattern("diffusion"),
                                  subdomain,
                                  subdomain,
                                  matrix_->components()[qq]);
        }
        copyLocalToGlobalVector(localSolvers_[subdomain]->systemVector("force")->fix(), subdomain, rhs_);
        // for the boundary contribution
        const typename LocalSolverType::AnsatzSpaceType& innerAnsatzSpace = *(localSolvers_[subdomain]->ansatzSpace());
        const typename LocalSolverType::TestSpaceType& innerTestSpace = *(localSolvers_[subdomain]->testSpace());
        if (msGrid_->boundary(subdomain)) {
          //   * initialize the boundary matrix and vector,
          typename std::map< unsigned int, Dune::shared_ptr< const PatternType > >::const_iterator result = boundaryPatternMap.find(subdomain);
          assert(result != boundaryPatternMap.end());
          const Dune::shared_ptr< const PatternType > boundaryPattern = result->second;
          Dune::shared_ptr< SeparableMatrixType > boundaryMatrix;
          if (!model_->diffusion()->parametric()) {
            boundaryMatrix
                = Dune::make_shared< SeparableMatrixType >(Dune::make_shared< MatrixType >(innerAnsatzSpace.map().size(),
                                                                                           innerTestSpace.map().size(),
                                                                                           *boundaryPattern));
          } else {
            // create one matrix for each component
            std::vector< Dune::shared_ptr< MatrixType > > boundaryMatrices;
            for (size_t qq = 0; qq < model_->diffusion()->numComponents(); ++qq)
              boundaryMatrices.push_back(Dune::make_shared< MatrixType >(innerAnsatzSpace.map().size(),
                                                                         innerTestSpace.map().size(),
                                                                         *boundaryPattern));
            boundaryMatrix = Dune::make_shared< SeparableMatrixType >(model_->diffusion()->paramSize(),
                                                                      boundaryMatrices,
                                                                      model_->diffusion()->coefficients());
          } // if (!model_->diffusion()->parametric())
          // we only work with trivial dirichlet at the moment
//          Dune::shared_ptr< VectorType > boundaryRhs = Dune::make_shared< VectorType >(innerTestSpace.map().size());
          //   * assemble them
          assembleBoundaryContribution(subdomain, boundaryMatrix/*, boundaryRhs*/);
          //   * and copy them into the global matrix and vector
          for (size_t qq = 0; qq < boundaryMatrix->numComponents(); ++qq) {
//            writeMatrixToDisc(*(boundaryMatrix->components()[qq]),
//                              *boundaryPattern,
//                              "boundaryMatrix_subdomain_"
//                              + Dune::Stuff::Common::toString(subdomain)
//                              + "_component_"
//                              + Dune::Stuff::Common::toString(qq));
//            std::cout << "== boundaryMatrix_subdomain_"
//                         + Dune::Stuff::Common::toString(subdomain)
//                         + "_component_"
//                         + Dune::Stuff::Common::toString(qq) << " =======================" << std::endl;
//            std::cout << localMatrix->components()[qq]->backend();
            copyLocalToGlobalMatrix(boundaryMatrix->components()[qq],
                                    boundaryPattern,
                                    subdomain,
                                    subdomain,
                                    matrix_->components()[qq]);
          }
//          copyLocalToGlobalVector(boundaryRhs, subdomain, rhs_);
        } // if (msGrid_->boundary(subdomain))
        // walk the neighbors
        std::map< unsigned int, std::map< std::string, Dune::shared_ptr< const PatternType > > >& couplingPatternMapMap
            = couplingPatternMapMaps[subdomain];
        for (size_t neighboringSubdomain : msGrid_->neighborsOf(subdomain)) {
          // visit each coupling only once (assemble primaly)
          if (subdomain < neighboringSubdomain) {
            // for the coupling contribution
            const typename LocalSolverType::TestSpaceType& outerAnsatzSpace = *(localSolvers_[neighboringSubdomain]->ansatzSpace());
            const typename LocalSolverType::TestSpaceType& outerTestSpace = *(localSolvers_[neighboringSubdomain]->testSpace());
            std::map< std::string, Dune::shared_ptr< const PatternType > >& couplingPatternMap
                = couplingPatternMapMap[neighboringSubdomain];
            const Dune::shared_ptr< const PatternType >& insideInsidePattern = couplingPatternMap["inside/inside"];
            const Dune::shared_ptr< const PatternType >& insideOutsidePattern = couplingPatternMap["inside/outside"];
            const Dune::shared_ptr< const PatternType >& outsideInsidePattern = couplingPatternMap["outside/inside"];
            const Dune::shared_ptr< const PatternType >& outsideOutsidePattern = couplingPatternMap["outside/outside"];
            // initialize the coupling matrices
            // * inside/inside
            Dune::shared_ptr< SeparableMatrixType > insideInsideMatrix;
            if (!model_->diffusion()->parametric()) {
              insideInsideMatrix
                  = Dune::make_shared< SeparableMatrixType >(Dune::make_shared< MatrixType >(innerAnsatzSpace.map().size(),
                                                                                             innerTestSpace.map().size(),
                                                                                             *insideInsidePattern));
            } else {
              // create one matrix for each component
              std::vector< Dune::shared_ptr< MatrixType > > insideInsideMatrices;
              for (size_t qq = 0; qq < model_->diffusion()->numComponents(); ++qq)
                insideInsideMatrices.push_back(Dune::make_shared< MatrixType >(innerAnsatzSpace.map().size(),
                                                                               innerTestSpace.map().size(),
                                                                               *insideInsidePattern));
              insideInsideMatrix = Dune::make_shared< SeparableMatrixType >(model_->diffusion()->paramSize(),
                                                                            insideInsideMatrices,
                                                                            model_->diffusion()->coefficients());
            } // if (!model_->diffusion()->parametric())
            // * inside/outside
            Dune::shared_ptr< SeparableMatrixType > insideOutsideMatrix;
            if (!model_->diffusion()->parametric()) {
              insideOutsideMatrix
                  = Dune::make_shared< SeparableMatrixType >(Dune::make_shared< MatrixType >(innerAnsatzSpace.map().size(),
                                                                                             outerTestSpace.map().size(),
                                                                                             *insideOutsidePattern));
            } else {
              // create one matrix for each component
              std::vector< Dune::shared_ptr< MatrixType > > insideOutsideMatrices;
              for (size_t qq = 0; qq < model_->diffusion()->numComponents(); ++qq)
                insideOutsideMatrices.push_back(Dune::make_shared< MatrixType >(innerAnsatzSpace.map().size(),
                                                                               outerTestSpace.map().size(),
                                                                               *insideOutsidePattern));
              insideOutsideMatrix = Dune::make_shared< SeparableMatrixType >(model_->diffusion()->paramSize(),
                                                                            insideOutsideMatrices,
                                                                            model_->diffusion()->coefficients());
            } // if (!model_->diffusion()->parametric())
            // * outside/inside
            Dune::shared_ptr< SeparableMatrixType > outsideInsideMatrix;
            if (!model_->diffusion()->parametric()) {
              outsideInsideMatrix
                  = Dune::make_shared< SeparableMatrixType >(Dune::make_shared< MatrixType >(outerAnsatzSpace.map().size(),
                                                                                             innerTestSpace.map().size(),
                                                                                             *outsideInsidePattern));
            } else {
              // create one matrix for each component
              std::vector< Dune::shared_ptr< MatrixType > > outsideInsideMatrices;
              for (size_t qq = 0; qq < model_->diffusion()->numComponents(); ++qq)
                outsideInsideMatrices.push_back(Dune::make_shared< MatrixType >(outerAnsatzSpace.map().size(),
                                                                               innerTestSpace.map().size(),
                                                                               *outsideInsidePattern));
              outsideInsideMatrix = Dune::make_shared< SeparableMatrixType >(model_->diffusion()->paramSize(),
                                                                            outsideInsideMatrices,
                                                                            model_->diffusion()->coefficients());
            } // if (!model_->diffusion()->parametric())
            // * outside/outside
            Dune::shared_ptr< SeparableMatrixType > outsideOutsideMatrix;
            if (!model_->diffusion()->parametric()) {
              outsideOutsideMatrix
                  = Dune::make_shared< SeparableMatrixType >(Dune::make_shared< MatrixType >(outerAnsatzSpace.map().size(),
                                                                                             outerTestSpace.map().size(),
                                                                                             *outsideOutsidePattern));
            } else {
              // create one matrix for each component
              std::vector< Dune::shared_ptr< MatrixType > > outsideOutsideMatrices;
              for (size_t qq = 0; qq < model_->diffusion()->numComponents(); ++qq)
                outsideOutsideMatrices.push_back(Dune::make_shared< MatrixType >(outerAnsatzSpace.map().size(),
                                                                               outerTestSpace.map().size(),
                                                                               *outsideOutsidePattern));
              outsideOutsideMatrix = Dune::make_shared< SeparableMatrixType >(model_->diffusion()->paramSize(),
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
      } // walk the subdomains for the third time
      out<< " done (took " << timer.elapsed() << " sek)" << std::endl;
      // done
      initialized_ = true;
    } // if (!initialized_)
  } // void init()

  Dune::shared_ptr< const LocalSolverType > localSolver(const size_t subdomain) const
  {
    if (!initialized_)
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " call init() before calling localSolver()!");
    assert(subdomain < localSolvers_.size());
    return localSolvers_[subdomain];
  }

  const AnsatzMapperType ansatzMapper() const
  {
    if (!initialized_)
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " call init() before calling ansatzMapper()!");
    return ansatzMapper_;
  }

  const AnsatzMapperType testMapper() const
  {
    if (!initialized_)
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " call init() before calling testMapper()!");
    return testMapper_;
  }

  std::vector< Dune::shared_ptr< VectorType > > createVector() const
  {
    assert(initialized_ && "Please call init() beafore calling createVector()!");
    std::vector< Dune::shared_ptr< VectorType > > ret(msGrid_->size());
    for (unsigned int subdomain = 0; subdomain < msGrid_->size(); ++subdomain) {
      const LocalSolverType& localSolver = *(localSolvers_[subdomain]);
      ret[subdomain] = localSolver.createAnsatzVector();
    }
    return ret;
  } // ... createVector() const

//  Dune::shared_ptr< DiscreteFunctionType > createDiscreteFunction(const std::string name = "discrete_function") const
//  {
//    assert(initialized_ && "Please call init() beafore calling createDiscreteFunction()!");
//    // create vector of local discrete functions
//    std::vector< Dune::shared_ptr< LocalDiscreteFunctionType > > localDiscreteFunctions;
//    for (unsigned int subdomain = 0; subdomain < msGrid_->size(); ++subdomain) {
//      localDiscreteFunctions.push_back(Dune::shared_ptr< LocalDiscreteFunctionType >(
//          new LocalDiscreteFunctionType(localSolvers_[subdomain]->ansatzSpace())));
//    }
//    // create multiscale discrete function
//    Dune::shared_ptr< DiscreteFunctionType > discreteFunction(new DiscreteFunctionType(*msGrid_,
//                                                                                       localDiscreteFunctions,
//                                                                                       name));
//    return discreteFunction;
//  }

  Dune::shared_ptr< DiscreteFunctionType > createDiscreteFunction(std::vector< Dune::shared_ptr< VectorType > >& vectors,
                                                                  const std::string name = "discrete_function") const
  {
    assert(initialized_ && "Please call init() beafore calling createDiscreteFunction()!");
    // create vector of local discrete functions
    std::vector< Dune::shared_ptr< LocalDiscreteFunctionType > > localDiscreteFunctions;
    for (unsigned int subdomain = 0; subdomain < msGrid_->size(); ++subdomain) {
      localDiscreteFunctions.push_back(Dune::shared_ptr< LocalDiscreteFunctionType >(
          new LocalDiscreteFunctionType(*(localSolvers_[subdomain]->ansatzSpace()),
                                        vectors[subdomain])));
    }
    // create multiscale discrete function
    Dune::shared_ptr< DiscreteFunctionType > discreteFunction(new DiscreteFunctionType(*msGrid_,
                                                                                       localDiscreteFunctions,
                                                                                       name));
    return discreteFunction;
  }

  std::vector< Dune::shared_ptr< MatrixType > > systemMatrices(const ParamType& mu) const
  {
    if (!initialized_)
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " call init() before calling systemMatrices()!");
    std::vector< Dune::shared_ptr< MatrixType > > matrices(msGrid_->size());
    for (size_t ii = 0; ii < msGrid_->size(); ++ii)
      matrices[ii] = localSolvers_[ii]->systemMatrix(mu);
    return matrices;
  }

  Dune::shared_ptr< MatrixType > systemMatrix(const ParamType& mu)
  {
    if (!initialized_)
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " call init() before calling systemMatrix()!");
    return matrix_->fix(model_->mapParam(mu, "diffusion"));
  }

  Dune::shared_ptr< const SeparableMatrixType > systemMatrix(const std::string type) const
  {
    if (!initialized_)
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " call init() before calling systemMatrix()!");
    assert(type == "diffusion");
    return matrix_;
  }

  Dune::shared_ptr< const VectorType > systemVector(const std::string type) const
  {
    if (!initialized_)
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " call init() before calling systemVector()!");
    assert(type == "force");
    return rhs_;
  }

  Dune::shared_ptr< VectorType > vector(const std::string type) const
  {
    if (!initialized_)
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " call init() before calling vector()!");
    assert(type == "force");
    return rhs_;
  }

  Dune::shared_ptr< SeparableMatrixType > matrix(const std::string type) const
  {
    if (!initialized_)
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " call init() before calling matrix()!");
    assert(type == "diffusion");
    return matrix_;
  }

  void solve(std::vector< Dune::shared_ptr< VectorType > >& solutionVector,
             const std::string linearSolverType = "bicgstab.diagonal",
             const size_t linearSolverMaxIter = 5000,
             const double linearSolverPrecision = 1e-12,
             const std::string prefix = "",
             std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
  {
    if (!initialized_)
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " call init() before calling solve()!");
    // check, that we are really in the nonparametric setting!
    if (parametric())
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " nonparametric solve() called for a parametric model!");
    generic_solve(solutionVector,
                  ParamType(),
                  linearSolverType, linearSolverMaxIter, linearSolverPrecision,
                  prefix, out);
  }

  void solve(std::vector< Dune::shared_ptr< VectorType > >& solutionVector,
             const ParamType& mu,
             const std::string linearSolverType = "bicgstab.diagonal",
             const size_t linearSolverMaxIter = 5000,
             const double linearSolverPrecision = 1e-12,
             const std::string prefix = "",
             std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
  {
    if (!initialized_)
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " call init() before calling solve()!");
    // check, that we are really in the nonparametric setting!
    if (!parametric())
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " parametric solve() called for a nonparametric model!");
    generic_solve(solutionVector,
                  mu,
                  linearSolverType, linearSolverMaxIter, linearSolverPrecision,
                  prefix, out);
  }

private:
  void generic_solve(std::vector< Dune::shared_ptr< VectorType > >& solutionVector,
                     const ParamType& mu,
                     const std::string& linearSolverType,
                     const size_t& linearSolverMaxIter,
                     const double linearSolverPrecision,
                     const std::string prefix,
                     std::ostream& out) const
  {
    // first of all, get the corect parameters (the model returns empty ones for nonparametric functions)
    const ParamType muDiffusion = model_->mapParam(mu, "diffusion");
//    const ParamType muForce = model_->mapParam(mu, "force");
//    const ParamType muDirichlet = model_->mapParam(mu, "dirichlet");
//    const ParamType muNeumann = model_->mapParam(mu, "neumann");
    Dune::Timer timer;
    out << prefix << "computing system matrix...   " << std::flush;
    Dune::shared_ptr< MatrixType > systemMatrix = matrix_->fix(muDiffusion);
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
    out << prefix << "solving linear system (of size " << systemMatrix->rows()
        << "x" << systemMatrix->cols() << ")" << std::endl;
    out << prefix << "  using '" << linearSolverType << "'... " << std::flush;
    timer.reset();
    // create global solution vector
    VectorType tmpSolutionVector(testMapper_.size());
    typedef typename Dune::Stuff::LA::Solver::Interface< MatrixType, VectorType > SolverType;
    const Dune::shared_ptr< const SolverType > solver(Dune::Stuff::LA::Solver::create< MatrixType, VectorType >(linearSolverType));
//    writeMatrixToDisc(*systemMatrix,
//                      *pattern_,
//                      "matrix.out");
//    writeVectorToDisc(*rhs_,
//                      "vector.out");
//    std::cout << "== system matrix ==============================" << std::endl;
//    std::cout << systemMatrix->backend() << std::endl;
//    std::cout << "== right hand side ==============================" << std::endl;
//    std::cout << rhs_->backend().transpose() << std::endl;
//    std::cout << "================================" << std::endl;
    const unsigned int failure = solver->apply(*systemMatrix,
                                               *rhs_,
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
      Dune::shared_ptr< VectorType >& localVector = solutionVector[subdomain];
      assert(localVector->size() == localSolvers_[subdomain]->ansatzSpace()->map().size());
      for (unsigned int localI = 0; localI < localVector->size(); ++localI) {
        const unsigned int globalI = testMapper_.toGlobal(subdomain, localI);
        localVector->set(localI, tmpSolutionVector.get(globalI));
      }
    } // copy global vector to local vector
    out << "done (took " << timer.elapsed() << " sec)" << std::endl;
  } // void solve(...)

public:
  void visualize(const std::vector< Dune::shared_ptr< VectorType > >& vector,
                 const std::string filename = "solution",
                 const std::string name = "solution",
                 const std::string prefix = "",
                 std::ostream& out = Dune::Stuff::Common::Logger().debug()) const
  {
    // preparations
    assert(initialized_ && "A vector can only be visualized after init() has been called! ");
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
    std::vector< Dune::shared_ptr< LocalDiscreteFunctionType > > localDiscreteFunctions;
    for (unsigned int subdomain = 0; subdomain < msGrid_->size(); ++subdomain) {
      localDiscreteFunctions.push_back(Dune::shared_ptr< LocalDiscreteFunctionType >(
          new LocalDiscreteFunctionType(*(localSolvers_[subdomain]->ansatzSpace()),
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

//  const Dune::shared_ptr< const LocalSolverType > localSolver(const unsigned int subdomain) const
//  {
//    assert(initialized_ && "Please call init() before calling localSolver()!");
//    assert(subdomain < msGrid_->size());
//    return localSolvers_[subdomain];
//  }

//  Dune::shared_ptr< LocalSolverType > localSolver(const unsigned int subdomain)
//  {
//    assert(initialized_ && "Please call init() before calling localSolver()!");
//    assert(subdomain < msGrid_->size());
//    return localSolvers_[subdomain];
//  }

//  const Dune::shared_ptr< const MatrixBackendType > systemMatrix() const
//  {
//    assert(initialized_);
//    return matrix_;
//  }

//  Dune::shared_ptr< MatrixBackendType > systemMatrix()
//  {
//    assert(initialized_);
//    return matrix_;
//  }

//  const Dune::shared_ptr< const VectorBackendType > rightHandSide() const
//  {
//    assert(initialized_);
//    return rhs_;
//  }

//  Dune::shared_ptr< VectorBackendType > rightHandSide()
//  {
//    assert(initialized_);
//    return rhs_;
//  }

//  const Dune::shared_ptr< const PatternType > pattern() const
//  {
//    assert(initialized_);
//    return pattern_;
//  }

  Dune::shared_ptr< VectorType > globalizeVectors(const std::vector< Dune::shared_ptr< VectorType > >& localVectors) const
  {
    if (!initialized_)
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " call init() before calling globalizeVectors()!");
    assert(localVectors.size() == msGrid_->size());
    Dune::shared_ptr< VectorType > globalVector = Dune::make_shared< VectorType >(ansatzMapper_.size());
    for (size_t ii = 0; ii < msGrid_->size(); ++ii) {
      const Dune::shared_ptr< const VectorType > localVector = localVectors[ii];
      assert(localVector->size() == localSolvers_[ii]->ansatzSpace()->map().size());
      copyLocalToGlobalVector(localVector, ii, globalVector);
    }
    return globalVector;
  } // ... globalizeVectors(...)

  void globalizeVector(const Dune::shared_ptr< const VectorType >& localVector,
                       const size_t subdomain,
                       Dune::shared_ptr< VectorType >& globalVector) const
  {
    if (!initialized_)
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " call init() before calling globalizeVectors()!");
    if (globalVector->size() != int(ansatzMapper_.size()))
      globalVector->backend().resize(ansatzMapper_.size(), 1);
    globalVector->backend().setZero();
    copyLocalToGlobalVector(localVector, subdomain, globalVector);
  } // ... globalizeVector(...)

  void localizeVector(const Dune::shared_ptr< const VectorType >& globalVector,
                      std::vector< Dune::shared_ptr< VectorType > >& localVectors) const
  {
    if (!initialized_)
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " call init() before calling globalizeVectors()!");
    assert(localVectors.size() == msGrid_->size());
    for (size_t ss = 0; ss < msGrid_->size(); ++ss) {
      VectorType& localVector = *(localVectors[ss]);
      if (localVector.size() != int(localSolvers_[ss]->ansatzSpace()->map().size()))
        localVector.backend().resize(localSolvers_[ss]->ansatzSpace()->map().size(), 1);
      for (size_t localII = 0; localII < localSolvers_[ss]->ansatzSpace()->map().size(); ++localII)
        localVector.set(localII, globalVector->get(ansatzMapper_.toGlobal(ss, localII)));
    }
  }

private:
  void addLocalToGlobalPattern(const Dune::shared_ptr< const PatternType >& local,
                               const unsigned int ansatzSubdomain,
                               const unsigned int testSubdomain,
                               Dune::shared_ptr< PatternType >& global) const
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

  void copyLocalToGlobalMatrix(const Dune::shared_ptr< const MatrixType >& localMatrix,
                               const Dune::shared_ptr< const PatternType >& localPattern,
                               const unsigned int ansatzSubdomain,
                               const unsigned int testSubdomain,
                               Dune::shared_ptr< MatrixType >& global) const
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

  void copyLocalToGlobalVector(const Dune::shared_ptr< const VectorType >& local,
                               const unsigned int subdomain,
                               Dune::shared_ptr< VectorType >& global) const
  {
    for (unsigned int localI = 0; localI < local->size(); ++localI) {
      const unsigned int globalI = testMapper_.toGlobal(subdomain, localI);
      global->add(globalI, local->get(localI));
    }
  } // void copyLocalToGlobalVector(...)

  void assembleBoundaryContribution(const unsigned int subdomain,
                                    Dune::shared_ptr< SeparableMatrixType >& boundaryMatrix/*,
                                    Dune::shared_ptr< VectorType >& boundaryRhs*/) const
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
    std::vector< Dune::shared_ptr< const LocalDirichletMatrixAssemblerType > > localDirichletMatrixAssembler;
    for (size_t qq = 0; qq < dirichletOperators.size(); ++qq)
      localDirichletMatrixAssembler.push_back(
            Dune::make_shared< LocalDirichletMatrixAssemblerType >(*(dirichletOperators[qq])));
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
                                    Dune::shared_ptr< SeparableMatrixType >& insideInsideMatrix,
                                    Dune::shared_ptr< SeparableMatrixType >& insideOutsideMatrix,
                                    Dune::shared_ptr< SeparableMatrixType >& outsideInsideMatrix,
                                    Dune::shared_ptr< SeparableMatrixType >& outsideOutsideMatrix) const
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
    std::vector< Dune::shared_ptr< const LocalCouplingMatrixAssemblerType > > localCouplingMatrixAssembler;
    for (size_t qq = 0; qq < ipdgOperators.size(); ++qq)
      localCouplingMatrixAssembler.push_back(
            Dune::make_shared< LocalCouplingMatrixAssemblerType >(*(ipdgOperators[qq])));
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

  void writeMatrixToDisc(const MatrixType& matrix,
                         const PatternType& pattern,
                         const std::string filename) const
  {
    std::ofstream file(filename);
    for (size_t ii = 0; ii < pattern.size(); ++ii)
      for (size_t jj : pattern.set(ii))
        file << ii << " " << jj << " " << matrix.get(ii,jj) << std::endl;
  }

  void writeVectorToDisc(const VectorType& vector,
                         const std::string filename) const
  {
    std::ofstream file(filename);
    for (int ii = 0; ii < vector.size(); ++ii)
      file << ii << " " << vector.get(ii) << std::endl;
  }

  const Dune::shared_ptr< const ModelType > model_;
  const Dune::shared_ptr< const MsGridType > msGrid_;
  const Dune::shared_ptr< const BoundaryInfoType > boundaryInfo_;
  const RangeFieldType penaltyFactor_;
  bool initialized_;
  AnsatzMapperType ansatzMapper_;
  TestMapperType testMapper_;
  std::vector< Dune::shared_ptr< LocalSolverType > > localSolvers_;
  Dune::shared_ptr< PatternType > pattern_;
  Dune::shared_ptr< SeparableMatrixType > matrix_;
  Dune::shared_ptr< VectorType > rhs_;
}; // class DetailedDiscretizations

} // namespace SemiCG
} // namespace MS
} // namespace Elliptic
} // namespace Linear
} // namespace Stationary
} // namespace Solvers
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MS_SEMICG_DETAILED_DISCRETIZATIONS_HH
