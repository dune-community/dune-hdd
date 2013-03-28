#include <dune/stuff/test/test_common.hh>

#include <utility>
#include <memory>

#include <dune/grid/sgrid.hh>

#include <dune/grid/part/leaf.hh>

#include <dune/stuff/grid/provider.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/common/logging.hh>

#include <dune/detailed-solvers/linearelliptic/model.hh>
#include <dune/detailed-solvers/linearelliptic/solver/interface.hh>
#include <dune/detailed-solvers/linearelliptic/solver/cg/detailed-discretizations.hh>
#include <dune/detailed-solvers/linearelliptic/solver/fv/pdelab.hh>

namespace Stuff = Dune::Stuff;
namespace Elliptic = Dune::DetailedSolvers::LinearElliptic;

static const int            dimDomain = 2;
typedef Dune::SGrid< 2, 2 > GridType;

typedef typename GridType::ctype  DomainFieldType;
typedef double                    RangeFieldType;
static const int                  dimRange = 1;

typedef Dune::grid::Part::Leaf::Const< GridType > GridPartType;

typedef testing::Types< std::pair<  Elliptic::SolverContinuousGalerkinDD< GridPartType, RangeFieldType, dimRange, 1 >,
                                    Elliptic::SolverContinuousGalerkinDDTraits< GridPartType, RangeFieldType, dimRange, 1 >
                                  >
//                      , std::pair<  Elliptic::SolverFiniteVolumePdelab< GridPartType, RangeFieldType, dimRange >,
//                                    Elliptic::SolverFiniteVolumePdelabTraits< GridPartType, RangeFieldType, dimRange >
//                                  >
                      > SolverTypes;

template< class T >
struct SolverCRTPtest
  : public ::testing::Test
{
  typedef typename T::first_type  SolverType;
  typedef typename T::second_type Traits;
  typedef Elliptic::SolverInterface< Traits > SolverInterfaceType;
  typedef Stuff::GridProviderInterface< GridType > GridProviderType;
  typedef typename GridPartType::GridViewType GridViewType;
  typedef Stuff::GridboundaryInterface< GridViewType > GridBoundaryinfoType;
  typedef Elliptic::ModelInterface< DomainFieldType, dimDomain, RangeFieldType, dimRange > ModelType;

  void check() const
  {
    typedef Stuff::GridProviders< GridType > GridProviders;
    const std::shared_ptr< const GridProviderType >
        gridProvider(GridProviders::create(GridProviders::available()[0],
                                           GridProviders::createSampleDescription(GridProviders::available()[0])));
    const std::shared_ptr< const GridPartType > gridPart(new GridPartType(*gridProvider->grid()));
    typedef Stuff::Gridboundaries< GridViewType > Gridboundaries;
    const std::shared_ptr< const GridBoundaryinfoType >
        boundaryInfo(Gridboundaries::create(Gridboundaries::available()[0],
                                            Gridboundaries::createSampleDescription(Gridboundaries::available()[0])));
    typedef Elliptic::Models< DomainFieldType, dimDomain, RangeFieldType, dimRange > EllipticModels;
    const std::shared_ptr< const ModelType >
        model(EllipticModels::create(EllipticModels::available()[0],
                                     EllipticModels::createSampleDescription(EllipticModels::available()[0])));
    SolverType solver(gridPart, boundaryInfo, model);
    SolverInterfaceType& solverAsInterface = static_cast< SolverInterfaceType& >(solver);
    // check for static information
    typedef typename SolverType::GridPartType TestGridPartType;
    const int DUNE_UNUSED(testPolOrder) = SolverType::polOrder;
    typedef typename SolverType::DomainFieldType TestDomainFieldType;
    const int DUNE_UNUSED(testDimDomain) = SolverType::dimDomain;
    typedef typename SolverType::RangeFieldType TestRangeFieldType;
    const int DUNE_UNUSED(testDimRange) = SolverType::dimRange;
    typedef typename SolverType::BoundaryInfoType TestBoundaryInfoType;
    typedef typename SolverType::ModelType TestModelType;
    typedef typename SolverType::VectorType VectorType;
    const std::string DUNE_UNUSED(id) = SolverType::id();
    // check for functionality
    const std::shared_ptr< const TestGridPartType > DUNE_UNUSED(testGridPart) = solverAsInterface.gridPart();
    const std::shared_ptr< const TestBoundaryInfoType > DUNE_UNUSED(testBoundaryInfo) = solverAsInterface.boundaryInfo();
    const std::shared_ptr< const TestModelType > DUNE_UNUSED(testModel) = solverAsInterface.model();
    std::shared_ptr< VectorType > vector = solverAsInterface.createVector();
    solverAsInterface.visualize(vector, "test_empty_vector", "empty vector", Dune::Stuff::Common::Logger().devnull(), "");
    solverAsInterface.init(Dune::Stuff::Common::Logger().devnull(), "");
    if (!solverAsInterface.initialized())
      DUNE_THROW(Dune::InvalidStateException, "This should not happen!");
    solverAsInterface.solve(vector, "bicgstab.ilut", 1e-6, 1000, Dune::Stuff::Common::Logger().devnull(), "");
    solverAsInterface.visualize(vector, "test_solution", "solution", Dune::Stuff::Common::Logger().devnull(), "");
  }
}; // struct SolverCRTPtest

TYPED_TEST_CASE(SolverCRTPtest, SolverTypes);
TYPED_TEST(SolverCRTPtest, SolverCRTP)
{
    this->check();
}


typedef testing::Types< std::pair<  Elliptic::SolverContinuousGalerkinDD< GridPartType, RangeFieldType, dimRange, 1 >,
                                    Elliptic::SolverContinuousGalerkinDDTraits< GridPartType, RangeFieldType, dimRange, 1 >
                                  >
//                      , std::pair<  Elliptic::SolverFiniteVolumePdelab< GridPartType, RangeFieldType, dimRange >,
//                                    Elliptic::SolverFiniteVolumePdelabTraits< GridPartType, RangeFieldType, dimRange >
//                                  >
                      > ParametricSolverTypes;

template< class T >
struct ParametricSolverCRTPtest
  : public ::testing::Test
{
  typedef typename T::first_type  SolverType;
  typedef typename T::second_type Traits;
  typedef Elliptic::SolverParametricInterface< Traits > SolverInterfaceType;
  typedef Stuff::GridProviderInterface< GridType > GridProviderType;
  typedef typename GridPartType::GridViewType GridViewType;
  typedef Stuff::GridboundaryInterface< GridViewType > GridBoundaryinfoType;
  typedef Elliptic::ModelInterface< DomainFieldType, dimDomain, RangeFieldType, dimRange > ModelType;

  void check() const
  {
    typedef Stuff::GridProviders< GridType > GridProviders;
    const std::shared_ptr< const GridProviderType >
        gridProvider(GridProviders::create(GridProviders::available()[0],
                                           GridProviders::createSampleDescription(GridProviders::available()[0])));
    const std::shared_ptr< const GridPartType > gridPart(new GridPartType(*gridProvider->grid()));
    typedef Stuff::Gridboundaries< GridViewType > Gridboundaries;
    const std::shared_ptr< const GridBoundaryinfoType >
        boundaryInfo(Gridboundaries::create(Gridboundaries::available()[0],
                                            Gridboundaries::createSampleDescription(Gridboundaries::available()[0])));
    typedef Elliptic::Models< DomainFieldType, dimDomain, RangeFieldType, dimRange > EllipticModels;
    const std::string modelType = "model.linearelliptic.affineparametric.thermalblock";
    const std::shared_ptr< const ModelType >
        model(EllipticModels::create(modelType,
                                     EllipticModels::createSampleDescription(modelType)));
    SolverType solver(gridPart, boundaryInfo, model);
    SolverInterfaceType& solverAsInterface = static_cast< SolverInterfaceType& >(solver);
    // check for static information
    typedef typename SolverType::GridPartType TestGridPartType;
    const int DUNE_UNUSED(testPolOrder) = SolverType::polOrder;
    typedef typename SolverType::DomainFieldType TestDomainFieldType;
    const int DUNE_UNUSED(testDimDomain) = SolverType::dimDomain;
    typedef typename SolverType::RangeFieldType TestRangeFieldType;
    const int DUNE_UNUSED(testDimRange) = SolverType::dimRange;
    typedef typename SolverType::BoundaryInfoType TestBoundaryInfoType;
    typedef typename SolverType::ModelType TestModelType;
    typedef typename SolverType::ParamType ParamType;
    typedef typename SolverType::VectorType VectorType;
    const std::string DUNE_UNUSED(id) = SolverType::id();
    // check for functionality
    const std::shared_ptr< const TestGridPartType > DUNE_UNUSED(testGridPart) = solverAsInterface.gridPart();
    const std::shared_ptr< const TestBoundaryInfoType > DUNE_UNUSED(testBoundaryInfo) = solverAsInterface.boundaryInfo();
    const std::shared_ptr< const TestModelType > DUNE_UNUSED(testModel) = solverAsInterface.model();
    std::shared_ptr< VectorType > vector = solverAsInterface.createVector();
    solverAsInterface.visualize(vector, "test_empty_vector", "empty vector", Dune::Stuff::Common::Logger().devnull(), "");
    solverAsInterface.init(Dune::Stuff::Common::Logger().devnull(), "");
    if (!solverAsInterface.initialized())
      DUNE_THROW(Dune::InvalidStateException, "This should not happen!");
    const ParamType mu = model->paramRange()[0];
    solverAsInterface.solve(vector, mu, "bicgstab.ilut", 1e-6, 1000, Dune::Stuff::Common::Logger().devnull(), "");
    solverAsInterface.visualize(vector, "test_solution", "solution", Dune::Stuff::Common::Logger().devnull(), "");
  }
}; // struct ParametricSolverCRTPtest

TYPED_TEST_CASE(ParametricSolverCRTPtest, ParametricSolverTypes);
TYPED_TEST(ParametricSolverCRTPtest, SolverCRTP)
{
    this->check();
}


int main(int argc, char** argv)
{
  test_init(argc, argv);
  return RUN_ALL_TESTS();
}
