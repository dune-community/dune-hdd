#include <dune/stuff/test/test_common.hh>

#include <memory>

#include <dune/grid/sgrid.hh>

#include <dune/grid/part/leaf.hh>

#include <dune/stuff/grid/provider.hh>
#include <dune/stuff/common/logging.hh>

#include <dune/hdd/linearelliptic/model.hh>
#include <dune/hdd/linearelliptic/model/default.hh>
#include <dune/hdd/linearelliptic/model/thermalblock.hh>
#include <dune/hdd/linearelliptic/model/affineparametric/default.hh>
#include <dune/hdd/linearelliptic/model/affineparametric/thermalblock.hh>
#include <dune/hdd/linearelliptic/model/affineparametric/twophase.hh>

namespace Stuff = Dune::Stuff;
namespace Elliptic = Dune::HDD::LinearElliptic;

static const int            dimDomain = 2;
typedef Dune::SGrid< 2, 2 > GridType;

typedef typename GridType::ctype  DomainFieldType;
typedef double                    RangeFieldType;
static const int                  dimRange = 1;

typedef Dune::grid::Part::Leaf::Const< GridType > GridPartType;

typedef Elliptic::ModelInterface< DomainFieldType, dimDomain, RangeFieldType, dimRange > ModelInterfaceType;


TEST(ModelTest, create)
{
  typedef Stuff::GridProviderInterface< GridType > GridProviderType;
  typedef Stuff::GridProviders< GridType > GridProviders;
  const std::shared_ptr< const GridProviderType >
      gridProvider(GridProviders::create(GridProviders::available()[0],
                                         GridProviders::createSampleDescription(GridProviders::available()[0])));
  const GridPartType gridPart(*gridProvider->grid());
  typedef Elliptic::Models< DomainFieldType, dimDomain, RangeFieldType, dimRange > EllipticModels;
  const std::vector< std::string > modelTypes = EllipticModels::available();
  for (auto modelType : modelTypes) {
    const Dune::ParameterTree modelDescription = EllipticModels::createSampleDescription(modelType);
    std::shared_ptr< ModelInterfaceType >
        DUNE_UNUSED(model) = std::shared_ptr< ModelInterfaceType >(EllipticModels::create(modelType, modelDescription));
    model->visualize(gridPart.gridView(), "test_model");
  }
}


typedef testing::Types< Elliptic::ModelDefault< DomainFieldType, dimDomain, RangeFieldType, dimRange >
                      , Elliptic::ModelThermalblock< DomainFieldType, dimDomain, RangeFieldType, dimRange >
                      > NonparametricModelTypes;

template< class T >
struct NonparametricModelTest
  : public ::testing::Test
{
  typedef T ModelType;
  void check() const
  {
    const Dune::ParameterTree modelDescription = ModelType::createSampleDescription();
    std::shared_ptr< ModelType > model = std::shared_ptr< ModelType >(ModelType::create(modelDescription));
    if (model->parametric())
      DUNE_THROW(Dune::InvalidStateException, "This should not happen!");
    if (model->affineparametric())
      DUNE_THROW(Dune::InvalidStateException, "This should not happen!");
    if (model->diffusion()->parametric())
      DUNE_THROW(Dune::InvalidStateException, "This should not happen!");
    if (model->diffusion()->affineparametric())
      DUNE_THROW(Dune::InvalidStateException, "This should not happen!");
    if (model->force()->parametric())
      DUNE_THROW(Dune::InvalidStateException, "This should not happen!");
    if (model->force()->affineparametric())
      DUNE_THROW(Dune::InvalidStateException, "This should not happen!");
    if (model->dirichlet()->parametric())
      DUNE_THROW(Dune::InvalidStateException, "This should not happen!");
    if (model->dirichlet()->affineparametric())
      DUNE_THROW(Dune::InvalidStateException, "This should not happen!");
    if (model->neumann()->parametric())
      DUNE_THROW(Dune::InvalidStateException, "This should not happen!");
    if (model->neumann()->affineparametric())
      DUNE_THROW(Dune::InvalidStateException, "This should not happen!");
    if (model->paramSize() != 0)
      DUNE_THROW(Dune::InvalidStateException, "This should not happen!");
  }
};

TYPED_TEST_CASE(NonparametricModelTest, NonparametricModelTypes);
TYPED_TEST(NonparametricModelTest, NonparametricModels)
{
    this->check();
}


typedef testing::Types< Elliptic::ModelAffineParametricDefault< DomainFieldType, dimDomain, RangeFieldType, dimRange >
                      , Elliptic::ModelAffineParametricThermalblock< DomainFieldType, dimDomain, RangeFieldType, dimRange >
                      , Elliptic::ModelAffineParametricTwoPhase< DomainFieldType, dimDomain, RangeFieldType, dimRange >
                      > AffineparametricModelTypes;


template< class T >
struct AffineparametricModelTest
  : public ::testing::Test
{
  typedef T ModelType;
  typedef typename ModelType::ParamType ParamType;

  void check() const
  {
    const Dune::ParameterTree modelDescription = ModelType::createSampleDescription();
    std::shared_ptr< ModelType > model = std::shared_ptr< ModelType >(ModelType::create(modelDescription));
    if (!model->parametric())
      DUNE_THROW(Dune::InvalidStateException, "This should not happen!");
    if (!model->affineparametric())
      DUNE_THROW(Dune::InvalidStateException, "This should not happen!");
    if (!model->diffusion()->parametric()
        && !model->force()->parametric()
        && !model->dirichlet()->parametric()
        && !model->neumann()->parametric())
      DUNE_THROW(Dune::InvalidStateException, "This should not happen!");
    if (model->diffusion()->parametric() && !model->diffusion()->affineparametric())
      DUNE_THROW(Dune::InvalidStateException, "This should not happen!");
    if (model->force()->parametric() && !model->force()->affineparametric())
      DUNE_THROW(Dune::InvalidStateException, "This should not happen!");
    if (model->dirichlet()->parametric() && !model->dirichlet()->affineparametric())
      DUNE_THROW(Dune::InvalidStateException, "This should not happen!");
    if (model->neumann()->parametric() && !model->neumann()->affineparametric())
      DUNE_THROW(Dune::InvalidStateException, "This should not happen!");
    if (model->paramSize() == 0)
      DUNE_THROW(Dune::InvalidStateException, "This should not happen!");
    const std::vector< ParamType >& paramRange = model->paramRange();
    if (paramRange.size() != 2)
      DUNE_THROW(Dune::InvalidStateException, "This should not happen!");
    if (paramRange[0].size() != model->paramSize())
      DUNE_THROW(Dune::InvalidStateException, "This should not happen!");
    if (paramRange[1].size() != model->paramSize())
      DUNE_THROW(Dune::InvalidStateException, "This should not happen!");
    for (size_t pp = 0; pp < model->paramSize(); ++pp)
      if (!(paramRange[0][pp] < paramRange[1][pp]))
        DUNE_THROW(Dune::InvalidStateException, "This should not happen!");
  }
};

TYPED_TEST_CASE(AffineparametricModelTest, AffineparametricModelTypes);
TYPED_TEST(AffineparametricModelTest, AffineparametricModels)
{
    this->check();
}


int main(int argc, char** argv)
{
  test_init(argc, argv);
  return RUN_ALL_TESTS();
}
