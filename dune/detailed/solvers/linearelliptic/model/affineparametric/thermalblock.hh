#ifndef DUNE_DETAILED_SOLVERS_LINEARELLIPTIC_MODEL_AFFINEPARAMETRIC_THERMALBLOCK_HH
#define DUNE_DETAILED_SOLVERS_LINEARELLIPTIC_MODEL_AFFINEPARAMETRIC_THERMALBLOCK_HH

#include <memory>
#include <vector>

#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/string.hh>
#include <dune/stuff/function/affineparametric/checkerboard.hh>
#include <dune/stuff/function.hh>

#include "../default.hh"
#include "default.hh"

namespace Dune {
namespace Detailed {
namespace Solvers {
namespace LinearElliptic {


// forward, to allow for specialization
template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim >
class ModelAffineParametricThermalblock;


template< class DomainFieldImp, int domainDim, class RangeFieldImp >
class ModelAffineParametricThermalblock< DomainFieldImp, domainDim, RangeFieldImp, 1 >
  : public ModelAffineParametricDefault< DomainFieldImp, domainDim, RangeFieldImp, 1 >
{
public:
  typedef ModelAffineParametricThermalblock< DomainFieldImp, domainDim, RangeFieldImp, 1 > ThisType;
  typedef ModelAffineParametricDefault< DomainFieldImp, domainDim, RangeFieldImp, 1 >      BaseType;

  typedef typename BaseType::DomainFieldType  DomainFieldType;
  static const int                            dimDomain = BaseType::dimDomain;
  typedef typename BaseType::RangeFieldType   RangeFieldType;
  static const int                            dimRange = BaseType::dimRange;
  typedef typename BaseType::ParamFieldType   ParamFieldType;
  static const int                            maxParamDim = BaseType::maxParamDim;

  typedef typename Dune::Stuff::FunctionAffineParametricCheckerboard< DomainFieldType, dimDomain,
                                                                      RangeFieldType, dimRange >
                                          CheckerboardFunctionType;
  typedef typename BaseType::FunctionType FunctionType;

  static const std::string id()
  {
    return BaseType::BaseType::id() + ".affineparametric.thermalblock";
  }

  ModelAffineParametricThermalblock(const Dune::shared_ptr< const CheckerboardFunctionType > _diffusion,
                                    const Dune::shared_ptr< const FunctionType > _force,
                                    const Dune::shared_ptr< const FunctionType > _dirichlet,
                                    const Dune::shared_ptr< const FunctionType > _neumann)
    : BaseType(_diffusion, _force, _dirichlet, _neumann)
  {}

  ModelAffineParametricThermalblock(const ThisType& other)
    : BaseType(other.diffusion(),
               other.force(),
               other.dirichlet(),
               other.neumann())
  {}

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      BaseType::operator=(other);
    }
    return *this;
  }

  static Dune::ParameterTree createSampleDescription(const std::string subName = "")
  {
    Dune::Stuff::Common::ExtendedParameterTree description;
    description.add(Dune::Stuff::Functions< DomainFieldType, dimDomain,
                                            RangeFieldType, dimRange >::createSampleDescription("function.affineparametric.checkerboard"),
                    "diffusion");
    description["diffusion.name"] = "diffusion";
    description.add(Dune::Stuff::Functions< DomainFieldType, dimDomain,
                                            RangeFieldType, dimRange >::createSampleDescription("function.expression"),
                    "force");
    description["force.name"] = "force";
    description["force.expression"] = "1.0";
    description["force.order"] = "0";
    description.add(Dune::Stuff::Functions< DomainFieldType, dimDomain,
                                            RangeFieldType, dimRange >::createSampleDescription("function.expression"),
                    "neumann");
    description["neumann.name"] = "neumann";
    description["neumann.expression"] = "0.1";
    description["neumann.order"] = "0";
    description.add(Dune::Stuff::Functions< DomainFieldType, dimDomain,
                                            RangeFieldType, dimRange >::createSampleDescription("function.expression"),
                    "dirichlet");
    description["dirichlet.name"] = "dirichlet";
    description["dirichlet.expression"] = "0.1*x[0]";
    description["dirichlet.order"] = "1";
    if (subName.empty())
      return description;
    else {
      Dune::Stuff::Common::ExtendedParameterTree extendedDescription;
      extendedDescription.add(description, subName);
      return extendedDescription;
    }
  }

  static ThisType* create(const Dune::ParameterTree& _description, const std::string _subName = id())
  {
    // get correct description
    typedef Stuff::Common::ExtendedParameterTree DescriptionType;
    DescriptionType description;
    if (_description.hasSub(_subName))
      description = _description.sub(_subName);
    else
      description = _description;
    // create the checkerboard diffusion
    DescriptionType diffusionDescription = description.sub("diffusion");
    diffusionDescription["diffusion.name"] = "checkerboard diffusion";
    const shared_ptr< const CheckerboardFunctionType >
        diffusion(CheckerboardFunctionType::create(diffusionDescription));
    // create the rest of the functions
    //   * therefore create a fake diffusion subdescription,
    description.sub("diffusion") = CheckerboardFunctionType::createSampleDescription();
    description["diffusion.type"] = "function.affineparametric.checkerboard";
    description["diffusion.name"] = "checkerboard diffusion";
    //   * create a default model
    const std::shared_ptr< const BaseType > baseModel(BaseType::create(description));
    // create and return
    return new ThisType(diffusion,
                        baseModel->force(),
                        baseModel->dirichlet(),
                        baseModel->neumann());
  } // static ThisType createFromParamTree(const Dune::ParameterTree& paramTree, const std::string subName = id())
}; // class ModelAffineParametricThermalblock


} // namespace LinearElliptic
} // namespace Solvers
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_LINEARELLIPTIC_MODEL_AFFINEPARAMETRIC_THERMALBLOCK_HH
