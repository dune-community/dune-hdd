#ifndef DUNE_DETAILED_SOLVERS_LINEARELLIPTIC_MODEL_DEFAULT_HH
#define DUNE_DETAILED_SOLVERS_LINEARELLIPTIC_MODEL_DEFAULT_HH

#include <memory>

#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/function.hh>

#include "interface.hh"

namespace Dune {
namespace Detailed {
namespace Solvers {
namespace LinearElliptic {


template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim >
class ModelDefault;


template< class DomainFieldImp, int domainDim, class RangeFieldImp >
class ModelDefault< DomainFieldImp, domainDim, RangeFieldImp, 1 >
  : public ModelInterface< DomainFieldImp, domainDim, RangeFieldImp, 1 >
{
public:
  typedef ModelDefault< DomainFieldImp, domainDim, RangeFieldImp, 1 >    ThisType;
  typedef ModelInterface<  DomainFieldImp, domainDim, RangeFieldImp, 1 > BaseType;

  typedef typename BaseType::DomainFieldType  DomainFieldType;
  static const int                            dimDomain = BaseType::dimDomain;

  typedef typename BaseType::RangeFieldType   RangeFieldType;
  static const int                            dimRange = BaseType::dimRange;

  typedef typename BaseType::FunctionType FunctionType;

  static const std::string id()
  {
    return BaseType::id() + ".default";
  }

  ModelDefault(const std::shared_ptr< const FunctionType > _diffusion,
               const std::shared_ptr< const FunctionType > _force,
               const std::shared_ptr< const FunctionType > _dirichlet,
               const std::shared_ptr< const FunctionType > _neumann)
    : diffusion_(_diffusion)
    , force_(_force)
    , dirichlet_(_dirichlet)
    , neumann_(_neumann)
  {
    if (BaseType::parametric())
      DUNE_THROW(Dune::RangeError,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " not implemented for parametric functions!");
  }

  ModelDefault(const ThisType& _other)
    : diffusion_(_other.diffusion_)
    , force_(_other.force_)
    , dirichlet_(_other.dirichlet_)
    , neumann_(_other.neumann_)
  {}

  ThisType& operator=(const ThisType& _other)
  {
    if (this != &_other) {
      diffusion_ = _other.diffusion();
      force_ = _other.force();
      dirichlet_ = _other.dirichlet();
      neumann_ = _other.neumann();
    }
    return this;
  }

  static Dune::ParameterTree createSampleDescription(const std::string subName = "")
  {
    Dune::Stuff::Common::ExtendedParameterTree description;
    description.add(Dune::Stuff::Functions< DomainFieldType, dimDomain,
                                            RangeFieldType, dimRange >::createSampleDescription("function.expression"),
                    "diffusion");
    description["diffusion.name"] = "diffusion";
    description["diffusion.expression"] = "1.0";
    description["diffusion.order"] = "0";
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
    Dune::Stuff::Common::ExtendedParameterTree description;
    if (_description.hasSub(_subName))
      description = _description.sub(_subName);
    else
      description = _description;
    return new ThisType(createFunction("diffusion", description),
                        createFunction("force", description),
                        createFunction("dirichlet", description),
                        createFunction("neumann", description));
  } // static ThisType createFromParamTree(const Dune::ParameterTree& paramTree)

  virtual std::shared_ptr< const FunctionType > diffusion() const
  {
    return diffusion_;
  }

  virtual std::shared_ptr< const FunctionType > force() const
  {
    return force_;
  }

  virtual std::shared_ptr< const FunctionType > dirichlet() const
  {
    return dirichlet_;
  }

  virtual std::shared_ptr< const FunctionType > neumann() const
  {
    return neumann_;
  }

private:
  static std::shared_ptr< const FunctionType > createFunction(const std::string& _id,
                                                               const Dune::Stuff::Common::ExtendedParameterTree& _description)
  {
    const Dune::Stuff::Common::ExtendedParameterTree functionDescription = _description.sub(_id);
    std::string type;
    // we do this instead of using get() with 'function.expression' as a default in order to prevent get() from
    // throwing a warning
    if (functionDescription.hasKey("type"))
      type = functionDescription.get< std::string >("type");
    else
      type = "function.expression";
    return std::shared_ptr< const FunctionType >(Dune::Stuff::Functions<  DomainFieldType, dimDomain,
                                                                          RangeFieldType, dimRange >::create(type,
                                                                                                             functionDescription));
  }

  std::shared_ptr< const FunctionType > diffusion_;
  std::shared_ptr< const FunctionType > force_;
  std::shared_ptr< const FunctionType > dirichlet_;
  std::shared_ptr< const FunctionType > neumann_;
}; // class LinearDefault

} // namespace LinearElliptic
} // namespace Solvers
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_LINEARELLIPTIC_MODEL_DEFAULT_HH
