#ifndef DUNE_DETAILED_SOLVERS_LINEARELLIPTIC_MODEL_INTERFACE_HH
#define DUNE_DETAILED_SOLVERS_LINEARELLIPTIC_MODEL_INTERFACE_HH

// we still need this for the vtk writer
#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include "config.h"
#endif

#include <vector>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/stuff/common/string.hh>
#include <dune/stuff/common/color.hh>
#include <dune/stuff/common/parameter.hh>
#include <dune/stuff/function/interface.hh>
#include <dune/stuff/function/fixed.hh>

namespace Dune {
namespace DetailedSolvers {
namespace LinearElliptic {


// forward of the nonparametric default
template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, bool scalarDiffusion >
class ModelDefault;


// forward to allow for specialization
template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, bool scalarDiffusion = true >
class ModelInterface
{
public:
  ModelInterface() = delete;
};


template< class DomainFieldImp, int domainDim, class RangeFieldImp >
class ModelInterface< DomainFieldImp, domainDim, RangeFieldImp, 1, true >
{
public:
  typedef ModelInterface< DomainFieldImp, domainDim, RangeFieldImp, 1, true > ThisType;

  typedef DomainFieldImp  DomainFieldType;
  static const int        dimDomain = domainDim;

  typedef RangeFieldImp   RangeFieldType;
  static const int        dimRange = 1;

  typedef typename Stuff::Common::Parameter::FieldType  ParamFieldType;
  static const int                                      maxParamDim = Stuff::Common::Parameter::maxDim;
  typedef typename Stuff::Common::Parameter::Type       ParamType;

  typedef Dune::Stuff::GenericStationaryFunctionInterface< DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 > DiffusionType;
  typedef Dune::Stuff::GenericStationaryFunctionInterface< DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 > ForceType;
  typedef Dune::Stuff::GenericStationaryFunctionInterface< DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 > DirichletType;
  typedef Dune::Stuff::GenericStationaryFunctionInterface< DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 > NeumannType;

  static const std::string id()
  {
    return "model.linearelliptic";
  }

  /** \defgroup type ´´These methods determine the type of the model (and also, which of the below methods have to be implemented).'' */
  /* @{ */
  virtual bool parametric() const
  {
    if (diffusion()->parametric() || force()->parametric() || dirichlet()->parametric() || neumann()->parametric())
      return true;
    else
      return false;
  }

  virtual bool affineparametric() const
  {
    if (parametric()) {
      if (diffusion()->parametric() && !diffusion()->affineparametric())
        return false;
      if (force()->parametric() && !force()->affineparametric())
        return false;
      if (dirichlet()->parametric() && !dirichlet()->affineparametric())
        return false;
      if (neumann()->parametric() && !neumann()->affineparametric())
        return false;
      return true;
    } else
      return false;
  } // ... affineparametric()
  /* @} */

  /** \defgroup purevirtual ´´These methods have to be implemented.'' */
  /* @{ */
  virtual std::shared_ptr< const DiffusionType > diffusion() const = 0;

  virtual std::shared_ptr< const ForceType > force() const = 0;

  virtual std::shared_ptr< const DirichletType > dirichlet() const = 0;

  virtual std::shared_ptr< const NeumannType > neumann() const = 0;
  /* @} */

  /** \defgroup parametric ´´These methods have to be implemented additionally, if parametric() == true.'' */
  /* @{ */
  virtual size_t paramSize() const
  {
    if (!parametric())
      return 0;
    else
      DUNE_THROW(Dune::NotImplemented,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " implement me if parametric() == true!");
  }

  virtual const std::vector< ParamType >& paramRange() const
  {
    DUNE_THROW(Dune::NotImplemented,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " implement me if parametric() == true!");
  }

  virtual const std::vector< std::string >& paramExplanation() const
  {
    DUNE_THROW(Dune::NotImplemented,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " implement me if parametric() == true!");
  }

  /**
   * \brief     Maps a global to a local parameter.
   * \param[in] ParamType _globalMu
   *            global parameter
   * \param[in] std::string _id
   *            To identify the function, wrt which the localization is carried out. Has to be one of diffusion,
   *            force, dirichlet or neumann
   * \return    ParamType
   *            local parameter
   */
  virtual ParamType mapParam(const ParamType& /*_globalMu*/, const std::string /*_id*/) const
  {
    if (!parametric())
      return ParamType();
    else
      DUNE_THROW(Dune::NotImplemented,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " implement me if parametric() == true!");
  }
  /* @} */

//  virtual std::shared_ptr< const ThisType > fix(const ParamType& mu) const
//  {
//    if (!parametric() && mu.size() != 0)
//      DUNE_THROW(Dune::InvalidStateException,
//                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
//                 << " do not call fix(mu) with a nonempty mu for a nonparametric model (check parametric() == true beforehand)!");
//    typedef Dune::Stuff::FunctionFixed< DomainFieldType, dimDomain, RangeFieldType, dimRange > FixedFunctionType;
//    typedef ModelDefault< DomainFieldType, dimDomain, RangeFieldType, dimRange > DefaultModelType;
//    return std::make_shared< DefaultModelType >(diffusion()->parametric() ? std::make_shared< const FixedFunctionType >(diffusion(), mapParam(mu, "diffusion")) : diffusion(),
//                                                force()->parametric() ? std::make_shared< const FixedFunctionType >(force(), mapParam(mu, "force")) : force(),
//                                                dirichlet()->parametric() ? std::make_shared< const FixedFunctionType >(dirichlet(), mapParam(mu, "dirichlet")) : dirichlet(),
//                                                neumann()->parametric() ? std::make_shared< const FixedFunctionType >(neumann(), mapParam(mu, "neumann")) : neumann());
//  }

  template< class GridViewType >
  void visualize(const GridViewType& gridView, std::string filename) const
  {
    // prepare data
    Dune::VTKWriter< GridViewType > vtkWriter(gridView);
    std::vector< std::pair< std::string, std::vector< RangeFieldType >* > > diffusionPlots
        = preparePlot(gridView, *(diffusion()), "diffusion");
    std::vector< std::pair< std::string, std::vector< RangeFieldType >* > > forcePlots
        = preparePlot(gridView, *(force()), "force");
    std::vector< std::pair< std::string, std::vector< RangeFieldType >* > > dirichletPlots
        = preparePlot(gridView, *(dirichlet()), "dirichlet");
    std::vector< std::pair< std::string, std::vector< RangeFieldType >* > > neumannPlots
        = preparePlot(gridView, *(neumann()), "neumann");
    // walk the grid view
    for (typename GridViewType::template Codim< 0 >::Iterator it = gridView.template begin< 0 >();
         it != gridView.template end< 0 >();
         ++it) {
      const typename GridViewType::template Codim< 0 >::Entity& entity = *it;
      const typename GridViewType::IndexSet::IndexType index = gridView.indexSet().index(entity);
      const auto center = entity.geometry().center();
      // do a piecewise constant projection of the data functions
      //   * diffusion
      if (diffusion()->parametric() && diffusion()->affineparametric())
        for (size_t qq = 0; qq < diffusion()->components().size(); ++qq)
          diffusionPlots[qq].second->operator[](index) = diffusion()->components()[qq]->evaluate(center);
      else if (!diffusion()->parametric())
        diffusionPlots[0].second->operator[](index) = diffusion()->evaluate(center);
      else
        DUNE_THROW(Dune::InvalidStateException,
                   "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                   << " visualize() not yet implemented for parametric but not affineparametric data functions!");
      //   * force
      if (force()->parametric() && force()->affineparametric())
        for (size_t qq = 0; qq < force()->components().size(); ++qq)
          forcePlots[qq].second->operator[](index) = force()->components()[qq]->evaluate(center);
      else if (!force()->parametric())
        forcePlots[0].second->operator[](index) = force()->evaluate(center);
      else
        DUNE_THROW(Dune::InvalidStateException,
                   "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                   << " visualize() not yet implemented for parametric but not affineparametric data functions!");
      //   * dirichlet
      if (dirichlet()->parametric() && dirichlet()->affineparametric())
        for (size_t qq = 0; qq < dirichlet()->components().size(); ++qq)
          dirichletPlots[qq].second->operator[](index) = dirichlet()->components()[qq]->evaluate(center);
      else if (!dirichlet()->parametric())
        dirichletPlots[0].second->operator[](index) = dirichlet()->evaluate(center);
      else
        DUNE_THROW(Dune::InvalidStateException,
                   "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                   << " visualize() not yet implemented for parametric but not affineparametric data functions!");
      //   * neumann
      if (neumann()->parametric() && neumann()->affineparametric())
        for (size_t qq = 0; qq < neumann()->components().size(); ++qq)
          neumannPlots[qq].second->operator[](index) = neumann()->components()[qq]->evaluate(center);
      else if (!neumann()->parametric())
        neumannPlots[0].second->operator[](index) = neumann()->evaluate(center);
      else
        DUNE_THROW(Dune::InvalidStateException,
                   "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                   << " visualize() not yet implemented for parametric but not affineparametric data functions!");
    } // walk the grid view
    // write
    //   * diffusion
    for (size_t qq = 0; qq < diffusionPlots.size(); ++qq)
      vtkWriter.addCellData(*(diffusionPlots[qq].second), diffusionPlots[qq].first);
    //   * force
    for (size_t qq = 0; qq < forcePlots.size(); ++qq)
      vtkWriter.addCellData(*(forcePlots[qq].second), forcePlots[qq].first);
    //   * dirichlet
    for (size_t qq = 0; qq < dirichletPlots.size(); ++qq)
      vtkWriter.addCellData(*(dirichletPlots[qq].second), dirichletPlots[qq].first);
    //   * neumann
    for (size_t qq = 0; qq < neumannPlots.size(); ++qq)
      vtkWriter.addCellData(*(neumannPlots[qq].second), neumannPlots[qq].first);
    vtkWriter.write(filename, Dune::VTK::ascii);
  } // void visualize(const GridViewType& gridView, std::string filename) const

private:
  template< class GridViewType, class FunctionType >
  std::vector< std::pair< std::string, std::vector< RangeFieldType >* > > preparePlot(const GridViewType& gridView,
                                                                                      const FunctionType& function,
                                                                                      const std::string name) const
  {
    std::vector< std::pair< std::string, std::vector< RangeFieldType >* > > ret;
    if (function.parametric()) {
      if (function.affineparametric()) {
        const std::vector< std::string >& paramExplanations = function.paramExplanation();
        for (size_t qq = 0; qq < function.coefficients().size(); ++qq)
          ret.push_back(std::pair< std::string, std::vector< RangeFieldType >* >(
                          name + "_component_" + Dune::Stuff::Common::toString(qq) + ": " + paramExplanations[qq],
                          new std::vector< RangeFieldType >(gridView.indexSet().size(0), RangeFieldType(0))));
        if (function.coefficients().size() > function.coefficients().size())
          ret.push_back(std::pair< std::string, std::vector< RangeFieldType >* >(
                          name + "_component_" + Dune::Stuff::Common::toString(function.coefficients().size()),
                          new std::vector< RangeFieldType >(gridView.indexSet().size(0), RangeFieldType(0))));
      } else
        DUNE_THROW(Dune::InvalidStateException,
                   "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                   << " visualize() not yet implemented for parametric but not affineparametric data functions!");
    } else {
      ret.push_back(std::pair< std::string, std::vector< RangeFieldType >* >(
                      name,
                      new std::vector< RangeFieldType >(gridView.indexSet().size(0), RangeFieldType(0))));
    }
    return ret;
  } // std::vector< std::pair< std::string, std::vector< RangeFieldType >* > > preparePlot(...) const
}; // class ModelInterface


} // namespace LinearElliptic
} // namespace DetailedSolvers
} // namespace Dune

#include "default.hh"

#endif // DUNE_DETAILED_SOLVERS_LINEARELLIPTIC_MODEL_INTERFACE_HH
