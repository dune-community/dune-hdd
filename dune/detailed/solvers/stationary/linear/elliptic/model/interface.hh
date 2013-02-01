#ifndef DUNE_RB_MODEL_STATIONARY_LINEAR_ELLIPTIC_INTERFACE_HH
#define DUNE_RB_MODEL_STATIONARY_LINEAR_ELLIPTIC_INTERFACE_HH

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG

#include <vector>

#include <dune/common/shared_ptr.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/stuff/common/string.hh>
#include <dune/stuff/common/color.hh>
#include <dune/stuff/function/interface.hh>

namespace Dune {
namespace Detailed {
namespace Solvers {
namespace Stationary {
namespace Linear {
namespace Elliptic {
namespace Model {


template< class DomainFieldImp, int domainDim,
          class RangeFieldImp, int rangeDim,
          class ParamFieldImp = double, int maxParamDim = 0 >
class Interface;



template< class DomainFieldImp, int domainDim,
          class RangeFieldImp,
          class ParamFieldImp, int maxParamDim >
class Interface< DomainFieldImp, domainDim, RangeFieldImp, 1, ParamFieldImp, maxParamDim >
{
public:
  typedef Interface<  DomainFieldImp, domainDim,
                      RangeFieldImp, 1,
                      ParamFieldImp, maxParamDim > ThisType;

  typedef DomainFieldImp  DomainFieldType;
  static const int        dimDomain = domainDim;

  typedef RangeFieldImp   RangeFieldType;
  static const int        dimRange = 1;

  typedef ParamFieldImp   ParamFieldType;
  static const int        maxDimParam = maxParamDim;

  typedef Dune::Stuff::Function::Interface< DomainFieldType, dimDomain,
                                            RangeFieldType, dimRange,
                                            ParamFieldImp, maxDimParam >  FunctionType;

  typedef typename FunctionType::ParamType                              ParamType;
  typedef typename FunctionType::size_type                              size_type;

  static const std::string id()
  {
    return "model.stationary.linear.elliptic";
  }

  /** \defgroup purevirtual-type ´´These methods have to be implemented and determine the type of the model (and also, which of the below methods have to be implemented).'' */
  /* @{ */
  virtual bool parametric() const = 0;
  /* @} */

  /** \defgroup purevirtual ´´These methods have to be implemented.'' */
  /* @{ */
  virtual Dune::shared_ptr< const FunctionType > diffusion() const = 0;

  virtual Dune::shared_ptr< const FunctionType > force() const = 0;

  virtual Dune::shared_ptr< const FunctionType > dirichlet() const = 0;

  virtual Dune::shared_ptr< const FunctionType > neumann() const = 0;
  /* @} */

  /** \defgroup parametric ´´These methods have to be implemented additionally, if parametric() == true.'' */
  /* @{ */
  virtual size_type paramSize() const
  {
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
    DUNE_THROW(Dune::NotImplemented,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " implement me if parametric() == true!");
  }
  /* @} */

//  virtual const Dune::shared_ptr< const NonparametricType > fix(const ParamType&) const = 0;

//  void report(std::ostream& out = std::cout, std::string prefix = "") const
//  {
//    out << prefix << "parameter explanation:" << std::endl;
//    assert(paramExplanation().size() == paramSize());
//    assert(paramRange().size() == 2);
//    assert(paramRange()[0].size() == paramSize());
//    assert(paramRange()[1].size() == paramSize());
//    for (unsigned int pp = 0; pp < paramSize(); ++pp)
//      out << prefix << "  " << paramExplanation()[pp] << ", between " << paramRange()[0](pp) << " and " << paramRange()[1](pp) << std::endl;
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
      const typename FunctionType::DomainType center = entity.geometry().center();
      // do a piecewise constant projection of the data functions
      //   * diffusion
      if (diffusion()->parametric() && diffusion()->separable())
        for (unsigned int qq = 0; qq < diffusion()->numComponents(); ++qq)
          diffusionPlots[qq].second->operator[](index) = diffusion()->components()[qq]->evaluate(center);
      else if (!diffusion()->parametric())
        diffusionPlots[0].second->operator[](index) = diffusion()->evaluate(center);
      else
        DUNE_THROW(Dune::InvalidStateException,
                   "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                   << " visualize() not yet implemented for parametric but not separable data functions!");
      //   * force
      if (force()->parametric() && force()->separable())
        for (unsigned int qq = 0; qq < force()->numComponents(); ++qq)
          forcePlots[qq].second->operator[](index) = force()->components()[qq]->evaluate(center);
      else if (!force()->parametric())
        forcePlots[0].second->operator[](index) = force()->evaluate(center);
      else
        DUNE_THROW(Dune::InvalidStateException,
                   "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                   << " visualize() not yet implemented for parametric but not separable data functions!");
      //   * dirichlet
      if (dirichlet()->parametric() && dirichlet()->separable())
        for (unsigned int qq = 0; qq < dirichlet()->numComponents(); ++qq)
          dirichletPlots[qq].second->operator[](index) = dirichlet()->components()[qq]->evaluate(center);
      else if (!dirichlet()->parametric())
        dirichletPlots[0].second->operator[](index) = dirichlet()->evaluate(center);
      else
        DUNE_THROW(Dune::InvalidStateException,
                   "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                   << " visualize() not yet implemented for parametric but not separable data functions!");
      //   * neumann
      if (neumann()->parametric() && neumann()->separable())
        for (unsigned int qq = 0; qq < neumann()->numComponents(); ++qq)
          neumannPlots[qq].second->operator[](index) = neumann()->components()[qq]->evaluate(center);
      else if (!neumann()->parametric())
        neumannPlots[0].second->operator[](index) = neumann()->evaluate(center);
      else
        DUNE_THROW(Dune::InvalidStateException,
                   "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                   << " visualize() not yet implemented for parametric but not separable data functions!");
    } // walk the grid view
    // write
    //   * diffusion
    for (unsigned int qq = 0; qq < diffusionPlots.size(); ++qq)
      vtkWriter.addCellData(*(diffusionPlots[qq].second), diffusionPlots[qq].first);
    //   * force
    for (unsigned int qq = 0; qq < forcePlots.size(); ++qq)
      vtkWriter.addCellData(*(forcePlots[qq].second), forcePlots[qq].first);
    //   * dirichlet
    for (unsigned int qq = 0; qq < dirichletPlots.size(); ++qq)
      vtkWriter.addCellData(*(dirichletPlots[qq].second), dirichletPlots[qq].first);
    //   * neumann
    for (unsigned int qq = 0; qq < neumannPlots.size(); ++qq)
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
      if (function.separable()) {
        const std::vector< std::string >& paramExplanations = function.paramExplanation();
        for (unsigned int qq = 0; qq < function.numComponents(); ++qq)
          ret.push_back(std::pair< std::string, std::vector< RangeFieldType >* >(
                          name + "_component_" + Dune::Stuff::Common::toString(qq) + ": " + paramExplanations[qq],
                          new std::vector< RangeFieldType >(gridView.indexSet().size(0), RangeFieldType(0))));
      } else
        DUNE_THROW(Dune::InvalidStateException,
                   "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                   << " visualize() not yet implemented for parametric but not separable data functions!");
    } else {
      ret.push_back(std::pair< std::string, std::vector< RangeFieldType >* >(
                      name,
                      new std::vector< RangeFieldType >(gridView.indexSet().size(0), RangeFieldType(0))));
    }
    return ret;
  } // std::vector< std::pair< std::string, std::vector< RangeFieldType >* > > preparePlot(...) const

}; // class Interface

} // namespace Model
} // namespace Elliptic
} // namespace Linear
} // namespace Stationary
} // namespace Solvers
} // namespace Detailed
} // namespace Dune

#endif // DUNE_RB_MODEL_STATIONARY_LINEAR_ELLIPTIC_INTERFACE_HH
