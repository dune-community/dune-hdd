// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_PROBLEM_INTERFACES_HH
#define DUNE_HDD_LINEARELLIPTIC_PROBLEM_INTERFACES_HH

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
#include <dune/stuff/common/parameter/tree.hh>

#include <dune/pymor/parameters/base.hh>
#include <dune/pymor/functions/interfaces.hh>

namespace Dune {
namespace HDD {
namespace LinearElliptic {


template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, bool scalarDiffusion = true >
class ProblemInterface;


template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim >
class ProblemInterface< DomainFieldImp, domainDim, RangeFieldImp, rangeDim, true >
  : public Pymor::Parametric
{
  typedef Pymor::Parametric BaseType;
public:

  typedef DomainFieldImp  DomainFieldType;
  static const int        dimDomain = domainDim;

  typedef RangeFieldImp   RangeFieldType;
  static const int        dimRange = rangeDim;

  ProblemInterface(const Pymor::ParameterType tt = Pymor::ParameterType())
    : BaseType(tt)
  {}

  ProblemInterface(const Pymor::Parametric& other)
    : BaseType(other)
  {}

  typedef Pymor::ParametricFunctionInterface< DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 > DiffusionType;
  typedef Pymor::ParametricFunctionInterface< DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 > ForceType;
  typedef Pymor::ParametricFunctionInterface< DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 > DirichletType;
  typedef Pymor::ParametricFunctionInterface< DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 > NeumannType;

  static const std::string static_id()
  {
    return "dune.hdd.linearelliptic.problem";
  }

  virtual const std::string id()
  {
    return "dune.hdd.linearelliptic.problem";
  }

  virtual std::shared_ptr< const DiffusionType > diffusion() const = 0;

  virtual std::shared_ptr< const ForceType > force() const = 0;

  virtual std::shared_ptr< const DirichletType > dirichlet() const = 0;

  virtual std::shared_ptr< const NeumannType > neumann() const = 0;

  template< class GridViewType >
  void visualize(const GridViewType& gridView, std::string filename) const
  {
    // prepare data
    Dune::VTKWriter< GridViewType > vtkWriter(gridView);
    auto diffusionPlots = preparePlot(gridView, *(diffusion()), "diffusion");
    auto forcePlots = preparePlot(gridView, *(force()), "force");
    auto dirichletPlots = preparePlot(gridView, *(dirichlet()), "dirichlet");
    auto neumannPlots = preparePlot(gridView, *(neumann()), "neumann");
    // walk the grid view
    for (typename GridViewType::template Codim< 0 >::Iterator it = gridView.template begin< 0 >();
         it != gridView.template end< 0 >();
         ++it) {
      const auto& entity = *it;
      const auto index = gridView.indexSet().index(entity);
      const auto center = entity.geometry().center();
      // do a piecewise constant projection of the data functions
      plot_local(*(diffusion()), diffusionPlots, index, center);
      plot_local(*(force()), forcePlots, index, center);
      plot_local(*(dirichlet()), dirichletPlots, index, center);
      plot_local(*(neumann()), neumannPlots, index, center);
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
    // clean up
    for (auto element : diffusionPlots)
      delete element.second;
    for (auto element : forcePlots)
      delete element.second;
    for (auto element : dirichletPlots)
      delete element.second;
    for (auto element : neumannPlots)
      delete element.second;
  } // ... visualize(...) const

private:
  template< class GridViewType, class FunctionType >
  std::vector< std::pair< std::string, std::vector< RangeFieldType >* > > preparePlot(const GridViewType& gridView,
                                                                                      const FunctionType& function,
                                                                                      const std::string name) const
  {
    std::vector< std::pair< std::string, std::vector< RangeFieldType >* > > ret;
    if (function.parametric()) {
      if (function.affinely_decomposable()) {
        for (size_t qq = 0; qq < function.num_components(); ++qq)
          ret.push_back(std::pair< std::string, std::vector< RangeFieldType >* >(
                          name + ", component " + Dune::Stuff::Common::toString(qq),
                          new std::vector< RangeFieldType >(gridView.indexSet().size(0), RangeFieldType(0))));
        if (function.has_affine_part())
          ret.push_back(std::pair< std::string, std::vector< RangeFieldType >* >(
                          name + ", affine part",
                          new std::vector< RangeFieldType >(gridView.indexSet().size(0), RangeFieldType(0))));
      } else
        DUNE_PYMOR_THROW(Pymor::Exception::requirements_not_met,
                         "not implemented for parametric functions which are not affinely decomposable!");
    } else {
      ret.push_back(std::pair< std::string, std::vector< RangeFieldType >* >(
                      name,
                      new std::vector< RangeFieldType >(gridView.indexSet().size(0), RangeFieldType(0))));
    }
    return ret;
  }

  template< class FunctionType, class IndexType, class DomainType >
  void plot_local(const FunctionType& function,
                  std::vector< std::pair< std::string, std::vector< RangeFieldType >* > >& plots,
                  const IndexType& index,
                  const DomainType& center) const
  {
    if (function.parametric() && function.affinely_decomposable()) {
      for (size_t qq = 0; qq < function.num_components(); ++qq)
        plots[qq].second->operator[](index) = function.component(qq)->evaluate(center);
      if (function.has_affine_part())
        plots[plots.size() - 1].second->operator[](index) = function.affine_part()->evaluate(center);
    } else if (!function.parametric())
      plots[0].second->operator[](index) = function.evaluate(center);
    else
      DUNE_PYMOR_THROW(Pymor::Exception::requirements_not_met,
                       "not implemented for parametric functions which are not affinely decomposable!");
  }
}; // ProblemInterface


} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

//#include "default.hh"

#endif // DUNE_HDD_LINEARELLIPTIC_PROBLEM_INTERFACES_HH
