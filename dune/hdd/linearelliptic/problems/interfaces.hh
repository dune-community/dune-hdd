// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_PROBLEMS_INTERFACES_HH
#define DUNE_HDD_LINEARELLIPTIC_PROBLEMS_INTERFACES_HH

#include <vector>
#include <ostream>

#include <dune/grid/common/gridview.hh>
#include <dune/grid/io/file/vtk.hh>

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/string.hh>
#include <dune/stuff/functions/default.hh>

#include <dune/pymor/parameters/base.hh>
#include <dune/pymor/functions/interfaces.hh>

namespace Dune {
namespace HDD {
namespace LinearElliptic {


template< class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim >
class ProblemInterface
  : public Pymor::Parametric
{
  typedef ProblemInterface< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim > ThisType;
public:
  typedef EntityImp EntityType;
  typedef DomainFieldImp    DomainFieldType;
  static const unsigned int dimDomain = domainDim;
  typedef RangeFieldImp     RangeFieldType;
  static const unsigned int dimRange = rangeDim;

  ProblemInterface(const Pymor::ParameterType tt = Pymor::ParameterType())
    : Pymor::Parametric(tt)
  {}

  ProblemInterface(const Pymor::Parametric& other)
    : Pymor::Parametric(other)
  {}

  typedef Pymor::AffinelyDecomposableFunctionInterface< EntityType, DomainFieldType, dimDomain
                                                      , RangeFieldType, 1, 1 >
    DiffusionFactorType;
  typedef Pymor::AffinelyDecomposableFunctionInterface< EntityType, DomainFieldType, dimDomain
                                                      , RangeFieldType, dimDomain, dimDomain >
    DiffusionTensorType;
  typedef Pymor::AffinelyDecomposableFunctionInterface< EntityType, DomainFieldType, dimDomain
                                                      , RangeFieldType, dimRange >
    FunctionType;

  static std::string static_id()
  {
    return "hdd.linearelliptic.problem";
  }

  virtual std::string type() const
  {
    return "hdd.linearelliptic.problem";
  }

  virtual const DiffusionFactorType& diffusion_factor() const = 0;

  virtual const DiffusionTensorType& diffusion_tensor() const = 0;

  virtual const FunctionType& force() const = 0;

  virtual const FunctionType& dirichlet() const = 0;

  virtual const FunctionType& neumann() const = 0;

  template< class G >
  void visualize(const GridView< G >& grid_view,
                 std::string filename,
                 const bool subsampling = true,
                 const VTK::OutputType vtk_output_type = VTK::appendedraw) const
  {
    if (subsampling) {
      SubsamplingVTKWriter< GridView< G > > vtk_writer(grid_view, VTK::nonconforming);
      add_visualizations_(grid_view, vtk_writer);
      vtk_writer.write(filename, vtk_output_type);
    } else {
      VTKWriter< GridView< G > > vtk_writer(grid_view, VTK::nonconforming);
      add_visualizations_(grid_view, vtk_writer);
      vtk_writer.write(filename, vtk_output_type);
    }
  } // ... visualize(...) const

  virtual void report(std::ostream& out, std::string prefix = "") const
  {
    out << prefix << "problem '" << type() << "':" << std::endl;
    out << prefix << "  diffusion_factor:" << std::endl;
    diffusion_factor().report(out, prefix + "    ");
    out << "\n" << prefix << "  diffusion_tensor:" << std::endl;
    diffusion_tensor().report(out, prefix + "    ");
    out << "\n" << prefix << "  force:" << std::endl;
    force().report(out, prefix + "    ");
    out << "\n"<< prefix << "  dirichlet:" << std::endl;
    dirichlet().report(out, prefix + "    ");
    out << "\n" << prefix << "  neumann:" << std::endl;
    neumann().report(out, prefix + "    ");
  } // ... report(...)

private:
  template< class GridViewType, class VTKWriterType >
  void add_visualizations_(const GridViewType& grid_view, VTKWriterType& vtk_writer) const
  {
    add_function_visualization_(grid_view, diffusion_factor(), vtk_writer);
//    add_function_visualization_(grid_view, diffusion_tensor(), vtk_writer);
    add_function_visualization_(grid_view, force(), vtk_writer);
    add_function_visualization_(grid_view, dirichlet(), vtk_writer);
    add_function_visualization_(grid_view, neumann(), vtk_writer);
  } // ... add_visualizations_(...)

  template< class GridViewType, class F, class VTKWriterType >
  void add_function_visualization_(const GridViewType& /*grid_view*/, const F& function, VTKWriterType& vtk_writer) const
  {
    typedef Stuff::Functions::VisualizationAdapter< GridViewType, F::dimRange > VisualizationAdapter;
    for (size_t qq = 0; qq < function.num_components(); ++qq)
      vtk_writer.addVertexData(std::make_shared< VisualizationAdapter >(*(function.component(qq))));
    if (function.has_affine_part())
      vtk_writer.addVertexData(std::make_shared< VisualizationAdapter >(*(function.affine_part())));
  } // ... add_function_visualization_(...)

private:
  template< class T >
  friend std::ostream& operator<<(std::ostream& /*out*/, const ThisType& /*problem*/);
}; // ProblemInterface


template< class E, class D, int d, class R, int r >
std::ostream& operator<<(std::ostream& out, const ProblemInterface< E, D, d, R, r >& problem)
{
  problem.report(out);
  return out;
} // ... operator<<(...)


} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_PROBLEMS_INTERFACES_HH
