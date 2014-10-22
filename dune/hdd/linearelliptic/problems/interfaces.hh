// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_PROBLEMS_INTERFACES_HH
#define DUNE_HDD_LINEARELLIPTIC_PROBLEMS_INTERFACES_HH

#include <vector>
#include <ostream>

#include <dune/stuff/common/disable_warnings.hh>
# include <dune/grid/common/gridview.hh>
# include <dune/grid/io/file/vtk.hh>
# if HAVE_ALUGRID
#   include <dune/grid/alugrid.hh>
# endif
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/string.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/stuff/functions/default.hh>

#include <dune/pymor/parameters/base.hh>
#include <dune/pymor/functions/interfaces.hh>

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Problems {


// forward, needed for with_mu()
template< class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim >
class Default;


} // namespace Problems


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

  typedef Problems::Default< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > NonparametricType;

  typedef Pymor::AffinelyDecomposableFunctionInterface
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, 1, 1 >                 DiffusionFactorType;
  typedef Pymor::AffinelyDecomposableFunctionInterface
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimDomain, dimDomain > DiffusionTensorType;
  typedef Pymor::AffinelyDecomposableFunctionInterface
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >             FunctionType;

  typedef typename FunctionType::DomainType DomainType;

  static const bool available = false;

  static std::string static_id()
  {
    return "hdd.linearelliptic.problem";
  }
  ProblemInterface(const Pymor::ParameterType tt = Pymor::ParameterType())
    : Pymor::Parametric(tt)
  {}

  ProblemInterface(const Pymor::Parametric& other)
    : Pymor::Parametric(other)
  {}

  virtual std::string type() const
  {
    return "hdd.linearelliptic.problem";
  }

  virtual const std::shared_ptr< const DiffusionFactorType >& diffusion_factor() const = 0;

  virtual const std::shared_ptr< const DiffusionTensorType >& diffusion_tensor() const = 0;

  virtual const std::shared_ptr< const FunctionType >& force() const = 0;

  virtual const std::shared_ptr< const FunctionType >& dirichlet() const = 0;

  virtual const std::shared_ptr< const FunctionType >& neumann() const = 0;

  template< class G >
  void visualize(const GridView< G >& grid_view,
                 std::string filename,
                 const bool subsampling = true,
                 const VTK::OutputType vtk_output_type = VTK::appendedraw) const
  {
    std::unique_ptr< VTKWriter< GridView< G > > >
        vtk_writer = subsampling ? DSC::make_unique< SubsamplingVTKWriter< GridView< G > > >(grid_view,
                                                                                             VTK::nonconforming)
                                 : DSC::make_unique< VTKWriter< GridView< G > > >(grid_view, VTK::nonconforming);
    add_visualizations_(grid_view, *vtk_writer);
    if (!diffusion_factor()->parametric() && !diffusion_tensor()->parametric()) {
      auto diffusion = Stuff::Functions::make_product(diffusion_factor()->affine_part(),
                                                      diffusion_tensor()->affine_part(),
                                                      "diffusion");
      auto diffusion_adapter = std::make_shared< Stuff::Functions::VisualizationAdapter
          < GridView< G >, dimDomain, dimDomain > >(*diffusion);
      vtk_writer->addVertexData(diffusion_adapter);
      vtk_writer->write(filename, vtk_output_type);
    } else
      vtk_writer->write(filename, vtk_output_type);
  } // ... visualize(...) const

  virtual void report(std::ostream& out, std::string prefix = "") const
  {
    out << prefix << "problem '" << type() << "':" << std::endl;
    out << prefix << "  diffusion_factor:" << std::endl;
    diffusion_factor()->report(out, prefix + "    ");
    out << "\n" << prefix << "  diffusion_tensor:" << std::endl;
    diffusion_tensor()->report(out, prefix + "    ");
    out << "\n" << prefix << "  force:" << std::endl;
    force()->report(out, prefix + "    ");
    out << "\n"<< prefix << "  dirichlet:" << std::endl;
    dirichlet()->report(out, prefix + "    ");
    out << "\n" << prefix << "  neumann:" << std::endl;
    neumann()->report(out, prefix + "    ");
  } // ... report(...)

  std::shared_ptr< NonparametricType > with_mu(const Pymor::Parameter mu = Pymor::Parameter()) const
  {
    if (mu.type() != this->parameter_type())
      DUNE_THROW(Pymor::Exceptions::wrong_parameter_type,
                 "mu is " << mu.type() << ", should be " << this->parameter_type() << "!");
    return std::make_shared< NonparametricType >(diffusion_factor()->with_mu(this->map_parameter(mu,
                                                                                                 "diffusion_factor")),
                                                 diffusion_tensor()->with_mu(this->map_parameter(mu,
                                                                                                 "diffusion_tensor")),
                                                 force()->with_mu(this->map_parameter(mu, "force")),
                                                 dirichlet()->with_mu(this->map_parameter(mu, "dirichlet")),
                                                 neumann()->with_mu(this->map_parameter(mu, "neumann")));
  } // ... with_mu(...)

private:
  template< class GridViewType, class VTKWriterType >
  void add_visualizations_(const GridViewType& grid_view, VTKWriterType& vtk_writer) const
  {
    add_function_visualization_(grid_view, *diffusion_factor(), vtk_writer);
    add_function_visualization_(grid_view, *diffusion_tensor(), vtk_writer);
    add_function_visualization_(grid_view, *force(), vtk_writer);
    add_function_visualization_(grid_view, *dirichlet(), vtk_writer);
    add_function_visualization_(grid_view, *neumann(), vtk_writer);
  } // ... add_visualizations_(...)

  template< class GridViewType, class F, class VTKWriterType >
  void add_function_visualization_(const GridViewType& /*grid_view*/, const F& function, VTKWriterType& vtk_writer) const
  {
    typedef Stuff::Functions::VisualizationAdapter< GridViewType, F::dimRange, F::dimRangeCols > VisualizationAdapter;
    for (DUNE_STUFF_SSIZE_T qq = 0; qq < function.num_components(); ++qq)
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

#include "default.hh"

#endif // DUNE_HDD_LINEARELLIPTIC_PROBLEMS_INTERFACES_HH
