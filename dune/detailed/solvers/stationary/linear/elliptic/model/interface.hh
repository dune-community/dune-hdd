#ifndef DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_INTERFACE_HH
#define DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_INTERFACE_HH

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#include <dune/common/shared_ptr.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/stuff/function/interface.hh>

namespace Dune {
namespace Detailed {
namespace Solvers {
namespace Stationary {
namespace Linear {
namespace Elliptic {
namespace Model {

template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim >
class Interface
{
public:
  typedef DomainFieldImp DomainFieldType;

  static const int dimDomain = domainDim;

  typedef RangeFieldImp RangeFieldType;

  static const int dimRange = rangeDim;

  typedef Dune::Stuff::Function::Interface< DomainFieldType, dimDomain, RangeFieldType, dimRange > DiffusionType;

  typedef Dune::Stuff::Function::Interface< DomainFieldType, dimDomain, RangeFieldType, dimRange > ForceType;

  typedef Dune::Stuff::Function::Interface< DomainFieldType, dimDomain, RangeFieldType, dimRange > DirichletType;

  typedef Dune::Stuff::Function::Interface< DomainFieldType, dimDomain, RangeFieldType, dimRange > NeumannType;

  static const std::string id()
  {
    return "detailed.solvers.stationary.linear.elliptic.model";
  }

  virtual const Dune::shared_ptr< const DiffusionType > diffusion() const = 0;

  virtual int diffusionOrder() const = 0;

  virtual const Dune::shared_ptr< const ForceType > force() const = 0;

  virtual int forceOrder() const = 0;

  virtual const Dune::shared_ptr< const DirichletType > dirichlet() const = 0;

  virtual int dirichletOrder() const = 0;

  virtual const Dune::shared_ptr< const NeumannType > neumann() const = 0;

  virtual int neumannOrder() const = 0;

  /**
   *  \attention  Does a piecewise constant projection, given integration orders are ignored!
   */
  template< class GridViewType >
  void visualize(const GridViewType& gridView, const std::string filename = id()) const
  {
    // vtk writer
    assert(dimRange == 1 && "Only implemented for scalar functions!");
    Dune::VTKWriter< GridViewType > vtkwriter(gridView);
    // data
    std::vector< double > _diffusion(gridView.indexSet().size(0));
    std::vector< double > _force(gridView.indexSet().size(0));
    std::vector< double > _dirichlet(gridView.indexSet().size(0));
    std::vector< double > _neumann(gridView.indexSet().size(0));
    // walk the grid
    typedef typename GridViewType::IndexSet::IndexType IndexType;
    typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
    for (typename GridViewType::template Codim< 0 >::Iterator it = gridView.template begin< 0 >();
         it != gridView.template end< 0 >();
         ++it)
    {
      const EntityType& entity = *it;
      typename DiffusionType::DomainType x = entity.geometry().center();
      const IndexType& index = gridView.indexSet().index(entity);
      // evaluate
      _diffusion[index] = diffusion()->evaluate(x);
      _force[index] = force()->evaluate(x);
      _dirichlet[index] = dirichlet()->evaluate(x);
      _neumann[index] = neumann()->evaluate(x);
    } // walk the grid
    // write
    vtkwriter.addCellData(_diffusion, "diffusion");
    vtkwriter.addCellData(_force, "force");
    vtkwriter.addCellData(_dirichlet, "dirichlet");
    vtkwriter.addCellData(_neumann, "neumann");
    vtkwriter.write(filename, Dune::VTK::ascii);
  } // void visualize(const std::string filename) const
}; // class Interface

} // namespace Model
} // namespace Elliptic
} // namespace Linear
} // namespace Stationary
} // namespace Solvers
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_SOLVERS_STATIONARY_LINEAR_ELLIPTIC_MODEL_INTERFACE_HH
