// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_MPI_CG_HH
#define DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_MPI_CG_HH

#include <dune/hdd/linearelliptic/discretizations/cg.hh>
#include <dune/hdd/linearelliptic/discretizations/mpi_cg.hh>
#include <dune/grid/spgrid.hh>

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Discretizations {

class MpiCG : public CG < Dune::SPGrid<double, 3>,
                          Stuff::Grid::ChooseLayer::leaf,
                          double,
                          1,
                          1,
                          Dune::GDT::ChooseSpaceBackend::pdelab,
      Dune::Stuff::LA::ChooseBackend::istl_sparse>
{
public:
  typedef Dune::SPGrid<double, 3> GridType;
  typedef CG < GridType,  Stuff::Grid::ChooseLayer::leaf, double, 1, 1, Dune::GDT::ChooseSpaceBackend::pdelab, Dune::Stuff::LA::ChooseBackend::istl_sparse> BaseType;
  typedef typename BaseType::SpaceProvider::Type SpaceType;


  MpiCG(GridProviderType& grid_provider,
     Stuff::Common::Configuration bound_inf_cfg,
     const ProblemType& prob)
    :BaseType(grid_provider, bound_inf_cfg, prob)
  {

  }
};

} // namespace Discretizations
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_MPI_CG_HH
