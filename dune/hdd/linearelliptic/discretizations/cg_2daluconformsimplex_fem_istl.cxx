#include "config.h"

#if HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_ISTL

# include "cg.hxx"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Discretizations {


template class CG< ALUGrid< 2, 2, simplex, conforming, No_Comm >,
                   Stuff::Grid::ChooseLayer::leaf,
                   double,
                   1,
                   1,
                   GDT::ChooseSpaceBackend::fem,
                   Stuff::LA::ChooseBackend::istl_sparse >;

# if HAVE_MPI

template class CG< ALUGrid< 2, 2, simplex, conforming, MPI_Comm >,
                   Stuff::Grid::ChooseLayer::leaf,
                   double,
                   1,
                   1,
                   GDT::ChooseSpaceBackend::fem,
                   Stuff::LA::ChooseBackend::istl_sparse >;

# endif // HAVE_MPI


} // namespace Discretizations
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_ISTL
