#include "config.h"

#if HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_ISTL

# include "default.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Discretizations {
namespace BlockSwipdg {


template class DefaultMultiscaleGrid< ALUGrid< 2, 2, simplex, conforming, No_Comm >,
                                      double,
                                      1,
                                      1,
                                      Stuff::LA::ChooseBackend::istl_sparse >;

# if HAVE_MPI

template class DefaultMultiscaleGrid< ALUGrid< 2, 2, simplex, conforming, MPI_Comm >,
                                      double,
                                      1,
                                      1,
                                      Stuff::LA::ChooseBackend::istl_sparse >;


# endif // HAVE_MPI

} // BlockSwipdg
} // namespace Discretizations
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_ISTL
