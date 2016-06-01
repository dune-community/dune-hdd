#include "config.h"

#if HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_EIGEN

# include "block-swipdg.hxx"
# include "cg.hxx"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Discretizations {


template class BlockSWIPDG< ALUGrid< 2, 2, simplex, conforming, No_Comm >,
                            double,
                            1,
                            1,
                            Stuff::LA::ChooseBackend::eigen_sparse >;

# if HAVE_MPI

template class BlockSWIPDG< ALUGrid< 2, 2, simplex, conforming, MPI_Comm >,
                            double,
                            1,
                            1,
                            Stuff::LA::ChooseBackend::eigen_sparse >;


# endif // HAVE_MPI

} // namespace Discretizations
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_EIGEN
