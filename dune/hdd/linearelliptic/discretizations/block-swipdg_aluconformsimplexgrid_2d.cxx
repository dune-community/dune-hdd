#ifdef ENABLE_MPI
# undef ENABLE_MPI
#endif
#define ENABLE_MPI 0

#ifdef ENABLE_PARMETIS
# undef ENABLE_PARMETIS
#endif
#define ENABLE_PARMETIS 0

#include "config.h"

#if HAVE_ALUGRID

# include "block-swipdg.hxx"
# include "cg.hxx"
# include "swipdg.hxx"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Discretizations {


template class BlockSWIPDG< ALUGrid< 2, 2, simplex, conforming >,
                            double,
                            1,
                            1,
                            Stuff::LA::ChooseBackend::istl_sparse >;


} // namespace Discretizations
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // HAVE_ALUGRID
