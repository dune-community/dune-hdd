#include "config.h"

#if HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_ISTL

# include <dune/grid/alugrid.hh>

# include "swipdg.hxx"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Discretizations {


template class SWIPDG< ALUGrid< 2, 2, simplex, conforming, No_Comm >,
                       Stuff::Grid::ChooseLayer::level,
                       double,
                       1,
                       1,
                       GDT::ChooseSpaceBackend::fem,
                       Stuff::LA::ChooseBackend::istl_sparse >;
template SWIPDG< ALUGrid< 2, 2, simplex, conforming, No_Comm >,
                 Stuff::Grid::ChooseLayer::level,
                 double,
                 1,
                 1,
                 GDT::ChooseSpaceBackend::fem,
                             Stuff::LA::ChooseBackend::istl_sparse >
    ::SWIPDG< Stuff::Grid::ChooseLayer::level >(GridProviderType&,
                                                const Stuff::Common::Configuration&,
                                                const ProblemType&,
                                                const int level_or_subdomain,
                                                const std::vector< std::string >&,
                                                void*);

template class SWIPDG< ALUGrid< 2, 2, simplex, conforming, No_Comm >,
                       Stuff::Grid::ChooseLayer::leaf,
                       double,
                       1,
                       1,
                       GDT::ChooseSpaceBackend::fem,
                       Stuff::LA::ChooseBackend::istl_sparse >;
template SWIPDG< ALUGrid< 2, 2, simplex, conforming, No_Comm >,
                 Stuff::Grid::ChooseLayer::leaf,
                 double,
                 1,
                 1,
                 GDT::ChooseSpaceBackend::fem,
                 Stuff::LA::ChooseBackend::istl_sparse >
    ::SWIPDG< Stuff::Grid::ChooseLayer::leaf >(GridProviderType&,
                                               const Stuff::Common::Configuration&,
                                               const ProblemType&,
                                               const int level_or_subdomain,
                                               const std::vector< std::string >&,
                                               void*);

#if HAVE_MPI

template class SWIPDG< ALUGrid< 2, 2, simplex, conforming, MPI_Comm  >,
                       Stuff::Grid::ChooseLayer::level,
                       double,
                       1,
                       1,
                       GDT::ChooseSpaceBackend::fem,
                       Stuff::LA::ChooseBackend::istl_sparse >;
template SWIPDG< ALUGrid< 2, 2, simplex, conforming, MPI_Comm  >,
                      Stuff::Grid::ChooseLayer::level,
                      double,
                      1,
                      1,
                      GDT::ChooseSpaceBackend::fem,
                      Stuff::LA::ChooseBackend::istl_sparse >
    ::SWIPDG< Stuff::Grid::ChooseLayer::level >(GridProviderType&,
                                                const Stuff::Common::Configuration&,
                                                const ProblemType&,
                                                const int level_or_subdomain,
                                                const std::vector< std::string >&,
                                                void*);

template class SWIPDG< ALUGrid< 2, 2, simplex, conforming, MPI_Comm  >,
                       Stuff::Grid::ChooseLayer::leaf,
                       double,
                       1,
                       1,
                       GDT::ChooseSpaceBackend::fem,
                       Stuff::LA::ChooseBackend::istl_sparse >;

template SWIPDG< ALUGrid< 2, 2, simplex, conforming, MPI_Comm  >,
                 Stuff::Grid::ChooseLayer::leaf,
                 double,
                 1,
                 1,
                 GDT::ChooseSpaceBackend::fem,
                 Stuff::LA::ChooseBackend::istl_sparse >
    ::SWIPDG< Stuff::Grid::ChooseLayer::leaf >(GridProviderType&,
                                               const Stuff::Common::Configuration&,
                                               const ProblemType&,
                                               const int level_or_subdomain,
                                               const std::vector< std::string >&,
                                               void*);

#endif // HAVE_MPI
#if HAVE_DUNE_GRID_MULTISCALE

template class SWIPDG< ALUGrid< 2, 2, simplex, conforming, No_Comm >,
                       Stuff::Grid::ChooseLayer::local,
                       double,
                       1,
                       1,
                       GDT::ChooseSpaceBackend::fem,
                       Stuff::LA::ChooseBackend::istl_sparse >;

template class SWIPDG< ALUGrid< 2, 2, simplex, conforming, No_Comm >,
                       Stuff::Grid::ChooseLayer::local_oversampled,
                       double,
                       1,
                       1,
                       GDT::ChooseSpaceBackend::fem,
                       Stuff::LA::ChooseBackend::istl_sparse >;

# if HAVE_MPI

template class SWIPDG< ALUGrid< 2, 2, simplex, conforming, MPI_Comm >,
                       Stuff::Grid::ChooseLayer::local,
                       double,
                       1,
                       1,
                       GDT::ChooseSpaceBackend::fem,
                       Stuff::LA::ChooseBackend::istl_sparse >;

template class SWIPDG< ALUGrid< 2, 2, simplex, conforming, MPI_Comm >,
                       Stuff::Grid::ChooseLayer::local_oversampled,
                       double,
                       1,
                       1,
                       GDT::ChooseSpaceBackend::fem,
                       Stuff::LA::ChooseBackend::istl_sparse >;

# endif // HAVE_MPI
#endif // HAVE_DUNE_GRID_MULTISCALE


} // namespace Discretizations
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_ISTL
