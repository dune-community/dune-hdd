#include "config.h"

#if HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_EIGEN

# include "generic_multiscale.hxx"


template class GenericLinearellipticMultiscaleExample< Dune::ALUGrid< 3, 3, Dune::simplex, Dune::conforming, Dune::No_Comm >,
                                                       Dune::GDT::ChooseSpaceBackend::fem,
                                                       Dune::Stuff::LA::ChooseBackend::eigen_sparse >;

# if HAVE_MPI

template class GenericLinearellipticMultiscaleExample< Dune::ALUGrid< 3, 3, Dune::simplex, Dune::conforming, MPI_Comm >,
                                                       Dune::GDT::ChooseSpaceBackend::fem,
                                                       Dune::Stuff::LA::ChooseBackend::eigen_sparse >;

# endif // HAVE_MPI
#endif // HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_EIGEN
