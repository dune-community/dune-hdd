#include "config.h"

#if HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_ISTL

# include "generic_multiscale.hxx"


template class GenericLinearellipticMultiscaleExample< Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming, Dune::No_Comm >,
                                                       Dune::GDT::ChooseSpaceBackend::fem,
                                                       Dune::Stuff::LA::ChooseBackend::istl_sparse >;

# if HAVE_MPI

template class GenericLinearellipticMultiscaleExample< Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming, MPI_Comm >,
                                                       Dune::GDT::ChooseSpaceBackend::fem,
                                                       Dune::Stuff::LA::ChooseBackend::istl_sparse >;

# endif // HAVE_MPI
#endif // HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_ISTL
