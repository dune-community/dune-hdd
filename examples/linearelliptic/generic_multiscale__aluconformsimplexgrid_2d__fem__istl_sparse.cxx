#ifdef ENABLE_MPI
# undef ENABLE_MPI
#endif
#define ENABLE_MPI 0

#ifdef ENABLE_PARMETIS
# undef ENABLE_PARMETIS
#endif
#define ENABLE_PARMETIS 0

#include "config.h"

#if HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_ISTL

#include "generic_multiscale.hxx"


template class GenericLinearellipticMultiscaleExample< Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming >,
                                                       Dune::GDT::ChooseSpaceBackend::fem,
                                                       Dune::Stuff::LA::ChooseBackend::istl_sparse >;


#endif // HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_ISTL
