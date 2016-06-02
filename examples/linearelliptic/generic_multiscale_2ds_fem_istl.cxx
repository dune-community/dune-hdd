#include "config.h"

#if HAVE_DUNE_FEM && HAVE_DUNE_ISTL

# include "generic_multiscale.hxx"


template class GenericLinearellipticMultiscaleExample< Dune::SGrid< 2, 2 >,
                                                       Dune::GDT::ChooseSpaceBackend::fem,
                                                       Dune::Stuff::LA::ChooseBackend::istl_sparse >;

#endif // HAVE_DUNE_FEM && HAVE_DUNE_ISTL
