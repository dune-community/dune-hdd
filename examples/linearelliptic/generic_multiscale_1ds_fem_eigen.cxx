#include "config.h"

#if HAVE_DUNE_FEM && HAVE_EIGEN

# include "generic_multiscale.hxx"


template class GenericLinearellipticMultiscaleExample< Dune::SGrid< 1, 1 >,
                                                       Dune::GDT::ChooseSpaceBackend::fem,
                                                       Dune::Stuff::LA::ChooseBackend::eigen_sparse >;

#endif // HAVE_DUNE_FEM && HAVE_EIGEN
