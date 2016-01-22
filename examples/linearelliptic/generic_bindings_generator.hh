#ifndef DUNE_HDD_EXAMPLES_LINEARELLIPTIC_GENERIC_BINDINGS_GENERATOR_HH
#define DUNE_HDD_EXAMPLES_LINEARELLIPTIC_GENERIC_BINDINGS_GENERATOR_HH

#ifdef ENABLE_MPI
# undef ENABLE_MPI
#endif
#define ENABLE_MPI 0

#ifdef ENABLE_PARMETIS
# undef ENABLE_PARMETIS
#endif
#define ENABLE_PARMETIS 0

#include <dune/stuff/common/disable_warnings.hh>
# include "config.h"

# if HAVE_ALUGRID
#   include <dune/grid/alugrid.hh>
# endif
# if HAVE_DUNE_SPGRID
# include <dune/grid/spgrid.hh>
# endif
# include <dune/grid/yaspgrid.hh>

# include <dune/pymor/bindings/pymor.hh>

# include "generic.hh"


template< class G, Dune::GDT::ChooseSpaceBackend s, Dune::Stuff::LA::ChooseBackend l >
class PbGenericLinearellipticExample
  : public GenericLinearellipticExample< G, s, l >
{
  typedef GenericLinearellipticExample< G, s, l > BaseType;
public:
  using typename BaseType::CgDiscretizationType;
  using typename BaseType::SwipdgDiscretizationType;

  template< class... Args >
  explicit PbGenericLinearellipticExample(Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
  {}

  CgDiscretizationType* pb_cg_discretization_and_return_ptr()
  {
    return new CgDiscretizationType(this->cg_discretization());
  }

  SwipdgDiscretizationType* pb_swipdg_discretization_and_return_ptr()
  {
    return new SwipdgDiscretizationType(this->swipdg_discretization());
  }
}; // class PbGenericLinearellipticExample

#endif // DUNE_HDD_EXAMPLES_LINEARELLIPTIC_GENERIC_BINDINGS_GENERATOR_HH
