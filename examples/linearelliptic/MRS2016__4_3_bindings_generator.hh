#ifndef DUNE_HDD_EXAMPLES_LINEARELLIPTIC_MRS2016__4_3_BINDINGS_GENERATOR_HH
#define DUNE_HDD_EXAMPLES_LINEARELLIPTIC_MRS2016__4_3_BINDINGS_GENERATOR_HH

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

# include "MRS2016__4_3.hh"


template< class G, Dune::GDT::ChooseSpaceBackend s, Dune::Stuff::LA::ChooseBackend l >
class PbSpe10Example
  : public Spe10Example< G, s, l >
{
  typedef Spe10Example< G, s, l > BaseType;
public:
  using typename BaseType::DiscretizationType;

  template< class... Args >
  explicit PbSpe10Example(Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
  {}

  DiscretizationType* pb_discretization_and_return_ptr()
  {
    return new DiscretizationType(this->discretization());
  }
}; // class PbGenericLinearellipticExample

#endif // DUNE_HDD_EXAMPLES_LINEARELLIPTIC_MRS2016__4_3_BINDINGS_GENERATOR_HH
