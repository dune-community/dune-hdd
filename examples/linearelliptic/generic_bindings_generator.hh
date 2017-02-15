#ifndef DUNE_HDD_EXAMPLES_LINEARELLIPTIC_GENERIC_BINDINGS_GENERATOR_HH
#define DUNE_HDD_EXAMPLES_LINEARELLIPTIC_GENERIC_BINDINGS_GENERATOR_HH

#include <dune/stuff/common/disable_warnings.hh>
# include "config.h"

# if HAVE_ALUGRID
#   include <dune/grid/alugrid.hh>
# endif
# if HAVE_DUNE_SPGRID
# include <dune/grid/spgrid.hh
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
  using typename BaseType::DiscretizationType;

  template< class... Args >
  explicit PbGenericLinearellipticExample(Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
  {}

  DiscretizationType* pb_discretization_and_return_ptr()
  {
    return new DiscretizationType(this->discretization());
  }
}; // class PbGenericLinearellipticExample

#endif // DUNE_HDD_EXAMPLES_LINEARELLIPTIC_GENERIC_BINDINGS_GENERATOR_HH
