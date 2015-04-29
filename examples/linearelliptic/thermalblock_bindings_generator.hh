// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_EXAMPLES_LINEARELLIPTIC_THERMALBLOCK_BINDINGS_GENERATOR_HH
#define DUNE_HDD_EXAMPLES_LINEARELLIPTIC_THERMALBLOCK_BINDINGS_GENERATOR_HH

#ifdef ENABLE_MPI
# undef ENABLE_MPI
#endif
#define ENABLE_MPI 0

#include <dune/stuff/common/disable_warnings.hh>
# include "config.h"

# include "thermalblock.hh"

# if HAVE_ALUGRID
#   include <dune/grid/alugrid.hh>
# endif

template< class G, Dune::GDT::ChooseSpaceBackend s, Dune::Stuff::LA::ChooseBackend l >
class PbCgExample
  : public CgExample< G, s, l >
{
  typedef CgExample< G, s, l > BaseType;

public:
  using typename BaseType::DiscretizationType;

  template< class... Args >
  explicit PbCgExample(Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
  {}

  DiscretizationType* pb_discretization_and_return_ptr()
  {
    return new DiscretizationType(this->discretization());
  }
}; // class PbCgExample

#endif // DUNE_HDD_EXAMPLES_LINEARELLIPTIC_THERMALBLOCK_BINDINGS_GENERATOR_HH
