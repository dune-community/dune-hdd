// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include <dune/stuff/common/string.hh>

#include "cg.hh"

template< class G, int p >
std::string LinearellipticExampleCG< G, p >::static_id()
{
  return "linearelliptic.cg-with-dune-gdt";
}

template< class G, int p >
void LinearellipticExampleCG< G, p >::write_settings_file(const std::string filename)
{
  DiscreteProblemType::writeSettingsFile(filename, static_id());
}

template< class G, int p >
LinearellipticExampleCG< G, p >::LinearellipticExampleCG()
  : initialized_(false)
{}

template< class G, int p >
void LinearellipticExampleCG< G, p >::initialize(const std::vector< std::string > arguments)
{
  discreteProblem_ = std::make_shared< const DiscreteProblemType >(static_id(), arguments);
  const bool debugLogging = discreteProblem_->debugLogging();
  auto& info = DSC_LOG_INFO;
  auto& debug = DSC_LOG_DEBUG;
  discretization_ = std::make_shared< DiscretizationType >(discreteProblem_->gridPart(),
                                                           discreteProblem_->boundaryInfo(),
                                                           discreteProblem_->problem());
  info << "initializing discretization";
  if (debugLogging)
    info << ":" << std::endl;
  else
    info << "... " << std::flush;
  Dune::Timer timer;
  discretization_->initialize(debug, "  ");
  if (!debugLogging)
    info << "done (took " << timer.elapsed() << "s)" << std::endl;
  initialized_ = true;
}

template< class G, int p >
bool LinearellipticExampleCG< G, p >::initialized() const
{
  return initialized_;
}

template< class G, int p >
typename LinearellipticExampleCG< G, p >::DiscreteProblemType LinearellipticExampleCG< G, p >::discrete_problem() const
{
  if (!initialized_)
    DUNE_PYMOR_THROW(Dune::Pymor::Exception::requirements_not_met,
                     "do not call discrete_problem() if initialized() == false!");
  return *discreteProblem_;
}

template< class G, int p >
typename LinearellipticExampleCG< G, p >::DiscretizationType LinearellipticExampleCG< G, p >::discretization() const
{
  if (!initialized_)
    DUNE_PYMOR_THROW(Dune::Pymor::Exception::requirements_not_met,
                     "do not call discretization() if initialized() == false!");
  return *discretization_;
}

template< class G, int p >
typename LinearellipticExampleCG< G, p >::DiscretizationType* LinearellipticExampleCG< G, p >::discretization_and_return_ptr() const
{
  return new DiscretizationType(discretization());
}


// =================
// ===== sgrid =====
// =================
#include <dune/stuff/common/disable_warnings.hh>
  #include <dune/grid/sgrid.hh>
#include <dune/stuff/common/reenable_warnings.hh>

template class LinearellipticExampleCG< Dune::SGrid< 2, 2 >, 1 >;
template class LinearellipticExampleCG< Dune::SGrid< 2, 2 >, 2 >;
template class LinearellipticExampleCG< Dune::SGrid< 2, 2 >, 3 >;
template class LinearellipticExampleCG< Dune::SGrid< 3, 3 >, 1 >;
template class LinearellipticExampleCG< Dune::SGrid< 3, 3 >, 2 >;
template class LinearellipticExampleCG< Dune::SGrid< 3, 3 >, 3 >;

// ====================
// ===== yaspgrid =====
// ====================
#include <dune/grid/yaspgrid.hh>

template class LinearellipticExampleCG< Dune::YaspGrid< 2 >, 1 >;
template class LinearellipticExampleCG< Dune::YaspGrid< 3 >, 1 >;


// ===================
// ===== alugrid =====
// ===================
#ifdef HAVE_ALUGRID
  #include <dune/grid/alugrid.hh>

  template class LinearellipticExampleCG< Dune::ALUGrid< 2, 2, Dune::cube, Dune::nonconforming >, 1 >;
  template class LinearellipticExampleCG< Dune::ALUGrid< 2, 2, Dune::cube, Dune::nonconforming >, 2 >;
  template class LinearellipticExampleCG< Dune::ALUGrid< 2, 2, Dune::cube, Dune::nonconforming >, 3 >;

  template class LinearellipticExampleCG< Dune::ALUGrid< 2, 2, Dune::simplex, Dune::nonconforming >, 1 >;
  template class LinearellipticExampleCG< Dune::ALUGrid< 2, 2, Dune::simplex, Dune::nonconforming >, 2 >;
  template class LinearellipticExampleCG< Dune::ALUGrid< 2, 2, Dune::simplex, Dune::nonconforming >, 3 >;

  template class LinearellipticExampleCG< Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming >, 1 >;
  template class LinearellipticExampleCG< Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming >, 2 >;
  template class LinearellipticExampleCG< Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming >, 3 >;

  template class LinearellipticExampleCG< Dune::ALUGrid< 3, 3, Dune::cube, Dune::nonconforming >, 1 >;
  template class LinearellipticExampleCG< Dune::ALUGrid< 3, 3, Dune::cube, Dune::nonconforming >, 2 >;
  template class LinearellipticExampleCG< Dune::ALUGrid< 3, 3, Dune::cube, Dune::nonconforming >, 3 >;

  template class LinearellipticExampleCG< Dune::ALUGrid< 3, 3, Dune::simplex, Dune::nonconforming >, 1 >;
  template class LinearellipticExampleCG< Dune::ALUGrid< 3, 3, Dune::simplex, Dune::nonconforming >, 2 >;
  template class LinearellipticExampleCG< Dune::ALUGrid< 3, 3, Dune::simplex, Dune::nonconforming >, 3 >;
#endif // HAVE_ALUGRID
