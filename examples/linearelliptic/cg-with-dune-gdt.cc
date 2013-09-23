// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include "config.h"
#endif

#if HAVE_ALUGRID
  #include <dune/grid/alugrid.hh>
#else
  #warning "AluGrid not found, defaulting to SGrid!"
  #include <dune/stuff/common/disable_warnings.hh>
  #include <dune/grid/sgrid.hh>
  #include <dune/stuff/common/reenable_warnings.hh>
  #include <dune/stuff/common/string.hh>
#endif

#include "cg-with-dune-gdt.hh"

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

#if HAVE_ALUGRID
  template class LinearellipticExampleCG< Dune::ALUSimplexGrid< 2, 2 >, 1 >;
  template class LinearellipticExampleCG< Dune::ALUConformGrid< 2, 2 >, 1 >;
#else
  template class LinearellipticExampleCG< Dune::SGrid< 2, 2 >, 1 >;
#endif
