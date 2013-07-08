// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "cg_gdt.hh"

const std::string LinearellipticExampleCG::static_id()
{
  return "example.linearelliptic.cg.gdt";
}

void LinearellipticExampleCG::writeSettingsFile(const std::string filename)
{
  LinearellipticExampleCG::ProblemType::writeSettingsFile(filename, LinearellipticExampleCG::static_id());

}

LinearellipticExampleCG::LinearellipticExampleCG(const std::vector< std::string > arguments)
  throw (Dune::Pymor::Exception::this_does_not_make_any_sense)
  : problem_(LinearellipticExampleCG::static_id(),
             arguments.size(),
             Dune::Stuff::Common::String::vectorToMainArgs(arguments))
{
  const bool debugLogging = problem_.debugLogging();
  Dune::Stuff::Common::LogStream& info  = Dune::Stuff::Common::Logger().info();
  Dune::Stuff::Common::LogStream& debug = Dune::Stuff::Common::Logger().debug();
  Dune::Timer timer;
  // check
  if (problem_.model()->parametric() && !problem_.model()->affineparametric())
    DUNE_PYMOR_THROW(Dune::Pymor::Exception::this_does_not_make_any_sense,
                     "only implemented for nonparametric or affineparametric models!");
  // grid part
  const std::shared_ptr< const GridPartType > gridPart(new GridPartType(*(problem_.grid())));
  info << "initializing solver";
  if (!debugLogging)
    info << "... " << std::flush;
  else
    info << ":" << std::endl;
  timer.reset();
  discretization_ = std::make_shared< DiscretizationType >(gridPart, problem_.boundaryInfo(), problem_.model());
  discretization_->init(debug, "  ");
  if (!debugLogging)
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;
}

bool LinearellipticExampleCG::parametric() const
{
  return problem_.model()->parametric();
}

Dune::Pymor::ParameterType LinearellipticExampleCG::parameter_type() const
{
  if (problem_.model()->parametric()) {
    std::vector< std::string > keys;
    std::vector< int > values;
    if (problem_.model()->diffusion()->parametric()) {
      keys.push_back("diffusion");
      values.push_back(problem_.model()->diffusion()->paramSize());
    }
    if (problem_.model()->force()->parametric()) {
      keys.push_back("force");
      values.push_back(problem_.model()->force()->paramSize());
    }
    if (problem_.model()->dirichlet()->parametric()) {
      keys.push_back("dirichlet");
      values.push_back(problem_.model()->dirichlet()->paramSize());
    }
    if (problem_.model()->neumann()->parametric()) {
      keys.push_back("neumann");
      values.push_back(problem_.model()->neumann()->paramSize());
    }
    return Dune::Pymor::ParameterType(keys, values);
  } else
    return Dune::Pymor::ParameterType();
}

Dune::Pymor::LA::EigenDenseVector* LinearellipticExampleCG::solve(const Dune::Pymor::Parameter mu) const
  throw (Dune::Pymor::Exception::wrong_parameter_type)
{
  if (mu.type() != parameter_type())
    DUNE_PYMOR_THROW(Dune::Pymor::Exception::wrong_parameter_type,
                     "type of mu (" << mu.type().report() << ") does not match the parameter_type of this ("
                     << parameter_type().report() << ")!");
  return nullptr;
}
