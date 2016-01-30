#ifndef DUNE_HDD_EXAMPLES_LINEARELLIPTIC_SPE10_HH
#define DUNE_HDD_EXAMPLES_LINEARELLIPTIC_SPE10_HH

#include <dune/stuff/common/fvector.hh>
#include <dune/stuff/common/string.hh>

#include "generic.hh"

template< class GridType, Dune::GDT::ChooseSpaceBackend space_backend, Dune::Stuff::LA::ChooseBackend la_backend >
class Spe10Example
  : GenericLinearellipticExample< GridType, space_backend, la_backend >
{
  typedef GenericLinearellipticExample< GridType, space_backend, la_backend > BaseType;

  static Dune::Stuff::Common::Configuration grid_config(const std::string& num_grid_elements)
  {
    auto cfg = BaseType::grid_options("stuff.grid.provider.cube");
    cfg["type"] = "stuff.grid.provider.cube";
    cfg["num_elements"] = num_grid_elements;
    cfg["upper_right"] = "[2 5 1]";
    return cfg;
  }

  static Dune::Stuff::Common::Configuration boundary_config()
  {
    Dune::Stuff::Common::Configuration cfg;
    cfg["type"] = "stuff.grid.boundaryinfo.normalbased";
    cfg["default"] = "neumann";
    cfg["dirichlet.0"] = "[0 1 0]";
    return cfg;
  }

  static Dune::Stuff::Common::Configuration problem_config(const std::string& num_grid_elements)
  {
    Dune::Stuff::Common::Configuration cfg;
    cfg["type"] = "hdd.linearelliptic.problem.default";
    cfg["diffusion_factor.type"] = "stuff.function.constant";
    cfg["diffusion_factor.name"] = "diffusion_factor";
    cfg["diffusion_factor.value"] = "1";
    cfg["diffusion_tensor.type"] = "pymor.function.spe10.model2";
    cfg["diffusion_tensor.upper_right"] = "[2 5 1]";
    cfg["diffusion_tensor.blockade_width"] = "0.2";
    cfg["force.type"] = "pymor.function.affinelydecomposabledefault";
    cfg["force.name"] = "force";
    cfg["force.component.0.type"] = "stuff.function.indicator";
    cfg["force.component.0.name"] = "sink";
    cfg["force.component.0.0.domain"] = "[0.167 0.5; 1.82 2.05; 0 1]";
    cfg["force.component.0.0.value"] = "-25";
    cfg["force.coefficient.0.sink"] = "1";
    cfg["force.coefficient.0.expression"] = "sink[0]";
    cfg["dirichlet.type"] = "stuff.function.constant";
    cfg["dirichlet.name"] = "dirichlet";
    cfg["dirichlet.value"] = "0";
    cfg["neumann.type"] = "stuff.function.indicator";
    cfg["neumann.name"] = "neumann";
    cfg["neumann.0.domain"]
        = std::string("[-999 999; -999 ")
          + DSC::toString(5./(2*DSC::fromString<Dune::Stuff::Common::FieldVector<size_t, 3>>(num_grid_elements)[2]))
          + "; -999 999]";
    cfg["neumann.0.value"] = "1";
    return cfg;
  }

  static Dune::Stuff::Common::Configuration logger_config(const bool logging_enabled)
  {
    auto cfg = BaseType::logger_options();
    cfg["info_color"] = "blue";
    cfg["info"] = logging_enabled ? "99" : "-1";
    return cfg;
  }

public:
  using typename BaseType::DiscretizationType;

  Spe10Example(const std::string& num_grid_elements = "[60 220 85]", const bool logging_enabled = true)
    : BaseType(logger_config(logging_enabled),
               grid_config(num_grid_elements),
               boundary_config(),
               problem_config(num_grid_elements))
  {}

  using BaseType::discretization;
  using BaseType::visualize_grid;
  using BaseType::visualize_problem;
  using BaseType::visualize_darcy_velocity;
}; // class Spe10Example


#endif // DUNE_HDD_EXAMPLES_LINEARELLIPTIC_SPE10_HH
