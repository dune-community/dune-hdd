// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_PROBLEMS_BATTERY_HH
#define DUNE_HDD_LINEARELLIPTIC_PROBLEMS_BATTERY_HH

#include <iostream>
#include <vector>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/stuff/functions/checkerboard.hh>
#include <dune/stuff/playground/functions/indicator.hh>

#include <dune/pymor/functions/default.hh>

#include "default.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Problems {


template< class E, class D, int d, class R, int r = 1 >
class Battery
  : public ProblemInterface< E, D, d, R, r >
{
  static_assert(AlwaysFalse< E >::value, "Not available for these dimensions!");
};


template< class E, class D, class R >
class Battery< E, D, 3, R, 1 >
  : public Problems::Default< E, D, 3, R, 1 >
{
  typedef Problems::Default< E, D, 3, R, 1 > BaseType;
  typedef Battery< E, D, 3, R, 1 >           ThisType;
public:
  typedef typename BaseType::EntityType      EntityType;
  typedef typename BaseType::DomainFieldType DomainFieldType;
  static const unsigned int                  dimDomain = BaseType::dimDomain;
  typedef typename BaseType::DomainType      DomainType;
  typedef typename BaseType::RangeFieldType  RangeFieldType;
  static const unsigned int                  dimRange = BaseType::dimRange;
  using typename BaseType::FunctionType;
  using typename BaseType::DiffusionFactorType;
  using typename BaseType::DiffusionTensorType;

  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::BaseType::static_id() + ".battery";
  }

  static Stuff::Common::Configuration default_config(const std::string sub_name = "")
  {
    auto config = BaseType::default_config();
    config["diffusion_factor."] = "";
    config["diffusion_factor.filename"]    = "geometry";
    config["diffusion_factor.lower_left"]  = "[0.0 0.0 0.0]";
    config["diffusion_factor.upper_right"] = "[0.0184 0.008 0.008]";
    config["diffusion_factor.num_elements"] = "[46 20 20]";
    config["diffusion_factor.name"] = "battery_geometry";
    config["diffusion_factor.separator"] = "[0.0084 0.01; 0 0.008; 0 0.008]";
    if (sub_name.empty())
      return config;
    else {
      Stuff::Common::Configuration tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  static std::unique_ptr< ThisType > create(const Stuff::Common::Configuration config = default_config(),
                                            const std::string sub_name = static_id())
  {
    const Stuff::Common::Configuration cfg = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Stuff::Common::Configuration def_cfg = default_config();
    return DSC::make_unique< ThisType >(cfg.get("diffusion_factor.filename", def_cfg.get<std::string>("diffusion_factor.filename")),
                                        cfg.get("diffusion_factor.lower_left",
                                                def_cfg.get<DomainType>("diffusion_factor.lower_left")),
                                        cfg.get("diffusion_factor.upper_right",
                                                def_cfg.get<DomainType>("diffusion_factor.upper_right")),
                                        cfg.get("diffusion_factor.num_elements",
                                                def_cfg.get<Stuff::Common::FieldVector< size_t, 3 >>("diffusion_factor.num_elements")),
                                        cfg.get("diffusion_factor.separator",
                                                def_cfg.get<std::string>("diffusion_factor.separator")),
                                        cfg.get("diffusion_factor.name",
                                                def_cfg.get<std::string>("diffusion_factor.name")),
                                        BaseType::create_matrix_function("diffusion_tensor", cfg),
                                        BaseType::create_vector_function("force", cfg),
                                        BaseType::create_vector_function("dirichlet", cfg),
                                        BaseType::create_vector_function("neumann", cfg));
  } // ... create(...)

  Battery(const std::string& filename,
          const DomainType& lower_left,
          const DomainType& upper_right,
          const Stuff::Common::FieldVector< size_t, 3 >& num_elements,
          const std::string& separator,
          const std::string& name,
          const std::shared_ptr< const DiffusionTensorType >& diff_ten,
          const std::shared_ptr< const FunctionType >& forc,
          const std::shared_ptr< const FunctionType >& dir,
          const std::shared_ptr< const FunctionType >& neum)
    : BaseType(make_diffusion_factor(filename, lower_left, upper_right, num_elements, separator, name),
               diff_ten, forc, dir, neum)
  {}

private:
  static std::shared_ptr< DiffusionFactorType > make_diffusion_factor(const std::string& filename,
                                                                      const DomainType& lower_left,
                                                                      const DomainType& upper_right,
                                                                      const Stuff::Common::FieldVector< size_t, 3 >& num_elements,
                                                                      const std::string& separator_domain,
                                                                      const std::string& name)
  {
    // Interpretation of the data as int:
    //    enum Subdomain
    //    {
    //      ANODE,      // 0
    //      CATHODE,    // 1
    //      CC_ANODE,   // 2
    //      CC_CATHODE, // 3
    //      ELECTROLYTE // 4
    //    };
    static_assert(sizeof(int) == 4, "This is not the correct architecture, sizeof(int) should be 4!");

    size_t num_entries = 1;
    for (size_t ii = 0; ii < 3; ++ii)
      num_entries *= num_elements[ii];
    std::vector< int > raw_data(num_entries);

    std::ifstream file(filename, std::ifstream::binary);
    if (!file)
      DUNE_THROW(IOError, "Could not open '" << filename << "'!");
    // check for the amount of data present in the file
    file.seekg(0, file.end);
    const auto length = file.tellg();
    if (length != boost::numeric_cast< decltype(length) >(raw_data.size()*sizeof(int)))
      DUNE_THROW(IOError,
                 "Given file '" << filename << "' has wrong size (should be " << raw_data.size()*sizeof(int) << ", is "
                 << length << ")!");
    // read the actual data
    file.seekg(0, file.beg);
    file.read((char *)(raw_data.data()), length);
    // check for the amount of read data
    const auto read_data = file.gcount();
    if (read_data != boost::numeric_cast< decltype(read_data) >(raw_data.size()*sizeof(int)))
      DUNE_THROW(IOError,
                 "Could not read the correct amount from given file '" << filename << "' (we would like to read "
                 << raw_data.size()*sizeof(int) << " and the file reported this many entries, but we could only read "
                 << read_data << "data)!");
    if (!file)
      DUNE_THROW(IOError, "Failed to read from file '" << filename << "'!");
    file.close();

    // create all but separator
    typedef Stuff::Functions::Checkerboard< E, D, 3, R, 1, 1 > PiecewiseConstantFunctionType;
    std::vector< typename PiecewiseConstantFunctionType::RangeType > data(num_entries);
    for (size_t ii = 0; ii < num_entries; ++ii)
      data[ii] = raw_data[ii];
    auto battery_function = std::make_shared< PiecewiseConstantFunctionType >(lower_left, upper_right, num_elements, data, name);

    // create separator
    auto separator = std::shared_ptr< Stuff::Functions::DomainIndicator< E, D, 3, R, 1 > >(
          Stuff::Functions::DomainIndicator< E, D, 3, R, 1 >::create(DSC::Configuration({"0.domain",        "0.value", "name"},
                                                                                        {separator_domain,  "1",       "SEPARATOR"})));

    // create parametric diffusion factor, one component per battery part
    typedef Stuff::Functions::Constant< E, D, 3, R, 1, 1 >    ConstantFunctionType;
    typedef Stuff::Functions::LevelIndicator< E, D, 3, R, 1 > LevelIndicator;
    auto diffusion_factor = std::make_shared< Pymor::Functions::AffinelyDecomposableDefault< E, D, 3, R, 1 > >();
    diffusion_factor->register_component(new LevelIndicator(battery_function, 0., 0., 1., "ANODE"),
                                         new Pymor::ParameterFunctional("ANODE", 1, "ANODE[0]"));
    diffusion_factor->register_component(new LevelIndicator(battery_function, 1., 1., 1., "CATHODE"),
                                         new Pymor::ParameterFunctional("CATHODE", 1, "CATHODE[0]"));
    diffusion_factor->register_component(new LevelIndicator(battery_function, 2., 2., 1., "CC_ANODE"),
                                         new Pymor::ParameterFunctional("CC_ANODE", 1, "CC_ANODE[0]"));
    diffusion_factor->register_component(new LevelIndicator(battery_function, 3., 3., 1., "CC_CATHODE"),
                                         new Pymor::ParameterFunctional("CC_CATHODE", 1, "CC_CATHODE[0]"));
    //   for the electrolyte, substract the seperator
    diffusion_factor->register_component(new LevelIndicator(
        Stuff::Functions::make_difference(battery_function,
                                          Stuff::Functions::make_product(std::make_shared< ConstantFunctionType >(4.),
                                                                         separator)),
        4., 4., 1., "ELECTROLYTE"),
                                         new Pymor::ParameterFunctional("ELECTROLYTE", 1, "ELECTROLYTE[0]"));
    diffusion_factor->register_component(separator, new Pymor::ParameterFunctional("SEPARATOR", 1, "SEPARATOR[0]"));

    return diffusion_factor;
  } // ... make_diffusion_factor(...)
}; // class Battery< ..., 3, ... 1 >


} // namespace Problems
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_PROBLEMS_SPE10_HH
