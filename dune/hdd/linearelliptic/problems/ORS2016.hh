// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_PROBLEMS_ORS2016_HH
#define DUNE_HDD_LINEARELLIPTIC_PROBLEMS_ORS2016_HH

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
namespace internal {


template< class E, class D, class R >
class BatteryGeometry
{
  typedef Problems::Default< E, D, 3, R, 1 >                   ProblemType;
  typedef typename ProblemType::DomainType                     DomainType;
  typedef Stuff::LocalizableFunctionInterface< E, D, 3, R, 1 > NonparametricFunctionType;
  typedef typename ProblemType::DiffusionFactorType            DiffusionFactorType;

public:
  BatteryGeometry(const std::string& filename,
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
    separator_ = std::shared_ptr< Stuff::Functions::DomainIndicator< E, D, 3, R, 1 > >(
          Stuff::Functions::DomainIndicator< E, D, 3, R, 1 >::create(DSC::Configuration({"0.domain",        "0.value", "name"},
                                                                                        {separator_domain,  "1",       "SEPARATOR"})));
    // create one component per battery part
    typedef Stuff::Functions::Constant< E, D, 3, R, 1, 1 >    ConstantFunctionType;
    typedef Stuff::Functions::LevelIndicator< E, D, 3, R, 1 > LevelIndicator;
    anode_      = std::make_shared< LevelIndicator >(battery_function, 0., 0., 1., "ANODE");
    cathode_    = std::make_shared< LevelIndicator >(battery_function, 1., 1., 1., "CATHODE");
    cc_anode_   = std::make_shared< LevelIndicator >(battery_function, 2., 2., 1., "CC_ANODE");
    cc_cathode_ = std::make_shared< LevelIndicator >(battery_function, 3., 3., 1., "CC_CATHODE");
    //   for the electrolyte, substract the seperator
    electrolyte_ = std::make_shared< LevelIndicator >(
                     Stuff::Functions::make_difference(battery_function,
                                                       Stuff::Functions::make_product(std::make_shared< ConstantFunctionType >(4.),
                                                                                      separator_)),
                     4., 4., 1., "ELECTROLYTE");
    // create parametric diffusion factor
    battery_geometry_ = std::make_shared< Pymor::Functions::AffinelyDecomposableDefault< E, D, 3, R, 1 > >();
//    battery_geometry_->register_component(anode_,       new Pymor::ParameterFunctional("ANODE",       1, "ANODE[0]"));
//    battery_geometry_->register_component(cathode_,     new Pymor::ParameterFunctional("CATHODE",     1, "CATHODE[0]"));
//    battery_geometry_->register_component(cc_anode_,    new Pymor::ParameterFunctional("CC_ANODE",    1, "CC_ANODE[0]"));
//    battery_geometry_->register_component(cc_cathode_,  new Pymor::ParameterFunctional("CC_CATHODE",  1, "CC_CATHODE[0]"));
    battery_geometry_->register_component(electrolyte_, new Pymor::ParameterFunctional("ELECTROLYTE", 1, "ELECTROLYTE[0]"));
//    battery_geometry_->register_component(separator_,   new Pymor::ParameterFunctional("SEPARATOR",   1, "SEPARATOR[0]"));
    battery_geometry_->register_affine_part(
          Stuff::Functions::make_sum(Stuff::Functions::make_product(std::make_shared< ConstantFunctionType >(1.04),   anode_),
          Stuff::Functions::make_sum(Stuff::Functions::make_product(std::make_shared< ConstantFunctionType >(1.58),   cathode_),
          Stuff::Functions::make_sum(Stuff::Functions::make_product(std::make_shared< ConstantFunctionType >(238.),   cc_anode_),
          Stuff::Functions::make_sum(Stuff::Functions::make_product(std::make_shared< ConstantFunctionType >(398.),   cc_cathode_),
                                     Stuff::Functions::make_product(std::make_shared< ConstantFunctionType >(0.3344), separator_)))),
                                     "NON_ELECTROLYTE"));
  } // BatteryGeometry(...)

protected:
  std::shared_ptr< NonparametricFunctionType > anode_;
  std::shared_ptr< NonparametricFunctionType > cathode_;
  std::shared_ptr< NonparametricFunctionType > cc_anode_;
  std::shared_ptr< NonparametricFunctionType > cc_cathode_;
  std::shared_ptr< NonparametricFunctionType > electrolyte_;
  std::shared_ptr< NonparametricFunctionType > separator_;
  std::shared_ptr< Pymor::Functions::AffinelyDecomposableDefault< E, D, 3, R, 1 > > battery_geometry_;
}; // class BatteryGeometry


} // namespace internal


template< class E, class D, int d, class R, int r = 1 >
class ORS2016
  : public ProblemInterface< E, D, d, R, r >
{
  ORS2016() { static_assert(AlwaysFalse< E >::value, "Not available for these dimensions!"); }
};


template< class E, class D, class R >
class ORS2016< E, D, 3, R, 1 >
  : internal::BatteryGeometry< E, D, R >
  , public Problems::Default< E, D, 3, R, 1 >
{
  typedef internal::BatteryGeometry< E, D, R > DataType;
  typedef Problems::Default< E, D, 3, R, 1 > BaseType;
  typedef ORS2016< E, D, 3, R, 1 >           ThisType;
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
    return BaseType::BaseType::static_id() + ".ORS2016";
  }

  static Stuff::Common::Configuration default_config(const std::string sub_name = "")
  {
    Stuff::Common::Configuration config;
    for (auto sub : {"diffusion_tensor", "dirichlet", "neumann"})
      config.add(BaseType::default_config().sub(sub), sub);
    config["type"] = static_id();
    config["diffusion_factor.filename"]    = "geometry";
    config["diffusion_factor.lower_left"]  = "[0.0 0.0 0.0]";
    config["diffusion_factor.upper_right"] = "[0.0184 0.008 0.008]";
    config["diffusion_factor.num_elements"] = "[46 20 20]";
    config["diffusion_factor.name"] = "battery_geometry";
    config["diffusion_factor.separator"] = "[0.0084 0.01; 0 0.008; 0 0.008]";
    config["force.value"] = "1";
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
                                        cfg.get("force.value",
                                                def_cfg.get<RangeFieldType>("force.value")),
                                        BaseType::create_matrix_function("diffusion_tensor", cfg),
                                        BaseType::create_vector_function("dirichlet", cfg),
                                        BaseType::create_vector_function("neumann", cfg));
  } // ... create(...)

  ORS2016(const std::string& filename,
          const DomainType& lower_left,
          const DomainType& upper_right,
          const Stuff::Common::FieldVector< size_t, 3 >& num_elements,
          const std::string& separator,
          const std::string& name,
          const RangeFieldType& force_value,
          const std::shared_ptr< const DiffusionTensorType >& diff_ten,
          const std::shared_ptr< const FunctionType >& dir,
          const std::shared_ptr< const FunctionType >& neum)
    : DataType(filename, lower_left, upper_right, num_elements, separator, name)
    , BaseType(DataType::battery_geometry_,
               diff_ten,
               std::make_shared< Pymor::Functions::AffinelyDecomposableDefault< E, D, 3, R, 1 > >(
                  Stuff::Functions::make_product(std::make_shared< Stuff::Functions::Constant< E, D, 3, R, 1 > >(force_value),
                                                 Stuff::Functions::make_sum(DataType::anode_, DataType::cathode_),
                                                 "force")),
               dir,
               neum)
  {}
}; // class ORS2016< ..., 3, ... 1 >


} // namespace Problems
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_PROBLEMS_SPE10_HH
