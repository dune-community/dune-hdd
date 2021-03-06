// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_PROBLEMS_OS2015_HH
#define DUNE_HDD_LINEARELLIPTIC_PROBLEMS_OS2015_HH

#include <memory>
#include <vector>
#include <utility>
#include <sstream>
#include <cmath>

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/common/fvector.hh>
#include <dune/stuff/common/float_cmp.hh>
#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/combined.hh>
#include <dune/stuff/functions/flattop.hh>
#include <dune/stuff/functions/spe10.hh>
#include <dune/stuff/functions/global.hh>
#include <dune/stuff/functions/expression.hh>
#include <dune/stuff/functions/ESV2007.hh>
#include <dune/stuff/playground/functions/indicator.hh>

#include <dune/pymor/functions/default.hh>

#include "default.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Problems {
namespace OS2015 {


template< class E, class D, int d, class R, int r = 1 >
class ParametricESV2007
  : public ProblemInterface< E, D, d, R, r >
{
  ParametricESV2007() { static_assert(AlwaysFalse< E >::value, "Not available for dimRange > 1!"); }
};


template< class EntityImp, class DomainFieldImp, class RangeFieldImp >
class ParametricESV2007< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 >
  : public Default< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 >
{
  typedef Default< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 >            BaseType;
  typedef ParametricESV2007< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 >  ThisType;

  typedef Stuff::Functions::Constant< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 >    ConstantScalarFunctionType;
  typedef Stuff::Functions::Constant< EntityImp, DomainFieldImp, 2, RangeFieldImp, 2, 2 > ConstantMatrixFunctionType;
  typedef Pymor::Functions::NonparametricDefault< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 >
      ParametricScalarFunctionType;
  typedef Pymor::Functions::NonparametricDefault< EntityImp, DomainFieldImp, 2, RangeFieldImp, 2, 2 >
      ParametricMatrixFunctionType;

  typedef Stuff::Functions::ESV2007::Testcase1Force< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 > ForceType;
  typedef Stuff::GlobalLambdaFunction< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 > LambdaFunctionType;
  typedef Pymor::Functions::AffinelyDecomposableDefault< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 >
      DefaultParametricFunctionType;

  typedef typename ConstantMatrixFunctionType::RangeType MatrixType;

  static MatrixType unit_matrix()
  {
    MatrixType matrix(RangeFieldImp(0));
    matrix[0][0] = RangeFieldImp(1);
    matrix[1][1] = RangeFieldImp(1);
    return matrix;
  } // ... unit_matrix(...)

  static std::shared_ptr< DefaultParametricFunctionType > create_diffusion_factor(const size_t integration_order)
  {
    auto ret = std::make_shared< DefaultParametricFunctionType >("diffusion_factor");
    ret->register_affine_part(new LambdaFunctionType(
                                [](const typename LambdaFunctionType::DomainType& xx) {
                                  return 1.0 + (std::cos(0.5*M_PIl*xx[0])*std::cos(0.5*M_PIl*xx[1]));
                                },
                                integration_order,
                                "affine_part"));
    ret->register_component(new LambdaFunctionType(
                              [](const typename LambdaFunctionType::DomainType& xx){
                                return -1.0*(std::cos(0.5*M_PIl*xx[0])*std::cos(0.5*M_PIl*xx[1]));
                              },
                              integration_order,
                              "component_0"),
                            new Pymor::ParameterFunctional("mu", 1, "mu[0]"));
    return ret;
  } // ... create_diffusion_factor(...)

public:
  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::BaseType::static_id() + ".OS2015.parametricESV2007";
  }

  static Stuff::Common::Configuration default_config(const std::string sub_name = "")
  {
    Stuff::Common::Configuration config("integration_order", "3");
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
    return Stuff::Common::make_unique< ThisType >(cfg.get("integration_order",
                                                          default_config().template get< size_t >("integration_order")));
  } // ... create(...)

  ParametricESV2007(const size_t integration_order = default_config().template get< size_t >("integration_order"))
    : BaseType(create_diffusion_factor(integration_order),
               std::make_shared< ParametricMatrixFunctionType >(new ConstantMatrixFunctionType(unit_matrix(),
                                                                                                  "diffusion_tensor")),
               std::make_shared< ParametricScalarFunctionType >(new ForceType(integration_order,  "force")),
               std::make_shared< ParametricScalarFunctionType >(new ConstantScalarFunctionType(0, "dirichlet")),
               std::make_shared< ParametricScalarFunctionType >(new ConstantScalarFunctionType(0, "neumann")))
  {}

  virtual std::string type() const override
  {
    return BaseType::BaseType::static_id() + ".OS2015.parametricESV2007";
  }
}; // class ParametricESV2007< ..., 2, ..., 1 >


template< class E, class D, int d, class R, int r = 1 >
class Spe10Model1
  : public ProblemInterface< E, D, d, R, r >
{
  Spe10Model1() { static_assert(AlwaysFalse< E >::value, "Not available for these dimensions!"); }
};


template< class E, class D, class R >
class Spe10Model1< E, D, 2, R, 1 >
  : public Problems::Default< E, D, 2, R, 1 >
{
  typedef Problems::Default< E, D, 2, R, 1 > BaseType;
  typedef Spe10Model1< E, D, 2, R, 1 >       ThisType;

  typedef Stuff::Functions::Constant< E, D, 2, R, 1 >         ConstantFunctionType;
  typedef Stuff::Functions::FlatTop< E, D, 2, R, 1 >          FlatTopFunctionType;
  typedef Stuff::Functions::DomainIndicator< E, D, 2, R, 1 >  IndicatorFunctionType;
  typedef Stuff::Functions::Spe10::Model1< E, D, 2, R, 2, 2 > Spe10FunctionType;

public:
  typedef typename BaseType::EntityType      EntityType;
  typedef typename BaseType::DomainFieldType DomainFieldType;
  static const unsigned int                  dimDomain = BaseType::dimDomain;
  typedef typename BaseType::DomainType      DomainType;
  typedef typename BaseType::RangeFieldType  RangeFieldType;
  static const unsigned int                  dimRange = BaseType::dimRange;

  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::BaseType::static_id() + ".OS2015.spe10model1";
  }

  static Stuff::Common::Configuration default_config(const std::string sub_name = "")
  {
    std::istringstream ss("# a definition of a channel would be analogue to the one of forces\n"
                          "forces.0.domain = [0.95 1.10; 0.30 0.45]\n"
                          "forces.0.value = 2000\n"
                          "forces.1.domain = [3.00 3.15; 0.75 0.90]\n"
                          "forces.1.value = -1000\n"
                          "forces.2.domain = [4.25 4.40; 0.25 0.40]\n"
                          "forces.2.value = -1000");
    Stuff::Common::Configuration config(ss);
    config["type"] = static_id();
    config["filename"]    = Stuff::Functions::Spe10::internal::model1_filename;
    config["lower_left"]  = "[0.0 0.0]";
    config["upper_right"] = "[5.0 1.0]";
    config["parametric_channel"] = "false";
    config.set("channel_boundary_layer", FlatTopFunctionType::default_config().template get< std::string >("boundary_layer"));
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
    return Stuff::Common::make_unique< ThisType >(
          cfg.get("filename",    def_cfg.get< std::string >("filename")),
          cfg.get("lower_left",  def_cfg.get< DomainType >("lower_left")),
          cfg.get("upper_right", def_cfg.get< DomainType >("upper_right")),
          get_values(cfg, "channel"),
          get_values(cfg, "forces"),
          cfg.get("channel_boundary_layer", def_cfg.get< DomainType >("channel_boundary_layer")),
          cfg.get("parametric_channel", def_cfg.get< bool >("parametric_channel")));
  } // ... create(...)

  Spe10Model1(const std::string filename,
              const DomainType& lower_left,
              const DomainType& upper_right,
              const std::vector< std::tuple< DomainType, DomainType, typename IndicatorFunctionType::RangeType > >& channel_values,
              const std::vector< std::tuple< DomainType, DomainType, typename IndicatorFunctionType::RangeType > >& force_values,
              const DomainType& channel_boundary_layer = default_config().template get< DomainType >("channel_boundary_layer"),
              const bool parametric_channel = default_config().template get< bool >("parametric_channel"))
    : BaseType(create_base(filename,
                           lower_left,
                           upper_right,
                           channel_values,
                           force_values,
                           channel_boundary_layer,
                           parametric_channel))
  {}

private:
  typedef std::vector< std::tuple< DomainType, DomainType, typename IndicatorFunctionType::RangeType > > Values;
  typedef typename BaseType::DiffusionFactorType::NonparametricType FlatTopIndicatorType;

  static BaseType create_base(const std::string filename,
                              const DomainType& lower_left,
                              const DomainType& upper_right,
                              const Values& channel_values,
                              const Values& force_values,
                              const DomainType& channel_boundary_layer,
                              const bool parametric_channel)
  {
    // build the channel as a sum of flattop functions
    std::shared_ptr< FlatTopIndicatorType > channel(nullptr);
    if (channel_values.empty())
      channel = std::make_shared< ConstantFunctionType >(0, "zero");
    else
      channel = create_indicator(channel_values[0], channel_boundary_layer);
    for (size_t ii = 1; ii < channel_values.size(); ++ii)
      channel = Stuff::Functions::make_sum(channel,
                                           create_indicator(channel_values[ii], channel_boundary_layer),
                                           "channel");
    // build the rest
    auto one = std::make_shared< ConstantFunctionType >(1, "one");
    auto diffusion_tensor = std::make_shared< Spe10FunctionType >(filename,
                                                                  lower_left,
                                                                  upper_right,
                                                                  Stuff::Functions::Spe10::internal::model1_min_value,
                                                                  Stuff::Functions::Spe10::internal::model1_max_value,
                                                                  "diffusion_tensor");
    auto force = std::make_shared< IndicatorFunctionType >(force_values, "force");
    auto dirichlet = std::make_shared< ConstantFunctionType >(0, "dirichlet");
    auto neumann = std::make_shared< ConstantFunctionType >(0, "neumann");
    if (parametric_channel) {
      typedef Pymor::Functions::NonparametricDefault< E, D, 2, R, 1 >        ScalarWrapper;
      typedef Pymor::Functions::NonparametricDefault< E, D, 2, R, 2, 2 >     MatrixWrapper;
      typedef Pymor::Functions::AffinelyDecomposableDefault< E, D, 2, R, 1 > ParametricFunctionType;
      auto diffusion_factor = std::make_shared< ParametricFunctionType >("diffusion_factor");
      auto minus_one = std::make_shared< ConstantFunctionType >(-1, "minus_one");
      diffusion_factor->register_affine_part(Stuff::Functions::make_sum(one, channel));
      diffusion_factor->register_component(Stuff::Functions::make_product(minus_one, channel),
                                           new Pymor::ParameterFunctional("mu", 1, "mu[0]"));
      return BaseType(diffusion_factor,
                      std::make_shared< MatrixWrapper >(diffusion_tensor),
                      std::make_shared< ScalarWrapper >(force),
                      std::make_shared< ScalarWrapper >(dirichlet),
                      std::make_shared< ScalarWrapper >(neumann));
    } else {
      auto zero_pt_nine = std::make_shared< ConstantFunctionType >(0.9, "0.9");
      return BaseType(Stuff::Functions::make_sum(one,
                                                 Stuff::Functions::make_product(zero_pt_nine,
                                                                                channel,
                                                                                "scaled_channel"),
                                                 "diffusion_factor"),
                      diffusion_tensor,
                      force,
                      dirichlet,
                      neumann);
    }
  } // ... create_base(...)

  static Values get_values(const Stuff::Common::Configuration& cfg, const std::string id)
  {
    Values values;
    DomainType tmp_lower;
    DomainType tmp_upper;
    if (cfg.has_sub(id)) {
      const Stuff::Common::Configuration sub_cfg = cfg.sub(id);
      size_t cc = 0;
      while (sub_cfg.has_sub(DSC::toString(cc))) {
        const Stuff::Common::Configuration local_cfg = sub_cfg.sub(DSC::toString(cc));
        if (local_cfg.has_key("domain") && local_cfg.has_key("value")) {
          auto domains = local_cfg.get< FieldMatrix< DomainFieldType, 2, 2 > >("domain");
          tmp_lower[0] = domains[0][0];
          tmp_lower[1] = domains[1][0];
          tmp_upper[0] = domains[0][1];
          tmp_upper[1] = domains[1][1];
          auto val = local_cfg.get< typename IndicatorFunctionType::RangeType >("value");
          values.emplace_back(tmp_lower, tmp_upper, val);
        } else
          break;
        ++cc;
      }
    }
    return values;
  } // ... get_values(...)

  static std::shared_ptr< FlatTopIndicatorType > create_indicator(const typename Values::value_type& value,
                                                                  const DomainType& channel_boundary_layer)
  {
    if (Stuff::Common::FloatCmp::eq(channel_boundary_layer, DomainType(0))) {
      return std::make_shared< IndicatorFunctionType >(Values(1, value), "channel");
    } else
      return std::make_shared< FlatTopFunctionType >(std::get< 0 >(value),
                                                     std::get< 1 >(value),
                                                     channel_boundary_layer,
                                                     std::get< 2 >(value),
                                                     "channel");
  } // ... create_indicator(...)
}; // class Spe10Model1< ..., 2, ... 1 >


} // namespace OS2015
} // namespace Problems
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_PROBLEMS_OS2015_HH
