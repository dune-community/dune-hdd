// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_PROBLEMS_SPE10MODEL2_HH
#define DUNE_HDD_LINEARELLIPTIC_PROBLEMS_SPE10MODEL2_HH

#include <memory>

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/common/type_utils.hh>

#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/combined.hh>
#include <dune/stuff/playground/functions/indicator.hh>
#include <dune/stuff/functions/spe10model2.hh>

#include <dune/pymor/functions/default.hh>

#include "default.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Problems {
namespace Spe10 {


template< class E, class D, int d, class R, int r = 1 >
class Model2
  : public ProblemInterface< E, D, d, R, r >
{
  Model2() { static_assert(AlwaysFalse< E >::value, "Not available for dimDomain != 3!"); }
};


template< class EntityImp, class DomainFieldImp, class RangeFieldImp >
class Model2< EntityImp, DomainFieldImp, 3, RangeFieldImp, 1 >
  : public Default< EntityImp, DomainFieldImp, 3, RangeFieldImp, 1 >
{
  typedef Default< EntityImp, DomainFieldImp, 3, RangeFieldImp, 1 > BaseType;
  typedef Model2< EntityImp, DomainFieldImp, 3, RangeFieldImp, 1 > ThisType;
  typedef Stuff::Functions::Constant< EntityImp, DomainFieldImp, 3, RangeFieldImp, 1 >  ScalarConstantFunctionType;
  typedef Stuff::Functions::DomainIndicator< EntityImp, DomainFieldImp, 3, RangeFieldImp, 1 > ScalarIndicatorFunctionType;
  typedef typename ScalarConstantFunctionType::DomainType DomainType;
  typedef Stuff::Functions::Spe10::Model2< EntityImp, DomainFieldImp, 3, RangeFieldImp, 3, 3 > Spe10FunctionType;
  using typename BaseType::DiffusionFactorWrapperType;
  using typename BaseType::DiffusionTensorWrapperType;
  using typename BaseType::FunctionWrapperType;
public:
  using typename BaseType::DomainFieldType;
  using typename BaseType::DiffusionTensorType;
  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::BaseType::static_id() + ".spe10.model2";
  }

  static Stuff::Common::Configuration default_config(const std::string sub_name = "")
  {
    Stuff::Common::Configuration config;
    config["type"] = static_id();
    config["filename"] = Spe10FunctionType::default_config().template get< std::string >("filename");
    config["upper_right"] = Spe10FunctionType::default_config().template get< std::string >("upper_right");
    config["channel_width"] = "0";
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
    return Stuff::Common::make_unique< ThisType >(
          cfg.get("filename", default_config().template get< std::string >("filename")),
          cfg.get("upper_right", default_config().template get< DomainType >("upper_right")),
          cfg.get("channel_width", default_config().template get< DomainFieldType >("channel_width")));
  } // ... create(...)

  Model2(const std::string& filename, const DomainType& upper_right, const DomainFieldType& channel_width)
    : BaseType(std::make_shared<DiffusionFactorWrapperType>(std::make_shared< ScalarConstantFunctionType >(1, "diffusion_factor")),
               make_diffusion_tensor(filename, upper_right, channel_width),
               std::make_shared<FunctionWrapperType>(std::make_shared< ScalarConstantFunctionType >(0, "force")),
               std::make_shared<FunctionWrapperType>(std::make_shared< ScalarConstantFunctionType >(0, "dirichlet")),
               std::make_shared<FunctionWrapperType>(make_neumann(upper_right)))
  {}

  virtual std::string type() const override
  {
    return BaseType::BaseType::static_id() + ".spe10.model2";
  }

private:
  static std::shared_ptr< DiffusionTensorType > make_diffusion_tensor(const std::string& filename,
                                                                      const DomainType& upper_right,
                                                                      const DomainFieldType& channel_width)
  {
    typedef Pymor::Functions::AffinelyDecomposableDefault< EntityImp, DomainFieldType, 3, RangeFieldImp, 3, 3 >
        ParametricTensorType;
    auto diffusion_tensor = std::make_shared<ParametricTensorType>("diffusion_tensor");
    auto spe10 = std::shared_ptr<Spe10FunctionType>(new Spe10FunctionType(
                                                      filename, "spe10", {0., 0., 0.}, upper_right));
    diffusion_tensor->register_affine_part(spe10);
    if (channel_width > 0) {
      typedef Stuff::Functions::DomainIndicator< EntityImp, DomainFieldImp, 3, RangeFieldImp, 3, 3 > TensorIndicatorFunctionType;
      auto one = std::make_shared<ScalarConstantFunctionType>(1, "one");

      auto channel_x = std::shared_ptr<ScalarIndicatorFunctionType>(new ScalarIndicatorFunctionType(
          {{{{0.,             upper_right[1]/2. - channel_width/2., 0.},
             {upper_right[0], upper_right[1]/2. + channel_width/2., upper_right[2]/2.}},
            1.}}));
      typename TensorIndicatorFunctionType::RangeType channel_x_scaled_value(0.);
      channel_x_scaled_value[0][0]
            = channel_x_scaled_value[1][1]
            = channel_x_scaled_value[2][2]
            = Spe10FunctionType::default_config().template get< RangeFieldImp >("max");
      auto channel_x_scaled = std::shared_ptr<TensorIndicatorFunctionType>(new TensorIndicatorFunctionType(
          {{{{0.,             upper_right[1]/2. - channel_width/2., 0.},
             {upper_right[0], upper_right[1]/2. + channel_width/2., upper_right[2]/2.}},
            channel_x_scaled_value}}));
      auto channel_x_remover = Stuff::Functions::make_difference(one, channel_x, "channel_x_remover");
      auto spe10_wo_channel_x = Stuff::Functions::make_product(channel_x_remover, spe10, "spe10_wo_channel_x");
      auto spe10_w_scaled_channel_x = Stuff::Functions::make_sum(spe10_wo_channel_x,
                                                                 channel_x_scaled,
                                                                 "spe10_w_scaled_channel_x");
      diffusion_tensor->register_component(spe10_w_scaled_channel_x,
                                           new Pymor::ParameterFunctional("channel", 2, "channel[0]"));

      auto channel_y = std::shared_ptr<ScalarIndicatorFunctionType>(new ScalarIndicatorFunctionType(
          {{{{upper_right[0]/2. - channel_width/2., 0.,             upper_right[2]/2.},
             {upper_right[0]/2. + channel_width/2., upper_right[1], upper_right[2]}},
            1.}}));
      typename TensorIndicatorFunctionType::RangeType channel_y_scaled_value(0.);
      channel_y_scaled_value[0][0]
            = channel_y_scaled_value[1][1]
            = channel_y_scaled_value[2][2]
            = Spe10FunctionType::default_config().template get< RangeFieldImp >("max");
      auto channel_y_scaled = std::shared_ptr<TensorIndicatorFunctionType>(new TensorIndicatorFunctionType(
          {{{{upper_right[0]/2. - channel_width/2., 0.,             upper_right[2]/2.},
             {upper_right[0]/2. + channel_width/2., upper_right[1], upper_right[2]}},
            channel_y_scaled_value}}));
      auto channel_y_remover = Stuff::Functions::make_difference(one, channel_y, "channel_y_remover");
      auto spe10_wo_channel_y = Stuff::Functions::make_product(channel_y_remover, spe10, "spe10_wo_channel_y");
      auto spe10_w_scaled_channel_y = Stuff::Functions::make_sum(spe10_wo_channel_y,
                                                                 channel_y_scaled,
                                                                 "spe10_w_scaled_channel_y");
      diffusion_tensor->register_component(spe10_w_scaled_channel_y,
                                           new Pymor::ParameterFunctional("channel", 2, "channel[1]"));
    }
    return diffusion_tensor;
  } // ... make_diffusion_tensor(...)

  static std::shared_ptr< ScalarIndicatorFunctionType > make_neumann(const DomainType& upper_right)
  {
    return std::shared_ptr<ScalarIndicatorFunctionType>(new ScalarIndicatorFunctionType(
                                                          {{{{0.,             0.,   0.},
                                                             {upper_right[0], upper_right[1]/440., upper_right[2]}},
                                                            1.}},
                                                          "neumann"));
  }
}; // class Model2< ..., 3, 3 >


} // namespace Spe10
} // namespace Problems
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_PROBLEMS_SPE10MODEL2_HH
