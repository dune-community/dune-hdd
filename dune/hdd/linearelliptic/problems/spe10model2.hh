// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_PROBLEMS_Spe10Model2_HH
#define DUNE_HDD_LINEARELLIPTIC_PROBLEMS_Spe10Model2_HH

#include <memory>

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/global.hh>
#include <dune/stuff/functions/spe10model2.hh>

#include <dune/pymor/functions/default.hh>

#include "default.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Problems {


template< class E, class D, int d, class R, int r = 1 >
class Spe10Model2
  : public ProblemInterface< E, D, d, R, r >
{
  static_assert(AlwaysFalse< E >::value, "Not available for dimDomain != 3!");
};


template< class EntityImp, class DomainFieldImp, class RangeFieldImp >
class Spe10Model2< EntityImp, DomainFieldImp, 3, RangeFieldImp, 1 >
  : public Default< EntityImp, DomainFieldImp, 3, RangeFieldImp, 1 >
{
  typedef Default< EntityImp, DomainFieldImp, 3, RangeFieldImp, 1 > BaseType;
  typedef Spe10Model2< EntityImp, DomainFieldImp, 3, RangeFieldImp, 1 > ThisType;
public:
  typedef Stuff::Functions::Constant< EntityImp, DomainFieldImp, 3, RangeFieldImp, 1 >    ScalarConstantFunctionType;
  typedef Stuff::Functions::Constant< EntityImp, DomainFieldImp, 3, RangeFieldImp, 3, 3 > ConstantFunctionType;
  typedef ScalarConstantFunctionType ForceType;
  typedef Pymor::Functions::NonparametricDefault< EntityImp, DomainFieldImp, 3, RangeFieldImp, 1 >     ScalarFunctionType;
  typedef Pymor::Functions::NonparametricDefault< EntityImp, DomainFieldImp, 3, RangeFieldImp, 3, 3 >  MatrixFunctionType;

  typedef Stuff::Functions::Spe10::Model2< EntityImp, DomainFieldImp, 3, RangeFieldImp, 3, 3> Spe10FunctionType;
  typedef Stuff::GlobalLambdaFunction< EntityImp, DomainFieldImp, 3, RangeFieldImp, 1> NeumannFunctionType;
public:
  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::BaseType::static_id() + ".Spe10Model2";
  }

  static Stuff::Common::Configuration default_config(const std::string sub_name = "")
  {
    Stuff::Common::Configuration config("integration_order", "3");
    config["lower_left"] = "[0 0 0]";
    config.set("upper_right", Spe10FunctionType::default_upper_right);
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
                                                          default_config().get< size_t >("integration_order")));
  } // ... create(...)

  Spe10Model2(const size_t integration_order = default_config().get< size_t >("integration_order"))
    : BaseType(std::make_shared< ScalarFunctionType >(new ScalarConstantFunctionType(1, "diffusion_factor")),
               std::make_shared< MatrixFunctionType >(new Spe10FunctionType()),
               std::make_shared< ScalarFunctionType >(new ScalarConstantFunctionType(0, "force")),
               std::make_shared< ScalarFunctionType >(new ScalarConstantFunctionType(0, "dirichlet")),
               std::make_shared< ScalarFunctionType >(make_neumann()))
  {}

  virtual std::string type() const override
  {
    return BaseType::BaseType::static_id() + ".Spe10Model2";
  }

  static NeumannFunctionType* make_neumann() {
    typedef typename std::remove_const<decltype(Spe10FunctionType::default_upper_right)>::type UR;
    typedef typename NeumannFunctionType::DomainType DomainType;
    auto lmb = [](DomainType x) { return default_config().get<UR>("upper_right")[2] == x[2];};
    return new NeumannFunctionType(lmb, 0 /*order = constant on element*/, "neumann");
  }

private:
  typedef typename ConstantFunctionType::RangeType MatrixType;

  static MatrixType unit_matrix()
  {
    MatrixType matrix(RangeFieldImp(0));
    matrix[0][0] = RangeFieldImp(1);
    matrix[1][1] = RangeFieldImp(1);
    matrix[2][2] = RangeFieldImp(1);
    return matrix;
  }

}; // class Spe10Model2< ..., 1 >


} // namespace Problems
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_PROBLEMS_Spe10Model2_HH
