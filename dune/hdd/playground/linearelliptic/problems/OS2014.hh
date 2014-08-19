// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_PROBLEMS_OS2014_HH
#define DUNE_HDD_LINEARELLIPTIC_PROBLEMS_OS2014_HH

#include <memory>

#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/expression.hh>
#include <dune/stuff/functions/ESV2007.hh>

#include <dune/pymor/functions/default.hh>

#include "../../../linearelliptic/problems/default.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Problems {


template< class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim = 1 >
class OS2014
{
  static_assert(AlwaysFalse< EntityImp >::value, "Not available for dimRange > 1!");
};


template< class EntityImp, class DomainFieldImp, class RangeFieldImp >
class OS2014< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 >
  : public Default< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 >
{
  typedef Default< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 > BaseType;
  typedef OS2014< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 >  ThisType;

  typedef Stuff::Functions::Constant< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 >    ConstantScalarFunctionType;
  typedef Stuff::Functions::Constant< EntityImp, DomainFieldImp, 2, RangeFieldImp, 2, 2 > ConstantMatrixFunctionType;
  typedef Pymor::Function::NonparametricDefault< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 >
      ParametricScalarFunctionType;
  typedef Pymor::Function::NonparametricDefault< EntityImp, DomainFieldImp, 2, RangeFieldImp, 2, 2 >
      ParametricMatrixFunctionType;

  typedef Stuff::Functions::ESV2007::Testcase1Force< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 > ForceType;
  typedef Stuff::Functions::Expression< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 > ExpressionFunctionType;
  typedef Pymor::Function::AffinelyDecomposableDefault< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 >
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
    auto ret = std::make_shared< DefaultParametricFunctionType >(new ExpressionFunctionType(
                                                                   "x",
                                                                   "1+0.75*(sin(4*pi*(x[0]+0.5*x[1])))",
                                                                   integration_order,
                                                                   "affine_part"));
    ret->register_component(new ExpressionFunctionType("x",
                                                       "-0.75*(sin(4*pi*(x[0]+0.5*x[1])))",
                                                       integration_order,
                                                       "component_0"),
                            new Pymor::ParameterFunctional("mu", 1, "mu"));
    return ret;
  } // ... create_diffusion_factor(...)

public:
  static std::string static_id()
  {
    return BaseType::BaseType::static_id() + ".OS2014";
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
                                                          default_config().get< size_t >("integration_order")));
  } // ... create(...)

  OS2014(const size_t integration_order = default_config().get< size_t >("integration_order"))
    : BaseType(create_diffusion_factor(integration_order),
               std::make_shared< ParametricMatrixFunctionType >(new ConstantMatrixFunctionType(unit_matrix(),
                                                                                                  "diffusion_tensor")),
               std::make_shared< ParametricScalarFunctionType >(new ForceType(integration_order,  "force")),
               std::make_shared< ParametricScalarFunctionType >(new ConstantScalarFunctionType(0, "dirichlet")),
               std::make_shared< ParametricScalarFunctionType >(new ConstantScalarFunctionType(0, "neumann")))
  {}

  virtual std::string type() const DS_OVERRIDE
  {
    return BaseType::BaseType::static_id() + ".OS2014";
  }
}; // class OS2014< ..., 2, ..., 1 >


} // namespace Problems
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_PROBLEMS_OS2014_HH
