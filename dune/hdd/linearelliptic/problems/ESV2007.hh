// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_PROBLEMS_ESV2007_HH
#define DUNE_HDD_LINEARELLIPTIC_PROBLEMS_ESV2007_HH

#include <memory>

#include <dune/common/static_assert.hh>

#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/ESV2007.hh>
#include <dune/stuff/common/memory.hh>

#include <dune/pymor/functions/default.hh>

#include "default.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Problems {


template< class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim = 1 >
class ESV2007
{
  static_assert(AlwaysFalse< EntityImp >::value, "Not available for dimRange > 1!");
};


template< class EntityImp, class DomainFieldImp, class RangeFieldImp >
class ESV2007< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 >
  : public Default< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 >
{
  typedef Default< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 > BaseType;
  typedef ESV2007< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 > ThisType;

  typedef Stuff::Functions::Constant< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 >    ScalarConstantFunctionType;
  typedef Stuff::Functions::Constant< EntityImp, DomainFieldImp, 2, RangeFieldImp, 2, 2 > ConstantFunctionType;
  typedef Stuff::Functions::ESV2007::Testcase1Force< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 > ForceType;
  typedef Pymor::Function::NonparametricDefault< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 >     ScalarFunctionType;
  typedef Pymor::Function::NonparametricDefault< EntityImp, DomainFieldImp, 2, RangeFieldImp, 2, 2 >  FunctionType;

public:
  static std::string static_id()
  {
    return BaseType::BaseType::static_id() + ".ESV2007";
  }

  static Stuff::Common::ConfigTree default_config(const std::string sub_name = "")
  {
    Stuff::Common::ConfigTree config("integration_order", "3");
    if (sub_name.empty())
      return config;
    else {
      Stuff::Common::ConfigTree tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  static std::unique_ptr< ThisType > create(const Stuff::Common::ConfigTree config = default_config(),
                                            const std::string sub_name = static_id())
  {
    const Stuff::Common::ConfigTree cfg = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    return Stuff::Common::make_unique< ThisType >(cfg.get("integration_order",
                                                          default_config().get< size_t >("integration_order")));
  } // ... create(...)

  ESV2007(const size_t integration_order = default_config().get< size_t >("integration_order"))
    : BaseType(std::make_shared< ScalarFunctionType >(new ScalarConstantFunctionType(1, "diffusion_factor")),
               std::make_shared< FunctionType >(new ConstantFunctionType(unit_matrix(), "diffusion_tensor")),
               std::make_shared< ScalarFunctionType >(new ForceType(integration_order,  "force")),
               std::make_shared< ScalarFunctionType >(new ScalarConstantFunctionType(0, "dirichlet")),
               std::make_shared< ScalarFunctionType >(new ScalarConstantFunctionType(0, "neumann")))
  {}

  virtual std::string type() const DS_OVERRIDE
  {
    return BaseType::BaseType::static_id() + ".ESV2007";
  }

private:
  typedef typename ConstantFunctionType::RangeType MatrixType;

  static MatrixType unit_matrix()
  {
    MatrixType matrix(RangeFieldImp(0));
    matrix[0][0] = RangeFieldImp(1);
    matrix[1][1] = RangeFieldImp(1);
    return matrix;
  }
}; // class ESV2007< ..., 1 >


} // namespace Problems
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_PROBLEMS_ESV2007_HH
