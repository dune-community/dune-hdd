// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_PROBLEMS_RANDOM_BLOCK_HH
#define DUNE_HDD_LINEARELLIPTIC_PROBLEMS_RANDOM_BLOCK_HH

#include <memory>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/timer.hh>
#include <dune/common/typetraits.hh>

#if HAVE_DUNE_GRID_MULTISCALE
# include <dune/grid/multiscale/provider/cube.hh>
#endif

#include <dune/stuff/common/string.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/expression.hh>
#include <dune/stuff/playground/functions/indicator.hh>
#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/grid/boundaryinfo.hh>

#include <dune/pymor/functions/default.hh>
#include <dune/pymor/functions/complement_pair_function.hh>

#include "default.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Problems {


template< class E, class D, int d, class R, int r = 1 >
class RandomBlockProblem
  : public ProblemInterface< E, D, d, R, r >
{
  RandomBlockProblem() { static_assert(AlwaysFalse< E >::value, "Not available for these dimensions!"); }
};


template< class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp >
class RandomBlockProblem< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 >
  : public Default< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 >
{
  typedef Default< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 > BaseType;
  typedef RandomBlockProblem< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 > ThisType;

public:
  typedef Pymor::Functions::ComplementPairPymorFunction< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 >
      ComplementPairPymorFunction;
  using typename BaseType::DiffusionTensorType;
  using typename BaseType::FunctionType;

  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::BaseType::static_id() + ".RandomBlockProblem";
  }

  static Stuff::Common::Configuration default_config(const std::string sub_name = "")
  {
    typedef Stuff::Functions::Constant< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1 >
        ConstantFunctionType;
    typedef Stuff::Functions::Constant< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, domainDim, domainDim >
        ConstantMatrixFunctionType;
    Stuff::Common::Configuration config;
    Stuff::Common::Configuration randomblock_config = ComplementPairPymorFunction::default_config();
    randomblock_config["name"] = "diffusion_factor";
    randomblock_config["type"] = ComplementPairPymorFunction::static_id();
    randomblock_config["parameter_name"] = "diffusion";
    randomblock_config["ellipsoids.count"] = "10";
    randomblock_config["ellipsoids.min_radius"] = "0.05";
    randomblock_config["ellipsoids.max_radius"] = "0.1";
    randomblock_config["ellipsoids.seed"] = "0";
    randomblock_config["ellipsoids.children"] = "1";
    randomblock_config["ellipsoids.recursion_depth"] = "1";
    randomblock_config["ellipsoids.recursion_scale"] = "0.5";
    config.add(randomblock_config, "diffusion_factor");
    Stuff::Common::Configuration diffusion_tensor_config = ConstantMatrixFunctionType::default_config();
    diffusion_tensor_config["name"] = "diffusion_tensor";
    diffusion_tensor_config["type"] = ConstantMatrixFunctionType::static_id();
    config.add(diffusion_tensor_config, "diffusion_tensor");
    Stuff::Common::Configuration constant_config = ConstantFunctionType::default_config();
    constant_config["type"] = ConstantFunctionType::static_id();
    constant_config["name"] = "force";
    constant_config["value"] = "1";
    config.add(constant_config, "force");
    constant_config["name"] = "dirichlet";
    constant_config["value"] = "0";
    config.add(constant_config, "dirichlet");
    constant_config["name"] = "neumann";
    config.add(constant_config, "neumann");
    if (sub_name.empty())
      return config;
    else {
      Stuff::Common::Configuration tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  static std::unique_ptr< ThisType > create(const Stuff::Common::Configuration config)
  {
    const Stuff::Common::Configuration cfg = config;
    std::shared_ptr< ComplementPairPymorFunction >
        randomblock_function(ComplementPairPymorFunction::create(cfg.sub("diffusion_factor")));
    return Stuff::Common::make_unique< ThisType >(randomblock_function,
                                                  BaseType::create_matrix_function("diffusion_tensor", cfg),
                                                  BaseType::create_vector_function("force", cfg),
                                                  BaseType::create_vector_function("dirichlet", cfg),
                                                  BaseType::create_vector_function("neumann", cfg));
  } // ... create(...)

  RandomBlockProblem(const std::shared_ptr< const ComplementPairPymorFunction >& checkerboard_function,
               const std::shared_ptr< const DiffusionTensorType >& diffusion_tensor,
               const std::shared_ptr< const FunctionType >& force,
               const std::shared_ptr< const FunctionType >& dirichlet,
               const std::shared_ptr< const FunctionType >& neumann)
    : BaseType(checkerboard_function, diffusion_tensor, force, dirichlet, neumann)
  {}

  virtual std::string type() const override
  {
    return BaseType::BaseType::static_id() + ".RandomBlockProblem";
  }
}; // class RandomBlockProblem< ..., 1 >

} // namespace Problems
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_PROBLEMS_RANDOM_BLOCK_HH
