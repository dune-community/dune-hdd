// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#if HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
# define ENABLE_ALUGRID 1
# include <dune/grid/alugrid.hh>
#else
# error This example requires alugrid!
#endif

#include <dune/common/timer.hh>

#include <dune/grid/io/file/dgfparser.hh>

#include <dune/grid/multiscale/provider/cube.hh>

#include <dune/stuff/functions/spe10.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/common/float_cmp.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/stuff/common/fvector.hh>
#include <dune/stuff/grid/boundaryinfo.hh>

#include <dune/hdd/linearelliptic/problems/default.hh>
#include <dune/hdd/linearelliptic/discretizations/cg.hh>

using namespace Dune;
using namespace HDD;
using namespace GDT;


template< class E, class D, int d, class R >
class IndicatorFunction
  : public Stuff::LocalizableFunctionInterface< E, D, d, R, 1 >
{
  typedef Stuff::LocalizableFunctionInterface< E, D, d, R, 1 > BaseType;
  typedef IndicatorFunction< E, D, d, R > ThisType;

  class Localfunction
    : public Stuff::LocalfunctionInterface< E, D, d, R, 1 >
  {
    typedef Stuff::LocalfunctionInterface< E, D, d, R, 1 > InterfaceType;
  public:
    using typename InterfaceType::EntityType;
    using typename InterfaceType::DomainType;
    using typename InterfaceType::RangeType;
    using typename InterfaceType::JacobianRangeType;

    Localfunction(const EntityType& entity, const RangeType& value)
      : InterfaceType(entity)
      , value_(value)
    {}

    virtual size_t order() const DS_OVERRIDE DS_FINAL
    {
      return 0;
    }

    virtual void evaluate(const DomainType& xx, RangeType& ret) const DS_OVERRIDE DS_FINAL
    {
      assert(this->is_a_valid_point(xx));
      ret = value_;
    }

    virtual void jacobian(const DomainType& xx, JacobianRangeType& ret) const DS_OVERRIDE DS_FINAL
    {
      assert(this->is_a_valid_point(xx));
      ret *= 0.0;
    }

  private:
    const RangeType value_;
  }; // class Localfunction

public:
  using typename BaseType::EntityType;
  typedef Stuff::Common::FieldVector< D, d > DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::LocalfunctionType;

  IndicatorFunction(std::vector< std::pair< std::pair< DomainType, DomainType >, R > > values,
                    const std::string name = "indicator")
    : values_(values)
    , name_(name)
  {}

  virtual ~IndicatorFunction() {}

  virtual ThisType* copy() const DS_OVERRIDE DS_FINAL
  {
    DUNE_THROW(NotImplemented, "");
  }

  virtual std::string name() const DS_OVERRIDE DS_FINAL
  {
    return name_;
  }

  virtual std::unique_ptr< LocalfunctionType > local_function(const EntityType& entity) const DS_OVERRIDE DS_FINAL
  {
    const auto center = entity.geometry().center();
    for (const auto& element : values_) {
      const auto& domain = element.first;
      const auto& lower_left = domain.first;
      const auto& upper_right = domain.second;
      if (Stuff::Common::FloatCmp::lt(lower_left, center) && Stuff::Common::FloatCmp::lt(center, upper_right)) {
        const auto& value = element.second;
        return Stuff::Common::make_unique< Localfunction >(entity, value);
      }
    }
    return Stuff::Common::make_unique< Localfunction >(entity, 0.0);
  } // ... local_function(...)

private:
  const std::vector< std::pair< std::pair< DomainType, DomainType >, R > > values_;
  const std::string name_;
}; // class IndicatorFunction


int main(int /*argc*/, char** /*argv*/)
{
  try {
    typedef ALUGrid< 2, 2, simplex, conforming > GridType;
    typedef GridType::Codim< 0 >::Entity         EntityType;
    typedef GridType::ctype   DomainFieldType;
    static const unsigned int dimDomain = GridType::dimension;
    typedef double            RangeFieldType;
    static const unsigned int dimRange = 1;

    typedef Stuff::Functions::Spe10Model1< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >
        Spe10Model1FunctionType;
    Stuff::Common::ConfigTree config = Spe10Model1FunctionType::default_config();
    config["upper_right"] = "[5 1]";
    config["num_elements"] = "[100 20]";
    config["num_partitions"] = "1";
    std::shared_ptr< Spe10Model1FunctionType > spe10_model1_function(Spe10Model1FunctionType::create(config));
    typedef Stuff::Functions::Constant< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >
        ConstantFunctionType;
    auto dirichlet = std::make_shared< ConstantFunctionType >(0.0, "dirichlet");
    auto neumann = std::make_shared< ConstantFunctionType >(0.0, "neumann");
    typedef Stuff::Functions::Constant< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimDomain, dimDomain >
        ConstantMatrixFunctionType;
    typename ConstantMatrixFunctionType::RangeType unit_matrix(0.0);
    for (size_t dd = 0; dd < dimDomain; ++dd)
      unit_matrix[dd][dd] = 1.0;
    auto unit_matrix_function = std::make_shared< ConstantMatrixFunctionType >(unit_matrix, "diffusion tensor");
    typedef IndicatorFunction< EntityType, DomainFieldType, dimDomain, RangeFieldType > IndicatorFunctionType;
    const RangeFieldType scale = 100.0;
    auto force = std::shared_ptr< IndicatorFunctionType >(
                   new IndicatorFunctionType({{{{0.70, 0.60}, {0.90, 0.80}},  1.0 * scale},
                                              {{{1.00, 0.30}, {1.20, 0.50}},  1.0 * scale},
                                              {{{4.15, 0.30}, {4.35, 0.50}}, -1.0 * scale},
                                              {{{3.30, 0.45}, {3.50, 0.65}}, -1.0 * scale}}, "force"));
    typedef LinearElliptic::Problems::Default< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >
        ProblemType;
    const ProblemType problem(spe10_model1_function, unit_matrix_function, force, dirichlet, neumann);

    typedef Stuff::Grid::Providers::Cube< GridType > GridProvider;
    auto grid_provider = GridProvider::create(config);
    grid_provider->grid()->globalRefine(1);
    const int level = grid_provider->grid()->maxLevel();
    grid_provider->grid()->globalRefine(2 * DGFGridInfo< GridType >::refineStepsForHalf());
//    const int reference_level = grid_provider->grid()->maxLevel();
    grid_provider->visualize("grid");

    problem.visualize(*grid_provider->leaf_view(), "problem");

    std::cout << "solving... " << std::flush;
    Dune::Timer timer;
    LinearElliptic::Discretizations::CG< GridType, Stuff::Grid::ChooseLayer::level, RangeFieldType, dimRange, 1 >
        discretization(*grid_provider,
                       Stuff::Grid::BoundaryInfoConfigs::AllDirichlet::default_config(),
                       problem,
                       level);
    discretization.init();
    auto solution = discretization.create_vector();
    discretization.solve(solution);
    std::cout << "done (took " << timer.elapsed() << "s)" << std::endl;
    discretization.visualize(solution, "solution", "solution");

    // if we came that far we can as well be happy about it
    return 0;
  } catch (Dune::Exception& e) {
    std::cerr << "\ndune reported error: " << e.what() << std::endl;
    std::abort();
  } catch (std::exception& e) {
    std::cerr << "\n" << e.what() << std::endl;
    std::abort();
  } catch (...) {
    std::cerr << "\nUnknown exception thrown!" << std::endl;
    std::abort();
  } // try
} // ... main(...)
