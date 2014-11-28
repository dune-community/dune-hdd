// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_TESTCASES_BASE_HH
#define DUNE_HDD_LINEARELLIPTIC_TESTCASES_BASE_HH

#include <limits>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/stuff/common/disable_warnings.hh>
# include <dune/grid/io/file/dgfparser.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/stuff/grid/layers.hh>
#include <dune/stuff/grid/provider/default.hh>
#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/common/color.hh>
#include <dune/stuff/common/exceptions.hh>

#if HAVE_DUNE_GRID_MULTISCALE
# include <dune/grid/multiscale/provider/cube.hh>
#endif

#include <dune/pymor/parameters/base.hh>
#include <dune/pymor/common/exceptions.hh>

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace TestCases {
namespace internal {


class ParametricBase
  : public Pymor::Parametric
{
public:
  typedef std::map< std::string, Pymor::ParameterType > ParameterTypesMapType;
  typedef std::map< std::string, Pymor::Parameter >     ParametersMapType;

  static ParameterTypesMapType required_parameters()
  {
    return ParameterTypesMapType();
  }

  const ParametersMapType& parameters() const
  {
    return empty_parameters_map_;
  }

protected:
  static void check_parameters(const ParameterTypesMapType& required_types, const ParametersMapType& actual_parameters)
  {
    for (auto parameter_type : required_types) {
      const auto search_result = actual_parameters.find(parameter_type.first);
      if (search_result == actual_parameters.end())
        DUNE_THROW(Stuff::Exceptions::wrong_input_given,
                   "Missing parameter '" << parameter_type.first << "' in given actual_parameters!");
      if (search_result->second.type() != parameter_type.second)
        DUNE_THROW(Pymor::Exceptions::wrong_parameter_type,
                   "Given parameter '" << search_result->first << "' is a " << search_result->second.type()
                   << " and should be a " << parameter_type.second << "!");
    }
  } // ... check_paramets(...)

private:
  const ParametersMapType empty_parameters_map_;
}; // class ParametricBase


} // namespace internal


/**
 *  The purpose of this class is to behave like a Stuff::Grid::ConstProviderInterface and at the same time to provide a
 *  means to obtain the real grid level corresponding to a refinement level.
 */
template< class GridType >
class Base
  : public Stuff::Grid::Providers::Default< GridType >
  , public internal::ParametricBase
{
  typedef Stuff::Grid::Providers::Default< GridType > BaseType;
public:
  typedef typename GridType::template Codim< 0 >::Entity EntityType;
  typedef typename GridType::ctype                       DomainFieldType;
  static const unsigned int                              dimDomain = GridType::dimension;

  Base(std::shared_ptr< GridType > grd, size_t num_refinements)
    : BaseType(grd)
  {
    levels_.push_back(this->grid().maxLevel());
    static const int refine_steps_for_half = DGFGridInfo< GridType >::refineStepsForHalf();
    for (size_t rr = 0; rr < num_refinements; ++rr) {
      this->grid().globalRefine(refine_steps_for_half);
      levels_.push_back(this->grid().maxLevel());
    }
    this->grid().globalRefine(refine_steps_for_half);
    reference_level_ = this->grid().maxLevel();
  } // Base(...)

  size_t num_refinements() const
  {
    assert(levels_.size() > 0);
    return levels_.size() - 1;
  }

  int level_of(const size_t refinement) const
  {
    assert(refinement <= num_refinements());
    return levels_[refinement];
  }

  int reference_level() const
  {
    return reference_level_;
  }

  typename BaseType::LevelGridViewType reference_grid_view() const
  {
    return this->level_view(reference_level_);
  }

private:
  std::vector< int > levels_;
  int reference_level_;
}; // class Base


#if HAVE_DUNE_GRID_MULTISCALE


template< class GridImp >
class MultiscaleCubeBase
  : public internal::ParametricBase
{
public:
  typedef GridImp GridType;
  typedef typename GridType::template Codim< 0 >::Entity EntityType;
  typedef typename GridType::ctype                       DomainFieldType;
  static const unsigned int                              dimDomain = GridType::dimension;
  typedef FieldVector< DomainFieldType, dimDomain >      DomainType;

  typedef Stuff::Grid::Providers::Cube< GridType >      GridProviderType;
  typedef grid::Multiscale::Providers::Cube< GridType > MsGridProviderType;

  MultiscaleCubeBase(const Stuff::Common::Configuration& grid_cfg,
                     const int initial_refinements,
                     const size_t num_refinements,
                     const bool H_with_h = false)
    : partitions_(H_with_h
                  ? grid_cfg.get< std::string >("num_partitions") + "_H_with_h"
                  : grid_cfg.get< std::string >("num_partitions"))
  {
#ifndef DUNE_HDD_LINEARELLIPTIC_TESTCASES_BASE_DISABLE_WARNING
    std::cerr << Stuff::Common::Colors::red
              << "warning: running a multiscale testcase!\n"
              << "      => the boundaryinfo is set to AllDirichlet and all boundary values are set to zero!\n"
              << "         please manually check the testcase for compliance!\n"
              << Stuff::Common::StreamModifiers::normal << std::endl;
#endif // DUNE_HDD_LINEARELLIPTIC_TESTCASES_BASE_DISABLE_WARNING
    const auto lower_left = grid_cfg.get< DomainType >("lower_left");
    const auto upper_right = grid_cfg.get< DomainType >("upper_right");
    const auto num_elements = grid_cfg.get< std::vector< unsigned int > >("num_elements", dimDomain);
    const auto num_partitions = grid_cfg.get< std::vector< size_t > >("num_partitions", dimDomain);
    const auto num_oversampling_layers = grid_cfg.get("oversampling_layers", size_t(0));
    static const int refine_steps_for_half = DGFGridInfo< GridType >::refineStepsForHalf();
    for (size_t rr = 0; rr <= num_refinements; ++rr) {
      auto grid_ptr = GridProviderType(lower_left, upper_right, num_elements).grid_ptr();
      grid_ptr->globalRefine(boost::numeric_cast< int >(initial_refinements + rr*refine_steps_for_half));
      std::vector< size_t > actual_partitions = num_partitions;
      if (H_with_h)
        for (auto& element : actual_partitions)
          element *= std::pow(2, rr);
      level_providers_.emplace_back(new MsGridProviderType(grid_ptr,
                                                           lower_left,
                                                           upper_right,
                                                           actual_partitions,
                                                           num_oversampling_layers));
    }
    auto grid_ptr = GridProviderType(lower_left, upper_right, num_elements).grid_ptr();
    grid_ptr->globalRefine(boost::numeric_cast< int >(initial_refinements + (num_refinements + 1)*refine_steps_for_half));
    reference_provider_ = Stuff::Common::make_unique< MsGridProviderType >(grid_ptr,
                                                                           lower_left,
                                                                           upper_right,
                                                                           num_partitions,
                                                                           num_oversampling_layers);
  } // MultiscaleCubeBase(...)

  std::string partitioning() const
  {
    return partitions_;
  }

  size_t num_refinements() const
  {
    assert(level_providers_.size() > 0);
    return level_providers_.size() - 1;
  }

  const std::unique_ptr< MsGridProviderType >& level_provider(const size_t level) const
  {
    assert(level < level_providers_.size());
    return level_providers_[level];
  }

  const std::unique_ptr< MsGridProviderType >& reference_provider() const
  {
    return reference_provider_;
  }

private:
  const std::string partitions_;
  std::vector< std::unique_ptr< MsGridProviderType > > level_providers_;
  std::unique_ptr< MsGridProviderType > reference_provider_;
}; // class MultiscaleCubeBase


#endif // HAVE_DUNE_GRID_MULTISCALE

} // namespace TestCases
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_TESTCASES_BASE_HH
