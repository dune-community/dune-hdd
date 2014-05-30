// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_TESTCASES_BASE_HH
#define DUNE_HDD_LINEARELLIPTIC_TESTCASES_BASE_HH

#include <limits>

#include <dune/grid/io/file/dgfparser.hh>

#include <dune/stuff/grid/layers.hh>
#include <dune/stuff/grid/provider/default.hh>
#include <dune/stuff/functions/interfaces.hh>

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace TestCases {


/**
 *  The purpose of this class is to behave like a Stuff::Grid::ConstProviderInterface and at the same time to provide a
 *  means to obtain the real grid level corresponding to a refinement level.
 */
template< class GridType >
class Base
  : public Stuff::Grid::Providers::Default< GridType >
{
  typedef Stuff::Grid::Providers::Default< GridType > BaseType;
public:
  typedef typename GridType::template Codim< 0 >::Entity EntityType;
  typedef typename GridType::ctype                       DomainFieldType;
  static const unsigned int                              dimDomain = GridType::dimension;

  Base(std::shared_ptr< GridType > grd, size_t num_refinements)
    : BaseType(grd)
  {
    levels_.push_back(this->grid()->maxLevel());
    static const int refine_steps_for_half = DGFGridInfo< GridType >::refineStepsForHalf();
    for (size_t rr = 0; rr < num_refinements; ++rr) {
      this->grid()->globalRefine(refine_steps_for_half);
      levels_.push_back(this->grid()->maxLevel());
    }
    this->grid()->globalRefine(refine_steps_for_half);
    reference_level_ = this->grid()->maxLevel();
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

  std::shared_ptr< const typename BaseType::LevelGridViewType > reference_grid_view() const
  {
    return this->level_view(reference_level_);
  }

private:
  std::vector< int > levels_;
  int reference_level_;
}; // class Base


} // namespace TestCases
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_TESTCASES_BASE_HH
