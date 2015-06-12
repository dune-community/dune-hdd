#ifndef DUNE_HDD_EXAMPLES_LINEARELLIPTIC_CG_MPI_HH
#define DUNE_HDD_EXAMPLES_LINEARELLIPTIC_CG_MPI_HH

#include <config.h> //mandatory for python bindings

#include <memory>
#include <dune/grid/spgrid.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/hdd/linearelliptic/discreteproblem.hh>
#include <dune/hdd/linearelliptic/discretizations/mpi_cg.hh>
#include <dune/hdd/linearelliptic/problems/ESV2007.hh>
#include <dune/hdd/linearelliptic/testcases/ESV2007.hh>
#include <dune/hdd/linearelliptic/problems/spe10model2.hh>
#include <dune/hdd/linearelliptic/testcases/spe10model2.hh>
#include <dune/hdd/linearelliptic/problems/spe10.hh>
#include <dune/hdd/linearelliptic/testcases/spe10.hh>
#include <dune/hdd/linearelliptic/problems/thermalblock.hh>
#include <dune/hdd/linearelliptic/testcases/thermalblock.hh>
#include <dune/hdd/linearelliptic/problems/random_block_problem.hh>
#include <dune/hdd/linearelliptic/testcases/random_block_testcase.hh>

class MpiCGExampleThermal
{
public:
  static constexpr int dimDomain = 2;
  typedef double RangeFieldType;
  typedef Dune::HDD::LinearElliptic::Discretizations::MpiCG<dimDomain>
       DiscretizationType;
  typedef Dune::SPGrid<RangeFieldType, dimDomain> GridType;
  typedef Dune::Stuff::Grid::Providers::Cube< GridType > GridProviderType;
  static constexpr unsigned int dimRange = 1;
  typedef typename DSG::Entity<typename GridType::LeafGridView>::Type EntityType;

  typedef Dune::HDD::LinearElliptic::TestCases::Thermalblock
      < GridType > TestcaseType;
  typedef Dune::HDD::LinearElliptic::TestCases::Spe10::ParametricModel1
      < GridType > ProblemType;

public:
  MpiCGExampleThermal(const std::size_t num_refinements,
               const unsigned int overlap_size,
               DSC::FieldVector< size_t, dimDomain > num_blocks
               = DSC::ValueInitFieldVector< size_t, dimDomain, 2u >())
    : testcase_( num_refinements, num_blocks, overlap_size)
    , discretization_(testcase_,
                      testcase_.boundary_info(),
                      testcase_.problem())
  {
    discretization_.init(DSC_LOG_DEBUG_0);
  } // MpiCGExample(...)

  static std::string static_id()
  {
    return "mpi_cg";
  }

  const DiscretizationType& discretization() const
  {
    return discretization_;
  }

  DiscretizationType* discretization_and_return_ptr() const
  {
    return new DiscretizationType(discretization_);
  }

private:
  TestcaseType testcase_;
  DiscretizationType discretization_;
}; // class LinearellipticExampleCG

class MpiCGExampleRandom
{
public:
  static constexpr int dimDomain = 2;
  typedef double RangeFieldType;
  typedef Dune::HDD::LinearElliptic::Discretizations::MpiCG<dimDomain>
       DiscretizationType;
  typedef Dune::SPGrid<RangeFieldType, dimDomain> GridType;
  typedef Dune::Stuff::Grid::Providers::Cube< GridType > GridProviderType;
  static constexpr unsigned int dimRange = 1;
  typedef typename DSG::Entity<typename GridType::LeafGridView>::Type EntityType;

  typedef Dune::HDD::LinearElliptic::TestCases::RandomBlockTestcase
      < GridType > TestcaseType;
  typedef Dune::HDD::LinearElliptic::Problems::RandomBlockProblem
      < EntityType, RangeFieldType, dimRange, RangeFieldType, dimRange> ProblemType;

public:
  MpiCGExampleRandom(const std::size_t num_refinements,
               const unsigned int overlap_size,
               DSC::FieldVector< size_t, dimDomain > num_blocks
               = DSC::ValueInitFieldVector< size_t, dimDomain, 2u >())
    : testcase_( num_refinements, num_blocks, overlap_size)
    , discretization_(testcase_,
                      testcase_.boundary_info(),
                      testcase_.problem())
  {
    discretization_.init(DSC_LOG_DEBUG_0);
  } // MpiCGExample(...)

  static std::string static_id()
  {
    return "mpi_cg";
  }

  const DiscretizationType& discretization() const
  {
    return discretization_;
  }

  DiscretizationType* discretization_and_return_ptr() const
  {
    return new DiscretizationType(discretization_);
  }

private:
  TestcaseType testcase_;
  DiscretizationType discretization_;
}; // class LinearellipticExampleCG

//typedef MpiCGExampleRandom MpiCGExample;
typedef MpiCGExampleThermal MpiCGExample;

#endif // CG_MPI_HH
