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

class MpiCGExample
{
public:
  typedef Dune::HDD::LinearElliptic::Discretizations::MpiCG
       DiscretizationType;
  typedef double RangeFieldType;
  typedef typename Dune::SPGrid<RangeFieldType, 2> GridType;
  typedef Dune::Stuff::Grid::Providers::Cube< GridType > GridProviderType;
  static constexpr unsigned int dimRange = 1;
  typedef typename DSG::Entity<typename GridType::LeafGridView>::Type EntityType;

//  typedef Dune::HDD::LinearElliptic::Problems::ESV2007
//      < EntityType, RangeFieldType, dimRange, RangeFieldType, dimRange> ProblemType;
//  typedef Dune::HDD::LinearElliptic::TestCases::ESV2007
//      < GridType > TestcaseType;

  typedef Dune::HDD::LinearElliptic::Problems::Spe10::Model1
      < EntityType, RangeFieldType, dimRange, RangeFieldType, dimRange> ProblemType;
  typedef Dune::HDD::LinearElliptic::TestCases::Spe10::ParametricModel1
      < GridType > TestcaseType;


public:
  MpiCGExample(const std::size_t num_refinements = 0,
              const DUNE_STUFF_SSIZE_T info_log_levels = 0,
              const DUNE_STUFF_SSIZE_T debug_log_levels = -1,
              const bool enable_warnings = true,
              const bool enable_colors = true,
               const std::string info_color  = DSC::TimedLogging::default_info_color(),
               const std::string debug_color = DSC::TimedLogging::default_debug_color(),
               const std::string warn_color  = DSC::TimedLogging::default_warning_color())
    : testcase_( {{"mu",     Dune::Pymor::Parameter("mu", 1)},
                  {"mu_hat", Dune::Pymor::Parameter("mu", 1)},
                  {"mu_bar", Dune::Pymor::Parameter("mu", 1)},
                  {"parameter_range_min", Dune::Pymor::Parameter("mu", 0.1)},
                  {"parameter_range_max", Dune::Pymor::Parameter("mu", 1.0)}}, num_refinements)
    , discretization_(testcase_,
                      testcase_.boundary_info(),
                      testcase_.problem())
  {
    DSC::TimedLogger().create(info_log_levels,
                                debug_log_levels,
                                enable_warnings,
                                enable_colors,
                                info_color,
                                debug_color,
                                warn_color);

    DSC::TimedLogger().get("cg.thermalblock.example").info() << "creating grid and problem... " << std::endl;
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


#endif // CG_MPI_HH
