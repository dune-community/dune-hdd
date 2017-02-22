// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#define DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS 1
#include <config.h>
#include <sstream>

#include <dune/grid/alugrid.hh>
#include <dune/stuff/common/disable_warnings.hh>
#include <dune/stuff/common/reenable_warnings.hh>
#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/filesystem.hh>
#include <dune/stuff/common/float_cmp.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/print.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/test/common.hh>
#include <dune/grid/spgrid.hh>
#include <dune/hdd/linearelliptic/testcases/OS2014.hh>

#include "linearelliptic-block-swipdg.hh"
#include "linearelliptic-swipdg.hh"

#include "benchmark_util.hh"


using namespace Dune;
using namespace HDD;

// typedef ALUGrid< 2, 2, simplex, conforming, MPI_Comm > GridType;
// typedef YaspGrid< 2 > GridType;


typedef SPGrid<double, 2, SPIsotropicRefinement, MPI_Comm> StudyGridType;

class OS2014_FVCA7__estimator_study {

  typedef LinearElliptic::TestCases::OS2014<StudyGridType> TestCaseType;
  typedef LinearElliptic::TestCases::OS2014Multiscale<StudyGridType>
      BlockTestCaseType;

  static const GDT::ChooseSpaceBackend space_backend =
      GDT::ChooseSpaceBackend::pdelab;
  static const Stuff::LA::ChooseBackend la_backend =
      Stuff::LA::ChooseBackend::istl_sparse;

  typedef LinearElliptic::Tests::SWIPDGStudy<TestCaseType, 1, space_backend,
                                             la_backend>
      EocStudyType;
  typedef LinearElliptic::Tests::BlockSWIPDGStudy<BlockTestCaseType, 1,
                                                  la_backend>
      BlockEocStudyType;

public:
  static void BlockSWIPDG_coarse_triangulation(const std::string partitioning) {
    BlockTestCaseType test_case(partitioning);
    BlockEocStudyType eoc_study(test_case, {"energy", "eta_OS2014"});
    Stuff::Test::check_eoc_study_for_success(eoc_study,
                                             eoc_study.run_eoc(DSC_LOG_INFO));
  } // ... BlockSWIPDG_coarse_triangulation(...)

  static void ESV2007_fine_triangulation() {
    TestCaseType test_case(DSC_CONFIG_GET("eoc_run.refinements", 4), DSC_CONFIG_GET("eoc_run.start_cell_count", 16));
    test_case.print_header(DSC_LOG_INFO_0);
    DSC_LOG_INFO_0 << std::endl;
    EocStudyType eoc_study(test_case, {"L2", "energy", "eta_ESV2007"}, std::vector<std::string>(), "OS2014");
    eoc_study.run_eoc(DSC_LOG_INFO_0);
  }
}; // class OS2014_FVCA7__estimator_study

void single_run_swipdg()
{
  static constexpr int polOrder{1};
  static constexpr GDT::ChooseSpaceBackend space_backend{GDT::ChooseSpaceBackend::pdelab};
  static constexpr Stuff::LA::ChooseBackend la_backend{Stuff::LA::ChooseBackend::istl_sparse};
  using GridType = SPGrid<double, 2, SPIsotropicRefinement, MPI_Comm>;
  using TestCaseType = LinearElliptic::TestCases::OS2014<GridType>;
  using DiscretizationType = typename LinearElliptic::Tests::internal::DiscretizationSWIPDG< TestCaseType, polOrder, space_backend, la_backend >::Type;

  TestCaseType test_case(0, DSC_CONFIG_GET("single_run.cells_per_dim", 16));
  {
    DSC::ScopedTiming timing("single_run.swipdg.solve");
    DiscretizationType current_discretization(test_case, test_case.boundary_info(), test_case.problem(),
                                              test_case.level_of(test_case.reference_level()));
    current_discretization.init();
    auto current_solution_vector_on_level_ = current_discretization.create_vector();
    const auto solver_options = DSC_CONFIG.sub("global_solver");
    current_discretization.solve(solver_options, current_solution_vector_on_level_);
  }
}

void single_run_blockswipdg()
{
  static constexpr int polOrder{1};
  static constexpr GDT::ChooseSpaceBackend space_backend{GDT::ChooseSpaceBackend::pdelab};
  static constexpr Stuff::LA::ChooseBackend la_backend{Stuff::LA::ChooseBackend::istl_sparse};
  using GridType = SPGrid<double, 2, SPIsotropicRefinement, MPI_Comm>;
  using TestCaseType = LinearElliptic::TestCases::OS2014Multiscale<GridType>;
  using DiscretizationType = typename LinearElliptic::Tests::internal::DiscretizationBlockSWIPDG< TestCaseType, polOrder, la_backend >::Type;

  TestCaseType test_case(DSC_CONFIG_GET("lrbms.partitioning", "[1 1 1]"), DSC_CONFIG_GET("lrbms.refinements", 0));
  {
    DSC::ScopedTiming timing("single_run.lrbms.solve");
    DiscretizationType current_discretization(*test_case.reference_provider(), test_case.boundary_info(), test_case.problem());
    auto grid_view = test_case.reference_provider()->grid().leafGridView();

    current_discretization.init();
    auto current_solution_vector_on_level_ = current_discretization.create_vector();
    const auto solver_options = DSC_CONFIG.sub("global_solver");
    current_discretization.solve(solver_options, current_solution_vector_on_level_);
  }
}

int main(int argc, char **argv) {
  try {
    DSC_CONFIG.read_command_line(argc, argv);
#if HAVE_DUNE_FEM
    Dune::Fem::MPIManager::initialize(argc, argv);
#else
    Dune::MPIHelper::instance(argc, argv);
#endif

    DSC::testCreateDirectory(DSC_CONFIG_GET("global.datadir", "data/"));

    // LOG_NONE = 1, LOG_ERROR = 2, LOG_INFO = 4,LOG_DEBUG = 8,LOG_CONSOLE =
    // 16,LOG_FILE = 32
    // --> LOG_ERROR | LOG_INFO | LOG_DEBUG | LOG_CONSOLE | LOG_FILE = 62
    DSC::Logger().create(
        DSC_CONFIG_GET("logging.level", 62),
        DSC_CONFIG_GET("logging.file", std::string(argv[0]) + ".log"),
        DSC_CONFIG_GET("global.datadir", "data"),
        DSC_CONFIG_GET("logging.dir", "log" /*path below datadir*/));
    DSC_PROFILER.set_outputdir(DSC_CONFIG_GET("global.datadir", "data"));

    DSC_LOG_INFO_0 << DSC_CONFIG.report_string() << std::endl;

    DSC::TimedLogger().create(-1, -1);
    DS::threadManager().set_max_threads(1);
    {
      DSC::OutputScopedTiming outs("all", DSC_LOG_INFO_0);
      DSC_CONFIG.set("grids.total_macro_cells", 256);
  //    single_run_swipdg();
  //    single_run_blockswipdg();
      OS2014_FVCA7__estimator_study::BlockSWIPDG_coarse_triangulation(DSC_CONFIG_GET("lrbms.partitioning", "[1 1 1]"));
  //    OS2014_FVCA7__estimator_study::ESV2007_fine_triangulation();
    }
    DSC_PROFILER.output_per_rank("profiler");
    DS::mem_usage();
    DS::dump_environment();

  } catch (Dune::Exception &e) {
    std::cerr << "\nDune reported error: " << e.what() << std::endl;
    std::abort();
  } catch (std::exception &e) {
    std::cerr << "\n" << e.what() << std::endl;
    std::abort();
  }
}
