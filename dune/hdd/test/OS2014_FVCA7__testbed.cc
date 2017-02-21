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

#include <dune/hdd/linearelliptic/testcases/OS2014.hh>

#include "linearelliptic-block-swipdg.hh"
#include "linearelliptic-swipdg.hh"

#include <dune/grid/spgrid.hh>

using namespace Dune;
using namespace HDD;

// typedef ALUGrid< 2, 2, simplex, conforming, MPI_Comm > GridType;
// typedef YaspGrid< 2 > GridType;
typedef SPGrid<double, 2, SPIsotropicRefinement, MPI_Comm> GridType;

class OS2014_FVCA7__estimator_study {
  typedef LinearElliptic::TestCases::OS2014<GridType> TestCaseType;
  typedef LinearElliptic::TestCases::OS2014Multiscale<GridType>
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
    TestCaseType test_case(DSC_CONFIG_GET("refinements", 4), DSC_CONFIG_GET("start_cell_count", 16));
    test_case.print_header(DSC_LOG_INFO_0);
    DSC_LOG_INFO_0 << std::endl;
    EocStudyType eoc_study(test_case, {"L2", "energy", "eta_ESV2007"}, std::vector<std::string>(), "OS2014");
    eoc_study.run_eoc(DSC_LOG_INFO_0);
  }
}; // class OS2014_FVCA7__estimator_study

extern template class LinearElliptic::Tests::SWIPDGStudyExpectations<
    LinearElliptic::TestCases::OS2014<GridType>, 1>;

extern template class LinearElliptic::Tests::BlockSWIPDGStudyExpectations<
    LinearElliptic::TestCases::OS2014Multiscale<GridType>, 1>;

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
    DSC_PROFILER.setOutputdir(DSC_CONFIG_GET("global.datadir", "data"));

    DSC_LOG_INFO_0 << DSC_CONFIG.report_string() << std::endl;

    DSC::TimedLogger().create(-1, -1);
    const size_t threads =
        DSC_CONFIG.has_key(
            "threading.max_count") // <- doing this so complicated to
            ? DSC_CONFIG.get<size_t>(
                  "threading.max_count") //    silence the WARNING: ...
#if HAVE_TBB
            : std::thread::hardware_concurrency();
#else
            : 1u;
#endif
    DS::threadManager().set_max_threads(threads);

    OS2014_FVCA7__estimator_study::ESV2007_fine_triangulation();

  } catch (Dune::Exception &e) {
    std::cerr << "\nDune reported error: " << e.what() << std::endl;
    std::abort();
  } catch (std::exception &e) {
    std::cerr << "\n" << e.what() << std::endl;
    std::abort();
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
    std::abort();
  } // try
}
