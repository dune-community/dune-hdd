// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include <dune/stuff/test/gtest/gtest.h>
#include <dune/stuff/test/common.hh>

#include "OS2014_nonparametric_convergence_study.hh"

#if HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_GRID_MULTISCALE


void nonparametric_convergence_study(const std::string visualization)
{
  const NonparametricEocTestCaseType test_case;
  test_case.print_header(DSC_LOG_INFO);
  DSC_LOG_INFO << std::endl;
  NonparametricEocStudyType study(test_case,
                                  {"energy", "eta_NC_ESV2007", "eta_R_ESV2007", "eta_DF_ESV2007", "eta_ESV2007",
                                   "eff_ESV2007"},
                                  {},
                                  visualization);
  Stuff::Test::check_eoc_study_for_success(study, study.run_eoc(DSC_LOG_INFO));
} // ... nonparametric_convergence_study(...)


void nonparametric_convergence_study_alternative_summation()
{
  const NonparametricEocTestCaseType test_case;
  NonparametricEocStudyType study(test_case,
                                  {"energy", "eta_ESV2007", "eff_ESV2007", "eta_ESV2007_alt", "eff_ESV2007_alt"});
  Stuff::Test::check_eoc_study_for_success(study, study.run_eoc(DSC_LOG_INFO));
} // ... nonparametric_convergence_study_alternative_summation(...)


void nonparametric_block_convergence_study(const std::string& partitioning)
{
  const NonparametricBlockEocTestCaseType test_case(partitioning);
  NonparametricBlockEocStudyType study(test_case,
                                       {"energy", "eta_NC_OS2014", "eta_R_OS2014", "eta_DF_OS2014", "eta_OS2014",
                                        "eff_OS2014"});
  Stuff::Test::check_eoc_study_for_success(study, study.run_eoc(DSC_LOG_INFO));
} // ... nonparametric_block_convergence_study(...)


#endif // HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_GRID_MULTISCALE
