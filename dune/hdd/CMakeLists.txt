set(lib_dune_hdd_sources
    linearelliptic/discretizations/block-swipdg/default__2daluconformsimplex_fem_eigen.cxx
    linearelliptic/discretizations/block-swipdg/default__2daluconformsimplex_fem_istl.cxx
    linearelliptic/discretizations/cg_2daluconformsimplex_fem_istl.cxx
    linearelliptic/discretizations/swipdg_2daluconformsimplex_fem_eigen.cxx
    linearelliptic/discretizations/swipdg_2daluconformsimplex_fem_istl.cxx
   )

dune_add_library("dunehdd" ${lib_dune_hdd_sources}
                 ADD_LIBS ${DUNE_LIBS})
#target_compile_definitions(dunehdd PUBLIC "-DDUNE_HDD_LIB=1")

add_dune_mpi_flags(dunehdd)
if(TBB_FOUND)
  add_dune_tbb_flags(dunehdd)
endif(TBB_FOUND)
add_dune_superlu_flags(dunehdd)
add_dune_alugrid_flags(dunehdd)

target_link_dune_default_libraries(dunehdd)

add_subdirectory(test EXCLUDE_FROM_ALL)
