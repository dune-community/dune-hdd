PROJECT(dune-detailed-solvers)

SET( CMAKE_CXX_COMPILER
  "g++-4.6" CACHE FORCE
    "" )

#set minimum cmake version
cmake_minimum_required(VERSION 2.6)
cmake_policy(SET CMP0009 NEW)
# tell cmake where to look for .cmake files
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )

# stuff needed for dune
include(DuneUtils)
include(UseDoxygen)
include (FindPkgConfig)

# General Settings
SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
SET( CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".lib" ".la" ".a" )
SET( CMAKE_COLOR_MAKEFILE ON )

# Options, that can be changed by the User in order to customise the module
SET( BLAS_LIB
  "blas" CACHE STRING
  "BLAS LIB" )

SET( GRIDTYPE
  "SGRID" CACHE STRING
  "GRIDTYPE" )

SET( GRIDDIM
  "2" CACHE STRING
  "GRIDDIM" )

SET( POLORDER
  "1" CACHE STRING
  "POLORDER" )

SET( ALUGRID_BASE_PATH
  "${CMAKE_CURRENT_SOURCE_DIR}/../local" CACHE STRING
  "ALUGRID_BASE_PATH" )

SET( CUSTOM_FLAGS
  "-Wall -O0 -DDEBUG -g -fdiagnostics-show-option" CACHE LIST
  "CUSTOM FLAGS" )

set( BOOST_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../local" )
find_package( Boost 1.48.0 REQUIRED )
include_directories(${Boost_INCLUDE_DIR})
link_directories( ${Boost_LIBRARY_DIRS} )

set( EIGEN3_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../local/include/eigen3")
find_package( Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package( SuperLU REQUIRED)
include_directories( ${SUPERLU_INCLUDES} )
link_directories( ${SUPERLU_LIBRARY_DIRS} )

# needs to come first to include correct config.h
include_dir( ${dune-detailed-solvers_SOURCE_DIR} )
add_cxx_flags( -DHAVE_CONFIG_H )
add_cxx_flags( -D${GRIDTYPE} )
add_cxx_flags( -DGRIDDIM=${GRIDDIM} )
add_cxx_flags( -DPOLORDER=${POLORDER} )
add_cxx_flags( -I${ALUGRID_BASE_PATH}/include
  -I${ALUGRID_BASE_PATH}/include/serial
  -I${ALUGRID_BASE_PATH}/include/duneinterface -DENABLE_ALUGRID )

add_cxx_flags(-DBOOST_DATE_TIME_NO_LOCALE)
link_directories( "${ALUGRID_BASE_PATH}/lib" )

add_cxx_flags(-std=c++0x)

add_cxx_flags(${CUSTOM_FLAGS})
# the other dune modules
ADD_DUNE_MODULES(common grid-multiscale grid fem geometry localfunctions stuff detailed-discretizations istl pdelab)

# Source listing, relative to the position of this file
FILE( GLOB_RECURSE detailed-discretizations "../dune-detailed-discretizations/*.hh" )
FILE( GLOB_RECURSE detailed-discretizations_cc "../dune-detailed-discretizations/*.cc" )
FILE( GLOB_RECURSE detailed-solvers "../dune-detailed-solvers/*.hh" )
FILE( GLOB_RECURSE detailed-solvers_cc "../dune-detailed-solvers/*.cc" )
FILE( GLOB_RECURSE common "../dune-common/*.hh" )
FILE( GLOB_RECURSE common_cc "../dune-common/*.cc" )
FILE( GLOB_RECURSE grid "../dune-grid/*.hh" )
FILE( GLOB_RECURSE grid_cc "../dune-grid/*.cc" )
FILE( GLOB_RECURSE grid_multiscale "../dune-grid-multiscale/*.hh" )
FILE( GLOB_RECURSE grid-multiscale_cc "../dune-grid-multiscale/*.cc" )
FILE( GLOB_RECURSE fem "../dune-fem/*.hh" )
FILE( GLOB_RECURSE fem_cc "../dune-fem/*.cc" )
FILE( GLOB_RECURSE geometry "../dune-geometry/*.hh" )
FILE( GLOB_RECURSE geometry_cc "../dune-geometry/*.cc" )
FILE( GLOB_RECURSE localfunctions "../dune-localfunctions/*.hh" )
FILE( GLOB_RECURSE localfunctions_cc "../dune-localfunctions/*.cc" )
FILE( GLOB_RECURSE stuff "../dune-stuff/*.hh" )
FILE( GLOB_RECURSE stuff_cc "../dune-stuff/*.cc" )
FILE( GLOB_RECURSE istl "../dune-istl/*.hh" )
FILE( GLOB_RECURSE istl_cc "../dune-istl/*.cc" )
FILE( GLOB_RECURSE pdelab "../dune-pdelab/*.hh" )
FILE( GLOB_RECURSE pdelab_cc "../dune-pdelab/*.cc" )

set_source_files_properties( ${common} ${grid} ${grid_multiscale} ${fem} ${localfunctions} ${geometry} ${stuff} ${detailed-discretizations} ${detailed-solvers} ${detailed-solvers_cc} ${istl_cc} ${pdelab} ${pdelab_cc}
${common_cc} ${grid_cc} ${grid_multiscale_cc} ${fem_cc} ${localfunctions_cc} ${geometry_cc} ${stuff_cc} ${detailed-discretizations_cc} ${istl} PROPERTIES HEADER_FILE_ONLY 1 )

set(LIST_LIBS ${common} ${grid} ${grid_multiscale} ${fem} ${localfunctions} ${geometry} ${stuff} ${detailed-discretizations} ${detailed-solvers} ${istl} ${pdelab} ${pdelab_cc}
${common_cc} ${grid_cc} ${fem_cc} ${localfunctions_cc} ${grid_multiscale_cc} ${geometry_cc} ${stuff_cc} ${detailed-discretizations_cc} ${istl_cc} )

set( COMMON_LIBS "fem" "grid" "dunegeometry" "common" ${BLAS_LIB} "alugrid" "lapack" "boost_filesystem" "boost_system" "superlu" )

ADD_EXECUTABLE( examples_stationary_linear_elliptic_continuousgalerkin_dune_detailed_discretizations "examples/stationary/linear/elliptic/continuousgalerkin/dune_detailed_discretizations.cc" ${LIST_LIBS} )
TARGET_LINK_LIBRARIES( examples_stationary_linear_elliptic_continuousgalerkin_dune_detailed_discretizations ${COMMON_LIBS} )

ADD_EXECUTABLE( examples_stationary_linear_elliptic_multiscale_semicontinuousgalerkin_dune_detailed_discretizations "examples/stationary/linear/elliptic/multiscale/semicontinuousgalerkin/dune_detailed_discretizations.cc" ${LIST_LIBS} )
TARGET_LINK_LIBRARIES( examples_stationary_linear_elliptic_multiscale_semicontinuousgalerkin_dune_detailed_discretizations ${COMMON_LIBS} )
