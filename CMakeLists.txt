# This file is part of the dune-hdd project:
#   http://users.dune-project.org/projects/dune-hdd
# Copyright holders: Felix Albrecht
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#
# Contributors: Rene Milk

if( NOT CMAKE_BUILD_TYPE )
  message( STATUS "No build type selected, setting default build type 'Debug'" )
  set( CMAKE_BUILD_TYPE "DEBUG" CACHE STRING "DEBUG/RELEASE")
endif( NOT CMAKE_BUILD_TYPE )

project( dune-hdd )

# cmake specific
cmake_minimum_required( VERSION 2.8 )
set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
set( CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".lib" ".la" ".a")
set( CMAKE_COLOR_MAKEFILE ON)
set( CMAKE_VERBOSE_MAKEFILE ON )

# local environment
set( ENV{PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/bin:$ENV{PATH}" )
set( ENV{LD_LIBRARY_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/lib:$ENV{LD_LIBRARY_PATH}" )
set( ENV{PKG_CONFIG_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}" )

# dune-stuff and dune-pymor cmake includes
set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_CURRENT_SOURCE_DIR}/../dune-stuff/cmake ${CMAKE_CURRENT_SOURCE_DIR}/../dune-pymor/cmake )
include(DuneUtils)
include(GridUtils)
include(PymorUtils)

# needs to come first to include correct config.h
include_dir( ${dune-hdd_SOURCE_DIR} )
include_dir( ${dune-hdd_BINARY_DIR} )

# user settings
set( ENABLE_MPI
     "0" CACHE STRING
     "set to 1 to enable mpi")

SET( POLORDER
  "1" CACHE STRING
  "POLORDER" )

# we need eigen
if( NOT EIGEN_FOUND )
  message( FATAL_ERROR "eigen3 not found, set PKG_CONFIGK_PATH to include eigen3.pc!" )
endif( NOT EIGEN_FOUND )

CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/cmake_config.h )

add_dune_modules(
  pymor
  stuff
  pdelab
  gdt
  fem
  fem-localfunctions
  localfunctions
  grid-multiscale
  grid
  istl
  geometry
  common
)

# header
file( GLOB_RECURSE hdd "${CMAKE_CURRENT_SOURCE_DIR}/dune/*.hh" )
set( COMMON_HEADER ${hdd} ${DUNE_HEADERS} )
set_source_files_properties( ${COMMON_HEADER} PROPERTIES HEADER_FILE_ONLY 1 )

# libs
set( COMMON_LIBS
     ${DUNE_LIBS}
     ${PARALIBS}
     ${CCGNU_LIBRARIES}
     ${BOOST_LIBS}
     ${GRIDLIBS}
)

add_subdirectory(doc)
add_subdirectory(examples)

