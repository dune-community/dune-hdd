# This file is part of the dune-hdd project:
#   http://users.dune-project.org/projects/dune-hdd
# Copyright holders: Felix Albrecht
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#
# Contributors: Rene Milk

if(NOT CMAKE_BUILD_TYPE)
  message( STATUS "No build type selected, setting default build type 'Debug'" )
  set( CMAKE_BUILD_TYPE "DEBUG" CACHE STRING "DEBUG/RELEASE")
endif(NOT CMAKE_BUILD_TYPE)

project(dune-hdd CXX)

# local environment
set( ENV{PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/bin:$ENV{PATH}" )
set( ENV{LD_LIBRARY_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/lib:$ENV{LD_LIBRARY_PATH}" )
set( ENV{PKG_CONFIG_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}" )

# cmake specific
cmake_minimum_required( VERSION 2.8 )

if(CMAKE_WITH_AUTOTOOLS)

  # cmake specific
  set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
  set( CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".lib" ".la" ".a")
  set( CMAKE_COLOR_MAKEFILE ON)
  set( CMAKE_VERBOSE_MAKEFILE ON )

  # local environment
  set( ENV{PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/bin:$ENV{PATH}" )
  set( ENV{LD_LIBRARY_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/lib:$ENV{LD_LIBRARY_PATH}" )
  set( ENV{PKG_CONFIG_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}" )

  # dune-stuff and dune-pymor cmake includes
  set(CMAKE_MODULE_PATH
      ${CMAKE_CURRENT_SOURCE_DIR}/cmake
      ${CMAKE_CURRENT_SOURCE_DIR}/../dune-stuff/cmake
      ${CMAKE_CURRENT_SOURCE_DIR}/../dune-stuff/cmake/modules
      ${CMAKE_CURRENT_SOURCE_DIR}/../dune-pymor/cmake
    )
  include(DuneUtils)
  include(GridUtils)
  include(PymorUtils)

  # needs to come first to include correct config.h
  include(${dune-hdd_SOURCE_DIR})
  include(${dune-hdd_BINARY_DIR})

  add_dune_modules(
    hdd
    pymor
    stuff
    pdelab
    typetree
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
  file( GLOB_RECURSE hdd "${CMAKE_CURRENT_SOURCE_DIR}/*.hh" )
  set( COMMON_HEADER ${hdd} ${DUNE_HEADERS} )
  set_source_files_properties( ${COMMON_HEADER} PROPERTIES HEADER_FILE_ONLY 1 )

  # libs
  set( COMMON_LIBS
      ${DUNE_LIBS}
      ${DUNE_DEFAULT_LIBS}
      ${PARALIBS}
      ${CCGNU_LIBRARIES}
      ${Boost_SYSTEM_LIBRARY}
      ${Boost_FILESYSTEM_LIBRARY}
      ${GRIDLIBS}
  )

  add_subdirectory(examples)
  add_subdirectory(test)

else(CMAKE_WITH_AUTOTOOLS)
  # this assumes that we run dunecontrol with the --use-cmake option

  set(CMAKE_COLOR_MAKEFILE ON)

  # guess dune-common build dir
  if(NOT (dune-common_DIR OR dune-common_ROOT OR "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR ${PROJECT_BINARY_DIR})
  endif()

  # find dune-common and set the module path
  find_package(dune-common REQUIRED)
  list(APPEND
       CMAKE_MODULE_PATH
       ${dune-common_MODULE_PATH}
       "${PROJECT_SOURCE_DIR}/cmake/modules")

  #include the dune macros
  include(DuneMacros)

  # start a dune project with information from dune.module
  dune_project()

  # dune-stuff and dune-pymor cmake includes
  list(APPEND CMAKE_MODULE_PATH
       ${CMAKE_CURRENT_SOURCE_DIR}/cmake
       ${CMAKE_CURRENT_SOURCE_DIR}/../dune-pymor/cmake
       ${CMAKE_CURRENT_SOURCE_DIR}/../dune-stuff/cmake)

  include(DuneUtils)
  include(GridUtils)
  include(PymorUtils)

  link_directories("/usr/lib" ${CCGNU_LIBRARY_DIRS})

  # header
  file(GLOB_RECURSE hdd "${CMAKE_CURRENT_SOURCE_DIR}/*.hh")
  set(COMMON_HEADER ${hdd} ${DUNE_HEADERS})
  # add header of dependent modules for header listing
  foreach(_mod ${ALL_DEPENDENCIES})
    file(GLOB_RECURSE HEADER_LIST "${CMAKE_CURRENT_SOURCE_DIR}/../${_mod}/*.hh")
    list(APPEND COMMON_HEADER ${HEADER_LIST})
  endforeach(_mod DEPENDENCIES)
  set_source_files_properties(${COMMON_HEADER} PROPERTIES HEADER_FILE_ONLY 1)
  
  #if (ALUGRID_FOUND)
  #  set(GRIDTYPE "ALUGRID_CONFORM" CACHE STRING "GRIDTYPE" FORCE)
  #endif (ALUGRID_FOUND)
  
  # libs
  set(COMMON_LIBS ${DUNE_DEFAULT_LIBS})

  #disable most warnings from dependent modules
  foreach(_mod ${ALL_DEPENDENCIES})
      dune_module_to_uppercase(_upper_case "${_mod}")
      if(${_mod}_INCLUDE_DIRS)
          foreach( _idir ${${_mod}_INCLUDE_DIRS} )
              add_definitions("-isystem ${_idir}")
          endforeach( _idir )
      endif(${_mod}_INCLUDE_DIRS)
  endforeach(_mod DEPENDENCIES)

  # search for spe10 data files
  find_file(SPE10MODEL1DATA
            NAMES perm_case1.dat
            PATHS "${dune-hdd_SOURCE_DIR}/../local/src/spe10/model1"
            DOC "Location of perm_case1.dat"
            NO_DEFAULT_PATH)

  add_subdirectory(examples EXCLUDE_FROM_ALL)
  add_subdirectory(test EXCLUDE_FROM_ALL)

  finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
endif(CMAKE_WITH_AUTOTOOLS)

