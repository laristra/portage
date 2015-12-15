#~----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

#------------------------------------------------------------------------------#
# If a C++11 compiler is available, then set the appropriate flags
#------------------------------------------------------------------------------#

include(cxx11)

check_for_cxx11_compiler(CXX11_COMPILER)

if(CXX11_COMPILER)
    enable_cxx11()
else()
    message(FATAL_ERROR "C++11 compatible compiler not found")
endif()

#------------------------------------------------------------------------------#
# Set up MPI builds
# (eventually most of this should be pushed down into cinch)
#------------------------------------------------------------------------------#

find_package(MPI REQUIRED)

# TODO:  Modify the below to use wrapper compilers instead of flags
#        (there isn't an obvious good way to do this)
add_definitions(${MPI_CXX_COMPILE_FLAGS})
include_directories(${MPI_CXX_INCLUDE_PATH})
link_directories(${MPI_CXX_LIBRARY_DIRS})

#------------------------------------------------------------------------------#
# Configure Jali
# (this includes the TPLs that Jali will need)
#------------------------------------------------------------------------------#

set(ARCHOS ${CMAKE_SYSTEM_PROCESSOR}_${CMAKE_SYSTEM_NAME})

set(Jali_DIR "$ENV{HOME}/devel/Jali" CACHE STRING
    "Set the locaiton of Jali")

if (NOT Jali_DIR)
  message(FATAL_ERROR "Error: Jali top level installation dir must be defined")
endif()

# Look for the Jali package

find_package(Jali REQUIRED)

message(STATUS "Located Jali")
message(STATUS "Jali_DIR=${Jali_DIR}")

# TODO:  do these steps inside of the Jali build?
list(REVERSE Jali_LIBRARIES)
# add full path to jali libs
unset(_LIBS)
foreach (_lib ${Jali_LIBRARIES})
    set(_LIBS ${_LIBS};${Jali_LIBRARY_DIRS}/lib${_lib}.a)
endforeach()
set(Jali_LIBRARIES ${_LIBS})

message(STATUS "Jali_INCLUDE_DIRS=${Jali_INCLUDE_DIRS}")
message(STATUS "Jali_LIBRARY_DIRS=${Jali_LIBRARY_DIRS}")
message(STATUS "Jali_LIBRARIES=${Jali_LIBRARIES}")
message(STATUS "Jali_TPL_INCLUDE_DIRS=${Jali_TPL_INCLUDE_DIRS}")
message(STATUS "Jali_TPL_LIBRARY_DIRS=${Jali_TPL_LIBRARY_DIRS}")
message(STATUS "Jali_TPL_LIBRARIES=${Jali_TPL_LIBRARIES}")

include_directories(${Jali_INCLUDE_DIRS} ${Jali_TPL_INCLUDE_DIRS})

#-----------------------------------------------------------------------------
# General NGC include directory information
#-----------------------------------------------------------------------------
set(NGC_INCLUDE_DIR "$ENV{NGC_INCLUDE_DIR}" CACHE PATH "NGC include directory")
if(NGC_INCLUDE_DIR)
  message(STATUS "Using NGC_INLCUDE_DIR=${NGC_INCLUDE_DIR}")
endif(NGC_INCLUDE_DIR)

#-----------------------------------------------------------------------------
# Thrust information
#-----------------------------------------------------------------------------
set(ENABLE_THRUST FALSE CACHE BOOL "Use Thrust")
if(ENABLE_THRUST)
  message(STATUS "Enabling compilation with Thrust")
  # allow the user to specify a THRUST_DIR, otherwise use ${NGC_INCLUDE_DIR}
  # NOTE: thrust internally uses include paths from the 'root' directory, e.g.
  #
  #       #include "thrust/device_vector.h"
  #
  #       so the path here should point to the directory that has thrust as
  #       a subdirectory.
  # Use THRUST_DIR directly if specified, otherwise try to build from NGC
  set(THRUST_DIR "${NGC_INCLUDE_DIR}" CACHE PATH "Thrust directory")
  message(STATUS "Using THRUST_DIR=${THRUST_DIR}")
  
  # Allow for swapping backends - should this be in CACHE?
  set(THRUST_BACKEND "THRUST_DEVICE_SYSTEM_OMP" CACHE STRING "Thrust backend")
  message(STATUS "Using ${THRUST_BACKEND} as Thrust backend.")
  include_directories(${THRUST_DIR})
  add_definitions(-DTHRUST)
  add_definitions(-DTHRUST_DEVICE_SYSTEM=${THRUST_BACKEND})

  if("${THRUST_BACKEND}" STREQUAL "THRUST_DEVICE_SYSTEM_OMP")
    FIND_PACKAGE( OpenMP REQUIRED)
    if(OPENMP_FOUND)
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
  endif ()

endif(ENABLE_THRUST)
