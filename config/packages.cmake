#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]
#------------------------------------------------------------------------------#
# If we are building with FleCSI, then we need a modern C++ compiler
#------------------------------------------------------------------------------#
if(ENABLE_FleCSI)
  # we already checked for CXX14 in project.cmake
  if(NOT CXX14_COMPILER)
    message(STATUS "C++14 compatible compiler not found")
  endif()
endif()

#------------------------------------------------------------------------------#
# Set up MPI builds
# (eventually most of this should be pushed down into cinch)
#------------------------------------------------------------------------------#
if (ENABLE_MPI)
  find_package(MPI REQUIRED)

# TODO:  Modify the below to use wrapper compilers instead of flags
#        (there isn't an obvious good way to do this)
  add_definitions(${MPI_CXX_COMPILE_FLAGS})
  include_directories(${MPI_CXX_INCLUDE_PATH})
  link_directories(${MPI_CXX_LIBRARY_DIRS})
endif ()


set(ARCHOS ${CMAKE_SYSTEM_PROCESSOR}_${CMAKE_SYSTEM_NAME})

#-----------------------------------------------------------------------------
# FleCSI and FleCSI-SP location
#-----------------------------------------------------------------------------

set(ENABLE_FleCSI FALSE CACHE BOOL "Use FleCSI")
if (ENABLE_FleCSI)
 
 find_package(FleCSI REQUIRED)
 message(STATUS "FleCSI_LIBRARIES=${FleCSI_LIBRARIES}" )
 include_directories(${FleCSI_INCLUDE_DIR})
 message(STATUS "FleCSI_INCLUDE_DIRS=${FleCSI_INCLUDE_DIR}")

 find_package(FleCSISP REQUIRED)
 message(STATUS "FleCSISP_LIBRARIES=${FleCSISP_LIBRARIES}" )
 include_directories(${FleCSISP_INCLUDE_DIR})
 message(STATUS "FleCSISP_INCLUDE_DIRS=${FleCSISP_INCLUDE_DIR}")

  ######################################################################
  # This is a placeholder for how we would do IO with FleCSI
  # There are still some issues with dumping the targetMesh data
  #
  # WARNING!!! THIS IS POTENTIALLY FRAGILE
  # it appears to work, but could cause conflicts with EXODUS and
  # other libraries used by Jali
  #
  # FOR NOW THIS IS DISABLED UNTIL WE CAN GET A PROPER WORKAROUND
  ######################################################################
  # STRING(REPLACE "flecsi" "flecsi-tpl" FLECSI_TPL_DIR ${FLECSI_INSTALL_DIR})
  # message(STATUS "FLECSI_TPL_DIR=${FLECSI_TPL_DIR}")
  # if(IS_DIRECTORY ${FLECSI_TPL_DIR})
  #   find_library(EXODUS_LIBRARY
  #     NAMES exodus
  #     PATHS ${FLECSI_TPL_DIR}
  #     PATH_SUFFIXES lib
  #     NO_DEFAULT_PATH)
  #   find_path(EXODUS_INCLUDE_DIR
  #     NAMES exodusII.h
  #     PATHS ${FLECSI_TPL_DIR}
  #     PATH_SUFFIXES include
  #     NO_DEFAULT_PATH)

  #   if(EXODUS_LIBRARY AND EXODUS_INCLUDE_DIR)
  #     set(FLECSI_LIBRARIES ${EXODUS_LIBRARY} ${FLECSI_LIBRARIES})
  #     include_directories(${EXODUS_INCLUDE_DIR})
  #     add_definitions(-DHAVE_EXODUS)
  #   endif(EXODUS_LIBRARY AND EXODUS_INCLUDE_DIR)

  # endif(IS_DIRECTORY ${FLECSI_TPL_DIR})
endif()


#------------------------------------------------------------------------------#
# Configure Jali
# (this includes the TPLs that Jali will need)
#------------------------------------------------------------------------------#

if (Jali_DIR)

   # Look for the Jali package

   find_package(Jali REQUIRED
                HINTS ${Jali_DIR}/lib)

   message(STATUS "Located Jali")
   message(STATUS "Jali_DIR=${Jali_DIR}")

   # add full path to jali libs
   unset(_LIBS)
   foreach (_lib ${Jali_LIBRARIES})
      set(_LIBS ${_LIBS};${Jali_LIBRARY_DIRS}/lib${_lib}.a)
   endforeach()
   set(Jali_LIBRARIES ${_LIBS})

   # message(STATUS "Jali_INCLUDE_DIRS=${Jali_INCLUDE_DIRS}")
   # message(STATUS "Jali_LIBRARY_DIRS=${Jali_LIBRARY_DIRS}")
   # message(STATUS "Jali_LIBRARIES=${Jali_LIBRARIES}")
   # message(STATUS "Jali_TPL_INCLUDE_DIRS=${Jali_TPL_INCLUDE_DIRS}")
   # message(STATUS "Jali_TPL_LIBRARY_DIRS=${Jali_TPL_LIBRARY_DIRS}")
   # message(STATUS "Jali_TPL_LIBRARIES=${Jali_TPL_LIBRARIES}")

   include_directories(${Jali_INCLUDE_DIRS} ${Jali_TPL_INCLUDE_DIRS})

endif (Jali_DIR)

#------------------------------------------------------------------------------#
# Configure LAPACKE
#------------------------------------------------------------------------------#

if (LAPACKE_DIR)
  # Directly look for cmake config file in LAPACKE_DIR
  find_package(LAPACKE REQUIRED NO_MODULE HINTS ${LAPACKE_DIR})
else (LAPACKE_DIR)
  # Use FindLAPACKE.cmake provided by cinch or cmake to find it
  # FindLAPACKE.cmake provided by cinch requires PC_LAPACKE_NCLUDE_DIRS and
  # PC_LAPACKE_LIBRARY to be able to find LAPACKE
  find_package(LAPACKE REQUIRED)
endif (LAPACKE_DIR)

if (LAPACKE_FOUND) 
   include_directories(${LAPACKE_INCLUDE_DIRS})
   add_definitions("-DHAVE_LAPACKE")
endif (LAPACKE_FOUND)

#-----------------------------------------------------------------------------
# General NGC include directory information
#-----------------------------------------------------------------------------
set(NGC_INCLUDE_DIR "$ENV{NGC_INCLUDE_DIR}" CACHE PATH "NGC include directory")
if(NGC_INCLUDE_DIR)
  message(STATUS "Using NGC_INCLUDE_DIR=${NGC_INCLUDE_DIR}")
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

  # Allow for swapping backends
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
    endif(OPENMP_FOUND)
  endif ()

  if("${THRUST_BACKEND}" STREQUAL "THRUST_DEVICE_SYSTEM_TBB")
    FIND_PACKAGE(TBB REQUIRED)
    if(TBB_FOUND)
      include_directories(${TBB_INCLUDE_DIRS})
      link_directories(${TBB_LIBRARY_DIRS})
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ltbb")
    endif(TBB_FOUND)
  endif()

endif(ENABLE_THRUST)


#-----------------------------------------------------------------------------
# Find Boost
#-----------------------------------------------------------------------------
if(NOT ENABLE_THRUST)
  find_package(Boost REQUIRED)
  if(Boost_FOUND)
    message(STATUS "Boost location: ${Boost_INCLUDE_DIRS}")
    include_directories( ${Boost_INCLUDE_DIRS} )
  endif(Boost_FOUND)
endif(NOT ENABLE_THRUST)
