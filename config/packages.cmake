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

