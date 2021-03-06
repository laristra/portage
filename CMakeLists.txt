#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]

cmake_minimum_required(VERSION 3.13)

project(portage CXX C)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (CMAKE_VERSION_MAJOR GREATER_EQUAL 3.13)
  CMAKE_POLICY(SET CMP0079 NEW)  # allow target_link_libraries to reference
                                 # targets from other directories
endif()

cmake_policy(SET CMP0074 NEW)  # Don't ignore Pkg_ROOT variables


# SEMANTIC VERSION NUMBERS - UPDATE DILIGENTLY
# As soon as a change with a new version number is merged into the master,
# tag the central repository.

set(PORTAGE_VERSION_MAJOR 2)
set(PORTAGE_VERSION_MINOR 3)
set(PORTAGE_VERSION_PATCH 0)


# Top level target
add_library(portage INTERFACE)

# Alias (Daniel Pfeiffer, Effective CMake) - this allows other
# projects that use Pkg as a subproject to find_package(Nmspc::Pkg)
# which does nothing because Pkg is already part of the project

add_library(portage::portage ALIAS portage)
set(PORTAGE_LIBRARIES portage::portage)

# Find our modules first
if (CMAKE_VERSION GREATER_EQUAL 3.15)
  list(PREPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
else ()
  set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")
endif ()

#------------------------------------------------------------------------------#
# Enable doxygen
#------------------------------------------------------------------------------#
option(ENABLE_DOXYGEN "Enable Doxygen documentation" OFF)
option(ENABLE_DOXYGEN_WARN "Enable Doxygen warnings" OFF)

if (ENABLE_DOXYGEN)
  include(cmake/doxygen.cmake)
  add_doxygen()
endif()

#------------------------------------------------------------------------------#
# Enable timings
#------------------------------------------------------------------------------#
option(ENABLE_APP_TIMINGS "Enable timings" OFF)

#------------------------------------------------------------------------------#
# Enable PORTAGE_DEBUG flag
#------------------------------------------------------------------------------#
option(PORTAGE_DEBUG "Additional checks will be performed and info printed" OFF)
if (PORTAGE_DEBUG)
  add_definitions(-DPORTAGE_DEBUG)
endif()

#----------------------------------------------------------------------------
# Find packages here and set CMake variables. Set link dependencies
# and compile definitions for portage_support target in
# portage/support and so on. Since the top level portage::portage
# target depends on portage_support, the transitive dependencies will
# be picked up by projects linking to Portage
# ----------------------------------------------------------------------------

#------------------------------------------------------------------------------#
# Find Tangram (interface reconstruction package)
# Will find Wonton support package as part of it
#------------------------------------------------------------------------------#

set(PORTAGE_HAS_TANGRAM False CACHE BOOL "Is Tangram support enabled?")
if (PORTAGE_ENABLE_TANGRAM)

  find_package(TANGRAM QUIET REQUIRED NAMES tangram)

  message(STATUS "TANGRAM FOUND? ${TANGRAM_FOUND}")

  target_link_libraries(portage INTERFACE ${TANGRAM_LIBRARIES})
  message(STATUS "TANGRAM_LIBRARIES ${TANGRAM_LIBRARIES}")
  message(STATUS "TANGRAM_INCLUDE_DIR ${TANGRAM_INCLUDE_DIR}")

  set(PORTAGE_HAS_TANGRAM True CACHE BOOL "Is Tangram support enabled?" FORCE)
  
  
else ()

  # No TANGRAM but we still need Wonton - find it
  
  # Link with an existing installation of Wonton, if provided. 
  find_package(WONTON QUIET REQUIRED NAMES wonton)
  
  target_link_libraries(portage INTERFACE ${WONTON_LIBRARIES})
  message(STATUS "WONTON_LIBRARIES=${WONTON_LIBRARIES}" )

endif() 


# Check that the Portage's options are compatible with Wonton (and Tangram)
if (PORTAGE_ENABLE_THRUST AND NOT WONTON_ENABLE_THRUST)
  message(FATAL_ERROR "Thrust enabled for Portage but Wonton are not built with Thrust")
endif ()
if (NOT PORTAGE_ENABLE_THRUST AND WONTON_ENABLE_THRUST)
  message(FATAL_ERROR "Thrust disabled for Portage but Wonton is built with Thrust")
endif ()

if (PORTAGE_ENABLE_MPI AND NOT WONTON_ENABLE_MPI)
  message(FATAL_ERROR "MPI enabled for Portage but Wonton is not compiled with MPI")
endif ()

# This combination will probably work but we are trying to enforce
# consistency throughout the software stack
if (NOT PORTAGE_ENABLE_MPI AND WONTON_ENABLE_MPI)
  message(FATAL_ERROR "MPI disabled for Portage but Wonton is compiled with MPI")
endif ()

if (PORTAGE_ENABLE_Jali AND NOT WONTON_ENABLE_Jali)
  message(FATAL_ERROR "Jali enabled for Portage but Wonton is not compiled with Jali")
endif ()

# This combination will probably work but we are trying to enforce
# consistency throughout the software stack
if (NOT PORTAGE_ENABLE_Jali AND WONTON_ENABLE_Jali)
  message(FATAL_ERROR "Jali not enabled for Portage but Wonton is compiled with Jali")
endif ()

if (PORTAGE_ENABLE_FleCSI AND NOT WONTON_ENABLE_FleCSI)
  message(FATAL_ERROR "FleCSI enabled for Portage but Wonton is not compiled with FleCSI")
endif ()

# This combination will probably work but we are trying to enforce
# consistency throughout the software stack
if (NOT PORTAGE_ENABLE_FleCSI AND WONTON_ENABLE_FleCSI)
  message(FATAL_ERROR "FleCSI not enabled for Portage but Wonton is compiled with FleCSI")
endif ()


if (ENABLE_UNIT_TESTS)

  # This needs to be set in the root directory for tests to be run by
  # 'make test' or ctest
  enable_testing()

  include(cmake/unittest.cmake)
  
endif ()


#-----------------------------------------------------------------------------
# Recurse down the source directories building up dependencies
#-----------------------------------------------------------------------------

add_subdirectory(portage)


# In addition to the include directories of the source, we need to
# include the build or directory to get the autogenerated
# portage-config.h (The first of these is needed if Wonton is included
# as a submodule, the second is needed for the auto-generated config
# file if Portage is included as a submodule, the third is to get the
# autogenerated config header if Tangram is being compiled separately
# and the last is for dependencies in installations)

target_include_directories(portage INTERFACE
  $<BUILD_INTERFACE:${portage_SOURCE_DIR}>
  $<BUILD_INTERFACE:${portage_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>)


# Tangram targets

install(TARGETS portage
  EXPORT portage_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include
  INCLUDES DESTINATION include
  )


#-----------------------------------------------------------------------------
# Add any applications built upon Portage
#-----------------------------------------------------------------------------

add_subdirectory(app)


#-----------------------------------------------------------------------------
# Prepare output for configuration files to be used by projects importing Portage
#-----------------------------------------------------------------------------

# Write a configuration file from template replacing only variables enclosed
# by the @ sign.
configure_file(${PROJECT_SOURCE_DIR}/cmake/portageConfig.cmake.in 
  portageConfig.cmake @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/portageConfig.cmake DESTINATION lib/cmake/portage)


# write out a version file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(portageConfigVersion.cmake
  VERSION "${PORTAGE_MAJOR_VERSION}.${PORTAGE_MINOR_VERSION}.${PORTAGE_PATCH_VERSION}"
  COMPATIBILITY SameMajorVersion)
install(FILES ${PROJECT_BINARY_DIR}/portageConfigVersion.cmake
  DESTINATION lib/cmake/portage)


# export targets

install(EXPORT portage_LIBRARIES
  FILE portageTargets.cmake
  NAMESPACE portage::
  EXPORT_LINK_INTERFACE_LIBRARIES
  DESTINATION lib/cmake/portage)


# Dynamically configured header files that contains defines like
# PORTAGE_HAS_TANGRAM etc. if enabled

configure_file(${PROJECT_SOURCE_DIR}/config/portage-config.h.in
  ${PROJECT_BINARY_DIR}/portage-config.h @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/portage-config.h
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include)

