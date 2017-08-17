



project(portage)

cinch_minimum_required(1.0)



# If a C++14 compiler is available, then set the appropriate flags
include(cxx14)
check_for_cxx14_compiler(CXX14_COMPILER)
if(CXX14_COMPILER)
  enable_cxx14()
else()
  message(STATUS "C++14 compatible compiler not found")
endif()

# If we couldn't find a C++14 compiler, try to see if a C++11 
# compiler is available, then set the appropriate flags
if (NOT CXX14_COMPILER)
  include(cxx11)
  check_for_cxx11_compiler(CXX11_COMPILER)
  if(CXX11_COMPILER)
    enable_cxx11()
  else()
    message(FATAL_ERROR "C++11 compatible compiler not found")
  endif()
endif()

cinch_add_application_directory(app)
cinch_add_library_target(portage src)


# Create some subcomponent libraries for the purposes of unit testing
#
# NOT ABLE TO DO EITHER OF THESE

if (ENABLE_UNIT_TESTS)
#  cinch_add_library_target(search src/search)
#  cinch_add_library_target(intersect src/intersect)
#  cinch_add_library_target(remap src/remap)
#  cinch_add_library_target(state src/state)

#  cinch_add_library_target(portage_search src/search)
#  cinch_add_library_target(portage_intersect src/intersect)
#  cinch_add_library_target(portage_remap src/remap)
#  cinch_add_library_target(portage_state src/state)
endif()

# Add application tests
# May pull this logic into cinch at some future point
option(ENABLE_APP_TESTS "Enable testing of full app" OFF)
if(ENABLE_APP_TESTS)
  enable_testing()
endif()

#------------------------------------------------------------------------------#
# 
#------------------------------------------------------------------------------#

set(CINCH_HEADER_SUFFIXES "\\.h")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")
