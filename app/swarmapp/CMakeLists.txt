#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]


#------------------------------------------------------------------------------#
# Add a rule to build the executable
#------------------------------------------------------------------------------#

# Swarm test
add_executable(swarmapp swarmapp.cc)
target_link_libraries(swarmapp PRIVATE portage)

# Compare
add_executable(swarmcompare swarmcompare.cc)

# Profiling option (to allow pausing and resuming profiling during execution)
set(ENABLE_PROFILE FALSE CACHE BOOL "Enable VTune Profiling")
if(ENABLE_PROFILE)
  set(VTUNE_INCLUDE_DIR "${HOME}" CACHE PATH "Set the locaiton of VTune include")
  set(VTUNE_LIB "${HOME}" CACHE STRING "Set the VTune library")

  target_include_directories(swarmcompare PRIVATE ${VTUNE_INCLUDE_DIR})

  target_link_libraries(swarmapp PRIVATE ${VTUNE_LIB})
endif(ENABLE_PROFILE)


if (ENABLE_APP_TESTS)
  add_subdirectory(test)
endif (ENABLE_APP_TESTS)
