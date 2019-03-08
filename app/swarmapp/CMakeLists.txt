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
target_link_libraries(swarmapp portage
  ${EXTRA_LIBS} ${LAPACKX_LIBRARIES} ${MPI_CXX_LIBRARIES})

# Compare
add_executable(compare compare.cc)

#TCMalloc option
set(ENABLE_TCMALLOC TRUE CACHE BOOL "Use TCMalloc")
if(ENABLE_TCMALLOC)
  set(TCMALLOC_LIB "${HOME}" CACHE STRING "Set the TCMalloc library")

  target_link_libraries(swarmapp ${TCMALLOC_LIB})
endif(ENABLE_TCMALLOC)

# Profiling option (to allow pausing and resuming profiling during execution)
set(ENABLE_PROFILE FALSE CACHE BOOL "Enable VTune Profiling")
if(ENABLE_PROFILE)
  set(VTUNE_INCLUDE_DIR "${HOME}" CACHE PATH "Set the locaiton of VTune include")
  set(VTUNE_LIB "${HOME}" CACHE STRING "Set the VTune library")

  include_directories(${VTUNE_INCLUDE_DIR})

  target_link_libraries(swarmapp ${VTUNE_LIB})
endif(ENABLE_PROFILE)


if (ENABLE_APP_TESTS)
  add_subdirectory(test)
endif (ENABLE_APP_TESTS)
