#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]


#------------------------------------------------------------------------------#
# Add a rule to build the executable
#------------------------------------------------------------------------------#

# Msm test
if (ENABLE_MPI)

add_executable(msmapp msmapp.cc)
target_link_libraries(msmapp portage
  ${EXTRA_LIBS} ${LAPACKE_LIBRARIES} ${MPI_CXX_LIBRARIES})

if (Jali_DIR)
    # MPI additions - need to push this down into cinch
    target_link_libraries(msmapp ${MPI_CXX_LIBRARIES})
    set_target_properties(msmapp PROPERTIES
                          LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

    # Jali additions
    include_directories(${Jali_INCLUDE_DIRS} ${Jali_TPL_INCLUDE_DIRS})
    target_link_libraries(msmapp ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES})
endif (Jali_DIR)

#TCMalloc option
set(ENABLE_TCMALLOC TRUE CACHE BOOL "Use TCMalloc")
if(ENABLE_TCMALLOC)
  set(TCMALLOC_LIB "${HOME}" CACHE STRING "Set the TCMalloc library")

  target_link_libraries(msmapp ${TCMALLOC_LIB})
endif(ENABLE_TCMALLOC)

# Profiling option (to allow pausing and resuming profiling during execution)
set(ENABLE_PROFILE FALSE CACHE BOOL "Enable VTune Profiling")
if(ENABLE_PROFILE)
  set(VTUNE_INCLUDE_DIR "${HOME}" CACHE PATH "Set the locaiton of VTune include")
  set(VTUNE_LIB "${HOME}" CACHE STRING "Set the VTune library")

  include_directories(${VTUNE_INCLUDE_DIR})

  target_link_libraries(msmapp ${VTUNE_LIB})
endif(ENABLE_PROFILE)


#if (ENABLE_APP_TESTS)
#  add_subdirectory(test)
#endif (ENABLE_APP_TESTS)

endif (ENABLE_MPI)
