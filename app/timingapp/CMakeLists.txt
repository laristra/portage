#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]





#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

if (Jali_DIR) 
# Add a rule to build the executable
# timing example

  add_executable(timingapp timingapp.cc)
  target_link_libraries(timingapp portage ${EXTRA_LIBS})
  if (ENABLE_MPI)
    target_link_libraries(timingapp ${MPI_CXX_LIBRARIES})
    set_target_properties(timingapp PROPERTIES
                              LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
  endif (ENABLE_MPI)
  target_link_libraries(timingapp ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES})


#TCMalloc option
  set(ENABLE_TCMALLOC TRUE CACHE BOOL "Use TCMalloc")
  if(ENABLE_TCMALLOC)
    set(TCMALLOC_LIB "${HOME}" CACHE STRING "Set the TCMalloc library")
    target_link_libraries(timingapp ${TCMALLOC_LIB})
  endif(ENABLE_TCMALLOC)

# Profiling option (to allow pausing and resuming profiling during execution)
  set(ENABLE_PROFILE FALSE CACHE BOOL "Enable VTune Profiling")
  if(ENABLE_PROFILE)
    set(VTUNE_INCLUDE_DIR "${HOME}" CACHE PATH "Set the locaiton of VTune include")
    set(VTUNE_LIB "${HOME}" CACHE STRING "Set the VTune library")

    include_directories(${VTUNE_INCLUDE_DIR})

    target_link_libraries(timingapp ${VTUNE_LIB})
  endif(ENABLE_PROFILE)
endif (Jali_DIR)
