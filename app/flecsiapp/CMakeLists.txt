#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]


if (ENABLE_FleCSI)
    add_executable(flecsiapp flecsiapp.cc)
    target_link_libraries(flecsiapp portage ${EXTRA_LIBS} ${LAPACKX_LIBRARIES})

    # FleCSI/FleCSISP additions
    include_directories(${FleCSI_INCLUDE_DIR} ${FleCSISP_INCLUDE_DIR})
    target_link_libraries(flecsiapp ${FleCSI_LIBRARIES} ${FleCSISP_LIBRARIES})
    
    # MPI additions - need to push this down into cinch
    if (ENABLE_MPI)
        target_link_libraries(flecsiapp ${MPI_CXX_LIBRARIES})
        set_target_properties(flecsiapp PROPERTIES
                              LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
    endif (ENABLE_MPI)

endif (ENABLE_FleCSI)


#TCMalloc option
set(ENABLE_TCMALLOC TRUE CACHE BOOL "Use TCMalloc")
if(ENABLE_TCMALLOC)
  set(TCMALLOC_LIB "${HOME}" CACHE STRING "Set the TCMalloc library")

  if (ENABLE_FleCSI)
    target_link_libraries(flecsiapp ${TCMALLOC_LIB})
  endif (ENABLE_FleCSI)
endif(ENABLE_TCMALLOC)

# Profiling option (to allow pausing and resuming profiling during execution)
set(ENABLE_PROFILE FALSE CACHE BOOL "Enable VTune Profiling")
if(ENABLE_PROFILE)
  set(VTUNE_INCLUDE_DIR "${HOME}" CACHE PATH "Set the locaiton of VTune include")
  set(VTUNE_LIB "${HOME}" CACHE STRING "Set the VTune library")

  include_directories(${VTUNE_INCLUDE_DIR})

  if (ENABLE_FleCSI)
    target_link_libraries(flecsiapp ${VTUNE_LIB})
  endif (ENABLE_FleCSI)
endif(ENABLE_PROFILE)

#if (ENABLE_APP_TESTS)
#  add_subdirectory(test)
#endif (ENABLE_APP_TESTS)
