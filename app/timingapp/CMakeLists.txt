#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

if (Jali_DIR AND ENABLE_MPI)

  # options
  set(ENABLE_TCMALLOC TRUE CACHE BOOL "Use TCMalloc")
  set(ENABLE_PROFILE FALSE CACHE BOOL "Enable VTune Profiling")

  add_executable(timingapp timingapp.cc)
  set_target_properties(timingapp
          PROPERTIES
          LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

  target_link_libraries(timingapp portage ${LAPACKX_LIBRARIES} ${EXTRA_LIBS})
  target_link_libraries(timingapp ${MPI_CXX_LIBRARIES})
  target_link_libraries(timingapp ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES})

  if (ENABLE_APP_TIMINGS)
    target_compile_definitions(timingapp PRIVATE ENABLE_TIMINGS=1)
  endif()

  if(ENABLE_TCMALLOC)
    set(TCMALLOC_LIB "${HOME}" CACHE STRING "Set the TCMalloc library")
    target_link_libraries(timingapp ${TCMALLOC_LIB})
  endif()

  if(ENABLE_PROFILE)
    set(VTUNE_INCLUDE_DIR "${HOME}" CACHE PATH "Set the locaiton of VTune include")
    set(VTUNE_LIB "${HOME}" CACHE STRING "Set the VTune library")

    include_directories(${VTUNE_INCLUDE_DIR})
    target_link_libraries(timingapp ${VTUNE_LIB})
  endif()
endif()
