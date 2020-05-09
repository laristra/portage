#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

if (WONTON_ENABLE_Jali AND WONTON_ENABLE_MPI)

  # options
  set(ENABLE_PROFILE FALSE CACHE BOOL "Enable VTune Profiling")

  add_executable(timingapp timingapp.cc)
  target_link_libraries(timingapp PRIVATE portage)

  if (ENABLE_APP_TIMINGS)
    target_compile_definitions(timingapp PRIVATE ENABLE_TIMINGS=1)
  endif()

  if(ENABLE_PROFILE)
    set(VTUNE_INCLUDE_DIR "${HOME}" CACHE PATH "Set the locaiton of VTune include")
    set(VTUNE_LIB "${HOME}" CACHE STRING "Set the VTune library")

    target_include_directories(timingapp PRIVATE ${VTUNE_INCLUDE_DIR})
    target_link_libraries(timingapp PRIVATE ${VTUNE_LIB})
  endif()
endif()
