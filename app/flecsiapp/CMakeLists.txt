#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]


if (WONTON_ENABLE_FleCSI)
    add_executable(flecsiapp flecsiapp.cc)
    target_link_libraries(flecsiapp portage)
endif ()


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
