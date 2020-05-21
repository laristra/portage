#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]

# put this with individual targets
include_directories(include)  # contains user_field.h for expression parsing

add_subdirectory(apptest_cmp)
add_subdirectory(multidef)
add_subdirectory(distributed_cmp)
add_subdirectory(simple_mesh_app)

if (WONTON_ENABLE_KOKKOS)
  add_subdirectory(swarmapp)
  add_subdirectory(msmapp)
endif()

add_subdirectory(momentumapp)

option(ENABLE_APP_TIMINGS "Enable timing in apps" OFF)

if (WONTON_ENABLE_Jali OR WONTON_ENABLE_FleCSI)
  add_subdirectory(portageapp)
endif ()

if (WONTON_ENABLE_Jali)
  add_subdirectory(timingapp)
  add_subdirectory(part-by-part)
  add_subdirectory(swept_face)
endif ()

if (WONTON_ENABLE_FleCSI)
  add_subdirectory(flecsiapp)
endif ()
