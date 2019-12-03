#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]

include_directories(include)  # contains user_field.h for expression parsing

add_subdirectory(apptest_cmp)
add_subdirectory(multidef)
add_subdirectory(distributed_cmp)
add_subdirectory(simple_mesh_app)
add_subdirectory(swarmapp)
add_subdirectory(msmapp)

option(ENABLE_APP_TIMINGS "Enable timing in apps" OFF)

if (Jali_DIR OR ENABLE_FleCSI)
  add_subdirectory(portageapp)
endif (Jali_DIR OR ENABLE_FleCSI)

if (Jali_DIR)
  add_subdirectory(timingapp)
  add_subdirectory(part-by-part)
  add_subdirectory(sweptfaceapp)
endif (Jali_DIR)

if (ENABLE_FleCSI)
  add_subdirectory(flecsiapp)
endif (ENABLE_FleCSI)
