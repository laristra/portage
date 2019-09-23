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

option(ENABLE_APP_TIMINGS "Enable timing in apps" ON)

if (Jali_DIR OR ENABLE_FleCSI)
  add_subdirectory(portageapp)
endif (Jali_DIR OR ENABLE_FleCSI)

if (Jali_DIR)
  add_subdirectory(timingapp)
  add_subdirectory(part-by-part)
endif (Jali_DIR)

if (ENABLE_FleCSI)
  add_subdirectory(flecsiapp)
endif (ENABLE_FleCSI)

if (Jali_DIR AND ENABLE_APP_TIMINGS)
  target_compile_definitions(portageapp_jali PRIVATE ENABLE_TIMINGS=1)
  target_compile_definitions(portageapp_t-junction_jali PRIVATE ENABLE_TIMINGS=1)
  target_compile_definitions(portageapp_rgmd_jali PRIVATE ENABLE_TIMINGS=1)
  target_compile_definitions(portageapp_multimat_jali PRIVATE ENABLE_TIMINGS=1)
  target_compile_definitions(timingapp PRIVATE ENABLE_TIMINGS=1)
endif()