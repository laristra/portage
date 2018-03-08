#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]

include_directories(include)  # contains user_field.h for expression parsing

add_subdirectory(apptest_cmp)
add_subdirectory(simple_mesh_app)
add_subdirectory(swarmapp)
add_subdirectory(msmapp)
if (Jali_DIR OR ENABLE_FleCSI)
  add_subdirectory(portageapp)
endif (Jali_DIR OR ENABLE_FleCSI)
#if (Jali_DIR)
#  add_subdirectory(timingapp)
#endif (Jali_DIR)
if(TANGRAM_FOUND AND XMOF2D_FOUND)
  add_subdirectory(ir_linetest_app)
endif(TANGRAM_FOUND AND XMOF2D_FOUND)

