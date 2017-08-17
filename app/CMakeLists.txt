#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]




add_subdirectory(apptest_cmp)
add_subdirectory(simple_mesh_app)
add_subdirectory(swarmapp)
if (Jali_DIR OR FLECSI_DIR)
  add_subdirectory(portageapp)
endif (Jali_DIR OR FLECSI_DIR)
if (Jali_DIR)
  add_subdirectory(shotshellapp)
  add_subdirectory(timingapp)
endif (Jali_DIR)

