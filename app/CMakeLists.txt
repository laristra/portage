# Copyright 2016, Los Alamos National Security, LLC.
# All rights reserved

add_subdirectory(apptest_cmp)
add_subdirectory(simple_mesh_app)
if (Jali_DIR OR FLECSI_DIR)
  add_subdirectory(portageapp)
endif (Jali_DIR OR FLECSI_DIR)
if (Jali_DIR)
  add_subdirectory(shotshellapp)
  add_subdirectory(timingapp)
endif (Jali_DIR)

