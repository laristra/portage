message(STATUS "Adding test directory")

add_subdirectory(apptest_cmp)
add_subdirectory(simple_mesh_app)
if (Jali_DIR)
  add_subdirectory(portageapp)
  add_subdirectory(shotshellapp)
endif (Jali_DIR)
