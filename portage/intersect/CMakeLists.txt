#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#

project(portage_intersect)

add_library(portage_intersect INTERFACE)

target_include_directories(portage_intersect INTERFACE
  $<BUILD_INTERFACE:${portage_BINARY_DIR}>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIRECTORY}>
  $<INSTALL_INTERFACE:include>
  )

# Add header files
set(portage_intersect_HEADERS
        intersect_boxes.h
        intersect_polys_r2d.h
        intersect_r2d.h
        intersect_polys_r3d.h
        intersect_r3d.h
        intersect_rNd.h
        intersect_swept_face.h
        dummy_interface_reconstructor.h
        )


# Not yet allowed for INTERFACE libraries
# 
# set_target_properties(portage_intersect PROPERTIES
#   PUBLIC_HEADER "${portage_intersect_HEADERS}")
#
# Directly install files instead
#
install(FILES ${portage_intersect_HEADERS} DESTINATION include/portage/intersect)

target_link_libraries(portage_intersect INTERFACE portage_support)

target_link_libraries(portage INTERFACE portage_intersect)

install(TARGETS portage_intersect
  EXPORT portage_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/portage/intersect
  INCLUDES DESTINATION include/portage/intersect
  )

      

# Add unit test
if (ENABLE_UNIT_TESTS)
  if (WONTON_ENABLE_Jali)
    portage_add_unittest(test_intersect_swept_face_2D
      SOURCES test/test_intersect_swept_face_2D.cc
      LIBRARIES portage_intersect
      POLICY MPI)

    portage_add_unittest(test_intersect_swept_face_3d
      SOURCES test/test_intersect_swept_face_3D.cc
      LIBRARIES portage_intersect
      POLICY MPI)
  endif ()

  portage_add_unittest(test_intersect_boxes
    SOURCES test/test_intersect_boxes.cc
    LIBRARIES portage_intersect
    POLICY SERIAL)
  
  portage_add_unittest(test_intersect_r2d
    SOURCES test/test_intersect_r2d.cc
    LIBRARIES portage_intersect
    POLICY SERIAL)
  
  portage_add_unittest(test_intersect_r3d
    SOURCES test/test_intersect_r3d.cc
    LIBRARIES portage_intersect
    POLICY SERIAL)
  
  if (PORTAGE_HAS_TANGRAM)
      portage_add_unittest(test_intersect_tangram_2d
        SOURCES test/test_tangram_intersect_2D.cc
        LIBRARIES portage_intersect
        POLICY SERIAL)
      
      portage_add_unittest(test_intersect_tangram_3d
        SOURCES test/test_tangram_intersect_3D.cc
        LIBRARIES portage_intersect
        POLICY SERIAL)
  endif ()

endif (ENABLE_UNIT_TESTS)

