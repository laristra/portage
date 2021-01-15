#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#

project(portage_search)

add_library(portage_search)

# Add header files
set(portage_search_HEADERS
    search_simple.h
    search_direct_product.h
    search_kdtree.h
    kdtree.h
    BoundBox.h
    pile.hh
    lretypes.hh
    pairs.hh
    search_swept_face.h
    search_simple_points.h
    search_points_bins.h
    search_points_by_cells.h)

set(portage_search_SOURCES pairs.cc)

target_include_directories(portage_search PUBLIC
    $<BUILD_INTERFACE:${portage_SOURCE_DIR}>
    $<BUILD_INTERFACE:${portage_BINARY_DIR}>
    $<INSTALL_INTERFACE:include>)

target_sources(portage_search PRIVATE "${portage_search_SOURCES}")

set_target_properties(portage_search PROPERTIES
    PUBLIC_HEADER "${portage_search_HEADERS}")

target_link_libraries(portage_search PUBLIC portage_support)
target_link_libraries(portage INTERFACE portage_search)

install(TARGETS portage_search
  EXPORT portage_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/portage/search
  INCLUDES DESTINATION include/portage/search
  )

# Add unit test
if(ENABLE_UNIT_TESTS)
   portage_add_unittest(test_search_direct_product
    SOURCES test/test_search_direct_product.cc
    LIBRARIES portage_search 
    POLICY SERIAL)
  
   portage_add_unittest(test_search_simple
    SOURCES test/test_search_simple.cc
    LIBRARIES portage_search  
    POLICY SERIAL)
  
  portage_add_unittest(test_search_kdtree2
    SOURCES test/test_search_kdtree2.cc
    LIBRARIES portage_search  
    POLICY SERIAL)

  portage_add_unittest(test_search_kdtree3
    SOURCES test/test_search_kdtree3.cc
    LIBRARIES portage_search 
    POLICY SERIAL)

  portage_add_unittest(test_search_swept_face
    SOURCES test/test_search_swept_face.cc
    LIBRARIES portage_search 
    POLICY SERIAL)

  portage_add_unittest(test_search_simple_points
    SOURCES test/test_search_simple_points.cc
    LIBRARIES portage_search
    POLICY SERIAL)

  portage_add_unittest(test_search_points_bins
    SOURCES test/test_search_points_bins.cc
    LIBRARIES portage_search
    POLICY SERIAL)

  portage_add_unittest(test_search_points_by_cells
    SOURCES pairs.cc test/test_search_points_by_cells.cc
    LIBRARIES portage_search
    POLICY SERIAL)

  if (nanoflann_FOUND)
    portage_add_unittest(test_search_kdtree_nanoflann
      SOURCES test/test_search_kdtree_nanoflann.cc
      LIBRARIES portage_search 
      POLICY SERIAL)
  endif (nanoflann_FOUND)


endif(ENABLE_UNIT_TESTS)

