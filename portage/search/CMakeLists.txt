#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#

project(portage_search)

# Add header files
set(search_HEADERS
    search_simple.h
    search_direct_product.h
    search_kdtree.h
    kdtree.h
    BoundBox.h
    search_swept_face.h)

if (WONTON_ENABLE_KOKKOS)
  add_library(portage_search)

  set(search_HEADERS
      ${search_HEADERS}
      search_simple_points.h
      search_points_by_cells.h
      pile.hh
      lretypes.hh
      pairs.hh)

  set(search_SOURCES pairs.cc)

  target_include_directories(portage_search PUBLIC
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIRECTORY}>
      $<BUILD_INTERFACE:${PROJECT_BINARY_DIRECTORY}>
      $<INSTALL_INTERFACE:include>)

  target_sources(portage_search PRIVATE "${search_SOURCES}")

else()
  add_library(portage_search INTERFACE)

  target_include_directories(portage_search INTERFACE
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIRECTORY}>
      $<BUILD_INTERFACE:${PROJECT_BINARY_DIRECTORY}>
      $<INSTALL_INTERFACE:include>)

  target_link_libraries(portage_search INTERFACE portage_support)
endif()

set_target_properties(portage_search PROPERTIES
  PUBLIC_HEADER "${portage_search_HEADERS}")

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
   cinch_add_unit(test_search_direct_product
    SOURCES test/test_search_direct_product.cc
    LIBRARIES portage_search 
    POLICY SERIAL)
  
   cinch_add_unit(test_search_simple
    SOURCES test/test_search_simple.cc
    LIBRARIES portage_search  
    POLICY SERIAL)
  
  cinch_add_unit(test_search_kdtree2
    SOURCES test/test_search_kdtree2.cc
    LIBRARIES portage_search  
    POLICY SERIAL)

  cinch_add_unit(test_search_kdtree3
    SOURCES test/test_search_kdtree3.cc
    LIBRARIES portage_search 
    POLICY SERIAL)

  cinch_add_unit(test_search_swept_face
    SOURCES test/test_search_swept_face.cc
    LIBRARIES portage_search 
    POLICY SERIAL)

  if (WONTON_ENABLE_KOKKOS)
    cinch_add_unit(test_search_simple_points
      SOURCES test/test_search_simple_points.cc
      LIBRARIES portage_search
      POLICY MPI
      THREADS 1)

    cinch_add_unit(test_search_points_by_cells
      SOURCES pairs.cc test/test_search_points_by_cells.cc
      LIBRARIES portage_search
      POLICY MPI
      THREADS 1)
  endif()

  if (nanoflann_FOUND)
    cinch_add_unit(test_search_kdtree_nanoflann
      SOURCES test/test_search_kdtree_nanoflann.cc
      LIBRARIES portage_search 
      POLICY MPI
      THREADS 1)
  endif (nanoflann_FOUND)


endif(ENABLE_UNIT_TESTS)

