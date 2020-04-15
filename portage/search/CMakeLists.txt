#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#


#-----------------------------------------------------------------------------~#

# Add header files
set(search_HEADERS
    search_simple.h
    search_direct_product.h
    search_kdtree.h
    search_simple_points.h
    search_points_by_cells.h
    kdtree.h
    pile.hh
    lretypes.hh
    pairs.hh
    BoundBox.h
    search_swept_face.h
    PARENT_SCOPE
)

if (nanoflann_FOUND)
  list(APPEND search_HEADERS search_kdtree_nanoflann.h)
endif (nanoflann_FOUND)

# Add source files
set(search_SOURCES
    pairs.cc
    PARENT_SCOPE
)

# Add unit test
if(ENABLE_UNIT_TESTS)
   cinch_add_unit(test_search_direct_product
    SOURCES test/test_search_direct_product.cc
    LIBRARIES portage  
    POLICY SERIAL)
  
   cinch_add_unit(test_search_simple
    SOURCES test/test_search_simple.cc
    LIBRARIES portage  
    POLICY SERIAL)
  
  cinch_add_unit(test_search_kdtree2
    SOURCES test/test_search_kdtree2.cc
    LIBRARIES portage  
    POLICY SERIAL)

  cinch_add_unit(test_search_kdtree3
    SOURCES test/test_search_kdtree3.cc
    LIBRARIES portage 
    POLICY SERIAL)

  cinch_add_unit(test_search_swept_face
    SOURCES test/test_search_swept_face.cc
    LIBRARIES portage 
    POLICY SERIAL)

  cinch_add_unit(test_search_simple_points
    SOURCES test/test_search_simple_points.cc
    LIBRARIES portage  
    POLICY MPI
    THREADS 1)

  cinch_add_unit(test_search_points_by_cells
    SOURCES pairs.cc test/test_search_points_by_cells.cc
    LIBRARIES portage 
    POLICY MPI
    THREADS 1)

  if (nanoflann_FOUND)
    cinch_add_unit(test_search_kdtree_nanoflann
      SOURCES test/test_search_kdtree_nanoflann.cc
      LIBRARIES portage 
      POLICY MPI
      THREADS 1)
  endif (nanoflann_FOUND)

endif(ENABLE_UNIT_TESTS)

