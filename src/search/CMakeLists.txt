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
    search_kdtree.h
    search_simple_points.h
    search_points_by_cells.h
    kdtree.h
    BoundBox.h
    PARENT_SCOPE
)

# Add source files
set(search_SOURCES
    pairs.cc
    PARENT_SCOPE
)

# Add unit test
if(ENABLE_UNIT_TESTS)
  if (Jali_DIR)
    cinch_add_unit(search_simple_test
      SOURCES search_simple_test.cc
      LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
      POLICY MPI
      THREADS 1)
    cinch_add_unit(search_kdtree2_test
      SOURCES search_kdtree2_test.cc
      LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
      POLICY MPI
      THREADS 1)
    cinch_add_unit(search_kdtree3_test
      SOURCES search_kdtree3_test.cc
      LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
      POLICY MPI
      THREADS 1)
  endif (Jali_DIR)
  cinch_add_unit(search_simple_points_test
    SOURCES search_simple_points_test.cc
    LIBRARIES portage 
    POLICY MPI
    THREADS 1)
  cinch_add_unit(search_points_by_cells_test
    SOURCES pairs.cc search_points_by_cells_test.cc 
    LIBRARIES portage
    POLICY MPI
    THREADS 1)
endif(ENABLE_UNIT_TESTS)

