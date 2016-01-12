#-----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#-----------------------------------------------------------------------------~#

# Add header files
set(search_HEADERS
    search_simple.h
    search_kdtree.h
    search_kdtree3.h
    kdtree.h
    BoundBox.h
    Point.h
    Vector.h
    PARENT_SCOPE
)

# Add source files
#set(search_SOURCES
#    search_simple.cc
#    PARENT_SCOPE
#)

# Add unit test
if(ENABLE_UNIT_TESTS)
  cinch_add_unit(search_simple_test
    SOURCES search_simple_test.cc
    LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
    POLICY MPI
    THREADS 1)
  cinch_add_unit(search_kdtree_test
    SOURCES search_kdtree_test.cc
    LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
    POLICY MPI
    THREADS 1)
  cinch_add_unit(search_kdtree3_test
    SOURCES search_kdtree3_test.cc
    LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
    POLICY MPI
    THREADS 1)
endif(ENABLE_UNIT_TESTS)

