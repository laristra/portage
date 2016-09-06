#-----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#-----------------------------------------------------------------------------~#

# Add header files
set(search_HEADERS
    search_simple.h
    search_kdtree.h
    kdtree.h
    BoundBox.h
    PARENT_SCOPE
)

# Add source files
#set(search_SOURCES
#    search_simple.cc
#    PARENT_SCOPE
#)

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
endif(ENABLE_UNIT_TESTS)

