#-----------------------------------------------------------------------------~#
# Copyright (c) 2015 Los Alamos National Security, LLC
# All rights reserved.
#-----------------------------------------------------------------------------~#

# Add header files
set(remap_HEADERS
    remap_1st_order.h
    remap_2nd_order.h
    gradient.h
    PARENT_SCOPE
)

# Add source files
#set(remap_SOURCES
#    remap_1st_order.cc
#    PARENT_SCOPE
#)

# Unit tests

if (ENABLE_UNIT_TESTS)

  cinch_add_unit(test_gradient
    SOURCES  test/test_gradient.cc
    LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
    POLICY MPI
    THREADS 1)

  cinch_add_unit(test_remap_1st_order
    SOURCES  test/test_remap_1st_order.cc
    LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
    POLICY MPI
    THREADS 1)

  cinch_add_unit(test_remap_2nd_order
    SOURCES test/test_remap_2nd_order.cc
    LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
    POLICY MPI
    THREADS 1)

endif(ENABLE_UNIT_TESTS)

