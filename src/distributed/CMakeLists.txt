#-----------------------------------------------------------------------------~#
# Copyright (c) 2015 Los Alamos National Security, LLC
# All rights reserved.
#-----------------------------------------------------------------------------~#

set(distributed_HEADERS
    mpi_bounding_boxes.h
    PARENT_SCOPE)

if (ENABLE_UNIT_TESTS)
  cinch_add_unit(test_mpi_bounding_boxes
                 SOURCES test/test_mpi_bounding_boxes.cc
                 LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
                 POLICY MPI
                 THREADS 1)
endif(ENABLE_UNIT_TESTS)

