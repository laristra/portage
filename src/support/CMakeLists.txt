#-----------------------------------------------------------------------------~#
# Copyright (c) 2015 Los Alamos National Security, LLC
# All rights reserved.
#-----------------------------------------------------------------------------~#

# Add header files
set(support_HEADERS
    portage.h
    matrix.h
    PARENT_SCOPE
)

# Add source files
set(support_SOURCES
    PARENT_SCOPE
)

# Unit tests

if (ENABLE_UNIT_TESTS) 

  cinch_add_unit(test_matfuncs
    SOURCES test/test_matfuncs.cc
    POLICY MPI
    THREADS 1)

endif(ENABLE_UNIT_TESTS)
