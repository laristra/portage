#-----------------------------------------------------------------------------~#
# Copyright (c) 2015 Los Alamos National Security, LLC
# All rights reserved.
#-----------------------------------------------------------------------------~#

# Add header files
set(remap_HEADERS
    remap_1st_order.h
    PARENT_SCOPE
)

# Add source files
#set(remap_SOURCES
#    remap_1st_order.cc
#    PARENT_SCOPE
#)

# Unit tests

if (ENABLE_UNIT_TESTS) 

  cinch_add_unit(test_remap_1st_order
    SOURCES  test/test_remap_1st_order.cc
    LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
    POLICY MPI
    THREADS 1)

endif(ENABLE_UNIT_TESTS)

#----------------------------------------------------------------------------~-#
# Formatting options for emacs and vim.
# vim: set tabstop=4 shiftwidth=4 expandtab :
#----------------------------------------------------------------------------~-#
