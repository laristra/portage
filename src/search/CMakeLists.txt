#-----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#-----------------------------------------------------------------------------~#

# Add header files
set(search_HEADERS
    search.h
    search_simple.h
    PARENT_SCOPE
)

# Add source files
set(search_SOURCES
    search.cc
    search_simple.cc
    PARENT_SCOPE
)

# Add unit test
if(ENABLE_UNIT_TESTS)
  cinch_add_unit(search
    SOURCES search_test.cc
    POLICY MPI
    THREADS 1)
  cinch_add_unit(search_simple_test
    SOURCES search_simple_test.cc
    LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
    POLICY MPI
    THREADS 1)
endif(ENABLE_UNIT_TESTS)

#----------------------------------------------------------------------------~-#
# Formatting options for emacs and vim.
# vim: set tabstop=4 shiftwidth=4 expandtab :
#----------------------------------------------------------------------------~-#
