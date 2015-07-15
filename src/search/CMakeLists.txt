#-----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#-----------------------------------------------------------------------------~#

# Add header files
set(search_HEADERS
    search.h
    PARENT_SCOPE
)

# Add source files
set(search_SOURCES
    search.cc
    PARENT_SCOPE
)

# Add unit test
if(ENABLE_UNIT_TESTS)
  cinch_add_unit(search search_test.cc)
endif(ENABLE_UNIT_TESTS)

#----------------------------------------------------------------------------~-#
# Formatting options for emacs and vim.
# vim: set tabstop=4 shiftwidth=4 expandtab :
#----------------------------------------------------------------------------~-#
