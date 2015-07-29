#-----------------------------------------------------------------------------~#
# Copyright (c) 2015 Los Alamos National Security, LLC
# All rights reserved.
#-----------------------------------------------------------------------------~#

# Add header files
set(state_HEADERS
    state_vector.h
    state.h
    PARENT_SCOPE
)

# Add source files
set(state_SOURCES
    state.cc
    PARENT_SCOPE
)


add_subdirectory(test)


#----------------------------------------------------------------------------~-#
# Formatting options for emacs and vim.
# vim: set tabstop=4 shiftwidth=4 expandtab :
#----------------------------------------------------------------------------~-#
