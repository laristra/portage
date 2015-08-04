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

if(ENABLE_UNIT_TESTS)
    cinch_add_unit(test_state_vector
        SOURCES test/test_state_vector.cc
        LIBRARIES ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
        POLICY MPI
        THREADS 2)
endif(ENABLE_UNIT_TESTS)


#----------------------------------------------------------------------------~-#
# Formatting options for emacs and vim.
# vim: set tabstop=4 shiftwidth=4 expandtab :
#----------------------------------------------------------------------------~-#
