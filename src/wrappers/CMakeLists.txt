#-----------------------------------------------------------------------------~#
# Copyright (c) 2015 Los Alamos National Security, LLC
# All rights reserved.
#-----------------------------------------------------------------------------~#

set(wrappers_HEADERS
    mesh/jali/jali_mesh_wrapper.h
    state/jali/jali_state_wrapper.h
    PARENT_SCOPE)

if (ENABLE_UNIT_TESTS)
  cinch_add_unit(test_jali_mesh_wrapper
                 SOURCES mesh/jali/test/test_jali_mesh_wrapper.cc
                 LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
                 POLICY MPI
                 THREADS 1)

  cinch_add_unit(test_jali_state_wrapper
                 SOURCES state/jali/test/test_jali_state_wrapper.cc
                 LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
                 POLICY MPI
                 THREADS 1)

endif(ENABLE_UNIT_TESTS)

#----------------------------------------------------------------------------~-#
# Formatting options for emacs and vim.
# vim: set tabstop=4 shiftwidth=4 expandtab :
#----------------------------------------------------------------------------~-#
