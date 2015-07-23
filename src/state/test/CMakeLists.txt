#-----------------------------------------------------------------------------~#
# Copyright (c) 2015 Los Alamos National Security, LLC
# All rights reserved.
#-----------------------------------------------------------------------------~#

# Add test
if (ENABLE_UNIT_TESTS)
  cinch_add_unit(test_state_vector test_state_vector.cc)

  message(STATUS "Jali_LIBRARIES ${Jali_LIBRARIES}")
  message(STATUS "Jali_TPL_LIBRARIES ${Jali_TPL_LIBRARIES}")
  
  target_link_libraries(test_state_vector ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES})
endif()

#----------------------------------------------------------------------------~-#
# Formatting options for emacs and vim.
# vim: set tabstop=4 shiftwidth=4 expandtab :
#----------------------------------------------------------------------------~-#
