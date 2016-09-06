#-----------------------------------------------------------------------------~#
# Copyright (c) 2015 Los Alamos National Security, LLC
# All rights reserved.
#-----------------------------------------------------------------------------~#

# Add header files
set(interpolate_HEADERS
    interpolate_1st_order.h
    interpolate_2nd_order.h
    gradient.h
    PARENT_SCOPE
)

# Unit tests

if (ENABLE_UNIT_TESTS)

  if (Jali_DIR)
    cinch_add_unit(test_gradient
      SOURCES  test/test_gradient.cc
      LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
      POLICY MPI
      THREADS 1)

    cinch_add_unit(test_interpolate_1st_order
      SOURCES  test/test_interp_1st_order.cc
      LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
      POLICY MPI
      THREADS 1)

    cinch_add_unit(test_interpolate_2nd_order
      SOURCES test/test_interp_2nd_order.cc
      LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
      POLICY MPI
      THREADS 1)
  endif (Jali_DIR)

endif(ENABLE_UNIT_TESTS)

