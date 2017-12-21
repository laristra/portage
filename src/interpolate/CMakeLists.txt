#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#


#-----------------------------------------------------------------------------~#

# Add header files
set(interpolate_HEADERS
    interpolate_1st_order.h
    interpolate_2nd_order.h
    interpolate_3rd_order.h
    gradient.h
    quadfit.h
    PARENT_SCOPE
)

# Unit tests

if (ENABLE_UNIT_TESTS)

  if (Jali_DIR)
    cinch_add_unit(test_gradient
      SOURCES  test/test_gradient.cc
      LIBRARIES portage ${LAPACKE_LIBRARIES} ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
      POLICY MPI
      THREADS 1)

    cinch_add_unit(test_quadfit
      SOURCES  test/test_quadfit.cc
      LIBRARIES portage  ${LAPACKE_LIBRARIES} ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
      POLICY MPI
      THREADS 1)

    cinch_add_unit(test_interpolate_1st_order
      SOURCES  test/test_interp_1st_order.cc
      LIBRARIES portage  ${LAPACKE_LIBRARIES} ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
      POLICY MPI
      THREADS 1)

    cinch_add_unit(test_interpolate_2nd_order
      SOURCES test/test_interp_2nd_order.cc
      LIBRARIES portage  ${LAPACKE_LIBRARIES} ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
      POLICY MPI
      THREADS 1)

    cinch_add_unit(test_interpolate_3rd_order
      SOURCES test/test_interp_3rd_order.cc
      LIBRARIES portage  ${LAPACKE_LIBRARIES} ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
      POLICY MPI
      THREADS 1)

  endif (Jali_DIR)

endif(ENABLE_UNIT_TESTS)

