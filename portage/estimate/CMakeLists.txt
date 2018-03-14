#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]


# Add header files
set (estimate_HEADERS
  estimate.h
  PARENT_SCOPE
  )

# Add source files
set (estimator_SOURCES
  PARENT_SCOPE
  )

# Unit tests

if (ENABLE_UNIT_TESTS)

  cinch_add_unit(test_estimate
    SOURCES test/test_estimate.cc
    LIBRARIES portage
    POLICY MPI
    THREADS 1
    )

endif (ENABLE_UNIT_TESTS)
