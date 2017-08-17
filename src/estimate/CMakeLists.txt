



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
    POLICY MPI
    THREADS 1
    )

endif (ENABLE_UNIT_TESTS)
