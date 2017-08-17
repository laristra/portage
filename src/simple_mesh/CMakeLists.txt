



# Add header files
set (simple_mesh_HEADERS
  simple_mesh.h
  simple_state.h
  PARENT_SCOPE
  )

# Add source files
set (simple_mesh_SOURCES
  PARENT_SCOPE
  )

# Unit tests

if (ENABLE_UNIT_TESTS)

  cinch_add_unit(test_simple_mesh
    SOURCES test/test_simple_mesh.cc
    POLICY MPI
    THREADS 1
    )

  cinch_add_unit(test_simple_state
    SOURCES test/test_simple_state.cc
    POLICY MPI
    THREADS 1
    )

endif (ENABLE_UNIT_TESTS)
