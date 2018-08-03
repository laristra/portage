#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]


# Add header files
set (simple_mesh_HEADERS
  simple_mesh.h
  simple_state.h
  SimpleStateVector.h
  SimpleState.h
  PARENT_SCOPE
  )

# Add source files
set (simple_mesh_SOURCES
  SimpleStateVector.cc
  SimpleState.cc
  PARENT_SCOPE
  )

# Unit tests

if (ENABLE_UNIT_TESTS)

  cinch_add_unit(test_simple_mesh
    SOURCES test/test_simple_mesh.cc
    POLICY SERIAL
    )

  cinch_add_unit(test_simple_state
    SOURCES test/test_simple_state.cc
    POLICY SERIAL
    )

endif (ENABLE_UNIT_TESTS)
