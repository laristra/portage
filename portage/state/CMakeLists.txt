#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]


# Add header files
set (state_HEADERS
  state_vector.h
  PARENT_SCOPE
  )

# Add source files
set (state_SOURCES
  PARENT_SCOPE
  )

# Unit tests

if (ENABLE_UNIT_TESTS)

  cinch_add_unit(test_state_vector
    SOURCES test/test_state_vector.cc
    POLICY SERIAL
    )

endif (ENABLE_UNIT_TESTS)
