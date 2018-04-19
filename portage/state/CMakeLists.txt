#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]


# Add header files
set (state_HEADERS
  state_vector_base.h
  state_vector_uni.h
  state_vector_multi.h
  state_manager.h
  PARENT_SCOPE
  )

# Add source files
set (state_SOURCES
  PARENT_SCOPE
  )

# Unit tests

if (ENABLE_UNIT_TESTS)

  cinch_add_unit(test_state_vector_uni
    SOURCES test/test_state_vector_uni.cc
    POLICY SERIAL
    )

  cinch_add_unit(test_state_vector_multi
    SOURCES test/test_state_vector_multi.cc
    POLICY SERIAL
    )

  cinch_add_unit(test_state_manager
    SOURCES test/test_state_manager.cc
    POLICY SERIAL
    )

endif (ENABLE_UNIT_TESTS)
