#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]


# Add header files
set (accumulate_HEADERS
  accumulate.h
  PARENT_SCOPE
  )

# Add source files
set (accumulate_SOURCES
  PARENT_SCOPE
  )

# Unit tests

if (ENABLE_UNIT_TESTS)

  cinch_add_unit(test_accumulate
    SOURCES test/test_accumulate.cc
    LIBRARIES portage 
    POLICY MPI
    THREADS 1
    )

endif (ENABLE_UNIT_TESTS)
