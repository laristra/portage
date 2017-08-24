#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]


# Add header files
set (swarm_HEADERS
  swarm.h
  swarm_state.h
  PARENT_SCOPE
  )

# Add source files
set (simple_mesh_SOURCES
  PARENT_SCOPE
  )

# Unit tests

if (ENABLE_UNIT_TESTS)

   if (Jali_DIR)
      cinch_add_unit(test_swarm
      SOURCES test/test_swarm.cc
      LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
      POLICY MPI
      THREADS 1
    )
   endif (Jali_DIR)

  cinch_add_unit(test_swarm_state
    SOURCES test/test_swarm_state.cc
    LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
    POLICY MPI
    THREADS 1
    )

endif (ENABLE_UNIT_TESTS)
