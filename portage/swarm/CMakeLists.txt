#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]

project(portage_swarm)

add_library(portage_swarm INTERFACE)

target_include_directories(portage_swarm INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIRECTORY}>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIRECTORY}>
  $<INSTALL_INTERFACE:include>
  )

# Add header files
set (portage_swarm_HEADERS
  swarm.h
  swarm_state.h
  )

# Not yet allowed for INTERFACE libraries
# 
# set_target_properties(portage_swarm PROPERTIES
#   PUBLIC_HEADER "${portage_swarm_HEADERS}")
#
# Directly install files instead
#
install(FILES ${portage_swarm_HEADERS} DESTINATION include/portage/swarm)

target_link_libraries(portage_swarm INTERFACE portage_support)

target_link_libraries(portage INTERFACE portage_swarm)

install(TARGETS portage_swarm
  EXPORT portage_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/portage/swarm
  INCLUDES DESTINATION include/portage/swarm
  )


# Unit tests

if (ENABLE_UNIT_TESTS)
  cinch_add_unit(test_swarm
    SOURCES test/test_swarm.cc
    LIBRARIES portage_swarm 
    POLICY MPI
    THREADS 1
    )

  cinch_add_unit(test_swarm_state
    SOURCES test/test_swarm_state.cc
    LIBRARIES portage_swarm 
    POLICY MPI
    THREADS 1
    )

  cinch_add_unit(test_swarm_state_multi
    SOURCES test/test_swarm_state_multi.cc
    LIBRARIES portage_swarm 
    POLICY MPI
    THREADS 1
    )
endif (ENABLE_UNIT_TESTS)
