#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#

project(portage_distribute)

add_library(portage_distribute INTERFACE)

target_include_directories(portage_distribute INTERFACE
  $<BUILD_INTERFACE:${portage_SOURCE_DIR}>
  $<BUILD_INTERFACE:${portage_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>
  )

# Add header files
set(portage_distribute_HEADERS
    mpi_bounding_boxes.h
    mpi_particle_distribute.h
    )

# Not yet allowed for INTERFACE libraries
# 
# set_target_properties(portage_distribute PROPERTIES
#   PUBLIC_HEADER "${portage_distribute_HEADERS}")
#
# Directly install files instead
#
install(FILES ${portage_distribute_HEADERS} DESTINATION include/portage/distributed)

target_link_libraries(portage_distribute INTERFACE portage_support)

target_link_libraries(portage INTERFACE portage_distribute)

install(TARGETS portage_distribute
  EXPORT portage_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/portage/distribute
  INCLUDES DESTINATION include/portage/distribute
  )



if (ENABLE_UNIT_TESTS)
  if (WONTON_ENABLE_Jali)
    portage_add_unittest(test_mpi_bounding_boxes
                   SOURCES test/test_mpi_bounding_boxes.cc
                   LIBRARIES portage_distributed
                   POLICY MPI
                   THREADS 4)

    portage_add_unittest(test_mpi_particle_distribute
                   SOURCES test/test_mpi_particle_distribute.cc
                   LIBRARIES portage_distributed
                   POLICY MPI
                   THREADS 4)
  endif ()
endif(ENABLE_UNIT_TESTS)

