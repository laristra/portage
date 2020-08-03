#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#

project(portage_distributed)

add_library(portage_distributed INTERFACE)

target_include_directories(portage_distributed INTERFACE
  $<BUILD_INTERFACE:${portage_SOURCE_DIR}>
  $<BUILD_INTERFACE:${portage_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>
  )

# Add header files
set(portage_distributed_HEADERS
    mpi_bounding_boxes.h
    mpi_particle_distribute.h
    mpi_ghost_manager.h)

# Not yet allowed for INTERFACE libraries
# 
# set_target_properties(portage_distributed PROPERTIES
#   PUBLIC_HEADER "${portage_distributed_HEADERS}")
#
# Directly install files instead
#
install(FILES ${portage_distributed_HEADERS} DESTINATION include/portage/distributed)

target_link_libraries(portage_distributed INTERFACE portage_support)

target_link_libraries(portage INTERFACE portage_distributed)

install(TARGETS portage_distributed
  EXPORT portage_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/portage/distributed
  INCLUDES DESTINATION include/portage/distributed
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

    portage_add_unittest(test_mpi_ghost_manager
               SOURCES test/test_mpi_ghost_manager.cc
               LIBRARIES portage_distribute
               POLICY MPI
               THREADS 4)
  endif ()
endif(ENABLE_UNIT_TESTS)

