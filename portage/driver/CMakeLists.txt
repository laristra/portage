#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#

project(portage_driver)

add_library(portage_driver INTERFACE)

target_include_directories(portage_driver INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIRECTORY}>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIRECTORY}>
  $<INSTALL_INTERFACE:include>
  )

set(portage_driver_HEADERS
  mmdriver.h
  driver_swarm.h
  driver_mesh_swarm_mesh.h
  fix_mismatch.h
  coredriver.h
  uberdriver.h
  parts.h)
if (PORTAGE_HAS_TANGRAM)
  list(APPEND portage_driver_headers
    write_to_gmv.h
    )
endif ()


# Not yet allowed for INTERFACE libraries
# 
# set_target_properties(portage_driver PROPERTIES
#   PUBLIC_HEADER "${portage_driver_HEADERS}")
#
# Directly install files instead
#
install(FILES ${portage_driver_HEADERS} DESTINATION include/portage/driver)

target_link_libraries(portage_driver INTERFACE portage_support)
target_link_libraries(portage_driver INTERFACE portage_search)
target_link_libraries(portage_driver INTERFACE portage_distribute)
target_link_libraries(portage_driver INTERFACE portage_accumulate)
target_link_libraries(portage_driver INTERFACE portage_estimate)
target_link_libraries(portage_driver INTERFACE portage_intersect)
target_link_libraries(portage_driver INTERFACE portage_interpolate)
if (PORTAGE_HAS_TANGRAM)
  target_link_libraries(portage_driver INTERFACE ${TANGRAM_LIBRARIES})
endif ()

target_link_libraries(portage INTERFACE portage_driver)

install(TARGETS portage_driver
  EXPORT portage_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/portage/driver
  INCLUDES DESTINATION include/portage/driver
  )



# Unit tests

if (ENABLE_UNIT_TESTS)

   if (Jali_DIR)
       cinch_add_unit(test_driver
       SOURCES test/test_driver.cc 
       LIBRARIES portage_driver
       POLICY MPI
       THREADS 1)

       cinch_add_unit(test_driver_multivar
       SOURCES test/test_driver_multivar.cc 
       LIBRARIES portage_driver
       POLICY MPI
       THREADS 1)

       cinch_add_unit(test_mismatch_fixup
       SOURCES test/test_driver_mismatch_fixup.cc 
       LIBRARIES portage_driver
       POLICY MPI
       THREADS 1)

       cinch_add_unit(test_mismatch_disjoint
       SOURCES test/test_driver_mismatch_disjoint.cc
       LIBRARIES portage_driver
       POLICY MPI
       THREADS 1)

       cinch_add_unit(test_driver_core_mismatch_fixup
       SOURCES test/test_driver_core_mismatch_fixup 
       LIBRARIES portage_driver
       POLICY MPI
       THREADS 1)

       cinch_add_unit(test_driver_swarm_distributed
       SOURCES test/test_driver_swarm_distributed.cc
       LIBRARIES portage_driver
       POLICY MPI
       THREADS 4)

     if (PORTAGE_HAS_TANGRAM)

	 cinch_add_unit(test_driver_multimat
	   SOURCES test/test_driver_multimat.cc
	   LIBRARIES portage_driver ${TANGRAM_LIBRARIES}
	   POLICY MPI
	   THREADS 1)
	 
	 cinch_add_unit(test_driver_multimat_second_order
	   SOURCES test/test_driver_multimat_second_order.cc
	   LIBRARIES portage_driver ${TANGRAM_LIBRARIES}
	   POLICY MPI
	   THREADS 1)
	 
	 cinch_add_unit(test_driver_multimat_uber
	   SOURCES test/test_driver_multimat_uber.cc
	   LIBRARIES portage_driver ${TANGRAM_LIBRARIES}
	   POLICY MPI
	   THREADS 1)
	 
	 cinch_add_unit(test_driver_multimat_uber_swept
	   SOURCES test/test_driver_multimat_uber_swept.cc
	   LIBRARIES portage_driver ${TANGRAM_LIBRARIES}
	   POLICY MPI
	   THREADS 1)

       cinch_add_unit(test_driver_multimat_distributed
       SOURCES test/test_driver_multimat_distributed.cc
       LIBRARIES portage_driver ${TANGRAM_LIBRARIES}
       POLICY MPI
       THREADS 4)

       cinch_add_unit(test_driver_multimat_uber_distributed
       SOURCES test/test_driver_multimat_uber_distributed
       LIBRARIES portage_driver ${TANGRAM_LIBRARIES}
       POLICY MPI
       THREADS 4)

     endif ()

     cinch_add_unit(test_driver_core_single_mat
       SOURCES test/test_driver_core_single_mat.cc
       LIBRARIES portage_driver
       POLICY MPI
       THREADS 1)

     cinch_add_unit(test_driver_part
       SOURCES test/test_driver_part.cc
       LIBRARIES portage_driver
       POLICY MPI
       THREADS 1)

     cinch_add_unit(test_driver_part_mismatch
       SOURCES test/test_driver_part_mismatch.cc
       LIBRARIES portage_driver
       POLICY MPI
       THREADS 1)
     
   endif ()

   cinch_add_unit(test_driver_swarm
     SOURCES test/test_driver_swarm.cc
     LIBRARIES portage_driver
     POLICY MPI
     THREADS 1)

   cinch_add_unit(test_driver_mesh_swarm_mesh
     SOURCES test/test_driver_mesh_swarm_mesh.cc
     LIBRARIES portage_driver
     POLICY MPI
     THREADS 1)

   cinch_add_unit(test_driver_part_by_particle
     SOURCES test/test_driver_part_by_particle
     LIBRARIES portage_driver
     POLICY MPI
     THREADS 1)

endif (ENABLE_UNIT_TESTS)
