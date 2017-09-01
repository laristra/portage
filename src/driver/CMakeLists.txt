#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#


#-----------------------------------------------------------------------------~#

# Add header files
set(driver_HEADERS
    driver.h driver_swarm.h driver_mesh_swarm_mesh.h
    PARENT_SCOPE
)

# Add source files
set(driver_SOURCES
    PARENT_SCOPE
)

# Unit tests

if (ENABLE_UNIT_TESTS)

   if (Jali_DIR)
       cinch_add_unit(test_driver_multivar
       SOURCES test/test_driver_multivar 
       LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
       POLICY MPI
       THREADS 1)
       cinch_add_unit(test_driver
       SOURCES test/test_driver 
       LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
       POLICY MPI
       THREADS 1)
       cinch_add_unit(test_driver_step
       SOURCES test/test_driver_step 
       LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
       POLICY MPI
       THREADS 1)
   endif (Jali_DIR)

   cinch_add_unit(test_driver_swarm
     SOURCES test/test_driver_swarm
     LIBRARIES portage ${LAPACKE_LIBRARIES}
     POLICY MPI
     THREADS 1)
   
   cinch_add_unit(test_driver_mesh_swarm_mesh
     SOURCES test/test_driver_mesh_swarm_mesh 
     LIBRARIES portage ${LAPACKE_LIBRARIES}
     POLICY MPI
     THREADS 1)

endif (ENABLE_UNIT_TESTS)
