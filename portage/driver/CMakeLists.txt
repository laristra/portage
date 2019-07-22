#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#


#-----------------------------------------------------------------------------~#

set(headers  mmdriver.h driver_swarm.h driver_mesh_swarm_mesh.h fix_mismatch.h)
if (TANGRAM_FOUND)
  list(APPEND headers write_to_gmv.h)
endif (TANGRAM_FOUND)

# Add header files
set(driver_HEADERS
   ${headers}
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

       cinch_add_unit(test_mismatch_fixup_serial
       SOURCES test/test_mismatch_fixup 
       LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
       POLICY MPI
       THREADS 1)

       cinch_add_unit(test_mismatch_disjoint
       SOURCES test/test_mismatch_disjoint
       LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
       POLICY MPI
       THREADS 1)

       cinch_add_unit(test_driver
       SOURCES test/test_driver 
       LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
       POLICY MPI
       THREADS 1)

       #cinch_add_unit(test_driver_step
       #SOURCES test/test_driver_step 
       #LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
       #POLICY MPI
       #THREADS 1)

       cinch_add_unit(test_driver_swarm_dist
       SOURCES test/test_driver_swarm_dist
       LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
       POLICY MPI
       THREADS 4)

     if (TANGRAM_FOUND AND XMOF2D_FOUND)

       cinch_add_unit(test_multimat_driver
       SOURCES test/test_multimat_driver
       LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES} ${XMOF2D_LIBRARIES}
       POLICY MPI
       THREADS 1)

       cinch_add_unit(test_multimat_driver_2ndOrder
       SOURCES test/test_multimat_driver_2ndOrder
       LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES} ${XMOF2D_LIBRARIES}
       POLICY MPI
       THREADS 1)

     endif (TANGRAM_FOUND AND XMOF2D_FOUND)

     if (TANGRAM_FOUND)

       cinch_add_unit(test_multimat_driver_distributed
       SOURCES test/test_multimat_driver_distributed
       LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES} 
       POLICY MPI
       THREADS 4)

     endif (TANGRAM_FOUND)

   endif (Jali_DIR)

   cinch_add_unit(test_driver_swarm
     SOURCES test/test_driver_swarm
     LIBRARIES portage
     POLICY MPI
     THREADS 1)

   cinch_add_unit(test_driver_mesh_swarm_mesh
     SOURCES test/test_driver_mesh_swarm_mesh
     LIBRARIES portage
     POLICY MPI
     THREADS 1)

endif (ENABLE_UNIT_TESTS)
