#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]
# Config file for a test directory.
# At some future point, much of what appears here will probably
# be pulled into a macro.

message(STATUS "Adding portageapp test")

macro(ADD_SPECIFIC_APPTEST TYPE SUFFIX)
  file(COPY portageapp_${TYPE}_${SUFFIX}.sh
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(GLOB GOLDFILES gold/GOLD_${TYPE}_${SUFFIX}*.txt*)
  file(COPY ${GOLDFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  add_test(NAME portageapp_${TYPE}_${SUFFIX}
    COMMAND ./portageapp_${TYPE}_${SUFFIX}.sh)
  set_property(TEST portageapp_${TYPE}_${SUFFIX}
    PROPERTY ENVIRONMENT
    TESTAPPDIR=${CMAKE_CURRENT_BINARY_DIR}/..
    CMPAPPDIR=${CMAKE_BINARY_DIR}/app/apptest_cmp)
endmacro(ADD_SPECIFIC_APPTEST)

macro(ADD_PORTAGEAPP_TEST DIM ENTITY FIELDORDER REMAPORDER)
  set(SUFFIX "rect_${DIM}d_${ENTITY}_f${FIELDORDER}_r${REMAPORDER}")

  if (WONTON_ENABLE_Jali)
    ADD_SPECIFIC_APPTEST(jali ${SUFFIX})
  endif ()

  if (WONTON_ENABLE_FleCSI)
    # For now, the flecsi app is hard-coded for 2d, cell-centered
    # remap with the field and remap order being the same
    if (${DIM} MATCHES "2" AND
	${ENTITY} MATCHES "cell" AND
	${FIELDORDER} MATCHES ${REMAPORDER})
      ADD_SPECIFIC_APPTEST(flecsi ${SUFFIX})
    endif()
  endif ()
endmacro(ADD_PORTAGEAPP_TEST)

# all of these scripts use mpirun to launch the job, so only include them
# if we know about MPI
if (WONTON_ENABLE_MPI)
  ADD_PORTAGEAPP_TEST(2 cell 0 1)
  ADD_PORTAGEAPP_TEST(2 cell 1 1)
  ADD_PORTAGEAPP_TEST(2 cell 1 2)
  ADD_PORTAGEAPP_TEST(2 cell 2 2)
  ADD_PORTAGEAPP_TEST(2 node 0 1)
  ADD_PORTAGEAPP_TEST(2 node 1 1)
  ADD_PORTAGEAPP_TEST(2 node 1 2)
  ADD_PORTAGEAPP_TEST(2 node 2 2)
  ADD_PORTAGEAPP_TEST(3 cell 0 1)
  ADD_PORTAGEAPP_TEST(3 cell 1 1)
  ADD_PORTAGEAPP_TEST(3 cell 1 2)
  ADD_PORTAGEAPP_TEST(3 cell 2 2)
  ADD_PORTAGEAPP_TEST(3 node 0 1)
  ADD_PORTAGEAPP_TEST(3 node 1 1)
  ADD_PORTAGEAPP_TEST(3 node 1 2)
  ADD_PORTAGEAPP_TEST(3 node 2 2)

  if (WONTON_ENABLE_Jali)
    ADD_SPECIFIC_APPTEST(jali_rect 3d_cell_gen_r2)
  endif ()
endif ()


macro(ADD_UNSTRUC_TEST DIM ENTITY FIELDORDER REMAPORDER)
  set(SUFFIX "${DIM}d_${ENTITY}_f${FIELDORDER}_r${REMAPORDER}")
  file(COPY portageapp_jali_unstruc_${SUFFIX}.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(GLOB GOLDFILES gold/GOLD_jali_unstruc_${SUFFIX}.txt*)
  file(COPY ${GOLDFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

  file(COPY ../../../test_data/shotshell.exo
            ../../../test_data/shotshell-v.exo
	    ../../../test_data/cube-poly1.exo
	    ../../../test_data/cube-poly2.exo
            ../../../test_data/coarse_poly2D.exo
            ../../../test_data/fine_poly2D.exo 
            ../../../test_data/3dvoro-voronoized-10k.exo
	    ../../../test_data/sector5.exo
	    ../../../test_data/sector6.exo
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  add_test(NAME portageapp_jali_unstruc_${SUFFIX}
    COMMAND ./portageapp_jali_unstruc_${SUFFIX}.sh)
  set_property(TEST portageapp_jali_unstruc_${SUFFIX}
    PROPERTY ENVIRONMENT
    TESTAPPDIR=${CMAKE_CURRENT_BINARY_DIR}/..
    CMPAPPDIR=${CMAKE_BINARY_DIR}/app/apptest_cmp)
endmacro(ADD_UNSTRUC_TEST)

if (WONTON_ENABLE_MPI AND WONTON_ENABLE_Jali)
  ADD_UNSTRUC_TEST(2 cell 0 1)
  ADD_UNSTRUC_TEST(2 cell 1 1)
  ADD_UNSTRUC_TEST(2 cell 1 2)
  ADD_UNSTRUC_TEST(2 cell 2 2)
  ADD_UNSTRUC_TEST(2 node 0 1)
  ADD_UNSTRUC_TEST(2 node 1 1)
  ADD_UNSTRUC_TEST(2 node 1 2)
  ADD_UNSTRUC_TEST(2 node 2 2)
  ADD_UNSTRUC_TEST(3 cell 0 1)
  ADD_UNSTRUC_TEST(3 cell 1 1)
  ADD_UNSTRUC_TEST(3 cell 1 2)
  ADD_UNSTRUC_TEST(3 cell 2 2)
  ADD_UNSTRUC_TEST(3 node 0 1)
  ADD_UNSTRUC_TEST(3 node 1 1)
  ADD_UNSTRUC_TEST(3 node 1 2)
  ADD_UNSTRUC_TEST(3 node 2 2)
  ADD_UNSTRUC_TEST(2 cell 1 8)
  ADD_UNSTRUC_TEST(3 cell 8 8)

  ADD_SPECIFIC_APPTEST(jali_unstruc 2d_cell_gen_r2)
  ADD_SPECIFIC_APPTEST(jali_unstruc 3d_cell_gen_r2)

  ADD_SPECIFIC_APPTEST(jali_radial 3d_node_f1_r2)
endif ()

#multi-material tests
if (WONTON_ENABLE_Jali AND PORTAGE_ENABLE_TANGRAM)
  file(COPY portageapp_multimat_jali_unstruc_2d_cell_f1_r2.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY gold/GOLD_multimat_jali_unstruc_2d_cell_f1_r2.gmv DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY ../../../test_data/reg10x10.exo
            ../../../test_data/reg10x10.bvf
            ../../../test_data/polymesh_10x10.exo
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  add_test(NAME portageapp_multimat_jali_unstruc_2d_cell_f1_r2
    COMMAND ./portageapp_multimat_jali_unstruc_2d_cell_f1_r2.sh)
  set_property(TEST portageapp_multimat_jali_unstruc_2d_cell_f1_r2
    PROPERTY ENVIRONMENT
    TESTAPPDIR=${CMAKE_CURRENT_BINARY_DIR}/..
    CMPAPPDIR=${CMAKE_BINARY_DIR}/app/apptest_cmp)


  file(COPY portageapp_multimat_jali_nonconvex_2d_cell_f1_r2.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY gold/GOLD_multimat_jali_nonconvex_2d_cell_f1_r2.gmv DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY ../../../test_data/non_convex.exo
            ../../../test_data/non_convex.bvf
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  add_test(NAME portageapp_multimat_jali_nonconvex_2d_cell_f1_r2
    COMMAND ./portageapp_multimat_jali_nonconvex_2d_cell_f1_r2.sh)
  set_property(TEST portageapp_multimat_jali_nonconvex_2d_cell_f1_r2
    PROPERTY ENVIRONMENT
    TESTAPPDIR=${CMAKE_CURRENT_BINARY_DIR}/..
    CMPAPPDIR=${CMAKE_BINARY_DIR}/app/apptest_cmp)


  file(COPY portageapp_rgmd_jali_tjunction_2d_swept.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY gold/GOLD_portageapp_rgmd_jali_tjunction_2d_swept.gmv DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  add_test(NAME portageapp_rgmd_jali_tjunction_2d_swept
    COMMAND ./portageapp_rgmd_jali_tjunction_2d_swept.sh)
  set_property(TEST portageapp_rgmd_jali_tjunction_2d_swept
    PROPERTY ENVIRONMENT
    TESTAPPDIR=${CMAKE_CURRENT_BINARY_DIR}/..
    CMPAPPDIR=${CMAKE_BINARY_DIR}/app/apptest_cmp)


  file(COPY portageapp_rgmd_jali_ball_limiter.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY gold/GOLD_portageapp_rgmd_jali_ball_limiter.gmv DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  add_test(NAME portageapp_rgmd_jali_ball_limiter
    COMMAND ./portageapp_rgmd_jali_ball_limiter.sh)
  set_property(TEST portageapp_rgmd_jali_ball_limiter
    PROPERTY ENVIRONMENT
    TESTAPPDIR=${CMAKE_CURRENT_BINARY_DIR}/..
    CMPAPPDIR=${CMAKE_BINARY_DIR}/app/apptest_cmp)

endif ()

####################################################################

message(STATUS "Adding portageapp_rgmd_jali tests")

# create a macro to run a distributed compare test
macro(ADD_DISTRIBUTED_APPTEST NP PROBLEM DIM NSOURCECELLS NTARGETCELLS FIELDS REMAP_ORDER )

  # copy the shell script
  file(
    COPY portageapp_rgmd_jali.sh 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
  )
  
  # add the test
  add_test(
    NAME portageapp_rgmd_jali_${NP}_${PROBLEM}_${DIM}_${NSOURCECELLS}_${NTARGETCELLS}_${FIELDS}_${REMAP_ORDER}
    COMMAND ./portageapp_rgmd_jali.sh ${NP} ${PROBLEM} ${DIM} ${NSOURCECELLS} ${NTARGETCELLS} ${FIELDS} ${REMAP_ORDER} 
  )
    
  # set environment variables for the script
  set_property(
    TEST portageapp_rgmd_jali_${NP}_${PROBLEM}_${DIM}_${NSOURCECELLS}_${NTARGETCELLS}_${FIELDS}_${REMAP_ORDER}
    PROPERTY ENVIRONMENT
    TESTAPPDIR=${CMAKE_CURRENT_BINARY_DIR}/..
    CMPAPPDIR=${CMAKE_BINARY_DIR}/app/distributed_cmp
  )

endmacro(ADD_DISTRIBUTED_APPTEST)

# add the tests
# These tests don't test whether the remap is correct or reproduces a certain
# function. They only test whether serial and parallel match. Hence higher 
# order fields or remap order changes don't really test anything new.
if (WONTON_ENABLE_MPI AND WONTON_ENABLE_Jali AND PORTAGE_ENABLE_TANGRAM)

  # 2D, t-junction, constant fields, vary number of processors with triple point
  ADD_DISTRIBUTED_APPTEST(2 tjunction 2 5 7 1,2,3 1)
  ADD_DISTRIBUTED_APPTEST(4 tjunction 2 5 7 1,2,3 1)
  ADD_DISTRIBUTED_APPTEST(2 tjunction 2 5 7 1,2,3 2)
  ADD_DISTRIBUTED_APPTEST(4 tjunction 2 5 7 1,2,3 2)
  
  # 2D, t-junction, linear fields, larger meshes 
  ADD_DISTRIBUTED_APPTEST(4 tjunction 2 51 71 x,y,x+.3*y 2)

  # 2D, pure source cells, constant fields, vary number of processors
  ADD_DISTRIBUTED_APPTEST(4 srcpurecells 2 4 8 1,2,3 1)
    
  # 3D, t-junction, linear fields, vary number of processors  
  ADD_DISTRIBUTED_APPTEST(4 tjunction 3 5 7 x,y,x+.3*y 1)

  # 3D, pure source cells, constant fields, vary number of processors
  ADD_DISTRIBUTED_APPTEST(4 srcpurecells 3 4 8 1,2,3 1)

endif ()


