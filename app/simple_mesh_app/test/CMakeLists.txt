#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]
message(STATUS "Adding simple_mesh_app test")

# this app can be run in serial or mpi mode - set an environment
# variable so that either will work
if(ENABLE_MPI)
  set(RUN_COMMAND "mpirun -np 1")
else()
  set(RUN_COMMAND "")
endif(ENABLE_MPI)

macro(ADD_SIMPLE_MESH_TEST TEST_NUM)
  file(COPY simple_mesh_test${TEST_NUM}.sh
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY field_gold${TEST_NUM}.txt
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  add_test(NAME simple_mesh_test${TEST_NUM}
    COMMAND ./simple_mesh_test${TEST_NUM}.sh)
  set_property(TEST simple_mesh_test${TEST_NUM}
    PROPERTY ENVIRONMENT
    TESTAPPDIR=${CMAKE_CURRENT_BINARY_DIR}/..
    CMPAPPDIR=${CMAKE_BINARY_DIR}/app/apptest_cmp
    RUN_COMMAND=${RUN_COMMAND})
  set_property(TEST simple_mesh_test${TEST_NUM}
    PROPERTY PROCESSORS 1)
endmacro(ADD_SIMPLE_MESH_TEST)

ADD_SIMPLE_MESH_TEST(0)
ADD_SIMPLE_MESH_TEST(1)
ADD_SIMPLE_MESH_TEST(2)
ADD_SIMPLE_MESH_TEST(3)
ADD_SIMPLE_MESH_TEST(4)
ADD_SIMPLE_MESH_TEST(5)
ADD_SIMPLE_MESH_TEST(6)
ADD_SIMPLE_MESH_TEST(7)
ADD_SIMPLE_MESH_TEST(8)
ADD_SIMPLE_MESH_TEST(9)
ADD_SIMPLE_MESH_TEST(10)
ADD_SIMPLE_MESH_TEST(11)

file(COPY triple_point_3x3_matdata.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if (PORTAGE_ENABLE_TANGRAM)
  # copy the shell script and, a file with material distribution and gold
  file(
    COPY simple_mesh_mm_app.sh
         ../../../test_data/non_convex.bvf
         GOLD_simple_mm_test.gmv
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
  )

  # add the test
  add_test(
    NAME simple_mesh_mm_app_2d_MOF
    COMMAND ./simple_mesh_mm_app.sh
  )

  # set environment variables for the script
  set_property(
    TEST simple_mesh_mm_app_2d_MOF
    PROPERTY ENVIRONMENT
    TESTAPPDIR=${CMAKE_CURRENT_BINARY_DIR}/..
    CMPAPPDIR=${CMAKE_BINARY_DIR}/app/apptest_cmp
  )

endif ()

