#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]
message(STATUS "Adding momentum remap test")

# this app can be run in serial or mpi mode - set an environment
# variable so that either will work
if(ENABLE_MPI)
  set(RUN_COMMAND "mpirun -np 1")
else()
  set(RUN_COMMAND "")
endif(ENABLE_MPI)

macro(ADD_MOMENTUM_TEST TEST_NUM)
  file(COPY momentumapp_test${TEST_NUM}.sh
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

  file(COPY errors_gold${TEST_NUM}.txt
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

  add_test(NAME momentumapp_test${TEST_NUM}
           COMMAND ./momentumapp_test${TEST_NUM}.sh)

  set_property(TEST momentumapp_test${TEST_NUM}
               PROPERTY ENVIRONMENT
               TESTAPPDIR=${CMAKE_CURRENT_BINARY_DIR}/..
               CMPAPPDIR=${CMAKE_BINARY_DIR}/app/apptest_cmp
               RUN_COMMAND=${RUN_COMMAND})

  set_property(TEST momentumapp_test${TEST_NUM}
               PROPERTY PROCESSORS 1)
endmacro(ADD_MOMENTUM_TEST)

ADD_MOMENTUM_TEST(0)
ADD_MOMENTUM_TEST(1)

