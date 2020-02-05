#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]
# Config file for a test directory.
# At some future point, much of what appears here will probably
# be pulled into a macro.


####################################################################

message(STATUS "Adding sweptfaceapp_jali_singlemat tests")

# create a macro to run a distributed compare test
macro(ADD_DISTRIBUTED_APPTEST NP DIM NCELLS FIELDS REMAP_ORDER )

  # copy the shell script
  file(
    COPY sweptfaceapp_jali_singlemat.sh 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
  )
  
  # add the test
  add_test(
    NAME sweptfaceapp_jali_singlemat_${NP}_${DIM}_${NCELLS}_${FIELDS}_${REMAP_ORDER}
    COMMAND ./sweptfaceapp_jali_singlemat.sh ${NP} ${DIM} ${NCELLS} ${FIELDS} ${REMAP_ORDER} 
  )
    
  # set environment variables for the script
  set_property(
    TEST sweptfaceapp_jali_singlemat_${NP}_${DIM}_${NCELLS}_${FIELDS}_${REMAP_ORDER}
    PROPERTY ENVIRONMENT
    TESTAPPDIR=${CMAKE_CURRENT_BINARY_DIR}/..
    CMPAPPDIR=${CMAKE_BINARY_DIR}/app/distributed_cmp
  )

endmacro(ADD_DISTRIBUTED_APPTEST)

# add the tests
# These tests don't test whether the remap is correct or reproduces a certain
# function. They only test whether serial and parallel match. Hence higher 
# order fields or remap order changes don't really test anything new.
if (ENABLE_MPI AND Jali_DIR AND TANGRAM_FOUND)

  # 2D, constant field, ro 1, vary number of processors 
  ADD_DISTRIBUTED_APPTEST(2 2 5 42 1)
  ADD_DISTRIBUTED_APPTEST(4 2 5 42 1)
  
  # 2D, constant field, ro 2, vary number of processors 
  ADD_DISTRIBUTED_APPTEST(2 2 5 42 2)
  ADD_DISTRIBUTED_APPTEST(4 2 5 42 2)

  # 2D, linear field, ro 1, vary number of processors 
  ADD_DISTRIBUTED_APPTEST(2 2 5 x+3*y 1)
  ADD_DISTRIBUTED_APPTEST(4 2 5 x+3*y 1)

  # 2D, linear field, ro 2, vary number of processors 
  ADD_DISTRIBUTED_APPTEST(2 2 5 x+3*y 2)
  ADD_DISTRIBUTED_APPTEST(4 2 5 x+3*y 2)

  # 3D, constant field, ro 1, vary number of processors 
  ADD_DISTRIBUTED_APPTEST(2 3 5 42 1)
  ADD_DISTRIBUTED_APPTEST(4 3 5 42 1)
  
  # 3D, constant field, ro 2, vary number of processors 
  ADD_DISTRIBUTED_APPTEST(2 3 5 42 2)
  ADD_DISTRIBUTED_APPTEST(4 3 5 42 2)

  # 3D, linear field, ro 1, vary number of processors 
  ADD_DISTRIBUTED_APPTEST(2 3 5 x+3*y+5*z 1)
  ADD_DISTRIBUTED_APPTEST(4 3 5 x+3*y+5*z 1)

  # 3D, linear field, ro 2, vary number of processors 
  ADD_DISTRIBUTED_APPTEST(2 3 5 x+3*y+5*z 2)
  ADD_DISTRIBUTED_APPTEST(4 3 5 x+3*y+5*z 2)

endif (ENABLE_MPI AND Jali_DIR AND TANGRAM_FOUND)


