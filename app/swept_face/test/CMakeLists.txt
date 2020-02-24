#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]
# Config file for a test directory.
# At some future point, much of what appears here will probably
# be pulled into a macro.


####################################################################

message(STATUS "Adding swept_face_demo tests")

# create a macro to run a distributed compare test
macro(ADD_DISTRIB_TEST NP DIM NCELLS FIELDS REMAP_ORDER)

  file(COPY swept_face_demo.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

  add_test(
    NAME swept_face_demo_${NP}_${DIM}_${NCELLS}_${FIELDS}_${REMAP_ORDER}
    COMMAND swept_face_demo.sh ${NP} ${DIM} ${NCELLS} ${FIELDS} ${REMAP_ORDER}
  )

  set_property(
    TEST swept_face_demo_${NP}_${DIM}_${NCELLS}_${FIELDS}_${REMAP_ORDER}
    PROPERTY ENVIRONMENT
    TESTAPPDIR="${CMAKE_CURRENT_BINARY_DIR}/.."
    CMPAPPDIR="${CMAKE_BINARY_DIR}/app/distributed_cmp"
  )

endmacro()

# add the tests
# These tests don't test whether the remap is correct or reproduces a certain
# function. They only test whether serial and parallel match. Hence higher 
# order fields or remap order changes don't really test anything new.
if (ENABLE_MPI AND Jali_DIR AND TANGRAM_FOUND)

  # 2D, constant field, order 2, vary number of processors
  ADD_DISTRIB_TEST(2 2 5 42 2)
  ADD_DISTRIB_TEST(4 2 5 42 2)
  ADD_DISTRIB_TEST(8 2 5 42 2)

  # 2D, linear field, order 2, vary number of processors
  ADD_DISTRIB_TEST(2 2 5 x+2*y 2)
  ADD_DISTRIB_TEST(4 2 5 x+2*y 2)
  ADD_DISTRIB_TEST(8 2 5 x+2*y 2)

  # 3D, constant field, order 2, vary number of processors
  ADD_DISTRIB_TEST(2 3 5 42 2)
  ADD_DISTRIB_TEST(4 3 5 42 2)
  ADD_DISTRIB_TEST(8 3 5 42 2)

  # 3D, linear field, order 2, vary number of processors
  ADD_DISTRIB_TEST(2 3 5 x+2*y+3*z 2)
  ADD_DISTRIB_TEST(4 3 5 x+2*y+3*z 2)
  ADD_DISTRIB_TEST(8 3 5 x+2*y+3*z 2)
endif()


