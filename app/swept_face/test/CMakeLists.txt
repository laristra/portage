#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]
# Config file for a test directory.
# At some future point, much of what appears here will probably
# be pulled into a macro.


####################################################################

if (Jali_DIR)

  message(STATUS "Adding swept-face remap tests")

  # create a macro to compare actual values with canonical ones
  macro(ADD_REMAP_TEST SUFFIX)

    file(COPY "swept_${SUFFIX}.sh" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    file(COPY "gold/gold_${SUFFIX}.txt" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

    add_test(NAME "test_swept-face_${SUFFIX}"
             COMMAND "./swept_${SUFFIX}.sh")

    set_property(TEST "test_swept-face_${SUFFIX}"
                 PROPERTY ENVIRONMENT
                 APPTEST=${CMAKE_BINARY_DIR}/app/swept_face/swept_face_app
                 COMPARE=${CMAKE_BINARY_DIR}/app/apptest_cmp/apptest_cmp
                 RUN_COMMAND=${RUN_COMMAND})
  endmacro()

  # create a macro to run a distributed compare test
  macro(ADD_DISTRIB_TEST NP DIM NCELLS FIELDS REMAP_ORDER FIELD_ALIAS)
    file(COPY swept_distrib.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

    add_test(NAME test_swept-face_distrib_rank_${NP}_dim_${DIM}_field_${FIELD_ALIAS}
             COMMAND swept_distrib.sh ${NP} ${DIM} ${NCELLS} ${FIELDS} ${REMAP_ORDER} ${FIELD_ALIAS})

    set_property(TEST test_swept-face_distrib_rank_${NP}_dim_${DIM}_field_${FIELD_ALIAS}
                 PROPERTY ENVIRONMENT
                 APPTEST=${CMAKE_BINARY_DIR}/app/swept_face/swept_face_app
                 COMPARE=${CMAKE_BINARY_DIR}/app/distributed_cmp/distributed_cmp)
  endmacro()

  # add the tests
  ADD_REMAP_TEST("dim_2_field_const_order_1_simple")
  ADD_REMAP_TEST("dim_2_field_const_order_1_vortex")
  ADD_REMAP_TEST("dim_2_field_linear_order_1_simple")
  ADD_REMAP_TEST("dim_2_field_linear_order_1_vortex")
  ADD_REMAP_TEST("dim_2_field_linear_order_2_simple")
  ADD_REMAP_TEST("dim_2_field_linear_order_2_vortex")
  ADD_REMAP_TEST("dim_3_field_const_order_1_simple")
  ADD_REMAP_TEST("dim_3_field_const_order_1_vortex")
  ADD_REMAP_TEST("dim_3_field_linear_order_1_simple")
  ADD_REMAP_TEST("dim_3_field_linear_order_1_vortex")
  ADD_REMAP_TEST("dim_3_field_linear_order_2_simple")
  ADD_REMAP_TEST("dim_3_field_linear_order_2_vortex")

  if (ENABLE_MPI)
    # 2D, constant field, order 2
    ADD_DISTRIB_TEST(2 2 5 "42" 2 const)
    ADD_DISTRIB_TEST(4 2 5 "42" 2 const)
    ADD_DISTRIB_TEST(8 2 5 "42" 2 const)

    # 2D, linear field, order 2
    ADD_DISTRIB_TEST(2 2 5 "x+2y" 2 linear)
    ADD_DISTRIB_TEST(4 2 5 "x+2y" 2 linear)
    ADD_DISTRIB_TEST(8 2 5 "x+2y" 2 linear)

    # 3D, constant field, order 2
    ADD_DISTRIB_TEST(2 3 5 "42" 2 const)
    ADD_DISTRIB_TEST(4 3 5 "42" 2 const)
    ADD_DISTRIB_TEST(8 3 5 "42" 2 const)

    # 3D, linear field, order 2
    ADD_DISTRIB_TEST(2 3 5 "x+2y+3z" 2 linear)
    ADD_DISTRIB_TEST(4 3 5 "x+2y+3z" 2 linear)
    ADD_DISTRIB_TEST(8 3 5 "x+2y+3z" 2 linear)
  endif()

endif(Jali_DIR)
