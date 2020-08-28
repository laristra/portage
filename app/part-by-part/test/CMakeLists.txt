#
# This file is part of the Ristra portage project.
# Please see the license file at the root of this repository, or at:
#     https://github.com/laristra/portage/blob/master/LICENSE
#

if (WONTON_ENABLE_Jali AND WONTON_ENABLE_MPI)

  message(STATUS "Adding part-by-part tests")

  # copy run script and exodus files in build directory
  file(COPY "run.sh" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY "../../../test_data/cast_source.exo" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY "../../../test_data/cast_target.exo" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY "../../../test_data/semi2-coarse.exo" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY "../../../test_data/semi3-coarse.exo" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY "../../../test_data/semi2-fine.exo" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY "../../../test_data/semi3-fine.exo" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

  # macro to compare remapped values with reference ones
  macro(ADD_REMAP_TEST SUFFIX)
    # copy gold files and input parameter file in build directory
    file(COPY "input_${SUFFIX}.json" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    file(COPY "gold/gold_${SUFFIX}_temperature.dat" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    if (NOT "${SUFFIX}" MATCHES "blocks$")
      file(COPY "gold/gold_${SUFFIX}_density.dat" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    endif()

    add_test(NAME "test_part-by-part_${SUFFIX}" COMMAND ./run.sh ${SUFFIX})

    set_property(TEST "test_part-by-part_${SUFFIX}"
        PROPERTY ENVIRONMENT
        ROOT_DIR=${CMAKE_BINARY_DIR}/app/part-by-part
        COMPARE=${CMAKE_BINARY_DIR}/app/apptest_cmp/apptest_cmp
        RUN_COMMAND=${RUN_COMMAND})
  endmacro()

# add tests then
ADD_REMAP_TEST("dim_2_analytic")
ADD_REMAP_TEST("dim_3_analytic")
ADD_REMAP_TEST("dim_3_blocks")

endif()