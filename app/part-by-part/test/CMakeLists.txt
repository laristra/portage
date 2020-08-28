#
# This file is part of the Ristra portage project.
# Please see the license file at the root of this repository, or at:
#     https://github.com/laristra/portage/blob/master/LICENSE
#

if (WONTON_ENABLE_Jali AND WONTON_ENABLE_MPI)

  # macro to compare remapped values with reference ones
  macro(ADD_REMAP_TEST SUFFIX)
    add_test(NAME "test_part-by-part_${SUFFIX}" COMMAND ./run.sh ${SUFFIX})

    set_property(TEST "test_part-by-part_${SUFFIX}"
        PROPERTY ENVIRONMENT
        ROOT_DIR=${CMAKE_BINARY_DIR}/app/part-by-part
        COMPARE=${CMAKE_BINARY_DIR}/app/apptest_cmp/apptest_cmp)
  endmacro()

  message(STATUS "Adding part-by-part tests")

  # copy run script, input meshes, reference values in build directory
  file(GLOB data
      "run.sh"
      "input_*.json"
      "gold/*.dat"
      "../../../test_data/cast*.exo"
      "../../../test_data/semi*.exo")

  file(COPY ${data} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

  # add tests then
  ADD_REMAP_TEST("dim_2_analytic")
  ADD_REMAP_TEST("dim_3_analytic")
  ADD_REMAP_TEST("dim_3_blocks")
endif()