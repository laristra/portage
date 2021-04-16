#
# This file is part of the Ristra portage project.
# Please see the license file at the root of this repository, or at:
#     https://github.com/laristra/portage/blob/master/LICENSE
#

if (WONTON_ENABLE_MPI)
  message(STATUS "Adding particle analysis test")

  # copy script, input, source, gold and notebook files
  file(GLOB data
      "run.sh"
      "*.json"
      "*.csv"
      "../analysis.ipynb"
      "../plot.py"
      "gold/*.dat")

  file(COPY ${data} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

  # add tests then
  add_test(NAME "test_particle-analysis" COMMAND ./run.sh)

  set_property(TEST "test_particle-analysis"
      PROPERTY ENVIRONMENT
      ROOT_DIR=${CMAKE_BINARY_DIR}/app/particle-analysis
      COMPARE=${CMAKE_BINARY_DIR}/app/apptest_cmp/apptest_cmp)
endif()