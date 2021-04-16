#
# This file is part of the Ristra portage project.
# Please see the license file at the root of this repository, or at:
#     https://github.com/laristra/portage/blob/master/LICENSE
#

if (WONTON_ENABLE_MPI)

  message(STATUS "Adding msmapp test")

  # copy run script, input parameters, and gold files
  file(GLOB data
       "run.sh"
       "example_input"
       "gold/gold_diagnostics.dat")
  
  file(COPY ${data} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

  # add tests
  add_test(NAME "test_msmapp" COMMAND ./run.sh)

  set_property(TEST "test_msmapp"
	       PROPERTY ENVIRONMENT
               ROOT_DIR=${CMAKE_BINARY_DIR}/app/msmapp
               COMPARE=${CMAKE_BINARY_DIR}/app/apptest_cmp/apptest_cmp)

endif()
