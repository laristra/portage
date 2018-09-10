# CMake generated Testfile for 
# Source directory: /home/shevitz/ngc/portage/portage/state
# Build directory: /home/shevitz/ngc/portage/build-flecsi/portage/state
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(test_state_vector_uni "/home/shevitz/ngc/portage/build-flecsi/test/state/test_state_vector_uni" "--gtest_color=no")
set_tests_properties(test_state_vector_uni PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/build-flecsi/test/state")
add_test(test_state_vector_multi "/home/shevitz/ngc/portage/build-flecsi/test/state/test_state_vector_multi" "--gtest_color=no")
set_tests_properties(test_state_vector_multi PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/build-flecsi/test/state")
add_test(test_state_manager "/home/shevitz/ngc/portage/build-flecsi/test/state/test_state_manager" "--gtest_color=no")
set_tests_properties(test_state_manager PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/build-flecsi/test/state")
add_test(test_state_vector_uni_raw "/home/shevitz/ngc/portage/build-flecsi/test/state/test_state_vector_uni_raw" "--gtest_color=no")
set_tests_properties(test_state_vector_uni_raw PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/build-flecsi/test/state")
add_test(test_state_vector_multi_raw "/home/shevitz/ngc/portage/build-flecsi/test/state/test_state_vector_multi_raw" "--gtest_color=no")
set_tests_properties(test_state_vector_multi_raw PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/build-flecsi/test/state")
add_test(test_state_manager_raw "/home/shevitz/ngc/portage/build-flecsi/test/state/test_state_manager_raw" "--gtest_color=no")
set_tests_properties(test_state_manager_raw PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/build-flecsi/test/state")
