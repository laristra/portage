# CMake generated Testfile for 
# Source directory: /home/shevitz/ngc/portage/portage/wonton
# Build directory: /home/shevitz/ngc/portage/app/portage/wonton
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(test_simple_mesh_wrapper "/home/shevitz/ngc/portage/app/test/wonton/test_simple_mesh_wrapper" "--gtest_color=no")
set_tests_properties(test_simple_mesh_wrapper PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/app/test/wonton")
add_test(test_simple_state_wrapper "/home/shevitz/ngc/portage/app/test/wonton/test_simple_state_wrapper" "--gtest_color=no")
set_tests_properties(test_simple_state_wrapper PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/app/test/wonton")
