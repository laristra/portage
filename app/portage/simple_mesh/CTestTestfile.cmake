# CMake generated Testfile for 
# Source directory: /home/shevitz/ngc/portage/portage/simple_mesh
# Build directory: /home/shevitz/ngc/portage/app/portage/simple_mesh
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(test_simple_mesh "/home/shevitz/ngc/portage/app/test/simple_mesh/test_simple_mesh" "--gtest_color=no")
set_tests_properties(test_simple_mesh PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/app/test/simple_mesh")
add_test(test_simple_state "/home/shevitz/ngc/portage/app/test/simple_mesh/test_simple_state" "--gtest_color=no")
set_tests_properties(test_simple_state PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/app/test/simple_mesh")
