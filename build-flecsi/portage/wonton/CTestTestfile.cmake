# CMake generated Testfile for 
# Source directory: /home/shevitz/ngc/portage/portage/wonton
# Build directory: /home/shevitz/ngc/portage/build-flecsi/portage/wonton
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(flecsi_2d "/opt/local/opnsoft/openmpi/2.1.2-gcc-6.4.0/bin/mpiexec" "-n" "1" "/home/shevitz/ngc/portage/build-flecsi/test/wonton/flecsi_2d" "--gtest_color=no")
set_tests_properties(flecsi_2d PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/build-flecsi/test/wonton")
add_test(flecsi_3d "/opt/local/opnsoft/openmpi/2.1.2-gcc-6.4.0/bin/mpiexec" "-n" "1" "/home/shevitz/ngc/portage/build-flecsi/test/wonton/flecsi_3d" "--gtest_color=no")
set_tests_properties(flecsi_3d PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/build-flecsi/test/wonton")
add_test(test_simple_mesh_wrapper "/home/shevitz/ngc/portage/build-flecsi/test/wonton/test_simple_mesh_wrapper" "--gtest_color=no")
set_tests_properties(test_simple_mesh_wrapper PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/build-flecsi/test/wonton")
add_test(test_simple_state_wrapper "/home/shevitz/ngc/portage/build-flecsi/test/wonton/test_simple_state_wrapper" "--gtest_color=no")
set_tests_properties(test_simple_state_wrapper PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/build-flecsi/test/wonton")
add_test(test_simple_state_mm_wrapper "/home/shevitz/ngc/portage/build-flecsi/test/wonton/test_simple_state_mm_wrapper" "--gtest_color=no")
set_tests_properties(test_simple_state_mm_wrapper PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/build-flecsi/test/wonton")
