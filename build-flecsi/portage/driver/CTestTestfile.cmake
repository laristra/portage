# CMake generated Testfile for 
# Source directory: /home/shevitz/ngc/portage/portage/driver
# Build directory: /home/shevitz/ngc/portage/build-flecsi/portage/driver
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(test_driver_swarm "/opt/local/opnsoft/openmpi/2.1.2-gcc-6.4.0/bin/mpiexec" "-n" "1" "/home/shevitz/ngc/portage/build-flecsi/test/driver/test_driver_swarm" "--gtest_color=no")
set_tests_properties(test_driver_swarm PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/build-flecsi/test/driver")
add_test(test_driver_mesh_swarm_mesh "/opt/local/opnsoft/openmpi/2.1.2-gcc-6.4.0/bin/mpiexec" "-n" "1" "/home/shevitz/ngc/portage/build-flecsi/test/driver/test_driver_mesh_swarm_mesh" "--gtest_color=no")
set_tests_properties(test_driver_mesh_swarm_mesh PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/build-flecsi/test/driver")
