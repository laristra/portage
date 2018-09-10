# CMake generated Testfile for 
# Source directory: /home/shevitz/ngc/portage/portage/swarm
# Build directory: /home/shevitz/ngc/portage/build-flecsi/portage/swarm
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(test_swarm "/opt/local/opnsoft/openmpi/2.1.2-gcc-6.4.0/bin/mpiexec" "-n" "1" "/home/shevitz/ngc/portage/build-flecsi/test/swarm/test_swarm" "--gtest_color=no")
set_tests_properties(test_swarm PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/build-flecsi/test/swarm")
add_test(test_swarm_state "/opt/local/opnsoft/openmpi/2.1.2-gcc-6.4.0/bin/mpiexec" "-n" "1" "/home/shevitz/ngc/portage/build-flecsi/test/swarm/test_swarm_state" "--gtest_color=no")
set_tests_properties(test_swarm_state PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/build-flecsi/test/swarm")
