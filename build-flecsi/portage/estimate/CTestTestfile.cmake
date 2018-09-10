# CMake generated Testfile for 
# Source directory: /home/shevitz/ngc/portage/portage/estimate
# Build directory: /home/shevitz/ngc/portage/build-flecsi/portage/estimate
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(test_estimate "/opt/local/opnsoft/openmpi/2.1.2-gcc-6.4.0/bin/mpiexec" "-n" "1" "/home/shevitz/ngc/portage/build-flecsi/test/estimate/test_estimate" "--gtest_color=no")
set_tests_properties(test_estimate PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/build-flecsi/test/estimate")
