# CMake generated Testfile for 
# Source directory: /home/shevitz/ngc/portage/portage/accumulate
# Build directory: /home/shevitz/ngc/portage/build-flecsi/portage/accumulate
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(test_accumulate "/opt/local/opnsoft/openmpi/2.1.2-gcc-6.4.0/bin/mpiexec" "-n" "1" "/home/shevitz/ngc/portage/build-flecsi/test/accumulate/test_accumulate" "--gtest_color=no")
set_tests_properties(test_accumulate PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/build-flecsi/test/accumulate")
