# CMake generated Testfile for 
# Source directory: /home/shevitz/ngc/portage/portage/support
# Build directory: /home/shevitz/ngc/portage/build-flecsi/portage/support
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(test_matfuncs "/opt/local/opnsoft/openmpi/2.1.2-gcc-6.4.0/bin/mpiexec" "-n" "1" "/home/shevitz/ngc/portage/build-flecsi/test/support/test_matfuncs" "--gtest_color=no")
set_tests_properties(test_matfuncs PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/build-flecsi/test/support")
add_test(test_weight "/opt/local/opnsoft/openmpi/2.1.2-gcc-6.4.0/bin/mpiexec" "-n" "1" "/home/shevitz/ngc/portage/build-flecsi/test/support/test_weight" "--gtest_color=no")
set_tests_properties(test_weight PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/build-flecsi/test/support")
add_test(test_basis "/opt/local/opnsoft/openmpi/2.1.2-gcc-6.4.0/bin/mpiexec" "-n" "1" "/home/shevitz/ngc/portage/build-flecsi/test/support/test_basis" "--gtest_color=no")
set_tests_properties(test_basis PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/build-flecsi/test/support")
add_test(test_operator "/opt/local/opnsoft/openmpi/2.1.2-gcc-6.4.0/bin/mpiexec" "-n" "1" "/home/shevitz/ngc/portage/build-flecsi/test/support/test_operator" "--gtest_color=no")
set_tests_properties(test_operator PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/build-flecsi/test/support")
