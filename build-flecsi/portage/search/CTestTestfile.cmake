# CMake generated Testfile for 
# Source directory: /home/shevitz/ngc/portage/portage/search
# Build directory: /home/shevitz/ngc/portage/build-flecsi/portage/search
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(search_simple_test "/home/shevitz/ngc/portage/build-flecsi/test/search/search_simple_test" "--gtest_color=no")
set_tests_properties(search_simple_test PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/build-flecsi/test/search")
add_test(search_kdtree2_test "/home/shevitz/ngc/portage/build-flecsi/test/search/search_kdtree2_test" "--gtest_color=no")
set_tests_properties(search_kdtree2_test PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/build-flecsi/test/search")
add_test(search_kdtree3_test "/home/shevitz/ngc/portage/build-flecsi/test/search/search_kdtree3_test" "--gtest_color=no")
set_tests_properties(search_kdtree3_test PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/build-flecsi/test/search")
add_test(search_simple_points_test "/opt/local/opnsoft/openmpi/2.1.2-gcc-6.4.0/bin/mpiexec" "-n" "1" "/home/shevitz/ngc/portage/build-flecsi/test/search/search_simple_points_test" "--gtest_color=no")
set_tests_properties(search_simple_points_test PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/build-flecsi/test/search")
add_test(search_points_by_cells_test "/opt/local/opnsoft/openmpi/2.1.2-gcc-6.4.0/bin/mpiexec" "-n" "1" "/home/shevitz/ngc/portage/build-flecsi/test/search/search_points_by_cells_test" "--gtest_color=no")
set_tests_properties(search_points_by_cells_test PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/build-flecsi/test/search")
