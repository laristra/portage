# CMake generated Testfile for 
# Source directory: /home/shevitz/ngc/portage/portage/interpolate
# Build directory: /home/shevitz/ngc/portage/build-flecsi/portage/interpolate
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(test_gradient "/home/shevitz/ngc/portage/build-flecsi/test/interpolate/test_gradient" "--gtest_color=no")
set_tests_properties(test_gradient PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/build-flecsi/test/interpolate")
add_test(test_interpolate_2nd_order "/home/shevitz/ngc/portage/build-flecsi/test/interpolate/test_interpolate_2nd_order" "--gtest_color=no")
set_tests_properties(test_interpolate_2nd_order PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/build-flecsi/test/interpolate")
add_test(test_interpolate_3rd_order "/home/shevitz/ngc/portage/build-flecsi/test/interpolate/test_interpolate_3rd_order" "--gtest_color=no")
set_tests_properties(test_interpolate_3rd_order PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/build-flecsi/test/interpolate")
add_test(test_quadfit "/home/shevitz/ngc/portage/build-flecsi/test/interpolate/test_quadfit" "--gtest_color=no")
set_tests_properties(test_quadfit PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/build-flecsi/test/interpolate")
