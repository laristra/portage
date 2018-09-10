# CMake generated Testfile for 
# Source directory: /home/shevitz/ngc/portage/portage/intersect
# Build directory: /home/shevitz/ngc/portage/build-flecsi/portage/intersect
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(intersect_r2d "/home/shevitz/ngc/portage/build-flecsi/test/intersect/intersect_r2d" "--gtest_color=no")
set_tests_properties(intersect_r2d PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/build-flecsi/test/intersect")
add_test(intersect_r3d "/home/shevitz/ngc/portage/build-flecsi/test/intersect/intersect_r3d" "--gtest_color=no")
set_tests_properties(intersect_r3d PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/build-flecsi/test/intersect")
add_test(test_tangram_intersect_2D "/home/shevitz/ngc/portage/build-flecsi/test/intersect/test_tangram_intersect_2D" "--gtest_color=no")
set_tests_properties(test_tangram_intersect_2D PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/build-flecsi/test/intersect")
add_test(test_tangram_intersect_3D "/home/shevitz/ngc/portage/build-flecsi/test/intersect/test_tangram_intersect_3D" "--gtest_color=no")
set_tests_properties(test_tangram_intersect_3D PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/build-flecsi/test/intersect")
