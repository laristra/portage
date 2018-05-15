# CMake generated Testfile for 
# Source directory: /home/shevitz/ngc/portage/portage/search
# Build directory: /home/shevitz/ngc/portage/app/portage/search
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(search_simple_test "/home/shevitz/ngc/portage/app/test/search/search_simple_test" "--gtest_color=no")
set_tests_properties(search_simple_test PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/app/test/search")
add_test(search_kdtree2_test "/home/shevitz/ngc/portage/app/test/search/search_kdtree2_test" "--gtest_color=no")
set_tests_properties(search_kdtree2_test PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/app/test/search")
add_test(search_kdtree3_test "/home/shevitz/ngc/portage/app/test/search/search_kdtree3_test" "--gtest_color=no")
set_tests_properties(search_kdtree3_test PROPERTIES  WORKING_DIRECTORY "/home/shevitz/ngc/portage/app/test/search")
