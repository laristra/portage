# CMake generated Testfile for 
# Source directory: /home/shevitz/ngc/portage/app/apptest_cmp/test
# Build directory: /home/shevitz/ngc/portage/build-flecsi/app/apptest_cmp/test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(apptest_cmp "./apptest_cmp.sh")
set_tests_properties(apptest_cmp PROPERTIES  ENVIRONMENT "APPDIR=/home/shevitz/ngc/portage/build-flecsi/app/apptest_cmp/test/..")
