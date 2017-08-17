#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]



file(COPY apptest_cmp.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY field1.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY field2.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY field3.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY field4.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME apptest_cmp
         COMMAND ./apptest_cmp.sh)
set_property(TEST apptest_cmp
         PROPERTY ENVIRONMENT APPDIR=${CMAKE_BINARY_DIR}/bin)


