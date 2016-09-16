# Config file for a test directory.
# At some future point, much of what appears here will probably
# be pulled into a macro.

message(STATUS "Adding portageapp test")

file(COPY runme.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME portageapp_test
         COMMAND ./runme.sh)
set_property(TEST portageapp_test 
         PROPERTY ENVIRONMENT APPDIR=${CMAKE_BINARY_DIR}/bin)


