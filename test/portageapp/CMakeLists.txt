# Config file for a test directory.
# At some future point, much of what appears here will probably
# be pulled into a macro.

message(STATUS "Adding portageapp test")

file(COPY runme.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME portageapp_test
         COMMAND ./runme.sh)
set_property(TEST portageapp_test
         PROPERTY ENVIRONMENT APPDIR=${CMAKE_BINARY_DIR}/bin)
set_property(TEST portageapp_test
         PROPERTY PROCESSORS 1)


file(COPY portageapp_test0.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY field_gold0.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME portageapp_test0
         COMMAND ./portageapp_test0.sh)
set_property(TEST portageapp_test0
         PROPERTY ENVIRONMENT APPDIR=${CMAKE_BINARY_DIR}/bin)
set_property(TEST portageapp_test0
         PROPERTY PROCESSORS 4)


file(COPY portageapp_test1.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY field_gold1.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME portageapp_test1
         COMMAND ./portageapp_test1.sh)
set_property(TEST portageapp_test1
         PROPERTY ENVIRONMENT APPDIR=${CMAKE_BINARY_DIR}/bin)
set_property(TEST portageapp_test1
         PROPERTY PROCESSORS 4)


file(COPY portageapp_test13.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY field_gold13.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME portageapp_test13
         COMMAND ./portageapp_test13.sh)
set_property(TEST portageapp_test13
         PROPERTY ENVIRONMENT APPDIR=${CMAKE_BINARY_DIR}/bin)
set_property(TEST portageapp_test13
         PROPERTY PROCESSORS 1)
