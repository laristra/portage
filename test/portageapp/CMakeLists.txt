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

macro(ADD_PORTAGE_TEST TEST_NUM NUM_PROC)
  file(COPY portageapp_test${TEST_NUM}.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY field_gold${TEST_NUM}.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  add_test(NAME portageapp_test${TEST_NUM}
    COMMAND ./portageapp_test${TEST_NUM}.sh)
  set_property(TEST portageapp_test${TEST_NUM}
    PROPERTY ENVIRONMENT APPDIR=${CMAKE_BINARY_DIR}/bin)
  set_property(TEST portageapp_test${TEST_NUM}
    PROPERTY PROCESSORS ${NUM_PROC})
endmacro(ADD_PORTAGE_TEST)

ADD_PORTAGE_TEST(0 4)
ADD_PORTAGE_TEST(1 4)
ADD_PORTAGE_TEST(8 8)
ADD_PORTAGE_TEST(9 8)
ADD_PORTAGE_TEST(10 8)
ADD_PORTAGE_TEST(11 8)
ADD_PORTAGE_TEST(13 1)
