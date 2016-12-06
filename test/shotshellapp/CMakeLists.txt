message(STATUS "Adding shotshellapp test")

macro(ADD_SHOTSHELLAPP_TEST TEST_NUM NUM_PROC)
  file(COPY shotshellapp_test${TEST_NUM}.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY field_gold${TEST_NUM}.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  add_test(NAME shotshellapp_test${TEST_NUM}
    COMMAND ./shotshellapp_test${TEST_NUM}.sh)
  set_property(TEST shotshellapp_test${TEST_NUM}
    PROPERTY ENVIRONMENT APPDIR=${CMAKE_BINARY_DIR}/bin)
  set_property(TEST shotshellapp_test${TEST_NUM}
    PROPERTY PROCESSORS ${NUM_PROC})
endmacro(ADD_SHOTSHELLAPP_TEST)

ADD_SHOTSHELLAPP_TEST(0 1)
ADD_SHOTSHELLAPP_TEST(1 1)
ADD_SHOTSHELLAPP_TEST(2 1)
ADD_SHOTSHELLAPP_TEST(3 1)
