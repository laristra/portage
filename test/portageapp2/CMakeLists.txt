file(COPY portageapp_test2.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY field_gold.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME portageapp_test2
         COMMAND ./portageapp_test2.sh)
set_property(TEST portageapp_test2
         PROPERTY ENVIRONMENT APPDIR=${CMAKE_BINARY_DIR}/bin)
set_property(TEST portageapp_test2
         PROPERTY PROCESSORS 4)
