file(COPY portageapp_test3.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY field_gold.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME portageapp_test3
         COMMAND ./portageapp_test3.sh)
set_property(TEST portageapp_test3
         PROPERTY ENVIRONMENT APPDIR=${CMAKE_BINARY_DIR}/bin)


