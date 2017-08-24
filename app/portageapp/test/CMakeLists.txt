#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]
# Config file for a test directory.
# At some future point, much of what appears here will probably
# be pulled into a macro.

message(STATUS "Adding portageapp test")

macro(ADD_PORTAGEAPP_TEST DIM ENTITY FIELDORDER REMAPORDER)
  set(SUFFIX "${DIM}d_${ENTITY}_f${FIELDORDER}_r${REMAPORDER}")

  file(COPY portageapp_${SUFFIX}.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(GLOB GOLDFILES gold/GOLD-field_${SUFFIX}*.txt*)
  file(COPY ${GOLDFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  add_test(NAME portageapp_${SUFFIX} COMMAND ./portageapp_${SUFFIX}.sh)
  set_property(TEST portageapp_${SUFFIX}
    PROPERTY ENVIRONMENT APPDIR=${CMAKE_BINARY_DIR}/bin)
endmacro(ADD_PORTAGEAPP_TEST)




ADD_PORTAGEAPP_TEST(2 cell 0 1)
ADD_PORTAGEAPP_TEST(2 cell 1 1)
ADD_PORTAGEAPP_TEST(2 cell 1 2)
ADD_PORTAGEAPP_TEST(2 cell 2 2)
ADD_PORTAGEAPP_TEST(2 node 0 1)
ADD_PORTAGEAPP_TEST(2 node 1 1)
ADD_PORTAGEAPP_TEST(2 node 1 2)
ADD_PORTAGEAPP_TEST(2 node 2 2)
ADD_PORTAGEAPP_TEST(3 cell 0 1)
ADD_PORTAGEAPP_TEST(3 cell 1 1)
ADD_PORTAGEAPP_TEST(3 cell 1 2)
ADD_PORTAGEAPP_TEST(3 cell 2 2)
ADD_PORTAGEAPP_TEST(3 node 0 1)
ADD_PORTAGEAPP_TEST(3 node 1 1)
ADD_PORTAGEAPP_TEST(3 node 1 2)
ADD_PORTAGEAPP_TEST(3 node 2 2)
