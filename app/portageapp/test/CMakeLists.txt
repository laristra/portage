#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]
# Config file for a test directory.
# At some future point, much of what appears here will probably
# be pulled into a macro.

message(STATUS "Adding portageapp test")

macro(ADD_SPECIFIC_APPTEST TYPE SUFFIX)
  file(COPY portageapp_${TYPE}_${SUFFIX}.sh
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(GLOB GOLDFILES gold/GOLD_${TYPE}_${SUFFIX}*.txt*)
  file(COPY ${GOLDFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  add_test(NAME portageapp_${TYPE}_${SUFFIX}
    COMMAND ./portageapp_${TYPE}_${SUFFIX}.sh)
  set_property(TEST portageapp_${TYPE}_${SUFFIX}
    PROPERTY ENVIRONMENT
    TESTAPPDIR=${CMAKE_CURRENT_BINARY_DIR}/..
    CMPAPPDIR=${CMAKE_BINARY_DIR}/app/apptest_cmp)
endmacro(ADD_SPECIFIC_APPTEST)

macro(ADD_PORTAGEAPP_TEST DIM ENTITY FIELDORDER REMAPORDER)
  set(SUFFIX "rect_${DIM}d_${ENTITY}_f${FIELDORDER}_r${REMAPORDER}")

  if (Jali_DIR)
    ADD_SPECIFIC_APPTEST(jali ${SUFFIX})
  endif (Jali_DIR)

  if (ENABLE_FleCSI)
    # For now, the flecsi app is hard-coded for 2d, cell-centered
    # remap with the field and remap order being the same
    if (${DIM} MATCHES "2" AND
	${ENTITY} MATCHES "cell" AND
	${FIELDORDER} MATCHES ${REMAPORDER})
      ADD_SPECIFIC_APPTEST(flecsi ${SUFFIX})
    endif()
  endif (ENABLE_FleCSI)
endmacro(ADD_PORTAGEAPP_TEST)

# all of these scripts use mpirun to launch the job, so only include them
# if we know about MPI
if (ENABLE_MPI)
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

  if (Jali_DIR)
    ADD_SPECIFIC_APPTEST(jali_rect 3d_cell_gen_r2)
  endif (Jali_DIR)
endif (ENABLE_MPI)


macro(ADD_UNSTRUC_TEST DIM ENTITY FIELDORDER REMAPORDER)
  set(SUFFIX "${DIM}d_${ENTITY}_f${FIELDORDER}_r${REMAPORDER}")
  file(COPY portageapp_jali_unstruc_${SUFFIX}.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(GLOB GOLDFILES gold/GOLD_jali_unstruc_${SUFFIX}.txt*)
  file(COPY ${GOLDFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

  file(COPY ../../../test_data/shotshell.exo
            ../../../test_data/shotshell-v.exo
	    ../../../test_data/cube-poly1.exo
	    ../../../test_data/cube-poly2.exo
            ../../../test_data/coarse_poly2D.exo
            ../../../test_data/fine_poly2D.exo 
            ../../../test_data/3dvoro-voronoized-10k.exo
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  add_test(NAME portageapp_jali_unstruc_${SUFFIX}
    COMMAND ./portageapp_jali_unstruc_${SUFFIX}.sh)
  set_property(TEST portageapp_jali_unstruc_${SUFFIX}
    PROPERTY ENVIRONMENT
    TESTAPPDIR=${CMAKE_CURRENT_BINARY_DIR}/..
    CMPAPPDIR=${CMAKE_BINARY_DIR}/app/apptest_cmp)
endmacro(ADD_UNSTRUC_TEST)

if (ENABLE_MPI AND Jali_DIR)
  ADD_UNSTRUC_TEST(2 cell 0 1)
  ADD_UNSTRUC_TEST(2 cell 1 1)
  ADD_UNSTRUC_TEST(2 cell 1 2)
  ADD_UNSTRUC_TEST(2 cell 2 2)
  ADD_UNSTRUC_TEST(2 node 0 1)
  ADD_UNSTRUC_TEST(2 node 1 1)
  ADD_UNSTRUC_TEST(2 node 1 2)
  ADD_UNSTRUC_TEST(2 node 2 2)
  ADD_UNSTRUC_TEST(3 cell 0 1)
  ADD_UNSTRUC_TEST(3 cell 1 1)
  ADD_UNSTRUC_TEST(3 cell 1 2)
  ADD_UNSTRUC_TEST(3 cell 2 2)
  ADD_UNSTRUC_TEST(3 node 0 1)
  ADD_UNSTRUC_TEST(3 node 1 1)
  ADD_UNSTRUC_TEST(3 node 1 2)
  ADD_UNSTRUC_TEST(3 node 2 2)
  ADD_UNSTRUC_TEST(2 cell 1 8)
  ADD_UNSTRUC_TEST(3 cell 8 8)

  ADD_SPECIFIC_APPTEST(jali_unstruc 2d_cell_gen_r2)
  ADD_SPECIFIC_APPTEST(jali_unstruc 3d_cell_gen_r2)
endif (ENABLE_MPI AND Jali_DIR)

#multimat test
if (Jali_DIR AND TANGRAM_FOUND)
  file(COPY portageapp_multimat_jali_unstruc_2d_cell_f1_r2.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY gold/GOLD_multimat_jali_unstruc_2d_cell_f1_r2.gmv DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY ../../../test_data/reg10x10.exo
            ../../../test_data/reg10x10.bvf
            ../../../test_data/polymesh_10x10.exo
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  add_test(NAME portageapp_multimat_jali_unstruc_2d_cell_f1_r2
    COMMAND ./portageapp_multimat_jali_unstruc_2d_cell_f1_r2.sh)
  set_property(TEST portageapp_multimat_jali_unstruc_2d_cell_f1_r2
    PROPERTY ENVIRONMENT
    TESTAPPDIR=${CMAKE_CURRENT_BINARY_DIR}/..
    CMPAPPDIR=${CMAKE_BINARY_DIR}/app/apptest_cmp)

endif (Jali_DIR AND TANGRAM_FOUND)

