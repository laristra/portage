

#-----------------------------------------------------------------------------~#


#-----------------------------------------------------------------------------~#

if (Jali_DIR)
    set(jali_wrapper_headers
        mesh/jali/jali_mesh_wrapper.h
        state/jali/jali_state_wrapper.h
        )
endif (Jali_DIR)

if (ENABLE_FleCSI)
    set(flecsi_wrapper_headers
        mesh/flecsi/flecsi_mesh_wrapper.h
        state/flecsi/flecsi_state_wrapper.h
        )
endif ()

set(wrappers_HEADERS
    ${jali_wrapper_headers}
    ${flecsi_wrapper_headers}
    mesh/AuxMeshTopology.h
    mesh/simple_mesh/simple_mesh_wrapper.h
    state/simple_state/simple_state_wrapper.h
    mesh/flat/flat_mesh_wrapper.h
    state/flat/flat_state_wrapper.h
    PARENT_SCOPE)

if (ENABLE_UNIT_TESTS)
   if (Jali_DIR)
      cinch_add_unit(test_jali_mesh_wrapper
                    SOURCES mesh/jali/test/test_jali_mesh_wrapper.cc
                    LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
                    # this is relative to the wrapper directory
                    INPUTS ../../test_data/shotshell.exo
                    ../../test_data/shotshell-v.exo
                    POLICY MPI
                    THREADS 1)

      cinch_add_unit(test_jali_state_wrapper
                     SOURCES state/jali/test/test_jali_state_wrapper.cc
                     LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
                     POLICY MPI
                     THREADS 1)

      # Flat mesh wrapper
      cinch_add_unit(test_flat_mesh_wrapper
                     SOURCES mesh/flat/test/test_flat_mesh_wrapper.cc
                     LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
                     # this is relative to the wrapper directory
                     INPUTS ../../test_data/shotshell.exo
                     ../../test_data/shotshell-v.exo
                     POLICY MPI
                     THREADS 1)

      # Flat state wrapper
      cinch_add_unit(test_flat_state_wrapper
                     SOURCES state/flat/test/test_flat_state_wrapper.cc
                     LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
                     POLICY MPI
                     THREADS 1)
   endif (Jali_DIR)

    if (ENABLE_FleCSI)
      cinch_add_unit(flecsi_2d
                     SOURCES mesh/flecsi/test/flecsi_2d.cc
                     LIBRARIES portage ${FleCSI_LIBRARIES} ${FleCSISP_LIBRARIES}
                     POLICY MPI
                     THREADS 1)
      
      cinch_add_unit(flecsi_3d
                     SOURCES mesh/flecsi/test/flecsi_3d.cc
                     LIBRARIES portage ${FleCSI_LIBRARIES} ${FleCSISP_LIBRARIES}
                     POLICY MPI
                     THREADS 1)

      #cinch_add_unit(test_flecsi_mesh_wrapper
      #               SOURCES mesh/flecsi/test/test_flecsi_mesh_wrapper.cc
      #               LIBRARIES portage ${FLECSI_LIBRARIES})
    endif ()

   # Simple Mesh wrapper
   cinch_add_unit(test_simple_mesh_wrapper
     SOURCES mesh/simple_mesh/test/test_simple_mesh_wrapper.cc
     LIBRARIES portage
     POLICY MPI
     THREADS 1)
   # Simple State wrapper
   cinch_add_unit(test_simple_state_wrapper
     SOURCES state/simple_state/test/test_simple_state_wrapper.cc
     LIBRARIES portage
     POLICY MPI
     THREADS 1)
endif(ENABLE_UNIT_TESTS)

