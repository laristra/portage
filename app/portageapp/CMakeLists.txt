#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]

# We can't compile any Jali apps without MPI
if (ENABLE_MPI)
  # options
  set(ENABLE_TCMALLOC TRUE CACHE BOOL "Use TCMalloc")
  set(ENABLE_PROFILE FALSE CACHE BOOL "Enable VTune Profiling")
  option(FIXED_JALI2FLECSI "Is jali to flecsi fixed" OFF)

  if (Jali_DIR)
    add_executable(portageapp_jali portageapp_jali.cc)
    include_directories(${Jali_INCLUDE_DIRS} ${Jali_TPL_INCLUDE_DIRS})
    set_target_properties(portageapp_jali
            PROPERTIES
            LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

    target_link_libraries(portageapp_jali portage ${EXTRA_LIBS})
    target_link_libraries(portageapp_jali ${LAPACKX_LIBRARIES})
    target_link_libraries(portageapp_jali ${MPI_CXX_LIBRARIES})
    target_link_libraries(portageapp_jali ${Jali_LIBRARIES})
    target_link_libraries(portageapp_jali ${Jali_TPL_LIBRARIES})

    if (ENABLE_APP_TIMINGS)
      target_compile_definitions(portageapp_jali PRIVATE ENABLE_TIMINGS=1)
    endif()

    if (TANGRAM_FOUND)
      add_executable(portageapp_rgmd_jali portageapp_rgmd_jali.cc)

      set_target_properties(portageapp_rgmd_jali
              PROPERTIES
              LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

      target_link_libraries(portageapp_rgmd_jali portage ${EXTRA_LIBS})
      target_link_libraries(portageapp_rgmd_jali ${LAPACKX_LIBRARIES})
      target_link_libraries(portageapp_rgmd_jali ${TANGRAM_LIBRARIES})
      target_link_libraries(portageapp_rgmd_jali ${MPI_CXX_LIBRARIES})
      target_link_libraries(portageapp_rgmd_jali ${Jali_LIBRARIES})
      target_link_libraries(portageapp_rgmd_jali ${Jali_TPL_LIBRARIES})

      if (ENABLE_APP_TIMINGS)
        target_compile_definitions(portageapp_rgmd_jali PRIVATE ENABLE_TIMINGS=1)
      endif()

      if (XMOF2D_FOUND)
        target_link_libraries(portageapp_rgmd_jali ${XMOF2D_LIBRARIES})

        add_executable(portageapp_multimat_jali portageapp_multimat_jali.cc)
        set_target_properties(portageapp_multimat_jali
                PROPERTIES
                LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

        target_link_libraries(portageapp_multimat_jali portage ${EXTRA_LIBS})
        target_link_libraries(portageapp_multimat_jali ${LAPACKX_LIBRARIES})
        target_link_libraries(portageapp_multimat_jali ${TANGRAM_LIBRARIES})
        target_link_libraries(portageapp_multimat_jali ${XMOF2D_LIBRARIES})
        target_link_libraries(portageapp_multimat_jali ${MPI_CXX_LIBRARIES})
        target_link_libraries(portageapp_multimat_jali ${Jali_LIBRARIES})
        target_link_libraries(portageapp_multimat_jali ${Jali_TPL_LIBRARIES})
        if (ENABLE_APP_TIMINGS)
          target_compile_definitions(portageapp_multimat_jali PRIVATE ENABLE_TIMINGS=1)
        endif()
      endif (XMOF2D_FOUND)
    endif (TANGRAM_FOUND)

    # HAS BEEN BROKEN FOR A WHILE UNBEKNOWNST TO US - FIX LATER
    if (ENABLE_FleCSI AND FIXED_JALI2FLECSI)
      add_executable(portageapp_jali2flecsi portageapp_jali2flecsi.cc)
      include_directories(${FleCSI_INCLUDE_DIR} ${FleCSISP_INCLUDE_DIR})
      set_target_properties(portageapp_jali2flecsi
              PROPERTIES
              LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

      target_link_libraries(portageapp_jali2flecsi portage ${EXTRA_LIBS})
      target_link_libraries(portageapp_jali2flecsi ${LAPACKX_LIBRARIES})
      target_link_libraries(portageapp_jali2flecsi ${FleCSI_LIBRARIES})
      target_link_libraries(portageapp_jali2flecsi ${FleCSISP_LIBRARIES})
      target_link_libraries(portageapp_jali2flecsi ${MPI_CXX_LIBRARIES})
      target_link_libraries(portageapp_jali2flecsi ${Jali_LIBRARIES})
      target_link_libraries(portageapp_jali2flecsi ${Jali_TPL_LIBRARIES})
    endif (ENABLE_FleCSI AND FIXED_JALI2FLECSI)
  endif (Jali_DIR)

  # FleCSI example
  if (ENABLE_FleCSI)
    add_executable(portageapp_flecsi portageapp_flecsi.cc)
    include_directories(${FleCSI_INCLUDE_DIR} ${FleCSISP_INCLUDE_DIR})
    set_target_properties(portageapp_flecsi
            PROPERTIES
            LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

    target_link_libraries(portageapp_flecsi portage ${EXTRA_LIBS})
    target_link_libraries(portageapp_flecsi ${LAPACKX_LIBRARIES})
    target_link_libraries(portageapp_flecsi ${FleCSI_LIBRARIES})
    target_link_libraries(portageapp_flecsi ${FleCSISP_LIBRARIES})
    target_link_libraries(portageapp_flecsi ${MPI_CXX_LIBRARIES})
  endif (ENABLE_FleCSI)


  #TCMalloc
  if(ENABLE_TCMALLOC)
    set(TCMALLOC_LIB "${HOME}" CACHE STRING "Set the TCMalloc library")

    if (Jali_DIR)
      target_link_libraries(portageapp_jali ${TCMALLOC_LIB})
      if (TANGRAM_FOUND)
        target_link_libraries(portageapp_rgmd_jali ${TCMALLOC_LIB})

        if (XMOF2D_FOUND)
          target_link_libraries(portageapp_multimat_jali ${TCMALLOC_LIB})
        endif (XMOF2D_FOUND)
      endif (TANGRAM_FOUND)
    endif (Jali_DIR)
    
    if (ENABLE_FleCSI)
      target_link_libraries(portageapp_flecsi ${TCMALLOC_LIB})
    endif (ENABLE_FleCSI)
  endif(ENABLE_TCMALLOC)

  # Profiling option (to allow pausing and resuming profiling during execution)
  if(ENABLE_PROFILE)
    set(VTUNE_INCLUDE_DIR "${HOME}" CACHE PATH "Set the locaiton of VTune include")
    set(VTUNE_LIB "${HOME}" CACHE STRING "Set the VTune library")

    include_directories(${VTUNE_INCLUDE_DIR})

    if (Jali_DIR)
      target_link_libraries(portageapp_jali ${VTUNE_LIB})
      if (TANGRAM_FOUND)
        target_link_libraries(portageapp_rgmd_jali ${VTUNE_LIB})
        if (XMOF2D_FOUND)
          target_link_libraries(portageapp_multimat_jali ${VTUNE_LIB})
        endif (XMOF2D_FOUND)
      endif(TANGRAM_FOUND)
    endif (Jali_DIR)
    if (ENABLE_FleCSI)
      target_link_libraries(portageapp_flecsi ${VTUNE_LIB})
    endif (ENABLE_FleCSI)
  endif(ENABLE_PROFILE)

  if (ENABLE_APP_TESTS)
    add_subdirectory(test)
  endif (ENABLE_APP_TESTS)

endif (ENABLE_MPI)
