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
    add_executable(sweptfaceapp_jali_singlemat sweptfaceapp_jali_singlemat.cc)
    add_executable(sweptfacedemo_jali_singlemat sweptfacedemo_jali_singlemat.cc)
    include_directories(${Jali_INCLUDE_DIRS} ${Jali_TPL_INCLUDE_DIRS})
    set_target_properties(sweptfaceapp_jali_singlemat
            PROPERTIES
            LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

    set_target_properties(sweptfacedemo_jali_singlemat
            PROPERTIES
            LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

    target_link_libraries(sweptfaceapp_jali_singlemat portage ${EXTRA_LIBS})
    target_link_libraries(sweptfaceapp_jali_singlemat ${LAPACKX_LIBRARIES})
    target_link_libraries(sweptfaceapp_jali_singlemat ${MPI_CXX_LIBRARIES})
    target_link_libraries(sweptfaceapp_jali_singlemat ${Jali_LIBRARIES})
    target_link_libraries(sweptfaceapp_jali_singlemat ${Jali_TPL_LIBRARIES})
    
    target_link_libraries(sweptfacedemo_jali_singlemat portage ${EXTRA_LIBS})
    target_link_libraries(sweptfacedemo_jali_singlemat ${LAPACKX_LIBRARIES})
    target_link_libraries(sweptfacedemo_jali_singlemat ${MPI_CXX_LIBRARIES})
    target_link_libraries(sweptfacedemo_jali_singlemat ${Jali_LIBRARIES})
    target_link_libraries(sweptfacedemo_jali_singlemat ${Jali_TPL_LIBRARIES})


    if (ENABLE_APP_TIMINGS)
      target_compile_definitions(sweptfaceapp_jali_singlemat PRIVATE ENABLE_TIMINGS=1)
      target_compile_definitions(sweptfacedemo_jali_singlemat PRIVATE ENABLE_TIMINGS=1)
    endif()

  endif (Jali_DIR)

  #TCMalloc
  if(ENABLE_TCMALLOC)
    set(TCMALLOC_LIB "${HOME}" CACHE STRING "Set the TCMalloc library")

    if (Jali_DIR)
      target_link_libraries(sweptfaceapp_jali_singlemat ${TCMALLOC_LIB})
      target_link_libraries(sweptfacedemo_jali_singlemat ${TCMALLOC_LIB})
    endif (Jali_DIR)
    
  endif(ENABLE_TCMALLOC)

  # Profiling option (to allow pausing and resuming profiling during execution)
  if(ENABLE_PROFILE)
    set(VTUNE_INCLUDE_DIR "${HOME}" CACHE PATH "Set the locaiton of VTune include")
    set(VTUNE_LIB "${HOME}" CACHE STRING "Set the VTune library")

    include_directories(${VTUNE_INCLUDE_DIR})

    if (Jali_DIR)
      target_link_libraries(sweptfaceapp_jali_singlemat ${VTUNE_LIB})
      target_link_libraries(sweptfacedemo_jali_singlemat ${VTUNE_LIB})
    endif (Jali_DIR)
  endif(ENABLE_PROFILE)

#  if (ENABLE_APP_TESTS)
#    add_subdirectory(test)
#  endif (ENABLE_APP_TESTS)

endif (ENABLE_MPI)
