#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]

  message(STATUS "TANGRAM_FOUND value in app/portageapp/CMakeLists.txt is ${TANGRAM_FOUND}")
  message(STATUS "TANGRAM_INCLUDE_DIRS value in app/portageapp/CMakeLists.txt is ${TANGRAM_INCLUDE_DIRS}")

if (Jali_DIR)
    add_executable(portageapp_jali portageapp_jali.cc)
    target_link_libraries(portageapp_jali portage ${EXTRA_LIBS})

    # MPI additions - need to push this down into cinch
    if (ENABLE_MPI)
        target_link_libraries(portageapp_jali ${MPI_CXX_LIBRARIES})
        set_target_properties(portageapp_jali PROPERTIES
                              LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
    endif (ENABLE_MPI)

    # Jali additions
    include_directories(${Jali_INCLUDE_DIRS} ${Jali_TPL_INCLUDE_DIRS})
    target_link_libraries(portageapp_jali ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES})

    add_executable(portageapp_jali_file portageapp_jali_file.cc)
    target_link_libraries(portageapp_jali_file portage ${EXTRA_LIBS})

    # MPI additions - need to push this down into cinch
    if (ENABLE_MPI)
        target_link_libraries(portageapp_jali_file ${MPI_CXX_LIBRARIES})
        set_target_properties(portageapp_jali_file PROPERTIES
                              LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
    endif (ENABLE_MPI)

    # Jali additions
    include_directories(${Jali_INCLUDE_DIRS} ${Jali_TPL_INCLUDE_DIRS})
    target_link_libraries(portageapp_jali_file ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES})

  if (TANGRAM_FOUND)    # tangram is header only package
    include_directories(${TANGRAM_INCLUDE_DIRS})
    if (XMOF2D_FOUND)   # Right now we are using tangram with XMOF2D
      include_directories(${XMOF2D_INCLUDE_DIRS})
      target_link_libraries(portageapp_jali ${XMOF2D_LIBRARIES})
      target_link_libraries(portageapp_jali_file ${XMOF2D_LIBRARIES})
    endif (XMOF2D_FOUND)
  endif (TANGRAM_FOUND)

endif (Jali_DIR)

# FleCSI example
if (ENABLE_FleCSI)
  add_executable(portageapp_flecsi portageapp_flecsi.cc)
  target_link_libraries(portageapp_flecsi portage ${EXTRA_LIBS})
  include_directories(${FleCSI_INCLUDE_DIR} ${FleCSISP_INCLUDE_DIR})
  target_link_libraries(portageapp_flecsi ${FleCSI_LIBRARIES} ${FleCSISP_LIBRARIES})
# these appear needed because of something we are pulling in within src
  if (ENABLE_MPI)
      target_link_libraries(portageapp_flecsi ${MPI_CXX_LIBRARIES})
      set_target_properties(portageapp_flecsi PROPERTIES
                            LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
  endif (ENABLE_MPI)

  if (TANGRAM_FOUND)    # tangram is header only package
    include_directories(${TANGRAM_INCLUDE_DIRS})
    if (XMOF2D_FOUND)   # Right now we are using tangram with XMOF2D
      target_link_libraries(portageapp_flecsi ${XMOF2D_LIBRARIES})
    endif (XMOF2D_FOUND)
  endif (TANGRAM_FOUND)
endif (ENABLE_FleCSI)

if (JALI_DIR AND ENABLE_FleCSI)
  add_executable(portageapp_jali2flecsi portageapp_jali2flecsi.cc)
  target_link_libraries(portageapp_jali2flecsi portageapp_jali2flecsi ${EXTRA_LIBS})
  include_directories(${FleCSI_INCLUDE_DIR} ${FleCSISP_INCLUDE_DIR})
  target_link_libraries(portageapp_jali2flecsi ${FleCSI_LIBRARIES} ${FleCSISP_LIBRARIES})
# these appear needed because of something we are pulling in within src
  if (ENABLE_MPI)
      target_link_libraries(portageapp_jali2flecsi ${MPI_CXX_LIBRARIES})
      set_target_properties(portageapp_jali2flecsi PROPERTIES
                            LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
  endif (ENABLE_MPI)
  target_link_libraries(portageapp_jali2flecsi ${Jali_LIBRARIES}
                        ${Jali_TPL_LIBRARIES})

  if (TANGRAM_FOUND)    # tangram is header only package
    include_directories(${TANGRAM_INCLUDE_DIRS})
    if (XMOF2D_FOUND)   # Right now we are using tangram with XMOF2D
      target_link_libraries(portageapp_jali2flecsi ${XMOF2D_LIBRARIES})
    endif (XMOF2D_FOUND)
  endif (TANGRAM_FOUND)
endif (JALI_DIR AND ENABLE_FleCSI)


#TCMalloc option
set(ENABLE_TCMALLOC TRUE CACHE BOOL "Use TCMalloc")
if(ENABLE_TCMALLOC)
  set(TCMALLOC_LIB "${HOME}" CACHE STRING "Set the TCMalloc library")

  if (Jali_DIR)
    target_link_libraries(portageapp_jali ${TCMALLOC_LIB})
  endif (Jali_DIR)
  if (ENABLE_FleCSI)
    target_link_libraries(portageapp_flecsi ${TCMALLOC_LIB})
  endif (ENABLE_FleCSI)
  if (JALI_DIR AND ENABLE_FleCSI)
    target_link_libraries(portageapp_jali2flecsi ${TCMALLOC_LIB})
  endif (JALI_DIR AND ENABLE_FleCSI)
endif(ENABLE_TCMALLOC)

# Profiling option (to allow pausing and resuming profiling during execution)
set(ENABLE_PROFILE FALSE CACHE BOOL "Enable VTune Profiling")
if(ENABLE_PROFILE)
  set(VTUNE_INCLUDE_DIR "${HOME}" CACHE PATH "Set the locaiton of VTune include")
  set(VTUNE_LIB "${HOME}" CACHE STRING "Set the VTune library")

  include_directories(${VTUNE_INCLUDE_DIR})

  if (Jali_DIR)
    target_link_libraries(portageapp_jali ${VTUNE_LIB})
  endif (Jali_DIR)
  if (ENABLE_FleCSI)
    target_link_libraries(portageapp_flecsi ${VTUNE_LIB})
  endif (ENABLE_FleCSI)
  if (JALI_DIR AND ENABLE_FleCSI)
    target_link_libraries(portageapp_jali2flecsi ${VTUNE_LIB})
  endif (JALI_DIR AND ENABLE_FleCSI)
endif(ENABLE_PROFILE)

if (ENABLE_APP_TESTS)
  add_subdirectory(test)
endif (ENABLE_APP_TESTS)
