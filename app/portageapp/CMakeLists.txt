#~----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

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
endif (Jali_DIR)

# FleCSI example
if (FLECSI_INSTALL_DIR)
  add_executable(portageapp_flecsi portageapp/portageapp_flecsi.cc)
  target_link_libraries(portageapp_flecsi portage ${EXTRA_LIBS})
  include_directories(${FLECSI_INCLUDE_DIRS})
  target_link_libraries(portageapp_flecsi ${FLECSI_LIBRARIES})
# these appear needed because of something we are pulling in within src
  if (ENABLE_MPI)
      target_link_libraries(portageapp_flecsi ${MPI_CXX_LIBRARIES})
      set_target_properties(portageapp_flecsi PROPERTIES
                            LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
  endif (ENABLE_MPI)
endif (FLECSI_INSTALL_DIR)

if (JALI_DIR AND FLECSI_INSTALL_DIR)
  add_executable(portageapp_jali2flecsi portageapp/portageapp_jali2flecsi.cc)
  target_link_libraries(portageapp_jali2flecsi portageapp_jali2flecsi ${EXTRA_LIBS})
  include_directories(${FLECSI_INCLUDE_DIRS})
  target_link_libraries(portageapp_jali2flecsi ${FLECSI_LIBRARIES})
# these appear needed because of something we are pulling in within src
  if (ENABLE_MPI)
      target_link_libraries(portageapp_jali2flecsi ${MPI_CXX_LIBRARIES})
      set_target_properties(portageapp_jali2flecsi PROPERTIES
                            LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
  endif (ENABLE_MPI)
  target_link_libraries(portageapp_jali2flecsi ${Jali_LIBRARIES}
                        ${Jali_TPL_LIBRARIES})

endif (JALI_DIR AND FLECSI_INSTALL_DIR)

#TCMalloc option
set(ENABLE_TCMALLOC TRUE CACHE BOOL "Use TCMalloc")
if(ENABLE_TCMALLOC)
  set(TCMALLOC_LIB "${HOME}" CACHE STRING "Set the TCMalloc library")

  if (Jali_DIR)
    target_link_libraries(portageapp_jali ${TCMALLOC_LIB})
  endif (Jali_DIR)
  if (FLECSI_INSTALL_DIR)
    target_link_libraries(portageapp_flecsi ${TCMALLOC_LIB})
  endif (FLECSI_INSTALL_DIR)
  if (JALI_DIR AND FLECSI_INSTALL_DIR)
    target_link_libraries(portageapp_jali2flecsi ${TCMALLOC_LIB})
  endif (JALI_DIR AND FLECSI_INSTALL_DIR)
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
  if (FLECSI_INSTALL_DIR)
    target_link_libraries(portageapp_flecsi ${VTUNE_LIB})
  endif (FLECSI_INSTALL_DIR)
  if (JALI_DIR AND FLECSI_INSTALL_DIR)
    target_link_libraries(portageapp_jali2flecsi ${VTUNE_LIB})
  endif (JALI_DIR AND FLECSI_INSTALL_DIR)
endif(ENABLE_PROFILE)


if (ENABLE_APP_TESTS)
  add_subdirectory(test)
endif (ENABLE_APP_TESTS)
