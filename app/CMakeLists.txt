#~----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

#------------------------------------------------------------------------------#
# Add a rule to build the executable
#------------------------------------------------------------------------------#

add_executable(portageapp main.cc)
target_link_libraries(portageapp portage ${EXTRA_LIBS})

# MPI additions - need to push this down into cinch
target_link_libraries(portageapp ${MPI_CXX_LIBRARIES})
set_target_properties(portageapp PROPERTIES
                      LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

# Jali additions
include_directories(${Jali_INCLUDE_DIRS} ${Jali_TPL_INCLUDE_DIRS})
target_link_libraries(portageapp ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES})

# Shotshell example
add_executable(shotshellapp shotshell.cc)
target_link_libraries(shotshellapp portage ${EXTRA_LIBS})
target_link_libraries(shotshellapp ${MPI_CXX_LIBRARIES})
set_target_properties(shotshellapp PROPERTIES
                      LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
target_link_libraries(shotshellapp ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES})

# timing example
add_executable(timingapp timing_test.cc)
target_link_libraries(timingapp portage ${EXTRA_LIBS})
target_link_libraries(timingapp ${MPI_CXX_LIBRARIES})
set_target_properties(timingapp PROPERTIES
                      LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
target_link_libraries(timingapp ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES})

# FleCSI example
if(FLECSI_INSTALL_DIR)
  add_executable(portage_flecsiapp main_flecsi.cc)
  target_link_libraries(portage_flecsiapp portage ${EXTRA_LIBS})
  include_directories(${FLECSI_INCLUDE_DIRS})
  target_link_libraries(portage_flecsiapp ${FLECSI_LIBRARIES})
# these appear needed because of something we are pulling in within src
  target_link_libraries(portage_flecsiapp ${MPI_CXX_LIBRARIES})
  set_target_properties(portage_flecsiapp PROPERTIES
                        LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

  add_executable(portage_jali2flecsiapp main_jali_to_flecsi.cc)
  target_link_libraries(portage_jali2flecsiapp portage ${EXTRA_LIBS})
  include_directories(${FLECSI_INCLUDE_DIRS})
  target_link_libraries(portage_jali2flecsiapp ${FLECSI_LIBRARIES})
# these appear needed because of something we are pulling in within src
  target_link_libraries(portage_jali2flecsiapp ${MPI_CXX_LIBRARIES})
  target_link_libraries(portage_jali2flecsiapp ${Jali_LIBRARIES}
                        ${Jali_TPL_LIBRARIES})
  set_target_properties(portage_jali2flecsiapp PROPERTIES
                        LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

endif(FLECSI_INSTALL_DIR)

#TCMalloc option
set(ENABLE_TCMALLOC TRUE CACHE BOOL "Use TCMalloc")
if(ENABLE_TCMALLOC)
  set(TCMALLOC_LIB "${HOME}" CACHE STRING "Set the TCMalloc library")
  target_link_libraries(portageapp ${TCMALLOC_LIB})
  target_link_libraries(shotshellapp ${TCMALLOC_LIB})
  target_link_libraries(timingapp ${TCMALLOC_LIB})
endif()

# Profiling option (to allow pausing and resuming profiling during execution)
set(ENABLE_PROFILE FALSE CACHE BOOL "Enable VTune Profiling")
if(ENABLE_PROFILE)
  set(VTUNE_INCLUDE_DIR "${HOME}" CACHE PATH "Set the locaiton of VTune include")
  set(VTUNE_LIB "${HOME}" CACHE STRING "Set the VTune library")

  include_directories(${VTUNE_INCLUDE_DIR})
  target_link_libraries(portageapp ${VTUNE_LIB})
  target_link_libraries(shotshellapp ${VTUNE_LIB})
endif()
