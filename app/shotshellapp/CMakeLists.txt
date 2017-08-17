#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]





if (Jali_DIR) 
    # Shotshell example
    add_executable(shotshellapp_jali shotshell_jali.cc)
    target_link_libraries(shotshellapp_jali portage ${EXTRA_LIBS})
    if (ENABLE_MPI)
        target_link_libraries(shotshellapp_jali ${MPI_CXX_LIBRARIES})
        set_target_properties(shotshellapp_jali PROPERTIES
                              LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
    endif (ENABLE_MPI)
    target_link_libraries(shotshellapp_jali ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES})

endif (Jali_DIR)

#TCMalloc option
set(ENABLE_TCMALLOC TRUE CACHE BOOL "Use TCMalloc")
if(ENABLE_TCMALLOC)
  set(TCMALLOC_LIB "${HOME}" CACHE STRING "Set the TCMalloc library")

  if (Jali_DIR)
    target_link_libraries(shotshellapp_jali ${TCMALLOC_LIB})
  endif (Jali_DIR)
endif(ENABLE_TCMALLOC)

# Profiling option (to allow pausing and resuming profiling during execution)
set(ENABLE_PROFILE FALSE CACHE BOOL "Enable VTune Profiling")
if(ENABLE_PROFILE)
  set(VTUNE_INCLUDE_DIR "${HOME}" CACHE PATH "Set the locaiton of VTune include")
  set(VTUNE_LIB "${HOME}" CACHE STRING "Set the VTune library")

  include_directories(${VTUNE_INCLUDE_DIR})

  if (Jali_DIR)
    target_link_libraries(shotshellapp_jali ${VTUNE_LIB})
  endif (Jali_DIR)
endif(ENABLE_PROFILE)


if (ENABLE_APP_TESTS)
  add_subdirectory(test)
endif (ENABLE_APP_TESTS)
