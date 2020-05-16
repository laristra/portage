#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]

# Profiling option (to allow pausing and resuming profiling during execution)
if (ENABLE_PROFILE)
  set(VTUNE_INCLUDE_DIR "${HOME}" CACHE PATH "Set the locaiton of VTune include")
  set(VTUNE_LIB "${HOME}" CACHE STRING "Set the VTune library")
endif ()  

# We can't compile any Jali apps without MPI
if (WONTON_ENABLE_Jali AND WONTON_ENABLE_MPI)
  # options
  set(ENABLE_PROFILE FALSE CACHE BOOL "Enable VTune Profiling")
  option(FIXED_JALI2FLECSI "Is jali to flecsi fixed" OFF) # where is this used?
  
  add_executable(portageapp_jali portageapp_jali.cc)
  target_link_libraries(portageapp_jali portage)

  if (ENABLE_APP_TIMINGS)
    target_compile_definitions(portageapp_jali PRIVATE ENABLE_TIMINGS=1)
  endif()
  
  if (ENABLE_PROFILE)
    target_link_libraries(portageapp_jali PRIVATE ${VTUNE_LIB})
  endif ()


  if (PORTAGE_ENABLE_TANGRAM)
    add_executable(portageapp_rgmd_jali portageapp_rgmd_jali.cc)    
    target_link_libraries(portageapp_rgmd_jali PRIVATE portage)
    
    if (ENABLE_APP_TIMINGS)
      target_compile_definitions(portageapp_rgmd_jali PRIVATE ENABLE_TIMINGS=1)
    endif()
    
    if (ENABLE_PROFILE)
      target_link_libraries(portageapp_rgmd_jali PRIVATE ${VTUNE_LIB})
    endif ()
    
    if (TANGRAM_ENABLE_XMOF2D)
      target_link_libraries(portageapp_rgmd_jali PRIVATE ${XMOF2D_LIBRARIES})

      add_executable(portageapp_multimat_jali portageapp_multimat_jali.cc)
      target_link_libraries(portageapp_multimat_jali PRIVATE portage)
      target_link_libraries(portageapp_multimat_jali PRIVATE ${XMOF2D_LIBRARIES})

      if (ENABLE_APP_TIMINGS)
        target_compile_definitions(portageapp_multimat_jali PRIVATE ENABLE_TIMINGS=1)
      endif()

      if (ENABLE_PROFILE)
        target_link_libraries(portageapp_multimat_jali PRIVATE ${VTUNE_LIB})
        target_include_directories(portageapp_multimat_jali PRIVATE ${VTUNE_INCLUDE_DIR})
      endif ()
    endif ()
  endif ()

  # HAS BEEN BROKEN FOR A WHILE UNBEKNOWNST TO US - FIX LATER
  if (WONTON_ENABLE_FleCSI AND FIXED_JALI2FLECSI)
    add_executable(portageapp_jali2flecsi portageapp_jali2flecsi.cc)
    target_link_libraries(portageapp_jali2flecsi PRIVATE portage)
  endif ()

endif ()

# FleCSI example
if (WONTON_ENABLE_FleCSI)
  add_executable(portageapp_flecsi portageapp_flecsi.cc)
  target_link_libraries(portageapp_flecsi PRIVATE portage)
endif ()


if (ENABLE_APP_TESTS)
  add_subdirectory(test)
endif (ENABLE_APP_TESTS)

