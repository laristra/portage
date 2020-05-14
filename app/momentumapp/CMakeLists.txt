# ------------------------------------------------------------------------------
# This file is part of the Ristra portage project.
# Please see the license file at the root of this repository, or at:
#    https://github.com/laristra/portage/blob/master/LICENSE
# ------------------------------------------------------------------------------

# cannot build jali apps without MPI.
if (WONTON_ENABLE_MPI AND WONTON_ENABLE_Jali)

  #
  # create 2D binary
  #
  add_executable(momentumapp2D momentumapp2D.cc)

  # available build options
  set(ENABLE_PROFILE FALSE CACHE BOOL "Enable VTune profiling")

  # link to portage and its dependencies
  target_link_libraries(momentumapp2D PRIVATE portage)

  # use vtune
  if(ENABLE_PROFILE)
    set(VTUNE_HEADERS "${HOME}" CACHE PATH   "Set VTune headers location")
    set(VTUNE_LIBS    "${HOME}" CACHE STRING "Set VTune libraries location")

    target_include_directories(momentumapp2D PRIVATE ${VTUNE_HEADERS})
    target_link_libraries(momentumapp2D PRIVATE ${VTUNE_LIBS})
  endif()

  #
  # create 3D binary
  #
  add_executable(momentumapp3D momentumapp3D.cc)

  # available build options
  set(ENABLE_PROFILE FALSE CACHE BOOL "Enable VTune profiling")

  # link to portage and its dependencies
  target_link_libraries(momentumapp3D PRIVATE portage)

  # use vtune
  if(ENABLE_PROFILE)
    set(VTUNE_HEADERS "${HOME}" CACHE PATH   "Set VTune headers location")
    set(VTUNE_LIBS    "${HOME}" CACHE STRING "Set VTune libraries location")

    target_include_directories(momentumapp3D PRIVATE ${VTUNE_HEADERS})
    target_link_libraries(momentumapp3D PRIVATE ${VTUNE_LIBS})
  endif()

  #
  # MM : create 2D binary
  #
  if (PORTAGE_ENABLE_TANGRAM)
    add_executable(momentumapp2D_mm momentumapp2D_mm.cc)

    # available build options
    set(ENABLE_PROFILE FALSE CACHE BOOL "Enable VTune profiling")

    # link to portage and its dependencies
    target_link_libraries(momentumapp2D_mm PRIVATE portage)

    # link to jali
    target_include_directories(momentumapp2D_mm PRIVATE ${Jali_INCLUDE_DIRS})
    target_include_directories(momentumapp2D_mm PRIVATE ${Jali_TPL_INCLUDE_DIRS})
    target_link_libraries(momentumapp2D_mm PRIVATE ${Jali_LIBRARIES})
    target_link_libraries(momentumapp2D_mm PRIVATE ${Jali_TPL_LIBRARIES})
  endif()

  if(ENABLE_APP_TESTS)
    add_subdirectory(test)
  endif()
endif()
