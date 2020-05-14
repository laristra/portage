# ------------------------------------------------------------------------------
# This file is part of the Ristra portage project.
# Please see the license file at the root of this repository, or at:
#    https://github.com/laristra/portage/blob/master/LICENSE
# ------------------------------------------------------------------------------

# cannot build jali apps without MPI.
if (WONTON_ENABLE_MPI AND WONTON_ENABLE_Jali)

  # create targets
  set(target swept_face_app)
  add_executable(${target} swept_face_app.cc)

  # available build options
  set(ENABLE_PROFILE FALSE CACHE BOOL "Enable VTune profiling")
  option(FIXED_JALI2FLECSI "Is jali to flecsi fixed" OFF)  # Where is this used?

  # link to portage and its dependencies
  target_link_libraries(${target} PRIVATE portage)

  # use vtune
  if(ENABLE_PROFILE)
    set(VTUNE_HEADERS "${HOME}" CACHE PATH   "Set VTune headers location")
    set(VTUNE_LIBS    "${HOME}" CACHE STRING "Set VTune libraries location")

    target_include_directories(${target} PRIVATE ${VTUNE_HEADERS})
    target_link_libraries(${target} PRIVATE ${VTUNE_LIBS})
  endif()

   # build tests as well
  if (ENABLE_APP_TESTS)
    add_subdirectory(test)
  endif()
endif()
