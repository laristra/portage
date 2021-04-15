# ------------------------------------------------------------------------------
# This file is part of the Ristra portage project.
# Please see the license file at the root of this repository, or at:
#    https://github.com/laristra/portage/blob/master/LICENSE
# ------------------------------------------------------------------------------

# cannot build jali apps without MPI.
if (WONTON_ENABLE_Jali AND WONTON_ENABLE_MPI)

    # create binary
    add_executable(part-remap part-remap.cc)

    # available build options
    set(ENABLE_PROFILE FALSE CACHE BOOL "Enable VTune profiling")

    # link to portage and its dependencies
    target_link_libraries(part-remap PRIVATE portage)

    # use vtune
    if(ENABLE_PROFILE)
        set(VTUNE_HEADERS "${HOME}" CACHE PATH   "Set VTune headers location")
        set(VTUNE_LIBS    "${HOME}" CACHE STRING "Set VTune libraries location")

        target_include_directories(part-remap PRIVATE ${VTUNE_HEADERS})
        target_link_libraries(part-remap PRIVATE ${VTUNE_LIBS})
    endif()

    # build tests as well
    if (ENABLE_APP_TESTS)
        add_subdirectory(test)
    endif()
endif()
