# ------------------------------------------------------------------------------
# This file is part of the Ristra portage project.
# Please see the license file at the root of this repository, or at:
#    https://github.com/laristra/portage/blob/master/LICENSE
# ------------------------------------------------------------------------------

# cannot build jali apps without MPI.
if (ENABLE_MPI AND Jali_DIR)

    # create binary
    add_executable(part-remap part-remap.cc)

    # available build options
    set(ENABLE_TCMALLOC TRUE CACHE BOOL "Use TCMalloc")
    set(ENABLE_PROFILE FALSE CACHE BOOL "Enable VTune profiling")

    # link to portage and its dependencies
    target_link_libraries(part-remap PRIVATE portage)
    target_link_libraries(part-remap PRIVATE ${EXTRA_LIBS})
    target_link_libraries(part-remap PRIVATE ${LAPACKX_LIBRARIES})
    target_link_libraries(part-remap PRIVATE ${MPI_CXX_LIBRARIES})
    set_target_properties(part-remap PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

    # link to jali
    target_include_directories(part-remap PRIVATE ${Jali_INCLUDE_DIRS})
    target_include_directories(part-remap PRIVATE ${Jali_TPL_INCLUDE_DIRS})
    target_link_libraries(part-remap PRIVATE ${Jali_LIBRARIES})
    target_link_libraries(part-remap PRIVATE ${Jali_TPL_LIBRARIES})

    # use tc-malloc
    if(ENABLE_TCMALLOC)
        set(TCMALLOC_LIB "${HOME}" CACHE STRING "Set the TCMalloc library")
        target_link_libraries(part-remap PRIVATE ${TCMALLOC_LIB})
    endif()

    # use vtune
    if(ENABLE_PROFILE)
        set(VTUNE_HEADERS "${HOME}" CACHE PATH   "Set VTune headers location")
        set(VTUNE_LIBS    "${HOME}" CACHE STRING "Set VTune libraries location")

        target_include_directories(part-remap PRIVATE ${VTUNE_HEADERS})
        target_link_libraries(part-remap PRIVATE ${VTUNE_LIBS})
    endif()

    # build tests as well
#    if (ENABLE_APP_TESTS)
#        add_subdirectory(test)
#    endif()
endif()