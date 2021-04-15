# ------------------------------------------------------------------------------
# This file is part of the Ristra portage project.
# Please see the license file at the root of this repository, or at:
#    https://github.com/laristra/portage/blob/master/LICENSE
# ------------------------------------------------------------------------------

if (WONTON_ENABLE_MPI)
  # create binary
  add_executable(particle-analysis remap.cc)

  # link to portage and its dependencies
  target_link_libraries(particle-analysis PRIVATE portage)

  # build tests as well
  if (ENABLE_APP_TESTS)
    add_subdirectory(test)
  endif()
endif()