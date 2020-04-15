#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]


#------------------------------------------------------------------------------#
# Add a rule to build the executable
#------------------------------------------------------------------------------#

# Msm test
if (WONTON_ENABLE_MPI AND WONTON_ENABLE_Jali)

  add_executable(msmapp msmapp.cc)
  target_link_libraries(msmapp portage)

  # Can we do a find_package(vtune)? and just link to an imported library?
  # It would make this a lot less verbose
  
  # Profiling option (to allow pausing and resuming profiling during execution)
  set(ENABLE_PROFILE FALSE CACHE BOOL "Enable VTune Profiling")
  if(ENABLE_PROFILE)
    set(VTUNE_INCLUDE_DIR "${HOME}" CACHE PATH "Set the locaiton of VTune include")
    set(VTUNE_LIB "${HOME}" CACHE STRING "Set the VTune library")

    include_directories(${VTUNE_INCLUDE_DIR})

    target_link_libraries(msmapp ${VTUNE_LIB})
  endif(ENABLE_PROFILE)


  #if (ENABLE_APP_TESTS)
  #  add_subdirectory(test)
  #endif (ENABLE_APP_TESTS)

endif ()
