#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]

#-------------------------------------------------------------------------------
# intersect_reconstruct_app - uses Jali and XMOF2D
#-------------------------------------------------------------------------------
add_executable(intersect_reconstruct_app intersect_reconstruct_app.cc)
target_link_libraries(intersect_reconstruct_app
                      portage ${LAPACKE_LIBRARIES} ${XMOF2D_LIBRARIES} 
                      ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
                      ${MPI_CXX_LIBRARIES} ${EXTRA_LIBS})
include_directories(${Jali_INCLUDE_DIRS} ${Jali_TPL_INCLUDE_DIRS})
#TCMalloc option
set(ENABLE_TCMALLOC TRUE CACHE BOOL "Use TCMalloc")
if(ENABLE_TCMALLOC)
  set(TCMALLOC_LIB "${HOME}" CACHE STRING "Set the TCMalloc library")

  target_link_libraries(intersect_reconstruct_app ${TCMALLOC_LIB})
endif(ENABLE_TCMALLOC)

# Profiling option (to allow pausing and resuming profiling during execution)
set(ENABLE_PROFILE FALSE CACHE BOOL "Enable VTune Profiling")
if(ENABLE_PROFILE)
  set(VTUNE_INCLUDE_DIR "${HOME}" CACHE PATH "Set the locaiton of VTune include")
  set(VTUNE_LIB "${HOME}" CACHE STRING "Set the VTune library")

  include_directories(${VTUNE_INCLUDE_DIR})

  target_link_libraries(intersect_reconstruct_app ${VTUNE_LIB})
endif(ENABLE_PROFILE)
