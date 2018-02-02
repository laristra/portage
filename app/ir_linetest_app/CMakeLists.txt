#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]

#-------------------------------------------------------------------------------
# ir_linetest_app - uses 2D mesh SimpleMesh
#-------------------------------------------------------------------------------
add_executable(ir_linetest_app ir_linetest_app.cc)
target_link_libraries(ir_linetest_app
                      portage ${LAPACKE_LIBRARIES} ${XMOF2D_LIBRARIES} ${EXTRA_LIBS})
#TCMalloc option
set(ENABLE_TCMALLOC TRUE CACHE BOOL "Use TCMalloc")
if(ENABLE_TCMALLOC)
  set(TCMALLOC_LIB "${HOME}" CACHE STRING "Set the TCMalloc library")

  target_link_libraries(ir_linetest_app ${TCMALLOC_LIB})
endif(ENABLE_TCMALLOC)

# Profiling option (to allow pausing and resuming profiling during execution)
set(ENABLE_PROFILE FALSE CACHE BOOL "Enable VTune Profiling")
if(ENABLE_PROFILE)
  set(VTUNE_INCLUDE_DIR "${HOME}" CACHE PATH "Set the locaiton of VTune include")
  set(VTUNE_LIB "${HOME}" CACHE STRING "Set the VTune library")

  include_directories(${VTUNE_INCLUDE_DIR})

  target_link_libraries(ir_linetest_app ${VTUNE_LIB})
endif(ENABLE_PROFILE)
