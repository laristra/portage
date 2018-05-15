#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]


#------------------------------------------------------------------------------#
# Add a rule to build the executable
#------------------------------------------------------------------------------#

# Simple_Mesh test
add_executable(simple_mesh_app simple_mesh_app.cc)
target_link_libraries(simple_mesh_app portage
  ${EXTRA_LIBS} ${LAPACKX_LIBRARIES} ${MPI_CXX_LIBRARIES})


# Simple_Mesh Multi Material test
if(TANGRAM_FOUND AND XMOF2D_FOUND)
	add_executable(simple_mesh_mm_app simple_mesh_mm_app.cc)
	target_link_libraries(simple_mesh_mm_app portage
  	${EXTRA_LIBS} ${LAPACKX_LIBRARIES} ${MPI_CXX_LIBRARIES})
endif(TANGRAM_FOUND AND XMOF2D_FOUND)


#TCMalloc option
set(ENABLE_TCMALLOC TRUE CACHE BOOL "Use TCMalloc")
if(ENABLE_TCMALLOC)
  set(TCMALLOC_LIB "${HOME}" CACHE STRING "Set the TCMalloc library")

  target_link_libraries(simple_mesh_app ${TCMALLOC_LIB})
endif(ENABLE_TCMALLOC)

# Profiling option (to allow pausing and resuming profiling during execution)
set(ENABLE_PROFILE FALSE CACHE BOOL "Enable VTune Profiling")
if(ENABLE_PROFILE)
  set(VTUNE_INCLUDE_DIR "${HOME}" CACHE PATH "Set the locaiton of VTune include")
  set(VTUNE_LIB "${HOME}" CACHE STRING "Set the VTune library")

  include_directories(${VTUNE_INCLUDE_DIR})

  target_link_libraries(simple_mesh_app ${VTUNE_LIB})
  target_link_libraries(simple_mesh_mm_app ${VTUNE_LIB})
endif(ENABLE_PROFILE)


if (ENABLE_APP_TESTS)
  add_subdirectory(test)
endif (ENABLE_APP_TESTS)
