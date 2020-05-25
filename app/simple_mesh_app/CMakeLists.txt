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
target_link_libraries(simple_mesh_app PRIVATE portage)

if (ENABLE_APP_TIMINGS)
  target_compile_definitions(simple_mesh_app PRIVATE ENABLE_TIMINGS=1)
endif ()

if (ENABLE_PROFILE)
  target_link_libraries(simple_mesh_app PRIVATE ${VTUNE_LIB})
  target_include_directories(simple_mesh_app PRIVATE ${VTUNE_INCLUDE_DIR})
endif ()


# Simple_Mesh Multi Material test
if (PORTAGE_ENABLE_TANGRAM)
  if (TANGRAM_HAS_XMOF2D)
    add_executable(simple_mesh_mm_app simple_mesh_mm_app.cc)
    target_link_libraries(simple_mesh_mm_app PRIVATE portage)

    if (ENABLE_APP_TIMINGS)
      target_compile_definitions(simple_mesh_mm_app PRIVATE ENABLE_TIMINGS=1)
    endif ()
    
    if (ENABLE_PROFILE)
      target_link_libraries(simple_mesh_mm_app PRIVATE ${VTUNE_LIB})
      target_include_directories(simple_mesh_mm_app PRIVATE ${VTUNE_INCLUDE_DIR})
    endif ()
  endif ()
endif()

if (ENABLE_APP_TESTS)
  add_subdirectory(test)
endif (ENABLE_APP_TESTS)
