#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#

project(portage_support)

add_library(portage_support INTERFACE)

target_include_directories(portage_support INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIRECTORY}>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIRECTORY}>
  $<INSTALL_INTERFACE:include>)

# Add header files
set(portage_support_HEADERS
    portage.h
    timer.h
    mpi_collate.h)

if (WONTON_ENABLE_KOKKOS)
  set(portage_support_HEADERS
      ${portage_support_HEADERS}
      weight.h
      basis.h
      operator.h
      operator_references.h
      faceted_setup.h)
endif ()

# Not yet allowed for INTERFACE libraries
# 
# set_target_properties(portage_support PROPERTIES
#   PUBLIC_HEADER "${portage_support_HEADERS}")
#
# Directly install files instead
#
install(FILES ${portage_support_HEADERS} DESTINATION include/portage/support)

target_link_libraries(portage_support INTERFACE ${WONTON_LIBRARIES})

if (WONTON_ENABLE_MPI)
  target_compile_definitions(portage_support INTERFACE WONTON_ENABLE_MPI)
endif ()

if (WONTON_ENABLE_THRUST)
  target_compile_definitions(portage_support INTERFACE WONTON_ENABLE_THRUST)
endif ()

if (PORTAGE_HAS_TANGRAM)
  # By making it a dependency of portage_support, we make sure that all other
  # modules (intersect, interpolate etc.) also see TANGRAM
  target_link_libraries(portage_support INTERFACE ${TANGRAM_LIBRARIES})
endif ()

target_link_libraries(portage INTERFACE portage_support)

install(TARGETS portage_support
  EXPORT portage_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/portage/support
  INCLUDES DESTINATION include/portage/support
  )


# Unit tests

if (WONTON_ENABLE_KOKKOS AND ENABLE_UNIT_TESTS)

  cinch_add_unit(test_weight
    SOURCES test/test_weight.cc
    LIBRARIES portage_support
    POLICY MPI
    THREADS 1
    )

  cinch_add_unit(test_basis
    SOURCES test/test_basis.cc
    LIBRARIES portage_support
    POLICY MPI
    THREADS 1
    )

  cinch_add_unit(test_operator
    SOURCES test/test_operator.cc
    LIBRARIES portage_support
    POLICY MPI
    THREADS 1
    )

  cinch_add_unit(test_faceted_setup
    SOURCES test/test_faceted_setup.cc
    LIBRARIES portage_support
    POLICY SERIAL
    POLICY MPI
    THREADS 1
    )

endif()
