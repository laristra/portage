#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]

project(portage_accumulate)

add_library(portage_accumulate INTERFACE)

target_include_directories(portage_accumulate INTERFACE
  $<BUILD_INTERFACE:${portage_SOURCE_DIR}>
  $<BUILD_INTERFACE:${portage_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>
  )

# Add header files
set (portage_accumulate_HEADERS
  accumulate.h
  )

# Not yet allowed for INTERFACE libraries
# 
# set_target_properties(portage_accumulate PROPERTIES
#   PUBLIC_HEADER "${portage_accumulate_HEADERS}")
#
# Directly install files instead
#
install(FILES ${portage_accumulate_HEADERS} DESTINATION include/portage/accumulate)

target_link_libraries(portage_accumulate INTERFACE portage_support)

target_link_libraries(portage INTERFACE portage_accumulate)

install(TARGETS portage_accumulate
  EXPORT portage_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/portage/accumulate
  INCLUDES DESTINATION include/portage/accumulate
  )




# Unit tests

if (ENABLE_UNIT_TESTS)

  portage_add_unittest(test_accumulate
    SOURCES test/test_accumulate.cc
    LIBRARIES portage_accumulate
    POLICY MPI
    THREADS 1
    )

endif (ENABLE_UNIT_TESTS)
