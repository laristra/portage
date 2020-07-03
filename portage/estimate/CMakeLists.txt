#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]

project(portage_estimate)    

add_library(portage_estimate INTERFACE)

target_include_directories(portage_accumulate INTERFACE
  $<BUILD_INTERFACE:${portage_SOURCE_DIR}>
  $<BUILD_INTERFACE:${portage_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>
  )

# Add header files
set (estimate_HEADERS
  estimate.h
  )


# Not yet allowed for INTERFACE libraries
# 
# set_target_properties(portage_estimate PROPERTIES
#   PUBLIC_HEADER "${portage_estimate_HEADERS}")
#
# Directly install files instead
#
install(FILES ${portage_estimate_HEADERS} DESTINATION include/portage/estimate)

target_link_libraries(portage_estimate INTERFACE portage_support)

target_link_libraries(portage INTERFACE portage_estimate)

install(TARGETS portage_estimate
  EXPORT portage_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/portage/estimate
  INCLUDES DESTINATION include/portage/estimate
  )



# Unit tests

if (ENABLE_UNIT_TESTS)

  portage_add_unittest(test_estimate
    SOURCES test/test_estimate.cc
    LIBRARIES portage_estimate
    POLICY MPI
    THREADS 1
    )

endif (ENABLE_UNIT_TESTS)
