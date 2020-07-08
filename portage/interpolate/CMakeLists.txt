#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#

project(portage_interpolate)

add_library(portage_interpolate INTERFACE)

target_include_directories(portage_interpolate INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIRECTORY}>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIRECTORY}>
  $<INSTALL_INTERFACE:include>
  )

# Add header files
set(portage_interpolate_HEADERS
    interpolate_1st_order.h
    interpolate_2nd_order.h
    interpolate_3rd_order.h
    interpolate_nth_order.h
    gradient.h
    quadfit.h
)

# Not yet allowed for INTERFACE libraries
# 
# set_target_properties(portage_interpolate PROPERTIES
#   PUBLIC_HEADER "${portage_interpolate_HEADERS}")
#
# Directly install files instead
#
install(FILES ${portage_interpolate_HEADERS} DESTINATION include/portage/interpolate)

target_link_libraries(portage_interpolate INTERFACE portage_support)

target_link_libraries(portage INTERFACE portage_interpolate)

install(TARGETS portage_interpolate
  EXPORT portage_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/portage/interpolate
  INCLUDES DESTINATION include/portage/interpolate
  )

      
# Unit tests

if (ENABLE_UNIT_TESTS)

    cinch_add_unit(test_gradient
      SOURCES  test/test_gradient.cc
      LIBRARIES portage_interpolate  
      POLICY SERIAL)

    cinch_add_unit(test_interpolate_first_order
      SOURCES test/test_interp_1st_order.cc
      LIBRARIES portage_interpolate  
      POLICY SERIAL)

    cinch_add_unit(test_interpolate_first_order_rz
      SOURCES test/test_interp_1st_order_rz.cc
      LIBRARIES portage_interpolate
      POLICY SERIAL)

    cinch_add_unit(test_interpolate_second_order
      SOURCES test/test_interp_2nd_order.cc
      LIBRARIES portage_interpolate  
      POLICY SERIAL)

    cinch_add_unit(test_interpolate_second_order_rz
      SOURCES test/test_interp_2nd_order_rz.cc
      LIBRARIES portage_interpolate
      POLICY SERIAL)

    cinch_add_unit(test_interpolate_third_order
      SOURCES test/test_interp_3rd_order.cc
      LIBRARIES portage_interpolate  
      POLICY SERIAL)

    cinch_add_unit(test_quadfit
      SOURCES  test/test_quadfit.cc
      LIBRARIES portage_interpolate  
      POLICY SERIAL)

    if (WONTON_ENABLE_Jali)
      cinch_add_unit(test_interpolate_first_order_gentype
        SOURCES test/test_interp_1st_order_gentype.cc
        LIBRARIES portage_interpolate  
        POLICY MPI
	THREADS 1)
    endif ()

endif(ENABLE_UNIT_TESTS)

