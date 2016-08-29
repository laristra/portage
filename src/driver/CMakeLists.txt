#-----------------------------------------------------------------------------~#
# Copyright (c) 2015 Los Alamos National Security, LLC
# All rights reserved.
#-----------------------------------------------------------------------------~#

# Add header files
set(driver_HEADERS
    driver.h
    PARENT_SCOPE
)

# Add source files
set(driver_SOURCES
    PARENT_SCOPE
)

# Unit tests

if (ENABLE_UNIT_TESTS)

   if (Jali_DIR)
       cinch_add_unit(test_driver_multivar
       SOURCES test/test_driver_multivar.cc
       LIBRARIES portage ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
       POLICY MPI
       THREADS 1)
   endif (Jali_DIR)

endif (ENABLE_UNIT_TESTS)
