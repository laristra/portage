#-----------------------------------------------------------------------------~#
# Copyright (c) 2015 Los Alamos National Security, LLC
# All rights reserved.
#-----------------------------------------------------------------------------~#

# Add header files
set(intersect_HEADERS
    clipper.hpp
    intersectClipper.h
    PARENT_SCOPE
)

# Add source files
set(intersect_SOURCES
    clipper.cpp
    PARENT_SCOPE
)

# Add unit test
if(ENABLE_UNIT_TESTS)
  cinch_add_unit(intersectClipper
    SOURCES intersectClipperTest.cc clipper.cpp
    LIBRARIES ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
    POLICY MPI
    )
endif(ENABLE_UNIT_TESTS)

