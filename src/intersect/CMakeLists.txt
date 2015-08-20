#-----------------------------------------------------------------------------~#
# Copyright (c) 2015 Los Alamos National Security, LLC
# All rights reserved.
#-----------------------------------------------------------------------------~#

# Add header files
set(intersect_HEADERS
    clipper.hpp
    intersect.h
    intersectClipper.h
    PARENT_SCOPE
)

# Add source files
set(intersect_SOURCES
    clipper.cpp
    intersect.cc
    PARENT_SCOPE
)

# Add unit test
if(ENABLE_UNIT_TESTS)
  cinch_add_unit(intersectClipper
    SOURCES intersectClipperTest.cc clipper.cpp
    )
endif(ENABLE_UNIT_TESTS)

#----------------------------------------------------------------------------~-#
# Formatting options for emacs and vim.
# vim: set tabstop=4 shiftwidth=4 expandtab :
#----------------------------------------------------------------------------~-#
