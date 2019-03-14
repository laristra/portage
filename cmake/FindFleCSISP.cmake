#------------------------------------------------------------------------------#
# Copyright (c) 2016 Triad National Security, LLC
# All rights reserved.
#------------------------------------------------------------------------------#

# - Find flecsi
# Find the native FleCSI headers and libraries.
#
#  FleCSI_INCLUDE_DIRS - where to find flecsi.h, etc.
#  FleCSI_LIBRARIES    - List of libraries when using flecsi.
#  FleCSI_FOUND        - True if flecsi found.

# Look for the header file.
FIND_PATH(FleCSISP_INCLUDE_DIR NAMES flecsi-sp.h)

# Look for the library.
FIND_LIBRARY(FleCSISP_LIBRARY NAMES flecsi-sp libflecsi-sp)

# handle the QUIETLY and REQUIRED arguments and set FleCSI_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(FleCSISP FleCSISP_LIBRARY FleCSISP_INCLUDE_DIR)

# Copy the results to the output variables.
SET(FleCSISP_LIBRARIES ${FleCSISP_LIBRARY})
SET(FleCSISP_INCLUDE_DIRS ${FleCSISP_INCLUDE_DIR})

MARK_AS_ADVANCED(FleCSISP_INCLUDE_DIR FleCSISP_LIBRARY)
