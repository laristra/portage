#------------------------------------------------------------------------------#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#------------------------------------------------------------------------------#

# - Find flecsi
# Find the native FleCSI headers and libraries.
#
#  FleCSI_INCLUDE_DIRS - where to find flecsi.h, etc.
#  FleCSI_LIBRARIES    - List of libraries when using flecsi.
#  FleCSI_FOUND        - True if flecsi found.

# Look for the header file.
FIND_PATH(WONTON_INCLUDE_DIR NAMES wonton/support/wonton.h)

# Look for the library.
FIND_LIBRARY(WONTON_LIBRARY NAMES wonton libwonton)

# handle the QUIETLY and REQUIRED arguments and set FleCSI_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(WONTON WONTON_LIBRARY WONTON_INCLUDE_DIR)

# Copy the results to the output variables.
SET(WONTON_LIBRARIES ${WONTON_LIBRARY})
SET(WONTON_INCLUDE_DIRS ${WONTON_INCLUDE_DIR})

MARK_AS_ADVANCED(WONTON_INCLUDE_DIR WONTON_LIBRARY)
