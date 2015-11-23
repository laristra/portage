#~----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

#------------------------------------------------------------------------------#
# Add a rule to build the executable
#------------------------------------------------------------------------------#

add_executable(portageapp main.cc)

#------------------------------------------------------------------------------#
# Add include and link dependencies
#------------------------------------------------------------------------------#

target_link_libraries(portageapp portage ${EXTRA_LIBS})

# MPI additions - need to push this down into cinch
target_link_libraries(portageapp ${MPI_CXX_LIBRARIES})
set_target_properties(portageapp PROPERTIES
                      LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

# Jali additions
include_directories(${Jali_INCLUDE_DIRS} ${Jali_TPL_INCLUDE_DIRS})
target_link_libraries(portageapp ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES})

