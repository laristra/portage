#[[
Copyright (c) 2016, Los Alamos National Security, LLC
All rights reserved.

Copyright 2016. Los Alamos National Security, LLC. This software was produced
under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National
Laboratory (LANL), which is operated by Los Alamos National Security, LLC for
the U.S. Department of Energy. The U.S. Government has rights to use,
reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS
NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
derivative works, such modified software should be clearly marked, so as not to
confuse it with the version available from LANL.

Additionally, redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. Neither the name of Los Alamos National Security, LLC, Los Alamos
   National Laboratory, LANL, the U.S. Government, nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL
SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
]]



if (Jali_DIR) 
    add_executable(portageapp_jali portageapp_jali.cc)
    target_link_libraries(portageapp_jali portage ${EXTRA_LIBS})

    # MPI additions - need to push this down into cinch
    if (ENABLE_MPI)
        target_link_libraries(portageapp_jali ${MPI_CXX_LIBRARIES})
        set_target_properties(portageapp_jali PROPERTIES
                              LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
    endif (ENABLE_MPI)

    # Jali additions
    include_directories(${Jali_INCLUDE_DIRS} ${Jali_TPL_INCLUDE_DIRS})
    target_link_libraries(portageapp_jali ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES})
endif (Jali_DIR)

# FleCSI example
if (FLECSI_INSTALL_DIR)
  add_executable(portageapp_flecsi portageapp_flecsi.cc)
  target_link_libraries(portageapp_flecsi portage ${EXTRA_LIBS})
  include_directories(${FLECSI_INCLUDE_DIRS})
  target_link_libraries(portageapp_flecsi ${FLECSI_LIBRARIES})
# these appear needed because of something we are pulling in within src
  if (ENABLE_MPI)
      target_link_libraries(portageapp_flecsi ${MPI_CXX_LIBRARIES})
      set_target_properties(portageapp_flecsi PROPERTIES
                            LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
  endif (ENABLE_MPI)
endif (FLECSI_INSTALL_DIR)

if (JALI_DIR AND FLECSI_INSTALL_DIR)
  add_executable(portageapp_jali2flecsi portageapp_jali2flecsi.cc)
  target_link_libraries(portageapp_jali2flecsi portageapp_jali2flecsi ${EXTRA_LIBS})
  include_directories(${FLECSI_INCLUDE_DIRS})
  target_link_libraries(portageapp_jali2flecsi ${FLECSI_LIBRARIES})
# these appear needed because of something we are pulling in within src
  if (ENABLE_MPI)
      target_link_libraries(portageapp_jali2flecsi ${MPI_CXX_LIBRARIES})
      set_target_properties(portageapp_jali2flecsi PROPERTIES
                            LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
  endif (ENABLE_MPI)
  target_link_libraries(portageapp_jali2flecsi ${Jali_LIBRARIES}
                        ${Jali_TPL_LIBRARIES})

endif (JALI_DIR AND FLECSI_INSTALL_DIR)

#TCMalloc option
set(ENABLE_TCMALLOC TRUE CACHE BOOL "Use TCMalloc")
if(ENABLE_TCMALLOC)
  set(TCMALLOC_LIB "${HOME}" CACHE STRING "Set the TCMalloc library")

  if (Jali_DIR)
    target_link_libraries(portageapp_jali ${TCMALLOC_LIB})
  endif (Jali_DIR)
  if (FLECSI_INSTALL_DIR)
    target_link_libraries(portageapp_flecsi ${TCMALLOC_LIB})
  endif (FLECSI_INSTALL_DIR)
  if (JALI_DIR AND FLECSI_INSTALL_DIR)
    target_link_libraries(portageapp_jali2flecsi ${TCMALLOC_LIB})
  endif (JALI_DIR AND FLECSI_INSTALL_DIR)
endif(ENABLE_TCMALLOC)

# Profiling option (to allow pausing and resuming profiling during execution)
set(ENABLE_PROFILE FALSE CACHE BOOL "Enable VTune Profiling")
if(ENABLE_PROFILE)
  set(VTUNE_INCLUDE_DIR "${HOME}" CACHE PATH "Set the locaiton of VTune include")
  set(VTUNE_LIB "${HOME}" CACHE STRING "Set the VTune library")

  include_directories(${VTUNE_INCLUDE_DIR})

  if (Jali_DIR)
    target_link_libraries(portageapp_jali ${VTUNE_LIB})
  endif (Jali_DIR)
  if (FLECSI_INSTALL_DIR)
    target_link_libraries(portageapp_flecsi ${VTUNE_LIB})
  endif (FLECSI_INSTALL_DIR)
  if (JALI_DIR AND FLECSI_INSTALL_DIR)
    target_link_libraries(portageapp_jali2flecsi ${VTUNE_LIB})
  endif (JALI_DIR AND FLECSI_INSTALL_DIR)
endif(ENABLE_PROFILE)


if (ENABLE_APP_TESTS)
  add_subdirectory(test)
endif (ENABLE_APP_TESTS)
