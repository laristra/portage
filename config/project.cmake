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



project(portage)

cinch_minimum_required(1.0)

cinch_add_application_directory(app)
cinch_add_library_target(portage src)


# While we wait for Portage to upgrade cinch lets add the ENABLE_MPI definition
# to the compile line if ENABLE_MPI is on in the CMAKE command

if (ENABLE_MPI)
  add_definitions(-DENABLE_MPI)
endif (ENABLE_MPI)


# Create some subcomponent libraries for the purposes of unit testing
#
# NOT ABLE TO DO EITHER OF THESE

if (ENABLE_UNIT_TESTS)
#  cinch_add_library_target(search src/search)
#  cinch_add_library_target(intersect src/intersect)
#  cinch_add_library_target(remap src/remap)
#  cinch_add_library_target(state src/state)

#  cinch_add_library_target(portage_search src/search)
#  cinch_add_library_target(portage_intersect src/intersect)
#  cinch_add_library_target(portage_remap src/remap)
#  cinch_add_library_target(portage_state src/state)
endif()

# Add application tests
# May pull this logic into cinch at some future point
option(ENABLE_APP_TESTS "Enable testing of full app" OFF)
if(ENABLE_APP_TESTS)
  enable_testing()
endif()

#------------------------------------------------------------------------------#
# 
#------------------------------------------------------------------------------#

set(CINCH_HEADER_SUFFIXES "\\.h")

