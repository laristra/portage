/*
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
*/





#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <cmath>

#include "mpi.h"

#ifdef ENABLE_PROFILE
#include "ittnotify.h"
#endif

#include "portage/driver/driver.h"
#include "portage/wonton/mesh/jali/jali_mesh_wrapper.h"
#include "portage/wonton/state/jali/jali_state_wrapper.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliStateVector.h"
#include "JaliState.h"

#define MSTK_HAVE_MPI
#include "Mesh_MSTK.hh"

using Portage::Jali_Mesh_Wrapper;
using Portage::Jali_State_Wrapper;

int usage() {
    std::printf("Usage: timingapp input-mesh output-mesh [order=1]");
    std::printf(" [output?]\n");
    std::printf("   order: integer order of remap, defaults to 1st\n");
    std::printf("   output?: to dump the meshes, set to 'y'\n");
    return 1;
}

int main(int argc, char** argv) {
  // Pause profiling until main loop
#ifdef ENABLE_PROFILE
  __itt_pause();
#endif

  // Get the example to run from command-line parameter
  if (argc < 3) return usage();

  const int order = (argc >= 4) ? atoi(argv[3]) : 1;
  const bool dumpMesh = (argc >= 5) ?
      ((std::string(argv[4]) == "y") ? true : false)
      : false;

  // Initialize MPI
  int mpi_init_flag;
  MPI_Initialized(&mpi_init_flag);
  if (!mpi_init_flag)
    MPI_Init(&argc, &argv);
  int numpe;
  MPI_Comm_size(MPI_COMM_WORLD, &numpe);
  if (numpe > 1) {
      std::printf("error - only 1 mpi rank is allowed\n");
      std::exit(1);
  }

  std::printf("starting timingapp...\n");

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::FACE,
                        Jali::Entity_kind::EDGE,
                        Jali::Entity_kind::WEDGE,
                        Jali::Entity_kind::CORNER});

  // Repeat the field names so that each gets remapped twice - this is
  // to match the number of variables (8) being remapped in FLAG

  std::vector<std::string> remap_fields;
  remap_fields.push_back("zone_pres");
  remap_fields.push_back("zone_dens");
  remap_fields.push_back("zone_ener");
  remap_fields.push_back("zone_temp");

  remap_fields.push_back("zone_pres");
  remap_fields.push_back("zone_dens");
  remap_fields.push_back("zone_ener");
  remap_fields.push_back("zone_temp");

  //*********amh***********placeholder  this should be adapted so that any dimension of mesh & order of interpolation can be used
  // overlay a 2x2x2 target mesh on a 3x3x3 source mesh
  // each target mesh cell gives eight candidate source cells
  const std::shared_ptr<Jali::Mesh> smesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 3, 3, 3);
  const std::shared_ptr<Jali::Mesh> tmesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2);
  const Portage::Jali_Mesh_Wrapper source_mesh_wrapper(*smesh);
  const Portage::Jali_Mesh_Wrapper target_mesh_wrapper(*tmesh);
  Jali::State source_state(smesh); //blank source state
  Jali::State target_state(tmesh);
  Portage::Jali_State_Wrapper sourceStateWrapper2(source_state);
  sourceStateWrapper2.init_from_mesh();
  Portage::Jali_State_Wrapper targetStateWrapper2(target_state);
  targetStateWrapper2.init_from_mesh();
  //*********amh***********placeholder

  Portage::Driver<
    Portage::SearchKDTree, 
    Portage::IntersectR3D, 
    Portage::Interpolate_2ndOrder,
    3, 
    Portage::Jali_Mesh_Wrapper, 
    Portage::Jali_State_Wrapper>
      d(source_mesh_wrapper, sourceStateWrapper2,
        target_mesh_wrapper, targetStateWrapper2);
  d.set_remap_var_names(remap_fields);  
  d.run(false);
  
  if (dumpMesh) {
    std::cerr << "Saving the source mesh" << std::endl;
    sourceStateWrapper2.export_to_mesh();
    dynamic_cast<Jali::Mesh_MSTK*>(smesh.get())->write_to_exodus_file("input.exo");

    std::cerr << "Saving the target mesh" << std::endl;
    targetStateWrapper2.export_to_mesh();
    dynamic_cast<Jali::Mesh_MSTK*>(tmesh.get())->write_to_exodus_file("output.exo");
  }

  std::printf("finishing timingapp...\n");

  MPI_Finalize();

  return 0;
}
