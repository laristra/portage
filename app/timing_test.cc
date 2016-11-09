/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

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
#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"
#include "portage/wrappers/state/jali/jali_state_wrapper.h"

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
  Portage::SearchKDTree<3, Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper> search(source_mesh_wrapper, target_mesh_wrapper);
  const Portage::IntersectR3D<Portage::Jali_Mesh_Wrapper> intersect{source_mesh_wrapper , target_mesh_wrapper};
  Jali::State source_state(smesh); //blank source state
  Jali::State target_state(tmesh);
  Portage::Jali_State_Wrapper sourceStateWrapper2(source_state);
  sourceStateWrapper2.init_from_mesh();
  Portage::Jali_State_Wrapper targetStateWrapper2(target_state);
  targetStateWrapper2.init_from_mesh();
  Portage::Interpolate_2ndOrder<Portage::Jali_Mesh_Wrapper,
                                Portage::Jali_Mesh_Wrapper,
                                Portage::Jali_State_Wrapper,
                                Portage::CELL, 3>
      interpolate(source_mesh_wrapper, target_mesh_wrapper, sourceStateWrapper2);
  //*********amh***********placeholder

  Portage::Driver<
    Portage::SearchKDTree<3, Portage::Jali_Mesh_Wrapper, 
                          Portage::Jali_Mesh_Wrapper>, 
    Portage::IntersectR3D<Portage::Jali_Mesh_Wrapper>, 
    Portage::Interpolate_2ndOrder<Portage::Jali_Mesh_Wrapper,
                                  Portage::Jali_Mesh_Wrapper,
                                  Portage::Jali_State_Wrapper,
                                  Portage::CELL, 3>, 
    Portage::Jali_Mesh_Wrapper, Portage::Jali_State_Wrapper>
      d(search, intersect, interpolate, 
        source_mesh_wrapper, sourceStateWrapper2,
        target_mesh_wrapper, targetStateWrapper2);
  d.set_remap_var_names(remap_fields);  
  d.run();
  
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
