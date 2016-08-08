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

#include <omp.h>

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

  # pragma omp parallel
  {
    std::printf(" threads: %d\n", omp_get_num_threads());
  }

  std::printf("starting timingapp...\n");

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::FACE,
                        Jali::Entity_kind::EDGE,
                        Jali::Entity_kind::WEDGE,
                        Jali::Entity_kind::CORNER});

  const std::shared_ptr<Jali::Mesh> inputMesh = mf(argv[1]);
  const Jali_Mesh_Wrapper inputMeshWrapper(*inputMesh);
  const int inputDim = inputMesh->space_dimension();

  const std::shared_ptr<Jali::Mesh> targetMesh = mf(argv[2]);
  const Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);
  const int targetDim = targetMesh->space_dimension();

  assert(inputDim == targetDim);

  std::cout << "Target mesh stats: " <<
      targetMeshWrapper.num_owned_cells() << " " <<
      targetMeshWrapper.num_owned_nodes() << std::endl;

  Jali::State sourceState(inputMesh);
  Jali_State_Wrapper sourceStateWrapper(sourceState);
  sourceStateWrapper.init_from_mesh();

  Jali::State targetState(targetMesh);
  Jali_State_Wrapper targetStateWrapper(targetState);
  targetStateWrapper.init_from_mesh();

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
  
  Portage::Driver<Jali_Mesh_Wrapper,
      Jali_State_Wrapper> d(inputMeshWrapper,
                            sourceStateWrapper,
                            targetMeshWrapper,
                            targetStateWrapper);
  d.set_remap_var_names(remap_fields);
  
  d.set_interpolation_order(order);
  
  d.run();
  
  if (dumpMesh) {
    std::cerr << "Saving the source mesh" << std::endl;
    sourceStateWrapper.export_to_mesh();
    dynamic_cast<Jali::Mesh_MSTK*>(inputMesh.get())->write_to_exodus_file("input.exo");

    std::cerr << "Saving the target mesh" << std::endl;
    targetStateWrapper.export_to_mesh();
    dynamic_cast<Jali::Mesh_MSTK*>(targetMesh.get())->write_to_exodus_file("output.exo");
  }

  std::printf("finishing timingapp...\n");

  MPI_Finalize();

  return 0;
}
