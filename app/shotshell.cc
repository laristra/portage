/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>

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

int main(int argc, char** argv) {

  // Pause profiling until main loop
  #ifdef ENABLE_PROFILE
    __itt_pause();
  #endif

  // Get the example to run from command-line parameter
  int example = 0;
  if (argc <= 3) {
    std::printf("Usage: shotshellapp example-number input-mesh output-mesh\n");
    std::printf("example 0: 2d 1st order cell-centered remap\n");
    std::printf("example 1: 2d 1st order node-centered remap\n");
    std::printf("example 2: 2d 2nd order cell-centered remap\n");
    return 0;
  }
  if (argc > 1) example = atoi(argv[1]);

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

  std::printf("starting shotshellapp...\n");
  std::printf("running example %d\n", example);

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::FACE,
                        Jali::Entity_kind::EDGE,
                        Jali::Entity_kind::WEDGE,
                        Jali::Entity_kind::CORNER});

  const std::shared_ptr<Jali::Mesh> inputMesh = mf(argv[2]);
  const Jali_Mesh_Wrapper inputMeshWrapper(*inputMesh);

  const std::shared_ptr<Jali::Mesh> targetMesh = mf(argv[3]);
  const Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);

  std::cout << "Target mesh stats: " <<
      targetMeshWrapper.num_owned_cells() << " " <<
      targetMeshWrapper.num_owned_nodes() << std::endl;

  Jali::State sourceState(inputMesh);
  std::vector<double> sourceData(inputMeshWrapper.num_owned_cells(), 0);

  #ifdef FIXED_SIZE_EXAMPLE
    for (int i=0; i < 3034; i++) {
      auto coord = inputMeshWrapper.cellToXY(i);
      double x = coord[0][0];
      double y = coord[0][1];
      sourceData[i] = std::sqrt(x*x+y*y);
    }
    for (int i=3034; i < 4646; i++) {
      auto coord = inputMeshWrapper.cellToXY(i);
      double x = coord[0][0];
      double y = coord[0][1];
      sourceData[i] = x*y;
    }
    for (int i=4646; i < 5238; i++) {
      auto coord = inputMeshWrapper.cellToXY(i);
      double x = coord[0][0];
      double y = coord[0][1];
      sourceData[i] = sin(x+y);
    }
  #else
    for (int i=0; i < sourceData.size(); i++) {
      auto coord = inputMeshWrapper.cellToXY(i);
      double x = coord[0][0];
      double y = coord[0][1];
      sourceData[i] = x*x;
    }
  #endif

  Jali::Entity_kind entityKind =
      (example == 0) || (example == 2) ?
      Jali::Entity_kind::CELL : Jali::Entity_kind::NODE;
  sourceState.add("celldata", inputMesh, entityKind,
                  Jali::Parallel_type::ALL, &(sourceData[0]));
  const Jali_State_Wrapper sourceStateWrapper(sourceState);

  Jali::State targetState(targetMesh);
  std::vector<double> targetData(targetMeshWrapper.num_owned_cells(), 0);
  const Jali::StateVector<double> & cellvecout =
      targetState.add("celldata", targetMesh, entityKind,
                      Jali::Parallel_type::ALL, &(targetData[0]));
  Jali_State_Wrapper targetStateWrapper(targetState);

  std::vector<std::string> remap_fields;
  remap_fields.push_back("celldata");

  // Directly run cell-centered examples
  if ((example == 0) || (example == 2)) {
    Portage::Driver<Jali_Mesh_Wrapper,
                    Jali_State_Wrapper> d(Portage::CELL,
                                          inputMeshWrapper,
                                          sourceStateWrapper,
                                          targetMeshWrapper,
                                          targetStateWrapper);
    d.set_remap_var_names(remap_fields);

    if (example == 2) d.set_interpolation_order(2);

    d.run();
  }

  // Create a dual mesh for node-centered examples
  else if ((example == 1) || (example == 3)) {
    Portage::Driver<Portage::Jali_Mesh_Wrapper,
                    Portage::Jali_State_Wrapper> d(Portage::NODE,
                                                   inputMeshWrapper,
                                                   sourceStateWrapper,
                                                   targetMeshWrapper,
                                                   targetStateWrapper);
    d.set_remap_var_names(remap_fields);

    if (example == 3) d.set_interpolation_order(2);

    d.run();
  }

  std::cerr << "Sizes: " << sourceData.size() << " " << targetData.size() <<
      std::endl;
  std::cerr << "Last result: " << cellvecout[cellvecout.size()-1] << std::endl;

  #ifdef OUTPUT_RESULTS
    std::cerr << "Saving the source mesh" << std::endl;
    sourceState.export_to_mesh();
    dynamic_cast<Jali::Mesh_MSTK*>(inputMesh.get())->write_to_exodus_file("input.exo");

    std::cerr << "Saving the target mesh" << std::endl;
    targetState.export_to_mesh();
    dynamic_cast<Jali::Mesh_MSTK*>(targetMesh.get())->write_to_exodus_file("output.exo");
  #endif

  std::printf("finishing shotshellapp...\n");

  MPI_Finalize();
}
