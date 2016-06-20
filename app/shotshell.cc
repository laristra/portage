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
    std::printf("Usage: shotshellapp example-number input-mesh output-mesh");
    std::printf(" [output?] [unit?]\n");
    std::printf("example 0: 1st order cell-centered remap\n");
    std::printf("example 1: 1st order node-centered remap\n");
    std::printf("example 2: 2nd order cell-centered remap\n");
    std::printf("example 3: 2nd order cell-centered remap of linear func\n");
    std::printf("   output?: to dump the meshes, set to 'y'\n");
    std::printf("   unit?:   for unit testing of example 3, set to 'y'\n");
    return 1;
}

int main(int argc, char** argv) {

  // Pause profiling until main loop
#ifdef ENABLE_PROFILE
  __itt_pause();
#endif

  // Get the example to run from command-line parameter
  if (argc < 4) return usage();

  int example = atoi(argv[1]);
  assert(example >= 0 && example < 4);

  const bool dumpMesh = (argc >= 5) ?
      ((std::string(argv[4]) == "y") ? true : false)
      : false;
  // unitTest will basically assert that the L2 norm for example 3 is roundoff
  const bool unitTest = (argc == 6) ?
      ((std::string(argv[5]) == "y") ? true : false)
      : false;

  const double TOL = 1e-4;

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
  const int inputDim = inputMesh->space_dimension();

  const std::shared_ptr<Jali::Mesh> targetMesh = mf(argv[3]);
  const Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);
  const int targetDim = targetMesh->space_dimension();

  assert(inputDim == targetDim);

  std::cout << "Target mesh stats: " <<
      targetMeshWrapper.num_owned_cells() << " " <<
      targetMeshWrapper.num_owned_nodes() << std::endl;

  Jali::State sourceState(inputMesh);
  std::vector<double> sourceData(inputMeshWrapper.num_owned_cells(), 0);

  std::vector<double> coord;

#ifdef FIXED_SIZE_EXAMPLE
  for (int i=0; i < 3034; i++) {
    inputMeshWrapper.cell_centroid(i, &coord);
    double x = coord[0];
    double y = coord[1];
    double z = (inputDim == 3) ? coord[2] : 0.0;
    sourceData[i] = std::sqrt(x*x+y*y+z*z);
  }
  for (int i=3034; i < 4646; i++) {
    inputMeshWrapper.cell_centroid(i, &coord);
    double x = coord[0];
    double y = coord[1];
    double z = (inputDim == 3) ? coord[2] : 1.0;
    sourceData[i] = x*y*z;
  }
  for (int i=4646; i < 5238; i++) {
    inputMeshWrapper.cell_centroid(i, &coord);
    double x = coord[0];
    double y = coord[1];
    double z = (inputDim == 3) ? coord[2] : 0.0;
    sourceData[i] = sin(x+y+z);
  }
#else
  for (int i=0; i < sourceData.size(); i++) {
    inputMeshWrapper.cell_centroid(i, &coord);
    double x = coord[0];
    double y = coord[1];
    double z = (inputDim > 2) ? coord[2] : 0.0;
    if (example == 3) {
      sourceData[i] = x+y+z;
    } else {
      sourceData[i] = x*x;
    }
  }
#endif

  Jali::Entity_kind entityKind =
      (example == 0) || (example == 2) || (example == 3) ?
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
  if ((example == 0) || (example == 2) || (example == 3)) {
    Portage::Driver<Jali_Mesh_Wrapper,
                    Jali_State_Wrapper> d(inputMeshWrapper,
                                          sourceStateWrapper,
                                          targetMeshWrapper,
                                          targetStateWrapper);
    d.set_remap_var_names(remap_fields);

    if ((example == 2) || (example == 3))
      d.set_interpolation_order(2);

    d.run();
  }

  // Create a dual mesh for node-centered examples
  else if (example == 1) {
    Portage::Driver<Portage::Jali_Mesh_Wrapper,
                    Portage::Jali_State_Wrapper> d(inputMeshWrapper,
                                                   sourceStateWrapper,
                                                   targetMeshWrapper,
                                                   targetStateWrapper);
    d.set_remap_var_names(remap_fields);

    d.run();
  }

  std::cerr << "Sizes: " << sourceData.size() << " " << targetData.size() <<
      std::endl;
  std::cerr << "Last result: " << cellvecout[cellvecout.size()-1] << std::endl;

  if ((example == 3) && unitTest) {
    double toterr = 0.0;
    double error;
    for (auto c = 0; c < targetData.size(); c++) {
      targetMeshWrapper.cell_centroid(c, &coord);
      error = coord[0] + coord[1] - cellvecout[c];
      if (inputDim > 2) error += coord[2];
      toterr += error*error;
    }
    double L2 = sqrt(toterr/targetMeshWrapper.num_owned_cells());
    std::printf("\n\nL2 NORM OF ERROR = %lf\n\n", L2);
    assert(L2 < TOL);
  }

  if (dumpMesh) {
    std::cerr << "Saving the source mesh" << std::endl;
    sourceState.export_to_mesh();
    dynamic_cast<Jali::Mesh_MSTK*>(inputMesh.get())->write_to_exodus_file("input.exo");

    std::cerr << "Saving the target mesh" << std::endl;
    targetState.export_to_mesh();
    dynamic_cast<Jali::Mesh_MSTK*>(targetMesh.get())->write_to_exodus_file("output.exo");
  }

  std::printf("finishing shotshellapp...\n");

  MPI_Finalize();

  return 0;
}
