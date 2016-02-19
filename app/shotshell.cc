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

int main(int argc, char** argv)
{
  // Pause profiling until main loop
  #ifdef ENABLE_PROFILE
    __itt_pause();
  #endif

  // Get the example to run from command-line parameter
  int example = 0;
  if (argc <= 3)
  {
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

  const std::unique_ptr<Jali::Mesh> inputMesh = mf(argv[2],
      NULL, true, true, true, true);
  const Jali_Mesh_Wrapper inputMeshWrapper(*inputMesh);

  const std::unique_ptr<Jali::Mesh> targetMesh = mf(argv[3],
      NULL, true, true, true, true);
  const Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);

  std::cout << "Target mesh stats: " << targetMeshWrapper.num_owned_cells() << " " << targetMeshWrapper.num_owned_nodes() << std::endl;

  Jali::State sourceState(inputMesh.get());
  std::vector<double> sourceData(inputMeshWrapper.num_owned_cells(), 0);

  #ifdef FIXED_SIZE_EXAMPLE
    for (int i=0; i < 3034; i++) {
      auto coord = inputMeshWrapper.cellToXY(i);
      double x = coord[0].first;
      double y = coord[0].second;
      sourceData[i] = std::sqrt(x*x+y*y);
    }
    for (int i=3034; i < 4646; i++) {
      auto coord = inputMeshWrapper.cellToXY(i);
      double x = coord[0].first;
      double y = coord[0].second;
      sourceData[i] = x*y;
    }
    for (int i=4646; i < 5238; i++) {
      auto coord = inputMeshWrapper.cellToXY(i);
      double x = coord[0].first;
      double y = coord[0].second;
      sourceData[i] = sin(x+y);
    }
  #else
    for (int i=0; i < sourceData.size(); i++) {
      auto coord = inputMeshWrapper.cellToXY(i);
      double x = coord[0].first;
      double y = coord[0].second;
      sourceData[i] = x*x;
    }
  #endif

  sourceState.add("celldata", (example == 0) || (example == 2) ? Jali::CELL : Jali::NODE, &(sourceData[0]));
  const Jali_State_Wrapper sourceStateWrapper(sourceState);

  Jali::State targetState(targetMesh.get());
  std::vector<double> targetData(targetMeshWrapper.num_owned_cells(), 0);
  const Jali::StateVector<double> & cellvecout = targetState.add("celldata", (example == 0) || (example == 2) ? Jali::CELL : Jali::NODE, &(targetData[0]));
  Jali_State_Wrapper targetStateWrapper(targetState);

  std::vector<std::string> remap_fields;
  remap_fields.push_back("celldata");

  // Directly run cell-centered examples
  if ((example == 0) || (example == 2))
  {
    Portage::Driver<Jali_Mesh_Wrapper> d(Portage::CELL, inputMeshWrapper, sourceStateWrapper,
                                                        targetMeshWrapper, targetStateWrapper);
    d.set_remap_var_names(remap_fields);

    if (example == 2) d.set_remap_order(2);

    d.run();
  }

  // Create a dual mesh for node-centered examples
  else if ((example == 1) || (example == 3))
  {
    const Portage::MeshWrapperDual sourceDualMeshWrapper(inputMeshWrapper);
    const Portage::MeshWrapperDual targetDualMeshWrapper(targetMeshWrapper);

    Portage::Driver<Portage::MeshWrapperDual> d(Portage::NODE, sourceDualMeshWrapper, sourceStateWrapper,
                                                               targetDualMeshWrapper, targetStateWrapper);
    d.set_remap_var_names(remap_fields);

    if (example == 3) d.set_remap_order(2);

    d.run();
  }

  // When done, the "remapped_data" vector on the target mesh cells should
  // be identical to the "celldata" vector on the destination mesh
  std::cerr << "Sizes: " << sourceData.size() << " " << targetData.size() << std::endl;
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


