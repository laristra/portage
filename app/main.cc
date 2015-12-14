#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>

#include "mpi.h"

#include "portage/support/portage.h"
#include "portage/driver/driver.h"
#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"
#include "portage/wrappers/state/jali/jali_state_wrapper.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliStateVector.h"
#include "JaliState.h"


int main(int argc, char** argv)
{
  // Get the example to run from command-line parameter
  int example = 0;
  if (argc <= 1)
  {
    std::printf("Usage: portageapp example-number\n");
    std::printf("example 0: 2d cell-centered remap\n");
    std::printf("example 1: 2d node-centered remap\n");
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

  std::printf("starting portageapp...\n");
  std::printf("running example %d\n", example);

  // Example 0 is a 2d cell-centered remap
  if (example == 0)
  {
    Jali::MeshFactory mf(MPI_COMM_WORLD);

    // Create a 2d quad input mesh from (0,0) to (1,1) with 3x3 zones
    Jali::Mesh* inputMesh = mf(0.0, 0.0, 1.0, 1.0, 3, 3);
    Portage::Jali_Mesh_Wrapper inputMeshWrapper(*inputMesh);

    // Create a 2d quad output mesh from (0,0) to (1,1) with 4x4 zones
    Jali::Mesh* targetMesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4);
    Portage::Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);

    Jali::State sourceState(inputMesh);
    std::vector<double> sourceData = {0.0,1.0,2.0,1.0,2.0,3.0,2.0,3.0,4.0};
    Jali::StateVector<double> & cellvecin = 
        sourceState.add("celldata", Jali::CELL, &(sourceData[0]));
    Portage::Jali_State_Wrapper sourceStateWrapper(sourceState);

    Jali::State targetState(targetMesh);
    std::vector<double> targetData(16,0.0);
    Jali::StateVector<double> & cellvecout = targetState.add("celldata", Jali::CELL, &(targetData[0]));
    Portage::Jali_State_Wrapper targetStateWrapper(targetState);

    Portage::Driver<Portage::Jali_Mesh_Wrapper> d(Portage::CELL, 
                                         inputMeshWrapper, sourceStateWrapper,
                                         targetMeshWrapper, targetStateWrapper);
    std::vector<std::string> remap_fields;
    remap_fields.push_back("celldata");
    d.set_remap_var_names(remap_fields);
    d.run();

    // When done, the "remapped_data" vector on the target mesh cells should
    // be identical to the "celldata" vector on the destination mesh

    std::cerr << "celldata vector on target mesh after remapping is:" << std::endl;
    std::cerr << cellvecout << std::endl;
  }

  // Example 1 is a 2d node-centered remap
  else if (example == 1)
  {
    Jali::MeshFactory mf(MPI_COMM_WORLD);

    // Create a 2d quad input mesh from (0,0) to (1,1) with 3x3 zones; 
    // The "true" arguments request that a dual mesh be constructed with wedges, corners, etc.
    Jali::Mesh* inputMesh = mf(0.0, 0.0, 1.0, 1.0, 3, 3, NULL, true, true, true, true);
    Portage::Jali_Mesh_Wrapper inputMeshWrapper(*inputMesh);

    // Create a 2d quad output mesh from (0,0) to (1,1) with 1x1 zones;
    // The "true" arguments request that a dual mesh be constructed with wedges, corners, etc.
    Jali::Mesh* targetMesh = mf(0.0, 0.0, 1.0, 1.0, 1, 1, NULL, true, true, true, true);
    Portage::Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);

    Jali::State sourceState(inputMesh);
    std::vector<double> sourceData = {1.5,1.5,1.5,1.5, 1.5,1.5,1.5,1.5, 1.5,1.5,1.5,1.5, 1.5,1.5,1.5,1.5};
    Jali::StateVector<double> & nodevecin = 
        sourceState.add("nodedata", Jali::NODE, &(sourceData[0]));
    Portage::Jali_State_Wrapper sourceStateWrapper(sourceState);

    Jali::State targetState(targetMesh);
    std::vector<double> targetData(4,0.0);
    Jali::StateVector<double> & nodevecout = targetState.add("nodedata", Jali::NODE, &(targetData[0]));
    Portage::Jali_State_Wrapper targetStateWrapper(targetState);

    // Since we are doing a node-centered remap, create the dual meshes
    Portage::MeshWrapperDual sourceDualMeshWrapper(inputMeshWrapper);
    Portage::MeshWrapperDual targetDualMeshWrapper(targetMeshWrapper);

    Portage::Driver<Portage::MeshWrapperDual> d(Portage::NODE, 
                                     sourceDualMeshWrapper, sourceStateWrapper,
                                     targetDualMeshWrapper, targetStateWrapper);
    std::vector<std::string> remap_fields;
    remap_fields.push_back("nodedata");
    d.set_remap_var_names(remap_fields);
    d.run();

    std::cerr << "nodedata vector on target mesh after remapping is:" << std::endl;
    std::cerr << nodevecout << std::endl;
  }

  std::printf("finishing portageapp...\n");

  MPI_Finalize();
}


