#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <sys/time.h>

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
  int n = 3;
  if (argc <= 2)
  {
    std::printf("Usage: portageapp example-number ncells\n");
    std::printf("example 0: 2d cell-centered remap\n");
    std::printf("example 1: 2d node-centered remap\n");
    return 0;
  }
  if (argc > 1) example = atoi(argv[1]);
  if (argc > 2) n = atoi(argv[2]);

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

    // Create a 2d quad input mesh from (0,0) to (1,1) with nxn zones
    Jali::Mesh* inputMesh = mf(0.0, 0.0, 1.0, 1.0, n, n);
    Portage::Jali_Mesh_Wrapper inputMeshWrapper(*inputMesh);

    // Create a 2d quad output mesh from (0,0) to (1,1) with (n+1)x(n+1) zones
    Jali::Mesh* targetMesh = mf(0.0, 0.0, 1.0, 1.0, n+1, n+1);
    Portage::Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);

    Jali::State sourceState(inputMesh);
    std::vector<double> sourceData(n*n);
    for (unsigned int i=0; i<n; i++) 
      for (unsigned int j=0; j<n; j++)
        sourceData[i*n+j] = 1.0f*i+j;  //{0.0,1.0,2.0,1.0,2.0,3.0,2.0,3.0,4.0};
    Jali::StateVector<double> & cellvecin = sourceState.add("celldata", Jali::CELL, &(sourceData[0]));
    Portage::Jali_State_Wrapper sourceStateWrapper(sourceState);

    Jali::State targetState(targetMesh);
    std::vector<double> targetData((n+1)*(n+1), 0.0);
    Jali::StateVector<double> & cellvecout = targetState.add("celldata", Jali::CELL, &(targetData[0]));
    Portage::Jali_State_Wrapper targetStateWrapper(targetState);

    Portage::Driver<Portage::Jali_Mesh_Wrapper> d(Portage::CELL, 
                                         inputMeshWrapper, sourceStateWrapper,
                                         targetMeshWrapper, targetStateWrapper);
    std::vector<std::string> remap_fields;
    remap_fields.push_back("celldata");
    d.set_remap_var_names(remap_fields);

    struct timeval begin, end, diff;
    gettimeofday(&begin, 0);

    d.run();

    gettimeofday(&end, 0);
    timersub(&end, &begin, &diff);
    float seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
    std::cout << "Time: " << seconds << std::endl; 

    // Output results for small test cases
    if (n < 10)
    {
      std::cerr << "celldata vector on target mesh after remapping is:" << std::endl;
      std::cerr << cellvecout;
    }
  }

  // Example 1 is a 2d node-centered remap
  else if (example == 1)
  {
    Jali::MeshFactory mf(MPI_COMM_WORLD);

    // Create a 2d quad input mesh from (0,0) to (1,1) with 3x3 zones; 
    // The "true" arguments request that a dual mesh be constructed with wedges, corners, etc.
    Jali::Mesh* inputMesh = mf(0.0, 0.0, 1.0, 1.0, n, n, NULL, true, true, true, true);
    Portage::Jali_Mesh_Wrapper inputMeshWrapper(*inputMesh);

    // Create a 2d quad output mesh from (0,0) to (1,1) with 1x1 zones;
    // The "true" arguments request that a dual mesh be constructed with wedges, corners, etc.
    Jali::Mesh* targetMesh = mf(0.0, 0.0, 1.0, 1.0, n-2, n-2, NULL, true, true, true, true);
    Portage::Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);

    Jali::State sourceState(inputMesh);
    std::vector<double> sourceData((n+1)*(n+1), 1.5); 
    Jali::StateVector<double> & nodevecin = sourceState.add("nodedata", Jali::NODE, &(sourceData[0]));
    Portage::Jali_State_Wrapper sourceStateWrapper(sourceState);

    Jali::State targetState(targetMesh);
    std::vector<double> targetData((n+1)*(n+1), 0.0);
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

    struct timeval begin, end, diff;
    gettimeofday(&begin, 0);

    d.run();

    gettimeofday(&end, 0);
    timersub(&end, &begin, &diff);
    float seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
    std::cout << "Time: " << seconds << std::endl;

    // Output results for small test cases
    if (n < 10)
    {
      std::cerr << "nodedata vector on target mesh after remapping is:" << std::endl;
      std::cerr << nodevecout;
    }
  }

  std::printf("finishing portageapp...\n");

  MPI_Finalize();
}


