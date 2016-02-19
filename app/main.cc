#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <sys/time.h>

#include "mpi.h"

#ifdef ENABLE_PROFILE
#include "ittnotify.h"
#endif

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
  // Pause profiling until main loop
  #ifdef ENABLE_PROFILE
    __itt_pause();
  #endif

  // Get the example to run from command-line parameter
  int example = 0;
  int n = 3;
  if (argc <= 2)
  {
    std::printf("Usage: portageapp example-number ncells\n");
    std::printf("example 0: 2d 1st order cell-centered remap of linear func\n");
    std::printf("example 1: 2d 1st order node-centered remap of linear func\n");
    std::printf("example 2: 2d 2nd order cell-centered remap of linear func\n");
    std::printf("example 3: 2d 1st order cell-centered remap of quadratic func\n");
    std::printf("example 4: 2d 2nd order cell-centered remap of quadratic func\n");
    std::printf("example 5: 2d 2nd order node-centered remap of linear func\n");
    std::printf("example 6: 3d 1st order cell-centered remap of linear func\n");
    std::printf("example 7: 3d 2nd order cell-centered remap of linear func\n");
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

  // Example 0,2,3,4,6 are cell-centered remaps
  if ((example != 1) && (example != 5)) 
  {
    Jali::MeshFactory mf(MPI_COMM_WORLD);

    std::unique_ptr<Jali::Mesh> inputMesh;
    std::unique_ptr<Jali::Mesh> targetMesh;

    if (example < 6) {
      // 2d quad input mesh from (0,0) to (1,1) with nxn zones
      inputMesh = mf(0.0, 0.0, 1.0, 1.0, n, n);
      // 2d quad output mesh from (0,0) to (1,1) with (n+1)x(n+1) zones
      targetMesh = mf(0.0, 0.0, 1.0, 1.0, n+1, n+1);
    }
    else {
      // 3d hex input mesh from (0,0,0) to (1,1,1) with nxn zones
      inputMesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, n, n, n,
		     NULL, true, true, true, false);
      // 3d hex output mesh from (0,0,0) to (1,1,1) with (n+1)x(n+1)x(n+1) zones
      targetMesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, n+1, n+1, n+1,
		      NULL, true, true, true, false);
    }

    Portage::Jali_Mesh_Wrapper inputMeshWrapper(*inputMesh);
    Portage::Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);

    int nsrccells = inputMeshWrapper.num_owned_cells();
    int ntarcells = targetMeshWrapper.num_owned_cells();

    Jali::State sourceState(inputMesh.get());
    std::vector<double> sourceData(nsrccells);

    if (example == 3 || example == 4) { // quadratic function      
      for (unsigned int c = 0; c < nsrccells; c++) {
        std::vector<double> cen;
        inputMeshWrapper.cell_centroid(c,&cen);
        sourceData[c] = cen[0]*cen[0]+cen[1]*cen[1];
      }
    }
    else { // linear function
      for (unsigned int c = 0; c < nsrccells; c++) {
        std::vector<double> cen;
        inputMeshWrapper.cell_centroid(c,&cen);
        sourceData[c] = cen[0]+cen[1];
	if (example > 5) sourceData[c] += cen[2];
      }
    }
    Jali::StateVector<double> & cellvecin = sourceState.add("celldata", Jali::CELL, &(sourceData[0]));
    Portage::Jali_State_Wrapper sourceStateWrapper(sourceState);

    Jali::State targetState(targetMesh.get());
    std::vector<double> targetData(ntarcells, 0.0);
    Jali::StateVector<double> & cellvecout = targetState.add("celldata", Jali::CELL, &(targetData[0]));
    Portage::Jali_State_Wrapper targetStateWrapper(targetState);

    Portage::Driver<Portage::Jali_Mesh_Wrapper> d(Portage::CELL, 
                                         inputMeshWrapper, sourceStateWrapper,
                                         targetMeshWrapper, targetStateWrapper);
    std::vector<std::string> remap_fields;
    remap_fields.push_back("celldata");
    d.set_remap_var_names(remap_fields);

    // Examples 2, 4 and 7 are 2nd order accurate remaps

    if (example == 2 || example == 4 || example == 7)
      d.set_remap_order(2);

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
      double toterr = 0.0;

      std::cerr << "celldata vector on target mesh after remapping is:" << std::endl;
      int ntargetcells = targetMeshWrapper.num_owned_cells();
      for (int c = 0; c < ntargetcells; c++) {
        std::vector<double> ccen;
        targetMeshWrapper.cell_centroid(c,&ccen);

        double error;
        if (example == 0 || example == 2)
          error = ccen[0]+ccen[1] - cellvecout[c];
        else if (example == 3 || example == 4)
          error = ccen[0]*ccen[0]+ccen[1]*ccen[1] - cellvecout[c];
	else if (example == 6 || example == 7)
	  error = ccen[0]+ccen[1]+ccen[2] - cellvecout[c];

	if (example < 6) {
	  std::printf("Cell=% 4d Centroid = (% 5.3lf,% 5.3lf)",c,
		      ccen[0],ccen[1]);
	}
	else {
	  std::printf("Cell=% 4d Centroid = (% 5.3lf,% 5.3lf,% 5.3lf)",c,
		      ccen[0],ccen[1],ccen[2]);
	}
        std::printf("  Value = % 10.6lf  Err = % lf\n",
                    cellvecout[c],error);        

        toterr += error*error;
      }
      std::printf("\n\nL2 NORM OF ERROR = %lf\n\n",sqrt(toterr));
    }
  }
  // Examples 1 and 5 are 2d node-centered remaps
  else
  {
    Jali::MeshFactory mf(MPI_COMM_WORLD);

    // Create a 2d quad input mesh from (0,0) to (1,1) with 3x3 zones; 
    // The "true" arguments request that a dual mesh be constructed with wedges, corners, etc.
    std::unique_ptr<Jali::Mesh> inputMesh = mf(0.0, 0.0, 1.0, 1.0, n, n, NULL, true, true, true, true);
    Portage::Jali_Mesh_Wrapper inputMeshWrapper(*inputMesh);

    // Create a 2d quad output mesh from (0,0) to (1,1) with 1x1 zones;
    // The "true" arguments request that a dual mesh be constructed with wedges, corners, etc.
    std::unique_ptr<Jali::Mesh> targetMesh = mf(0.0, 0.0, 1.0, 1.0, n-2, n-2, NULL, true, true, true, true);
    Portage::Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);

    Jali::State sourceState(inputMesh.get());
    std::vector<double> sourceData((n+1)*(n+1), 1.5); 
    Jali::StateVector<double> & nodevecin = sourceState.add("nodedata", Jali::NODE, &(sourceData[0]));
    Portage::Jali_State_Wrapper sourceStateWrapper(sourceState);

    Jali::State targetState(targetMesh.get());
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

    if (example == 5)
      d.set_remap_order(2);

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


