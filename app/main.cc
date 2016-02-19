/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include <sys/time.h>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>

#include <mpi.h>

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


int main(int argc, char** argv) {

  // Pause profiling until main loop
  #ifdef ENABLE_PROFILE
    __itt_pause();
  #endif

  // Get the example to run from command-line parameter
  int example = 0;
  int n = 3;
  if (argc <= 2) {
    std::printf("Usage: portageapp example-number ncells\n");
    std::printf("example 0: 2d 1st order cell-centered remap of linear func\n");
    std::printf("example 1: 2d 2nd order cell-centered remap of linear func\n");
    std::printf("example 2: 2d 1st order cell-centered remap of quadratic func\n");
    std::printf("example 3: 2d 2nd order cell-centered remap of quadratic func\n");
    std::printf("example 4: 3d 1st order cell-centered remap of quadratic func\n");
    std::printf("example 5: 3d 2nd order cell-centered remap of quadratic func\n");

    std::printf("\n\n");
    std::printf("example 6: 2d 1st order node-centered remap of quadratic func\n");
    std::printf("example 7: 2d 2nd order node-centered remap of quadratic func\n");
    std::printf("example 8: 3d 1st order node-centered remap of quadratic func\n");
    std::printf("example 9: 3d 2nd order node-centered remap of quadratic func\n");
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

  // Example 0-5 are cell-centered remaps
  if (example <= 5) {
    Jali::MeshFactory mf(MPI_COMM_WORLD);

    std::unique_ptr<Jali::Mesh> inputMesh;
    std::unique_ptr<Jali::Mesh> targetMesh;

    if (example <= 3) {
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

    if (example < 2) {  // linear function in 2D
      for (unsigned int c = 0; c < nsrccells; c++) {
        std::vector<double> cen;
        inputMeshWrapper.cell_centroid(c, &cen);
        sourceData[c] = cen[0]+cen[1];
      }
    }
    else {  // quadratic function in 2d and 3d
      for (unsigned int c = 0; c < nsrccells; c++) {
        std::vector<double> cen;
        inputMeshWrapper.cell_centroid(c, &cen);
        sourceData[c] = cen[0]*cen[0]+cen[1]*cen[1];
        if (example > 3)
          sourceData[c] = cen[2]*cen[2];
      }
    }
    Jali::StateVector<double> & cellvecin =
        sourceState.add("celldata", Jali::CELL, &(sourceData[0]));
    Portage::Jali_State_Wrapper sourceStateWrapper(sourceState);

    Jali::State targetState(targetMesh.get());
    std::vector<double> targetData(ntarcells, 0.0);
    Jali::StateVector<double> & cellvecout =
        targetState.add("celldata", Jali::CELL, &(targetData[0]));
    Portage::Jali_State_Wrapper targetStateWrapper(targetState);

    Portage::Driver<Portage::Jali_Mesh_Wrapper> d(Portage::CELL,
                                                  inputMeshWrapper,
                                                  sourceStateWrapper,
                                                  targetMeshWrapper,
                                                  targetStateWrapper);
    std::vector<std::string> remap_fields;
    remap_fields.push_back("celldata");
    d.set_remap_var_names(remap_fields);

    // Examples 1, 2 and 5 are 2nd order accurate remaps
    if (example%2)
      d.set_remap_order(2);

    struct timeval begin, end, diff;
    gettimeofday(&begin, 0);

    d.run();

    gettimeofday(&end, 0);
    timersub(&end, &begin, &diff);
    float seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
    std::cout << "Time: " << seconds << std::endl;

    // Output results for small test cases
    if (n < 10) {
      double toterr = 0.0;

      std::cerr << "celldata vector on target mesh after remapping is:" <<
          std::endl;
      int ntargetcells = targetMeshWrapper.num_owned_cells();
      for (int c = 0; c < ntargetcells; c++) {
        std::vector<double> ccen;
        targetMeshWrapper.cell_centroid(c, &ccen);

        double error;
        if (example == 0 || example == 1)
          error = ccen[0]+ccen[1] - cellvecout[c];
        else if (example == 2 || example == 3)
          error = ccen[0]*ccen[0]+ccen[1]*ccen[1] - cellvecout[c];
        else
          error = ccen[0]*ccen[0]+ccen[1]*ccen[1]+ccen[2]*ccen[2] -
              cellvecout[c];

        if (example <= 3) {
          std::printf("Cell=% 4d Centroid = (% 5.3lf,% 5.3lf)", c,
                      ccen[0], ccen[1]);
        }
        else {
          std::printf("Cell=% 4d Centroid = (% 5.3lf,% 5.3lf,% 5.3lf)", c,
                      ccen[0], ccen[1], ccen[2]);
        }
        std::printf("  Value = % 10.6lf  Err = % lf\n",
                    cellvecout[c], error);

        toterr += error*error;
      }
      std::printf("\n\nL2 NORM OF ERROR = %lf\n\n", sqrt(toterr));
    }
  }
  // Examples 6-9 are node-centered remaps in 2d and 3d
  else
  {
    Jali::MeshFactory mf(MPI_COMM_WORLD);

    std::unique_ptr<Jali::Mesh> inputMesh;
    std::unique_ptr<Jali::Mesh> targetMesh;

    // Create a 2d quad input mesh from (0,0) to (1,1) with nxn zones or
    // a 3d hex input mesh from (0,0,0) to (1,1,1) with nxnxn zones.
    // The "true" arguments request that a dual mesh be constructed with
    // edges, faces, wedges, corners
    
    if (example == 6 || example == 7)
      inputMesh = mf(0.0, 0.0, 1.0, 1.0, n, n, NULL,
                     true, true, true, true);
    else
      inputMesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, n, n, n, NULL,
                     true, true, true);

    // Create a 2d quad output mesh from (0,0) to (1,1) with
    // (n-2)x(n-2) zones or a 3d hex input mesh from (0,0,0) to
    // (1,1,1) with (n-2)x(n-2)x(n-2) zones. The "true" arguments
    // request that a dual mesh be constructed with edges, faces,
    // wedges and corners

    if (example == 6 || example == 7)
      targetMesh = mf(0.0, 0.0, 1.0, 1.0, n-2, n-2, NULL,
                      true, true, true, true);
    else
      targetMesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, n-2, n-2, n-2, NULL,
                      true, true, true, true);

    Portage::Jali_Mesh_Wrapper inputMeshWrapper(*inputMesh);
    Portage::Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);

    int nn = inputMeshWrapper.num_owned_nodes();
    int dim = inputMeshWrapper.space_dimension();

    Jali::State sourceState(inputMesh.get());
    std::vector<double> sourceData(nn);

    // populate input field with quadratic function

    if (dim == 2) {
      for (int i = 0; i < nn; ++i) {
        std::pair<double,double> nodexy;
        inputMeshWrapper.node_get_coordinates(i, &nodexy);
        sourceData[i] = nodexy.first*nodexy.first + nodexy.second*nodexy.second;
      }
    }
    else if (dim == 3) {
      for (int i = 0; i < nn; ++i) {
        std::tuple<double,double,double> nodexyz;
        inputMeshWrapper.node_get_coordinates(i, &nodexyz);
        sourceData[i] = std::get<0>(nodexyz)*std::get<0>(nodexyz) +
            std::get<1>(nodexyz)*std::get<1>(nodexyz) +
            std::get<2>(nodexyz)*std::get<2>(nodexyz);
      }
    }

 
    Jali::StateVector<double> & nodevecin =
        sourceState.add("nodedata", Jali::NODE, &(sourceData[0]));
    Portage::Jali_State_Wrapper sourceStateWrapper(sourceState);

    Jali::State targetState(targetMesh.get());
    nn = targetMeshWrapper.num_owned_nodes();
    std::vector<double> targetData(nn, 0.0);
    Jali::StateVector<double> & nodevecout =
        targetState.add("nodedata", Jali::NODE, &(targetData[0]));
    Portage::Jali_State_Wrapper targetStateWrapper(targetState);

    Portage::Driver<Portage::Jali_Mesh_Wrapper> d(Portage::NODE,
                                                  inputMeshWrapper,
                                                  sourceStateWrapper,
                                                  targetMeshWrapper,
                                                  targetStateWrapper);
    std::vector<std::string> remap_fields;
    remap_fields.push_back("nodedata");
    d.set_remap_var_names(remap_fields);

    if (example%2)
      d.set_remap_order(2);

    struct timeval begin, end, diff;
    gettimeofday(&begin, 0);

    d.run();

    gettimeofday(&end, 0);
    timersub(&end, &begin, &diff);
    float seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
    std::cout << "Time: " << seconds << std::endl;

    // Output results for small test cases
    if (n < 10) {
      nn = targetMeshWrapper.num_owned_nodes();
      double toterr = 0.0;
      if (example == 6 || example == 7) {
        std::pair<double,double> nodexy;
        for (int i = 0; i < nn; i++) {
          targetMeshWrapper.node_get_coordinates(i, &nodexy);
          double stdval = nodexy.first*nodexy.first +
              nodexy.second*nodexy.second;
          double err = fabs(stdval-nodevecout[i]);
          std::printf("Node=% 4d Coords = (% 5.3lf,% 5.3lf) ", i,
                      nodexy.first, nodexy.second);
          std::printf("Value = %10.6lf Err = % lf\n", nodevecout[i], err);
          toterr += err;
        }
      }
      else {
        std::tuple<double,double,double> nodexyz;
        for (int i = 0; i < nn; i++) {
          targetMeshWrapper.node_get_coordinates(i, &nodexyz);
          double stdval = std::get<0>(nodexyz)*std::get<0>(nodexyz) +
              std::get<1>(nodexyz)*std::get<1>(nodexyz) +
              std::get<2>(nodexyz)*std::get<2>(nodexyz);
          
          double err = fabs(stdval-nodevecout[i]);
          std::printf("Node=% 4d Coords = (% 5.3lf,% 5.3lf,% 5.3lf) ", i,
                      std::get<0>(nodexyz), std::get<1>(nodexyz),
                      std::get<2>(nodexyz));
          std::printf("Value = %10.6lf Err = % lf\n", nodevecout[i], err);
          toterr += err;
        }
      }
      std::printf("\n\nL2 NORM OF ERROR = %lf\n\n", sqrt(toterr));
    }
  }

  std::printf("finishing portageapp...\n");

  MPI_Finalize();
}


