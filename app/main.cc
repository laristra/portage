/*---------------------------------------------------------------------------~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/
#include <sys/time.h>

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <memory>

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

/*!
  @file main.cc
  @brief A simple application that drives our remap routines.

  This program is used to showcase our capabilities with various types of
  remap operations (e.g. remap order) on various types of meshes (2d or 3d;
  node-centered or zone-centered) of some simple linear or quadratic data.
  For the cases of remapping linear data with a second-order remap function,
  the L2 norm output at the end should be identically zero.
 */

//////////////////////////////////////////////////////////////////////
// Helper routines and data structures

struct example_properties {
  example_properties(const int dim, const int order, const bool cell_centered,
                     const bool linear) : dim(dim), order(order),
                                          cell_centered(cell_centered),
                                          linear(linear) { }

  int dim;             // dimensionality of meshes in example
  int order;           // remap order in example
  bool cell_centered;  // is this example a cell-centered remap?
  bool linear;         // is this example a remap of linear data?
};

// Use this to add new problems.  If needed, we can extend the
// example_properties struct to contain more information.
std::vector<example_properties> setup_examples() {
  std::vector<example_properties> examples;

  // 0: 2d 1st order cell-centered remap of linear func
  examples.emplace_back(2, 1, true, true);

  // 1: 2d 1st order node-centered remap of linear func
  examples.emplace_back(2, 1, false, true);

  // 2: 2d 2nd order cell-centered remap of linear func
  examples.emplace_back(2, 2, true, true);

  // 3: 2d 1st order cell-centered remap of quadratic func
  examples.emplace_back(2, 1, true, false);

  // 4: 2d 2nd order cell-centered remap of quadratic func
  examples.emplace_back(2, 2, true, false);

  // 5: 2d 2nd order node-centered remap of linear func
  examples.emplace_back(2, 2, false, true);

  // 6: 3d 1st order cell-centered remap of linear func
  examples.emplace_back(3, 1, true, true);

  // 7: 3d 2nd order cell-centered remap of linear func
  examples.emplace_back(3, 2, true, true);

  return examples;
}

// Dump the usage with example information based on the registered
// examples from setup_examples()
void print_usage() {
  auto examples = setup_examples();
  std::printf("Usage: portageapp example-number ncells\n");
  std::printf("List of example numbers:\n");
  int i = 0;
  for (const auto &example : examples) {
    std::printf("  %d: %dd %s order %s-centered remap of %s func\n",
                i, example.dim,
                (example.order == 1) ? "1st" : "2nd",
                example.cell_centered ? "cell" : "node",
                example.linear ? "linear" : "quadratic");
    i++;
  }
}
//////////////////////////////////////////////////////////////////////


int main(int argc, char** argv) {
  // Pause profiling until main loop
#ifdef ENABLE_PROFILE
  __itt_pause();
#endif

  // Get the example to run from command-line parameter
  int example_num, n;
  // Must specify both the example number and size
  if (argc <= 2) {
    print_usage();
    return 0;
  }
  example_num = atoi(argv[1]);
  n = atoi(argv[2]);

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
  std::printf("running example %d\n", example_num);

  // Get the properties for the requested example problem
  example_properties example = setup_examples()[example_num];

  // The mesh factory and mesh pointers.
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  std::unique_ptr<Jali::Mesh> inputMesh(nullptr);
  std::unique_ptr<Jali::Mesh> targetMesh(nullptr);

  // Cell-centered remaps
  if (example.cell_centered) {
    // Construct the meshes
    if (example.dim == 2) {
      // 2d quad input mesh from (0,0) to (1,1) with nxn zones
      inputMesh = std::unique_ptr<Jali::Mesh>(mf(0.0, 0.0, 1.0, 1.0,
                                                 n, n));
      // 2d quad output mesh from (0,0) to (1,1) with (n+1)x(n+1) zones
      targetMesh = std::unique_ptr<Jali::Mesh>(mf(0.0, 0.0, 1.0, 1.0,
                                                  n+1, n+1));
    } else {  // 3d
      // 3d hex input mesh from (0,0,0) to (1,1,1) with nxnxn zones
      inputMesh = std::unique_ptr<Jali::Mesh>(mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                                 n, n, n,
                                                 NULL,
                                                 true, true, true, false));
      // 3d hex output mesh from (0,0,0) to (1,1,1) with (n+1)x(n+1)x(n+1) zones
      targetMesh = std::unique_ptr<Jali::Mesh>(mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                                  n+1, n+1, n+1,
                                                  NULL,
                                                  true, true, true, false));
    }

    // Wrappers for interfacing with the underlying mesh data structures
    Portage::Jali_Mesh_Wrapper inputMeshWrapper(*inputMesh);
    Portage::Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);

    const int nsrccells = inputMeshWrapper.num_owned_cells();
    const int ntarcells = targetMeshWrapper.num_owned_cells();

    // Fill the source state data with the specified profile
    Jali::State sourceState(inputMesh.get());
    std::vector<double> sourceData(nsrccells);

    std::vector<double> cen;
    if (example.linear) {
      for (unsigned int c = 0; c < nsrccells; c++) {
        inputMeshWrapper.cell_centroid(c, &cen);
        sourceData[c] = cen[0]+cen[1];
        if (example.dim == 3)
          sourceData[c] += cen[2];
      }
    } else {  // quadratic function
      for (unsigned int c = 0; c < nsrccells; c++) {
        inputMeshWrapper.cell_centroid(c, &cen);
        sourceData[c] = cen[0]*cen[0]+cen[1]*cen[1];
        if (example.dim == 3)
          sourceData[c] += cen[2]*cen[2];
      }
    }

    sourceState.add("celldata", Jali::CELL, &(sourceData[0]));
    Portage::Jali_State_Wrapper sourceStateWrapper(sourceState);

    // Build the target state storage
    Jali::State targetState(targetMesh.get());
    std::vector<double> targetData(ntarcells, 0.0);
    auto& cellvecout = targetState.add("celldata",
                                       Jali::CELL,
                                       &(targetData[0]));
    Portage::Jali_State_Wrapper targetStateWrapper(targetState);

    // Build the main driver data for this mesh type
    Portage::Driver<Portage::Jali_Mesh_Wrapper> d(Portage::CELL,
                                                  inputMeshWrapper,
                                                  sourceStateWrapper,
                                                  targetMeshWrapper,
                                                  targetStateWrapper);
    // Register the variable name and remap order with the driver
    std::vector<std::string> remap_fields;
    remap_fields.push_back("celldata");
    d.set_remap_var_names(remap_fields);

    d.set_remap_order(example.order);

    struct timeval begin, end, diff;
    gettimeofday(&begin, 0);

    // Do the remap
    d.run();

    // Dump some timing information
    gettimeofday(&end, 0);
    timersub(&end, &begin, &diff);
    const float seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
    std::cout << "Time: " << seconds << std::endl;

    // Output results for small test cases
    if (n < 10) {
      double toterr = 0.0;

      std::cout << "celldata vector on target mesh after remapping is:"
                << std::endl;

      for (int c = 0; c < ntarcells; c++) {
        std::vector<double> ccen;
        targetMeshWrapper.cell_centroid(c, &ccen);

        double error;
        if (example.linear) {
          error = ccen[0]+ccen[1] - cellvecout[c];
          if (example.dim == 3)
            error += ccen[2];
        } else {  // quadratic
          error = ccen[0]*ccen[0]+ccen[1]*ccen[1] - cellvecout[c];
          if (example.dim == 3)
            error += ccen[2]*ccen[2];
        }

        // dump diagnostics for each cell
        if (example.dim == 2) {
          std::printf("Cell=% 4d Centroid = (% 5.3lf,% 5.3lf)", c,
                      ccen[0], ccen[1]);
        } else {
          std::printf("Cell=% 4d Centroid = (% 5.3lf,% 5.3lf,% 5.3lf)", c,
                      ccen[0], ccen[1], ccen[2]);
        }
        std::printf("  Value = % 10.6lf  Err = % lf\n",
                    cellvecout[c], error);

        toterr += error*error;
      }
      // total L2 norm
      std::printf("\n\nL2 NORM OF ERROR = %lf\n\n", sqrt(toterr));
    }
  } else {  // node-centered remaps
    // Create a 2d quad input mesh from (0,0) to (1,1) with nxn zones;
    // The "true" arguments request that a dual mesh be constructed with
    // wedges, corners, etc.
    inputMesh = std::unique_ptr<Jali::Mesh>(mf(0.0, 0.0, 1.0, 1.0,
                                               n, n,
                                               NULL,
                                               true, true, true, true));
    Portage::Jali_Mesh_Wrapper inputMeshWrapper(*inputMesh);

    // Create a 2d quad output mesh from (0,0) to (1,1) with (n-2)x(n-2) zones;
    // The "true" arguments request that a dual mesh be constructed with
    // wedges, corners, etc.
    targetMesh = std::unique_ptr<Jali::Mesh>(mf(0.0, 0.0, 1.0, 1.0,
                                                n-2, n-2,
                                                NULL,
                                                true, true, true, true));
    Portage::Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);

    // Fill the source state datat with the specified profile
    Jali::State sourceState(inputMesh.get());
    std::vector<double> sourceData((n+1)*(n+1), 1.5);

    sourceState.add("nodedata", Jali::NODE, &(sourceData[0]));
    Portage::Jali_State_Wrapper sourceStateWrapper(sourceState);

    // Build the target state storage
    Jali::State targetState(targetMesh.get());
    std::vector<double> targetData((n+1)*(n+1), 0.0);
    Jali::StateVector<double> & nodevecout = targetState.add("nodedata",
                                                             Jali::NODE,
                                                             &(targetData[0]));
    Portage::Jali_State_Wrapper targetStateWrapper(targetState);

    // Since we are doing a node-centered remap, create the dual meshes
    Portage::MeshWrapperDual sourceDualMeshWrapper(inputMeshWrapper);
    Portage::MeshWrapperDual targetDualMeshWrapper(targetMeshWrapper);

    // Build the main driver data for this mesh type
    Portage::Driver<Portage::MeshWrapperDual> d(Portage::NODE,
                                                sourceDualMeshWrapper,
                                                sourceStateWrapper,
                                                targetDualMeshWrapper,
                                                targetStateWrapper);

    // Register the variable name and remap order with the driver
    std::vector<std::string> remap_fields;
    remap_fields.push_back("nodedata");
    d.set_remap_var_names(remap_fields);
    d.set_remap_order(example.order);

    struct timeval begin, end, diff;
    gettimeofday(&begin, 0);

    // Do the remap
    d.run();

    // Dump some timing information
    gettimeofday(&end, 0);
    timersub(&end, &begin, &diff);
    const float seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
    std::cout << "Time: " << seconds << std::endl;

    // Output results for small test cases
    if (n < 10) {
      std::cout << "nodedata vector on target mesh after remapping is:"
                << std::endl;
      std::cout << nodevecout;
    }
  }

  std::printf("finishing portageapp...\n");

  MPI_Finalize();
}


