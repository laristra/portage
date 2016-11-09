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
#include <utility>
#include <cmath>

#include <mpi.h>

#ifdef ENABLE_PROFILE
#include "ittnotify.h"
#endif

#include "portage/support/portage.h"
#include "portage/support/Point.h"
#include "portage/driver/driver.h"
#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"
#include "portage/wrappers/state/jali/jali_state_wrapper.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliStateVector.h"
#include "JaliState.h"

using Portage::Jali_Mesh_Wrapper;
/*!
  @file main.cc
  @brief A simple application that drives our remap routines.

  This program is used to showcase our capabilities with various types
  of remap operations (e.g. interpolation order) on various types of
  meshes (2d or 3d; node-centered or zone-centered) of some simple
  linear or quadratic data.  For the cases of remapping linear data
  with a second-order interpolator, the L2 norm output at the end
  should be identically zero.
 */

//////////////////////////////////////////////////////////////////////
// Helper routines and data structures

struct example_properties {
  example_properties(const int dim, const int order, const bool cell_centered,
                     const bool linear, const bool conformal = true)
      : dim(dim), order(order), cell_centered(cell_centered), linear(linear),
        conformal(conformal) { }

  int dim;             // dimensionality of meshes in example
  int order;           // interpolation order in example
  bool cell_centered;  // is this example a cell-centered remap?
  bool linear;         // is this example a remap of linear data?
  bool conformal;      // are the two meshes boundary-conformal?
};

// Use this to add new problems.  If needed, we can extend the
// example_properties struct to contain more information.
std::vector<example_properties> setup_examples() {
  std::vector<example_properties> examples;

  // Cell-centered remaps:

  // 2d 1st order cell-centered remap of linear func
  examples.emplace_back(2, 1, true, true);

  // 2d 2nd order cell-centered remap of linear func
  examples.emplace_back(2, 2, true, true);

  // 2d 1st order cell-centered remap of linear func on non-conformal meshes
  examples.emplace_back(2, 1, true, true, false);

  // 2d 2nd order cell-centered remap of linear func on non-conformal meshes
  examples.emplace_back(2, 2, true, true, false);

  // 2d 1st order cell-centered remap of quadratic func
  examples.emplace_back(2, 1, true, false);

  // 2d 2nd order cell-centered remap of quadratic func
  examples.emplace_back(2, 2, true, false);

  // 2d 1st order cell-centered remap of quadratic func on non-conformal meshes
  examples.emplace_back(2, 1, true, false, false);

  // 2d 2nd order cell-centered remap of quadratice func on non-conformal meshes
  examples.emplace_back(2, 2, true, false, false);

  // 3d 1st order cell-centered remap of linear func
  examples.emplace_back(3, 1, true, false);

  // 3d 2nd order cell-centered remap of linear func
  examples.emplace_back(3, 2, true, false);

  // 3d 1st order cell-centered remap of linear func on non-conformal meshes
  examples.emplace_back(3, 1, true, false, false);

  // 3d 2nd order cell-centered remap of linear func on non-conformal meshes
  examples.emplace_back(3, 2, true, false, false);

  // Node-centered remaps:

  // 2d 1st order node-centered remap of quadratic func
  examples.emplace_back(2, 1, false, false);

  // 2d 2nd order node-centered remap of quadratic func
  examples.emplace_back(2, 2, false, false);

  // 2d 1st order node-centered remap of quadratic func on non-conformal meshes
  examples.emplace_back(2, 1, false, false, false);

  // 2d 2nd order node-centered remap of quadratic func on non-conformal meshes
  examples.emplace_back(2, 2, false, false, false);

  // 3d 1st order node-centered remap of quadratic func
  examples.emplace_back(3, 1, false, false);

  // 3d 2nd order node-centered remap of quadratic func
  examples.emplace_back(3, 2, false, false);

  // 3d 1st order node-centered remap of quadratic func on non-conformal meshes
  examples.emplace_back(3, 1, false, false, false);

  // 3d 2nd order node-centered remap of quadratic func on non-conformal meshes
  examples.emplace_back(3, 2, false, false, false);

  return examples;
}

// Dump the usage with example information based on the registered
// examples from setup_examples()
void print_usage() {
  auto examples = setup_examples();
  std::printf("Usage: portageapp example-number nsourcecells ntargetcells [y] [y]\n");
  std::printf("If 'y' specified, dump data to input.exo and output.exo\n");
  std::printf("If second 'y' specified, reverse source ranks in distributed examples\n");
  std::printf("List of example numbers:\n");
  int i = 0;
  bool separated = false;
  std::printf("CELL-CENTERED EXAMPLES\n");
  for (const auto &example : examples) {
    if (!separated && !example.cell_centered) {
      std::printf("\n NODE-CENTERED EXAMPLES:\n");
      separated = true;
    }
    std::printf("  %d: %dd %s order %s-centered remap of %s func %s\n",
                i, example.dim,
                (example.order == 1) ? "1st" : "2nd",
                example.cell_centered ? "cell" : "node",
                example.linear ? "linear" : "quadratic",
                !example.conformal ? "on non-conformal mesh" : "");
    ++i;
  }
}
//////////////////////////////////////////////////////////////////////


int main(int argc, char** argv) {
  // Pause profiling until main loop
#ifdef ENABLE_PROFILE
  __itt_pause();
#endif

  // Get the example to run from command-line parameter
  int example_num, n_source, n_target;
  // Should we dump data?
  bool dump_output, reverse_source_ranks;
  // Must specify both the example number and size
  if (argc <= 3) {
    print_usage();
    return 0;
  }
  example_num = atoi(argv[1]);
  n_source = atoi(argv[2]);
  n_target = atoi(argv[3]);
  dump_output = (argc == 5) ?
      ((std::string(argv[4]) == "y") ? true : false)
      : false;
  reverse_source_ranks = (argc >= 6) ? ((std::string(argv[5]) == "y") ? true : false) : false;

  // Initialize MPI
  int mpi_init_flag;
  MPI_Initialized(&mpi_init_flag);
  if (!mpi_init_flag)
    MPI_Init(&argc, &argv);
  int numpe, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &numpe);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0)
  {
    std::printf("starting portageapp...\n");
    std::printf("running example %d\n", example_num);
  }

  // Get the properties for the requested example problem
  example_properties example = setup_examples()[example_num];

  // The mesh factory and mesh pointers.
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  std::shared_ptr<Jali::Mesh> inputMesh;
  std::shared_ptr<Jali::Mesh> targetMesh;

  // Set up a local communicator so that we can define mesh partitions
  // explicitly on each rank without Jali distributing it for us
  MPI_Group world_group, local_group;
  MPI_Comm_group(MPI_COMM_WORLD, &world_group);
  int ranks[1];  ranks[0] = rank;
  MPI_Group_incl(world_group, 1, ranks, &local_group);
  MPI_Comm local_comm;
  MPI_Comm_create(MPI_COMM_WORLD, local_group, &local_comm);
  Jali::MeshFactory mf_local(local_comm);

  // Cell-centered remaps
  if (example.cell_centered) {
    // Construct the meshes
    if (example.dim == 2) {
      mf.included_entities({Jali::Entity_kind::FACE});
      // 2d quad input mesh from (0,0) to (1,1) with nxn zones
      inputMesh = mf(0.0, 0.0, 1.0, 1.0, n_source, n_source);
      if (example.conformal) {
        // 2d quad output mesh from (0,0) to (1,1) with (n+1)x(n+1) zones
        targetMesh = mf(0.0, 0.0, 1.0, 1.0, n_target, n_target);
      } else {
        // 2d quad output mesh from (0,0) to (1+1.5dx,1) with (n+1)x(n+1)
        // zones and dx equal to the inputMesh grid spacing
        double dx = 1.0/static_cast<double>(n_target);
        targetMesh = mf(0.0, 0.0, 1.0+1.5*dx, 1.0, n_target, n_target);
      }
    } else {  // 3d
      mf.included_entities({Jali::Entity_kind::FACE,
                            Jali::Entity_kind::EDGE,
                            Jali::Entity_kind::WEDGE});
      // generate the input and target meshes for the non-distributed case
      if (numpe == 1) {
        // 3d hex input mesh from (0,0,0) to (1,1,1) with n_source x n_source x n_source zones
        inputMesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, n_source, n_source, n_source);
        if (example.conformal) {
           // 3d hex output mesh from (0,0,0) to (1,1,1) with n_target x n_target x n_target zones
          targetMesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, n_target, n_target, n_target);
        } else {
          // 3d hex output mesh from (0,0,0) to (1+1.5dx,1+1.5dx,1+1.5dx) with
          // (n_target)x(n_target)x(n_target) zones and dx equal to the inputMesh grid spacing
          double dx = 1.0/static_cast<double>(n_target);
          targetMesh = mf(0.0, 0.0, 0.0, 1.0+1.5*dx, 1.0+1.5*dx, 1.0+1.5*dx,
                          n_target, n_target, n_target);
        }
      // generate the input and target meshes for the distributed case
      } else {

        int source_dim = cbrt(1.0f*numpe) + 0.01f;
        mf_local.included_entities({Jali::Entity_kind::FACE,
                                    Jali::Entity_kind::EDGE,
                                    Jali::Entity_kind::WEDGE});
#ifdef MANUAL_SOURCE_DECOMPOSITION
        // compute the local partition of the source mesh based on the rank;
        // n_source is the number of cells in each dimension in each partition;
        // the number of ranks must be a perfect cube (1, 8, 27, etc.)
        double source_step = 1.0f / source_dim;
        int rrank = reverse_source_ranks ? numpe - rank - 1 : rank;
        int source_x = rrank % source_dim;
        int source_y = (rrank / source_dim) % source_dim;
        int source_z = rrank / (source_dim*source_dim);

        inputMesh = mf_local(source_step*source_x, source_step*source_y, source_step*source_z,
                             source_step*(source_x+1), source_step*(source_y+1), source_step*(source_z+1),
                             n_source, n_source, n_source);

#else
        mf.boundary_ghosts_requested(false);
        mf.num_ghost_layers_distmesh(1);
        inputMesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, n_source*source_dim, n_source*source_dim, n_source*source_dim);
#endif
                             
        // compute the local partition of the target mesh based on the rank;
        // n_target is the number of cells in each dimension in each partition;
        // the number of ranks must be a perfect cube (1, 8, 27, etc.)
        int target_dim = cbrt(1.0f*numpe) + 0.01f;
        double target_step = 1.0f / target_dim;
        if (!example.conformal)
        {
          double dx = 1.0/static_cast<double>(n_target*target_dim);
          target_step = (1.0f + 1.5*dx) / target_dim;
        }
        int target_x = rank % target_dim;
        int target_y = (rank / target_dim) % target_dim;
        int target_z = rank / (target_dim*target_dim);
        targetMesh = mf_local(target_step*target_x, target_step*target_y, target_step*target_z,
                              target_step*(target_x+1), target_step*(target_y+1), target_step*(target_z+1),
                              n_target, n_target, n_target);
      }    
    }

    // Wrappers for interfacing with the underlying mesh data structures
    Portage::Jali_Mesh_Wrapper inputMeshWrapper(*inputMesh);
    Portage::Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);

    const int nsrccells = inputMeshWrapper.num_owned_cells() + inputMeshWrapper.num_ghost_cells();
    const int ntarcells = targetMeshWrapper.num_owned_cells();
    
    // Fill the source state data with the specified profile
    Jali::State sourceState(inputMesh);
    std::vector<double> sourceData(nsrccells);

    if (example.linear) {
      for (unsigned int c = 0; c < nsrccells; ++c) {
        JaliGeometry::Point cen = inputMesh->cell_centroid(c);
        sourceData[c] = cen[0]+cen[1];
        if (example.dim == 3)
          sourceData[c] += cen[2];
      }
    } else {  // quadratic function
      for (unsigned int c = 0; c < nsrccells; ++c) {
        JaliGeometry::Point cen = inputMesh->cell_centroid(c);
        sourceData[c] = cen[0]*cen[0]+cen[1]*cen[1];
        if (example.dim == 3)
          sourceData[c] += cen[2]*cen[2];
      }
    }

    sourceState.add("celldata", inputMesh, Jali::Entity_kind::CELL,
                    Jali::Entity_type::ALL, &(sourceData[0]));
    Portage::Jali_State_Wrapper sourceStateWrapper(sourceState);
    
    // Build the target state storage
    Jali::State targetState(targetMesh);
    std::vector<double> targetData(ntarcells, 0.0);
    auto& cellvecout = targetState.add("celldata", targetMesh,
                                       Jali::Entity_kind::CELL,
                                       Jali::Entity_type::ALL,
                                       &(targetData[0]));
    Portage::Jali_State_Wrapper targetStateWrapper(targetState);

    // Build the main driver data for this mesh type

    // Register the variable name and interpolation order with the driver
    std::vector<std::string> remap_fields;
    remap_fields.push_back("celldata");

    if(example.dim == 2 && example.order == 2){
      Portage::SearchKDTree<2, Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper> 
               search(inputMeshWrapper, targetMeshWrapper);
      Portage::IntersectR2D<Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper> 
               intersect(inputMeshWrapper, targetMeshWrapper);
      Portage::Interpolate_2ndOrder<Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper, 
          Portage::Jali_State_Wrapper, Portage::CELL, 2> 
          interpolate(inputMeshWrapper, targetMeshWrapper, sourceStateWrapper);
    
      Portage::Driver<Portage::SearchKDTree<2, Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper>, 
          Portage::IntersectR2D<Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper>, 
          Portage::Interpolate_2ndOrder<Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper, 
          Portage::Jali_State_Wrapper, Portage::CELL, 2>,
          Portage::Jali_Mesh_Wrapper, Portage::Jali_State_Wrapper,
          Portage::Jali_Mesh_Wrapper, Portage::Jali_State_Wrapper>  
          d(search, intersect, interpolate, inputMeshWrapper, sourceStateWrapper,targetMeshWrapper,
            targetStateWrapper);
      d.set_remap_var_names(remap_fields);    
      d.run();
    }

    if(example.dim == 2 && example.order == 1){
      Portage::SearchKDTree<2, Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper> 
               search(inputMeshWrapper, targetMeshWrapper);
      Portage::IntersectR2D<Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper> 
               intersect(inputMeshWrapper, targetMeshWrapper);
      Portage::Interpolate_1stOrder<Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper, 
          Portage::Jali_State_Wrapper, Portage::CELL, 2> 
          interpolate(inputMeshWrapper, targetMeshWrapper, sourceStateWrapper);
    
      Portage::Driver<Portage::SearchKDTree<2, Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper>, 
          Portage::IntersectR2D<Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper>, 
          Portage::Interpolate_1stOrder<Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper, 
          Portage::Jali_State_Wrapper, Portage::CELL, 2>,
          Portage::Jali_Mesh_Wrapper, Portage::Jali_State_Wrapper,
          Portage::Jali_Mesh_Wrapper, Portage::Jali_State_Wrapper>  
          d(search, intersect, interpolate, inputMeshWrapper, sourceStateWrapper,targetMeshWrapper,
            targetStateWrapper);
      d.set_remap_var_names(remap_fields);    
      d.run();
    }


    if(example.dim == 3 && example.order == 1){
      Portage::SearchKDTree<3, Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper> 
               search(inputMeshWrapper, targetMeshWrapper);
      Portage::IntersectR3D<Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper> 
               intersect(inputMeshWrapper, targetMeshWrapper);
      Portage::Interpolate_1stOrder<Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper, 
          Portage::Jali_State_Wrapper, Portage::CELL, 3> 
          interpolate(inputMeshWrapper, targetMeshWrapper, sourceStateWrapper);
    
      Portage::Driver<Portage::SearchKDTree<3, Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper>, 
          Portage::IntersectR3D<Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper>, 
          Portage::Interpolate_1stOrder<Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper, 
          Portage::Jali_State_Wrapper, Portage::CELL, 3>,
          Portage::Jali_Mesh_Wrapper, Portage::Jali_State_Wrapper,
          Portage::Jali_Mesh_Wrapper, Portage::Jali_State_Wrapper>  
          d(search, intersect, interpolate, inputMeshWrapper, sourceStateWrapper,targetMeshWrapper,
            targetStateWrapper);
      d.set_remap_var_names(remap_fields);    
      d.run();
    }

    if(example.dim == 3 && example.order == 2){
      Portage::SearchKDTree<3, Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper> 
               search(inputMeshWrapper, targetMeshWrapper);
      Portage::IntersectR3D<Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper> 
               intersect(inputMeshWrapper, targetMeshWrapper);
      Portage::Interpolate_2ndOrder<Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper, 
          Portage::Jali_State_Wrapper, Portage::CELL, 3> 
          interpolate(inputMeshWrapper, targetMeshWrapper, sourceStateWrapper);
    
      Portage::Driver<Portage::SearchKDTree<3, Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper>, 
          Portage::IntersectR3D<Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper>, 
          Portage::Interpolate_2ndOrder<Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper, 
          Portage::Jali_State_Wrapper, Portage::CELL, 3>,
          Portage::Jali_Mesh_Wrapper, Portage::Jali_State_Wrapper,
          Portage::Jali_Mesh_Wrapper, Portage::Jali_State_Wrapper>  
          d(search, intersect, interpolate, inputMeshWrapper, sourceStateWrapper,targetMeshWrapper,
            targetStateWrapper);
      d.set_remap_var_names(remap_fields);    
      d.run();
    }

    struct timeval begin, end, diff;
    gettimeofday(&begin, 0);

    // Dump some timing information
    if (numpe > 1) MPI_Barrier(MPI_COMM_WORLD);
    gettimeofday(&end, 0);
    timersub(&end, &begin, &diff);
    const float seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
    if (rank == 0) std::cout << "Time: " << seconds << std::endl;

    // Output results for small test cases
    double toterr = 0.0;

    if (n_target < 10)
      std::cout << "celldata vector on target mesh after remapping is:"
                << std::endl;

    for (int c = 0; c < ntarcells; ++c) {
      JaliGeometry::Point ccen = targetMesh->cell_centroid(c);

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

      if (n_target < 10) {
        // dump diagnostics for each cell
        if (example.dim == 2) {
          std::printf("Cell=% 4d Centroid = (% 5.3lf,% 5.3lf)", c,
                      ccen[0], ccen[1]);
        } else {
          std::printf("%d Cell=% 4d Centroid = (% 5.3lf,% 5.3lf,% 5.3lf)", rank, c,
                      ccen[0], ccen[1], ccen[2]);
        }
        std::printf("  Value = % 10.6lf  Err = % lf\n",
                    cellvecout[c], error);
      }
      toterr += error*error;
    }

    if (numpe == 1) {
      // total L2 norm
      std::printf("\n\nL2 NORM OF ERROR = %lf\n\n", sqrt(toterr));
    } else {
      std::cout << std::flush << std::endl;
      MPI_Barrier(MPI_COMM_WORLD);
      double globalerr;
      MPI_Reduce(&toterr, &globalerr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

      if (rank == 0) std::printf("\n\nL2 NORM OF ERROR = %lf\n\n", sqrt(globalerr));
    }

    // Dump output, if requested
    if (dump_output) {
      std::cout << "Dumping data to Exodus files..." << std::endl;
      sourceState.export_to_mesh();
      targetState.export_to_mesh();
      inputMesh->write_to_exodus_file("input.exo");
      targetMesh->write_to_exodus_file("output.exo");
      std::cout << "...done." << std::endl;
    }

  } 

else {  // node-centered remaps
    mf.included_entities({Jali::Entity_kind::FACE,
                          Jali::Entity_kind::EDGE,
                          Jali::Entity_kind::WEDGE,
                          Jali::Entity_kind::CORNER});
    // Create the meshes
    if (example.dim == 2) {
      // Create a 2d quad input mesh from (0,0) to (1,1) with nxn zones;
      // The "true" arguments request that a dual mesh be constructed with
      // wedges, corners, etc.
      inputMesh = mf(0.0, 0.0, 1.0, 1.0, n_source, n_source);
      if (example.conformal) {
        // Create a 2d quad output mesh from (0,0) to (1,1) with (n-2)x(n-2)
        // zones.  The "true" arguments request that a dual mesh be constructed
        // with wedges, corners, etc.
        targetMesh = mf(0.0, 0.0, 1.0, 1.0, n_target, n_target);
      } else {
        // 2d quad output mesh from (0,0) to (1+1.5dx,1) with (n-2)x(n-2)
        // zones and dx equal to the inputMesh grid spacing
        double dx = 1.0/static_cast<double>(n_target);
        targetMesh = mf(0.0, 0.0, 1.0+1.5*dx, 1.0, n_target, n_target);
      }
    } else {  // 3d
      // Create a 3d hex input mesh from (0,0,0) to (1,1,1) with nxnxn zones;
      // The "true" arguments request that a dual mesh be constructed with
      // wedges, corners, etc.
      inputMesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, n_source, n_source, n_source);
      if (example.conformal) {
        // Create a 3d hex output mesh from (0,0,0) to (1,1,1) with
        // (n-2)x(n-2)x(n-2) zones.  The "true" arguments request that a dual
        // mesh be constructed with wedges, corners, etc.
        targetMesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, n_target, n_target, n_target);
      } else {
        // 3d hex output mesh from (0,0,0) to (1+1.5dx,1+1.5dx,1+1.5dx) with
        // (n-2)x(n-2)x(n-2) zones and dx equal to the inputMesh grid spacing
        double dx = 1.0/static_cast<double>(n_target);
        targetMesh = mf(0.0, 0.0, 0.0, 1.0+1.5*dx, 1.0+1.5*dx, 1.0+1.5*dx,
                        n_target, n_target, n_target);
      }
    }

    // Wrappers for interfacing with the underlying mesh data structures.
    Portage::Jali_Mesh_Wrapper inputMeshWrapper(*inputMesh);
    Portage::Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);

    const int nsrcnodes = inputMeshWrapper.num_owned_nodes();
    const int ntarnodes = targetMeshWrapper.num_owned_nodes();

    // Fill the source state datat with the specified profile.
    Jali::State sourceState(inputMesh);
    std::vector<double> sourceData(nsrcnodes);

    /*!
      @todo make node_get_coordinates be consistent in data type with
      cell_centroid?
    */
    // populate input field with quadratic function
    if (example.dim == 2) {
      Portage::Point<2> nodexy;
      for (int i = 0; i < nsrcnodes; ++i) {
        inputMeshWrapper.node_get_coordinates(i, &nodexy);
        sourceData[i] = nodexy[0]*nodexy[0] + nodexy[1]*nodexy[1];
      }
    } else {  // 3d
      Portage::Point<3> nodexyz;
      for (int i = 0; i < nsrcnodes; ++i) {
        inputMeshWrapper.node_get_coordinates(i, &nodexyz);
        sourceData[i] = nodexyz[0]*nodexyz[0] + nodexyz[1]*nodexyz[1]
            + nodexyz[2]*nodexyz[2];
      }
    }

    sourceState.add("nodedata", inputMesh, Jali::Entity_kind::NODE,
                    Jali::Entity_type::ALL, &(sourceData[0]));
    Portage::Jali_State_Wrapper sourceStateWrapper(sourceState);

    // Build the target state storage
    Jali::State targetState(targetMesh);
    auto& nodevecout = targetState.add("nodedata", targetMesh,
                                       Jali::Entity_kind::NODE,
                                       Jali::Entity_type::ALL,
                                       0.0);
    Portage::Jali_State_Wrapper targetStateWrapper(targetState);


    // Register the variable name and remap order with the driver
    std::vector<std::string> remap_fields;
    remap_fields.push_back("nodedata");

    // Build the main driver data for this mesh type

    if(example.dim == 2 && example.order == 1){

      Portage::MeshWrapperDual<Portage::Jali_Mesh_Wrapper> sourceDualWrapper(inputMeshWrapper);
      Portage::MeshWrapperDual<Portage::Jali_Mesh_Wrapper> targetDualWrapper(targetMeshWrapper);

      //Create the search, intersect functors

      Portage::SearchKDTree<2, Portage::MeshWrapperDual<Jali_Mesh_Wrapper>, Portage::MeshWrapperDual<Jali_Mesh_Wrapper>> search(sourceDualWrapper, targetDualWrapper);
      Portage::IntersectR2D<Portage::MeshWrapperDual<Jali_Mesh_Wrapper>, Portage::MeshWrapperDual<Jali_Mesh_Wrapper>> intersect(sourceDualWrapper, targetDualWrapper);

      Portage::Interpolate_1stOrder<Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper, 
          Portage::Jali_State_Wrapper, Portage::NODE, 2> 
          interpolate(inputMeshWrapper, targetMeshWrapper, sourceStateWrapper);

      Portage::Driver<Portage::SearchKDTree<2, Portage::MeshWrapperDual<Jali_Mesh_Wrapper>, Portage::MeshWrapperDual<Jali_Mesh_Wrapper>>, 
          Portage::IntersectR2D<Portage::MeshWrapperDual<Jali_Mesh_Wrapper>, Portage::MeshWrapperDual<Jali_Mesh_Wrapper>>, 
          Portage::Interpolate_1stOrder<Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper, 
          Portage::Jali_State_Wrapper, Portage::NODE, 2>,
          Jali_Mesh_Wrapper,Portage::Jali_State_Wrapper,  
          Jali_Mesh_Wrapper,Portage::Jali_State_Wrapper> 
          d(search, intersect, interpolate, inputMeshWrapper, sourceStateWrapper,targetMeshWrapper,
            targetStateWrapper);
      d.set_remap_var_names(remap_fields);    
      d.run();
    }

    if(example.dim == 2 && example.order == 2){
      Portage::MeshWrapperDual<Portage::Jali_Mesh_Wrapper> sourceDualWrapper(inputMeshWrapper);
      Portage::MeshWrapperDual<Portage::Jali_Mesh_Wrapper> targetDualWrapper(targetMeshWrapper);

      //Create the search, intersect functors

      Portage::SearchKDTree<2, Portage::MeshWrapperDual<Jali_Mesh_Wrapper>, Portage::MeshWrapperDual<Jali_Mesh_Wrapper>> search(sourceDualWrapper, targetDualWrapper);
      Portage::IntersectR2D<Portage::MeshWrapperDual<Jali_Mesh_Wrapper>, Portage::MeshWrapperDual<Jali_Mesh_Wrapper>> intersect(sourceDualWrapper, targetDualWrapper);

      Portage::Interpolate_2ndOrder<Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper, 
          Portage::Jali_State_Wrapper, Portage::NODE, 2> 
          interpolate(inputMeshWrapper, targetMeshWrapper, sourceStateWrapper);

      Portage::Driver<Portage::SearchKDTree<2, Portage::MeshWrapperDual<Jali_Mesh_Wrapper>, Portage::MeshWrapperDual<Jali_Mesh_Wrapper>>, 
          Portage::IntersectR2D<Portage::MeshWrapperDual<Jali_Mesh_Wrapper>, Portage::MeshWrapperDual<Jali_Mesh_Wrapper>>, 
          Portage::Interpolate_2ndOrder<Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper, 
          Portage::Jali_State_Wrapper, Portage::NODE, 2>,
          Jali_Mesh_Wrapper,Portage::Jali_State_Wrapper,  
          Jali_Mesh_Wrapper,Portage::Jali_State_Wrapper> 
          d(search, intersect, interpolate, inputMeshWrapper, sourceStateWrapper,targetMeshWrapper,
            targetStateWrapper);
      d.set_remap_var_names(remap_fields);    
      d.run();
    }


    if(example.dim == 3 && example.order == 1){
      Portage::MeshWrapperDual<Portage::Jali_Mesh_Wrapper> sourceDualWrapper(inputMeshWrapper);
      Portage::MeshWrapperDual<Portage::Jali_Mesh_Wrapper> targetDualWrapper(targetMeshWrapper);

      //Create the search, intersect functors

      Portage::SearchKDTree<3, Portage::MeshWrapperDual<Jali_Mesh_Wrapper>, Portage::MeshWrapperDual<Jali_Mesh_Wrapper>> search(sourceDualWrapper, targetDualWrapper);
      Portage::IntersectR3D<Portage::MeshWrapperDual<Jali_Mesh_Wrapper>, Portage::MeshWrapperDual<Jali_Mesh_Wrapper>> intersect(sourceDualWrapper, targetDualWrapper);

      Portage::Interpolate_1stOrder<Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper, 
          Portage::Jali_State_Wrapper, Portage::NODE, 3> 
          interpolate(inputMeshWrapper, targetMeshWrapper, sourceStateWrapper);

      Portage::Driver<Portage::SearchKDTree<3, Portage::MeshWrapperDual<Jali_Mesh_Wrapper>, Portage::MeshWrapperDual<Jali_Mesh_Wrapper>>, 
          Portage::IntersectR3D<Portage::MeshWrapperDual<Jali_Mesh_Wrapper>, Portage::MeshWrapperDual<Jali_Mesh_Wrapper>>, 
          Portage::Interpolate_1stOrder<Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper, 
          Portage::Jali_State_Wrapper, Portage::NODE, 3>,
          Jali_Mesh_Wrapper,Portage::Jali_State_Wrapper,  
          Jali_Mesh_Wrapper,Portage::Jali_State_Wrapper> 
          d(search, intersect, interpolate, inputMeshWrapper, sourceStateWrapper,targetMeshWrapper,
            targetStateWrapper);
          d.set_remap_var_names(remap_fields);    
          d.run();
    }

    if(example.dim == 3 && example.order == 2){
      Portage::MeshWrapperDual<Portage::Jali_Mesh_Wrapper> sourceDualWrapper(inputMeshWrapper);
      Portage::MeshWrapperDual<Portage::Jali_Mesh_Wrapper> targetDualWrapper(targetMeshWrapper);

      //Create the search, intersect functors

      Portage::SearchKDTree<3, Portage::MeshWrapperDual<Jali_Mesh_Wrapper>, Portage::MeshWrapperDual<Jali_Mesh_Wrapper>> search(sourceDualWrapper, targetDualWrapper);
      Portage::IntersectR3D<Portage::MeshWrapperDual<Jali_Mesh_Wrapper>, Portage::MeshWrapperDual<Jali_Mesh_Wrapper>> intersect(sourceDualWrapper, targetDualWrapper);

      Portage::Interpolate_2ndOrder<Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper, 
          Portage::Jali_State_Wrapper, Portage::NODE, 3> 
          interpolate(inputMeshWrapper, targetMeshWrapper, sourceStateWrapper);

      Portage::Driver<Portage::SearchKDTree<3, Portage::MeshWrapperDual<Jali_Mesh_Wrapper>, Portage::MeshWrapperDual<Jali_Mesh_Wrapper>>, 
          Portage::IntersectR3D<Portage::MeshWrapperDual<Jali_Mesh_Wrapper>, Portage::MeshWrapperDual<Jali_Mesh_Wrapper>>, 
          Portage::Interpolate_2ndOrder<Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper, 
          Portage::Jali_State_Wrapper, Portage::NODE, 3>,
          Jali_Mesh_Wrapper,Portage::Jali_State_Wrapper,  
          Jali_Mesh_Wrapper,Portage::Jali_State_Wrapper> 
          d(search, intersect, interpolate, inputMeshWrapper, sourceStateWrapper,targetMeshWrapper,
            targetStateWrapper);
          d.set_remap_var_names(remap_fields);    
          d.run();
    }

  //FIXME: amh: timing issues
    struct timeval begin, end, diff;
    gettimeofday(&begin, 0);

    // Dump some timing information
    gettimeofday(&end, 0);
    timersub(&end, &begin, &diff);
    const float seconds = diff.tv_sec + 1.0E-6*diff.tv_usec;
    if (rank == 0) std::cout << "Time: " << seconds << std::endl;

    // Output results for small test cases
    if (n_target < 10) {
      double toterr = 0.0;
      double stdval, err;
      if (example.dim == 2) {
        Portage::Point<2> nodexy;
        for (int i = 0; i < ntarnodes; ++i) {
          targetMeshWrapper.node_get_coordinates(i, &nodexy);
          stdval = nodexy[0]*nodexy[0] + nodexy[1]*nodexy[1];
          err = fabs(stdval-nodevecout[i]);
          std::printf("Node=% 4d Coords = (% 5.3lf,% 5.3lf) ", i,
                      nodexy[0], nodexy[1]);
          std::printf("Value = %10.6lf Err = % lf\n", nodevecout[i], err);
          toterr += err;
        }
      } else {  // 3d
        Portage::Point<3> nodexyz;
        for (int i = 0; i < ntarnodes; ++i) {
          targetMeshWrapper.node_get_coordinates(i, &nodexyz);
          stdval = nodexyz[0]*nodexyz[0] + nodexyz[1]*nodexyz[1]
              + nodexyz[2]*nodexyz[2];

          err = fabs(stdval-nodevecout[i]);
          std::printf("Node=% 4d Coords = (% 5.3lf,% 5.3lf,% 5.3lf) ", i,
                      nodexyz[0], nodexyz[1], nodexyz[2]);
          std::printf("Value = %10.6lf Err = % lf\n", nodevecout[i], err);
          toterr += err;
        }
      }
      std::printf("\n\nL2 NORM OF ERROR = %lf\n\n", sqrt(toterr));
    }

    // Dump output, if requested
    if (dump_output) {
      std::cout << "Dumping data to Exodus files..." << std::endl;
      sourceState.export_to_mesh();
      targetState.export_to_mesh();
      inputMesh->write_to_exodus_file("input.exo");
      targetMesh->write_to_exodus_file("output.exo");
      std::cout << "...done." << std::endl;
    }

  }
  

  std::printf("finishing portageapp...\n");

  MPI_Finalize();
}
