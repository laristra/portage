/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>

#include "mpi.h"

#ifdef ENABLE_PROFILE
#include "ittnotify.h"
#endif

#include "portage/distributed/mpi_collate.h"
#include "portage/driver/driver.h"
#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"
#include "portage/wrappers/state/jali/jali_state_wrapper.h"
#include "portage/intersect/intersect_r2d.h"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliStateVector.h"
#include "JaliState.h"

#define MSTK_HAVE_MPI
#include "Mesh_MSTK.hh"

using Portage::Jali_Mesh_Wrapper;
using Portage::Jali_State_Wrapper;
using Portage::collate;
using Portage::argsort;
using Portage::reorder;

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

// This is a 2-D test!  Find the example data in
// test_data/shotshell.exo, shotshell-v.exo.

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

  // CELL-centered or NODE-centered remapping?

  Jali::Entity_kind entityKind =
      (example == 0) || (example == 2) || (example == 3) ?
      Jali::Entity_kind::CELL : Jali::Entity_kind::NODE;



  const double TOL = 1e-4;


  // Initialize MPI
  int mpi_init_flag;
  MPI_Initialized(&mpi_init_flag);
  if (!mpi_init_flag)
    MPI_Init(&argc, &argv);
  int numpe, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &numpe);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::printf("starting shotshellapp...\n");
  std::printf("running example %d\n", example);

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::FACE,
                        Jali::Entity_kind::EDGE,
                        Jali::Entity_kind::WEDGE,
                        Jali::Entity_kind::CORNER});

  const std::shared_ptr<Jali::Mesh> sourceMesh = mf(argv[2]);
  const Jali_Mesh_Wrapper sourceMeshWrapper(*sourceMesh);
  const int inputDim = sourceMesh->space_dimension();

  const std::shared_ptr<Jali::Mesh> targetMesh = mf(argv[3]);
  const Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);
  const int targetDim = targetMesh->space_dimension();

  assert(inputDim == targetDim);

  std::cout << "Target mesh stats: " <<
      targetMeshWrapper.num_owned_cells() << " " <<
      targetMeshWrapper.num_owned_nodes() << std::endl;

  Jali::State sourceState(sourceMesh);
  Jali::State targetState(targetMesh);


  JaliGeometry::Point coord;

  if (entityKind == Jali::Entity_kind::CELL) {
    
    int ncall = sourceMeshWrapper.num_owned_cells() +
        sourceMeshWrapper.num_ghost_cells();

    std::vector<double> sourceData(ncall, 0);

#ifdef FIXED_SIZE_EXAMPLE
    for (int i = 0; i < 3034; i++) {
      coord = sourceMesh->cell_centroid(i);
      double x = coord[0];
      double y = coord[1];
      double z = (inputDim == 3) ? coord[2] : 0.0;
      sourceData[i] = std::sqrt(x*x+y*y+z*z);
    }
    for (int i = 3034; i < 4646; i++) {
      coord = sourceMesh->cell_centroid(i);
      double x = coord[0];
      double y = coord[1];
      double z = (inputDim == 3) ? coord[2] : 1.0;
      sourceData[i] = x*y*z;
    }
    for (int i = 4646; i < 5238; i++) {
      coord = sourceMesh->cell_centroid(i);
      double x = coord[0];
      double y = coord[1];
      double z = (inputDim == 3) ? coord[2] : 0.0;
      sourceData[i] = sin(x+y+z);
    }
#else
    for (int i = 0; i < ncall; i++) {
      coord = sourceMesh->cell_centroid(i);
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

    sourceState.add("celldata", sourceMesh, entityKind,
                    Jali::Entity_type::ALL, &(sourceData[0]));
    const Jali_State_Wrapper sourceStateWrapper(sourceState);

    const Jali::StateVector<double> & cellvecout =
        targetState.add("celldata", targetMesh, entityKind,
                        Jali::Entity_type::ALL, 0.0);
    Jali_State_Wrapper targetStateWrapper(targetState);

    std::vector<std::string> remap_fields;
    remap_fields.push_back("celldata");


    // Directly run cell-centered examples
    if ((example == 0) || (example == 2) || (example == 3)) {

      if (example == 0) {
        Portage::Driver<Portage::SearchKDTree,
                        Portage::IntersectR2D,
                        Portage::Interpolate_1stOrder,
                        2,
                        Portage::Jali_Mesh_Wrapper,
                        Portage::Jali_State_Wrapper>
            d(sourceMeshWrapper, sourceStateWrapper, targetMeshWrapper,
              targetStateWrapper);
        d.set_remap_var_names(remap_fields);
        d.run(false);
      }
      if ((example == 2) || (example == 3)) {

        Portage::Driver<Portage::SearchKDTree,
                        Portage::IntersectR2D,
                        Portage::Interpolate_2ndOrder,
                        2,
                        Portage::Jali_Mesh_Wrapper,
                        Portage::Jali_State_Wrapper>
            d(sourceMeshWrapper, sourceStateWrapper, targetMeshWrapper,
              targetStateWrapper);
        d.set_remap_var_names(remap_fields);
        d.run(false);
      }
    }

    std::cerr << "Last result: " << cellvecout[cellvecout.size()-1] <<
        std::endl;

    if ((example == 3) && unitTest) {
      double toterr = 0.0;
      double error;
      for (auto c = 0; c < cellvecout.size(); c++) {
        coord = targetMesh->cell_centroid(c);
        error = coord[0] + coord[1] - cellvecout[c];
        if (inputDim > 2) error += coord[2];
        toterr += error*error;
      }
      double L2 = sqrt(toterr/targetMeshWrapper.num_owned_cells());
      std::printf("\n\nL2 NORM OF ERROR = %lf\n\n", L2);
      assert(L2 < TOL);
    }

    if (dumpMesh) {
      // The `static_cast` is a workaround for an Intel compiler's header
      // files, which are missing a `std::to_string` function for ints.
      std::string example_num = std::to_string(static_cast<long long>(example));
      std::cerr << "Saving the source mesh" << std::endl;
      sourceState.export_to_mesh();
      dynamic_cast<Jali::Mesh_MSTK*>(sourceMesh.get())->
          write_to_exodus_file("input" + example_num + ".exo");
      
      std::cerr << "Saving the target mesh" << std::endl;
      targetState.export_to_mesh();
      dynamic_cast<Jali::Mesh_MSTK*>(targetMesh.get())->
          write_to_exodus_file("output" + example_num + ".exo");
      
      int field_len = targetMeshWrapper.num_owned_cells();

      // We concatenate the global IDs and field values on the rank 0 processor,
      // sort it by the global ID and save into a file.
      std::vector<int> lgid(field_len), gid;
      std::vector<double> lvalues(field_len), values;
      for (int i = 0; i < field_len; i++) {
        lgid[i] = targetMesh->GID(i, entityKind);
        lvalues[i] = cellvecout[i];
      }
      collate(MPI_COMM_WORLD, rank, numpe, lgid, gid);
      collate(MPI_COMM_WORLD, rank, numpe, lvalues, values);
      if (rank == 0) {
        std::vector<int> idx;
        argsort(gid, idx);
        reorder(gid, idx);
        reorder(values, idx);
        std::ofstream fout("field" + example_num + ".txt");
        fout << std::scientific;
        fout.precision(17);
        for (int i=0; i < gid.size(); i++) {
          fout << gid[i] << " " << values[i] << std::endl;
        }
      }
    }
    
  } else {

    int nnall = sourceMeshWrapper.num_owned_nodes() + 
        sourceMeshWrapper.num_ghost_nodes();
    std::vector<double> sourceData(nnall, 0);

    for (int i = 0; i < nnall; i++) {
      sourceMesh->node_get_coordinates(i, &coord);
      double x = coord[0];
      double y = coord[1];
      double z = (inputDim > 2) ? coord[2] : 0.0;
      if (example == 1)
        sourceData[i] = x+y+z;
      else
        sourceData[i] = x*x;
    }

    sourceState.add("nodedata", sourceMesh, entityKind,
                    Jali::Entity_type::ALL, &(sourceData[0]));
    const Jali_State_Wrapper sourceStateWrapper(sourceState);

    const Jali::StateVector<double> & nodevecout =
        targetState.add("nodedata", targetMesh, entityKind,
                        Jali::Entity_type::ALL, 0.0);
    Jali_State_Wrapper targetStateWrapper(targetState);

    std::vector<std::string> remap_fields;
    remap_fields.push_back("nodedata");

    // Create a dual mesh for node-centered examples
    if (example == 1) {

      Portage::Driver<Portage::SearchKDTree,
                      Portage::IntersectR2D,
                      Portage::Interpolate_2ndOrder,
                      2,
                      Portage::Jali_Mesh_Wrapper,
                      Portage::Jali_State_Wrapper>
          d(sourceMeshWrapper, sourceStateWrapper, targetMeshWrapper,
            targetStateWrapper);
      d.set_remap_var_names(remap_fields);
      d.run(false);
    }

    std::cerr << "Last result: " << nodevecout[nodevecout.size()-1] <<
        std::endl;

    if (dumpMesh) {
      // The `static_cast` is a workaround for an Intel compiler's header
      // files, which are missing a `std::to_string` function for ints.
      std::string example_num = std::to_string(static_cast<long long>(example));
      std::cerr << "Saving the source mesh" << std::endl;
      sourceState.export_to_mesh();
      dynamic_cast<Jali::Mesh_MSTK*>(sourceMesh.get())->
          write_to_exodus_file("input" + example_num + ".exo");
      
      std::cerr << "Saving the target mesh" << std::endl;
      targetState.export_to_mesh();
      dynamic_cast<Jali::Mesh_MSTK*>(targetMesh.get())->
          write_to_exodus_file("output" + example_num + ".exo");
      
      int field_len = targetMeshWrapper.num_owned_nodes();
      
      // We concatenate the global IDs and field values on the rank 0 processor,
      // sort it by the global ID and save into a file.
      std::vector<int> lgid(field_len), gid;
      std::vector<double> lvalues(field_len), values;
      for (int i = 0; i < field_len; i++) {
        lgid[i] = targetMesh->GID(i, entityKind);
        lvalues[i] = nodevecout[i];
      }
      collate(MPI_COMM_WORLD, rank, numpe, lgid, gid);
      collate(MPI_COMM_WORLD, rank, numpe, lvalues, values);
      if (rank == 0) {
        std::vector<int> idx;
        argsort(gid, idx);
        reorder(gid, idx);
        reorder(values, idx);
        std::ofstream fout("field" + example_num + ".txt");
        fout << std::scientific;
        fout.precision(17);
        for (int i=0; i < gid.size(); i++) {
          fout << gid[i] << " " << values[i] << std::endl;
        }
      }
    }

  }



  std::printf("finishing shotshellapp...\n");

  MPI_Finalize();

  return 0;
}
