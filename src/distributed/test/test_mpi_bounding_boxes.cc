/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include "portage/support/Point.h"
#include "portage/distributed/mpi_bounding_boxes.h"
#include "portage/wrappers/state/jali/jali_state_wrapper.h"
#include "portage/wrappers/state/flat/flat_state_wrapper.h"
#include "portage/wrappers/mesh/flat/flat_mesh_wrapper.h"

#include <iostream>

#include "gtest/gtest.h"
#include "mpi.h"

#include "Mesh.hh"
#include "MeshFactory.hh"

#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"


TEST(MPI_Bounding_Boxes, SimpleTest3D) {

  int commRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &commRank);

  // Add multiple state vector types
  double base1 = 10. + commRank * 2.;
  double dtest1[8] = {base1, base1 + 1., 0., 0., 0., 0., 0., 0.};
  double base2 = 100. + commRank * 2.;
  double dtest2[8] = {base2, base2 + 1., 0., 0., 0., 0., 0., 0.};

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  // Source mesh
  std::shared_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
        2, 2, 2);
  Portage::Jali_Mesh_Wrapper inputMeshWrapper(*source_mesh);
  Jali::State state(source_mesh);
  Portage::Jali_State_Wrapper wrapper(state);

  state.add("d1", source_mesh, Jali::Entity_kind::CELL,
            Jali::Entity_type::ALL, dtest1);
  state.add("d2", source_mesh, Jali::Entity_kind::CELL,
            Jali::Entity_type::ALL, dtest2);

  Portage::Flat_Mesh_Wrapper<> source_mesh_flat(8, inputMeshWrapper);
  Portage::Flat_State_Wrapper<> source_state_flat(wrapper, {"d1", "d2"});

  // Target mesh
  std::shared_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
        3, 3, 3);
  Portage::Jali_Mesh_Wrapper target_mesh_(*target_mesh);
  Jali::State target_state(target_mesh);
  Portage::Jali_State_Wrapper target_state_(target_state);

  // Use a bounding box distributor to send the source cells to the target
  // partitions where they are needed
  Portage::MPI_Bounding_Boxes distributor;
  distributor.distribute(source_mesh_flat, source_state_flat, target_mesh_,
                         target_state_);

  // Check number of cells and nodes received
  int exp_num_owned_cells = (commRank == 0 ? 2 :
                             commRank <= 2 ? 4 : 8);
  int num_owned_cells = source_mesh_flat.num_owned_cells();
  ASSERT_EQ(exp_num_owned_cells, num_owned_cells);
  int exp_num_owned_nodes = exp_num_owned_cells * 8;
  int num_owned_nodes = source_mesh_flat.num_owned_nodes();
  ASSERT_EQ(exp_num_owned_nodes, num_owned_nodes);
  int exp_num_cells = exp_num_owned_cells * 4;
  int num_cells = num_owned_cells + source_mesh_flat.num_ghost_cells();
  ASSERT_EQ(exp_num_cells, num_cells);
  int exp_num_nodes = exp_num_cells * 8;
  int num_nodes = num_owned_nodes + source_mesh_flat.num_ghost_nodes();
  ASSERT_EQ(exp_num_nodes, num_nodes);

  // Check coordinates
  // List coordinates of cell 0 - others are equal to this
  // with a shift
  std::vector<Portage::Point<3>> cell0Coords =
    {{0.0, 0.0, 0.0},  {0.5, 0.0, 0.0},  {0.5, 0.5, 0.0},  {0.0, 0.5, 0.0},
     {0.0, 0.0, 0.5},  {0.5, 0.0, 0.5},  {0.5, 0.5, 0.5},  {0.0, 0.5, 0.5}};
                                                              
  // List owned cells that should have been sent to each rank
  std::vector<int> expOwnedGids;
  switch (commRank) {
    case 0:
      expOwnedGids = {0, 1};
      break;
    case 1:
      expOwnedGids = {0, 1, 2, 3};
      break;
    case 2:
      expOwnedGids = {0, 1, 4, 5};
      break;
    case 3:
      expOwnedGids = {0, 1, 2, 3, 4, 5, 6, 7};
      break;
  }

  for (int c=0; c<num_owned_cells; ++c) {
    std::vector<Portage::Point<3>> myCoords;
    source_mesh_flat.cell_get_coordinates(c, &myCoords);
    ASSERT_EQ(8, myCoords.size());
    int gid = expOwnedGids[c];
    double dx = (gid & 4 ? 0.5 : 0.0);
    double dy = (gid & 2 ? 0.5 : 0.0);
    double dz = (gid & 1 ? 0.5 : 0.0);
    for (int n=0; n<8; ++n) {
      ASSERT_EQ(cell0Coords[n][0] + dx, myCoords[n][0]);
      ASSERT_EQ(cell0Coords[n][1] + dy, myCoords[n][1]);
      ASSERT_EQ(cell0Coords[n][2] + dz, myCoords[n][2]);
    }
  }

  // Check global IDs
  std::vector<int>& gids = source_mesh_flat.get_global_cell_ids();
  ASSERT_EQ(num_cells, gids.size());
  for (int c=0; c<num_owned_cells; ++c)
    ASSERT_EQ(expOwnedGids[c], gids[c]);

  // Check neighbors
  // Each cell should have all of the 7 other cells as neighbors
  for (int c=0; c<num_cells; ++c) {
    // Get my 7 neighbors
    std::vector<int> myNeighbors;
    source_mesh_flat.cell_get_node_adj_cells(c, Portage::ALL, &myNeighbors);
    ASSERT_EQ(7, myNeighbors.size());
    // Convert to global IDs
    for (int n=0; n<7; ++n) myNeighbors[n] = gids[myNeighbors[n]];
    // Add my own ID
    myNeighbors.push_back(gids[c]);
    // Now make sure all 8 neighbors are present, in any order
    std::sort(myNeighbors.begin(), myNeighbors.end());
    for (int n=0; n<8; ++n)
      ASSERT_EQ(n, myNeighbors[n]);
  }

  // Check field values
  double* ddata1 = nullptr;
  source_state_flat.get_data(Portage::CELL, "d1", &ddata1);
  double* ddata2 = nullptr;
  source_state_flat.get_data(Portage::CELL, "d2", &ddata2);
  for (int c=0; c<num_owned_cells; ++c) {
    int gid = gids[c];
    int expValue1 = 10. + gid;
    ASSERT_EQ(expValue1, ddata1[c]);
    int expValue2 = 100. + gid;
    ASSERT_EQ(expValue2, ddata2[c]);
  }

}


TEST(MPI_Bounding_Boxes, SimpleTest2D) {

  int commRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &commRank);

  // Add multiple state vector types
  double base1 = 10. + commRank * 4.;
  double dtest1[9] = {base1, base1 + 1., base1 + 2., base1 + 3., 0., 0., 0., 0., 0.};
  double base2 = 100. + commRank * 4.;
  double dtest2[9] = {base2, base2 + 1., base2 + 2., base2 + 3., 0., 0., 0., 0., 0.};

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  // Source mesh
  // NOTE:  The gid's of the source mesh are as follows:
  //   5  7 13 15
  //   4  6 12 14
  //   1  3  9 11
  //   0  2  8 10
  // It's tricky to reason from gid's to geometry, so in several
  // places I'll use the "cartesian id" or "cid" instead:
  //   3  7 11 15
  //   2  6 10 14
  //   1  5  9 13
  //   0  4  8 12
  // This is computed by swapping the middle 2 bits of the 4-bit gid.
  std::shared_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4);
  Portage::Jali_Mesh_Wrapper inputMeshWrapper(*source_mesh);
  Jali::State state(source_mesh);
  Portage::Jali_State_Wrapper wrapper(state);

  state.add("d1", source_mesh, Jali::Entity_kind::CELL,
            Jali::Entity_type::ALL, dtest1);
  state.add("d2", source_mesh, Jali::Entity_kind::CELL,
            Jali::Entity_type::ALL, dtest2);

  Portage::Flat_Mesh_Wrapper<> source_mesh_flat(4, inputMeshWrapper);
  Portage::Flat_State_Wrapper<> source_state_flat(wrapper, {"d1", "d2"});

  // Target mesh
  std::shared_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 1.0, 1.0, 3, 3);
  Portage::Jali_Mesh_Wrapper target_mesh_(*target_mesh);
  Jali::State target_state(target_mesh);
  Portage::Jali_State_Wrapper target_state_(target_state);

  // Use a bounding box distributor to send the source cells to the target
  // partitions where they are needed
  Portage::MPI_Bounding_Boxes distributor;
  distributor.distribute(source_mesh_flat, source_state_flat, target_mesh_,
                         target_state_);

  // Check number of cells and nodes received
  int exp_num_owned_cells = (commRank == 0 ? 4 :
                             commRank <= 2 ? 8 : 16);
  int num_owned_cells = source_mesh_flat.num_owned_cells();
  ASSERT_EQ(exp_num_owned_cells, num_owned_cells);
  int exp_num_owned_nodes = exp_num_owned_cells * 4;
  int num_owned_nodes = source_mesh_flat.num_owned_nodes();
  ASSERT_EQ(exp_num_owned_nodes, num_owned_nodes);
  int exp_num_cells = exp_num_owned_cells * 9 / 4;
  int num_cells = num_owned_cells + source_mesh_flat.num_ghost_cells();
  ASSERT_EQ(exp_num_cells, num_cells);
  int exp_num_nodes = exp_num_cells * 4;
  int num_nodes = num_owned_nodes + source_mesh_flat.num_ghost_nodes();
  ASSERT_EQ(exp_num_nodes, num_nodes);

  // Check coordinates
  // List coordinates of cell 0 - others are equal to this
  // with a shift
  std::vector<Portage::Point<2>> cell0Coords =
    {{0.0, 0.0},  {0.25, 0.0},  {0.25, 0.25},  {0.0, 0.25}};
  // List owned cells that should have been sent to each rank
  std::vector<int> expOwnedGids;
  switch (commRank) {
    case 0:
      expOwnedGids = {0, 1, 2, 3};
      break;
    case 1:
      expOwnedGids = {0, 1, 2, 3, 4, 5, 6, 7};
      break;
    case 2:
      expOwnedGids = {0, 1, 2, 3, 8, 9, 10, 11};
      break;
    case 3:
      expOwnedGids = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
      break;
  }

  for (int c=0; c<num_owned_cells; ++c) {
    std::vector<Portage::Point<2>> myCoords;
    source_mesh_flat.cell_get_coordinates(c, &myCoords);
    ASSERT_EQ(4, myCoords.size());
    int gid = expOwnedGids[c];
    int cid = (gid & 9) | ((gid & 4) >> 1) | ((gid & 2) << 1);
    double dx = (cid / 4) * 0.25;
    double dy = (cid % 4) * 0.25;
    for (int n=0; n<4; ++n) {
      ASSERT_EQ(cell0Coords[n][0] + dx, myCoords[n][0]);
      ASSERT_EQ(cell0Coords[n][1] + dy, myCoords[n][1]);
    }
  }

  // Check global IDs
  std::vector<int>& gids = source_mesh_flat.get_global_cell_ids();
  ASSERT_EQ(num_cells, gids.size());
  for (int c=0; c<num_owned_cells; ++c)
    ASSERT_EQ(expOwnedGids[c], gids[c]);

  // Check neighbors
  for (int c=0; c<num_owned_cells; ++c) {
    // Get my neighbors
    std::vector<int> myNeighbors;
    source_mesh_flat.cell_get_node_adj_cells(c, Portage::ALL, &myNeighbors);
    int count = myNeighbors.size();
    // Convert to global IDs
    for (int n=0; n<count; ++n) myNeighbors[n] = gids[myNeighbors[n]];
    // Add my own ID
    myNeighbors.push_back(gids[c]);
    count += 1;
    // Convert to cartesian IDs
    for (int n=0; n<count; ++n) {
      int ngid = myNeighbors[n];
      myNeighbors[n] = (ngid & 9) | ((ngid & 4) >> 1) | ((ngid & 2) << 1);
    }
    // Which neighbors do I expect?
    std::vector<int> expNeighbors;
    int gid = gids[c];
    int cid = (gid & 9) | ((gid & 4) >> 1) | ((gid & 2) << 1);
    // Loop through all possible neighbors, skipping if on boundary
    for (int dx=-1; dx<=+1; ++dx) {
      if (dx == -1 && cid / 4 == 0) continue;
      if (dx == +1 && cid / 4 == 3) continue;
      for (int dy=-1; dy<=+1; ++dy) {
        if (dy == -1 && cid % 4 == 0) continue;
        if (dy == +1 && cid % 4 == 3) continue;
        expNeighbors.push_back(cid + dx * 4 + dy);
      }
    }
    // Now make sure all expected neighbors are present, in any order
    ASSERT_EQ(expNeighbors.size(), count);
    std::sort(myNeighbors.begin(), myNeighbors.end());
    for (int n=0; n<count; ++n)
      ASSERT_EQ(expNeighbors[n], myNeighbors[n]);
  }

  // Check field values
  double* ddata1 = nullptr;
  source_state_flat.get_data(Portage::CELL, "d1", &ddata1);
  double* ddata2 = nullptr;
  source_state_flat.get_data(Portage::CELL, "d2", &ddata2);
  for (int c=0; c<num_owned_cells; ++c) {
    int gid = gids[c];
    int expValue1 = 10. + gid;
    ASSERT_EQ(expValue1, ddata1[c]);
    int expValue2 = 100. + gid;
    ASSERT_EQ(expValue2, ddata2[c]);
  }

}

