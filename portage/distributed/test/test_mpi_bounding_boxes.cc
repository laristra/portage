/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <iostream>
#include <memory>
#include <vector>
#include <algorithm>

#include "gtest/gtest.h"

#include "mpi.h"

// portage includes
#include "portage/support/portage.h"
#include "portage/distributed/mpi_bounding_boxes.h"

// Jali includes
#include "Mesh.hh"
#include "MeshFactory.hh"

// wonton includes
#include "wonton/state/jali/jali_state_wrapper.h"
#include "wonton/state/flat/flat_state_mm_wrapper.h"
#include "wonton/mesh/flat/flat_mesh_wrapper.h"
#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"

TEST(MPI_Bounding_Boxes, SimpleTest3D) {

  int commRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &commRank);

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  // Source mesh
  std::shared_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
        2, 2, 2);
  Wonton::Jali_Mesh_Wrapper inputMeshWrapper(*source_mesh);
  Wonton::Flat_Mesh_Wrapper<> source_mesh_flat;
  source_mesh_flat.initialize(inputMeshWrapper);
  
  // Add multiple state vector types
  std::vector<int>& gids = source_mesh_flat.get_global_cell_ids();
  
  // states are merged based on the gid of each cell, so we need to be consistent
  // across ranks. the easiest way to generate fields like this is to make the
  // field a function of gid, that way we are guaranteed to be consistent
  double dtest1[8];
  for (int i=0; i<8; ++i) dtest1[i]=gids[i]+10.;
  double dtest2[8];
  for (int i=0; i<8; ++i) dtest2[i]=gids[i]*gids[i]+100.;

  std::shared_ptr<Jali::State> state(Jali::State::create(source_mesh));
  Wonton::Jali_State_Wrapper wrapper(*state);

  state->add("d1", source_mesh, Jali::Entity_kind::CELL,
             Jali::Entity_type::ALL, dtest1);
  state->add("d2", source_mesh, Jali::Entity_kind::CELL,
             Jali::Entity_type::ALL, dtest2);
  
  Wonton::Flat_State_Wrapper<Wonton::Flat_Mesh_Wrapper<>> source_state_flat(source_mesh_flat);
  
  source_state_flat.initialize(wrapper, {"d1", "d2"});

  // Target mesh
  std::shared_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
        3, 3, 3);
  Wonton::Jali_Mesh_Wrapper target_mesh_(*target_mesh);
  std::shared_ptr<Jali::State> target_state(Jali::State::create(target_mesh));
  Wonton::Jali_State_Wrapper target_state_(*target_state);

  // Use a bounding box distributor to send the source cells to the target
  // partitions where they are needed
  Portage::MPI_Bounding_Boxes distributor;
  distributor.distribute(source_mesh_flat, source_state_flat, target_mesh_,
                         target_state_);

  // After distributing, the meshes are merged and all entities are considered
  // owned. There is no reason to make a distinction between owned and ghost 
  // entities.
  
  // all target ranks get all 8 cells
  int num_owned_cells = source_mesh_flat.num_owned_cells();
  ASSERT_EQ(8, num_owned_cells);

  // all target ranks get all 27 cells  
  int num_owned_nodes = source_mesh_flat.num_owned_nodes();
  ASSERT_EQ(27, num_owned_nodes);
  
  // all target ranks get all 36 faces 
  int num_owned_faces = source_mesh_flat.num_owned_faces();
  ASSERT_EQ(36, num_owned_faces);
  
  // there are no ghosts, all cells are considered owned
  int num_cells = num_owned_cells + source_mesh_flat.num_ghost_cells();
  ASSERT_EQ(8, num_cells);
  
  // there are no ghosts, all nodes are considered owned
  int num_nodes = num_owned_nodes + source_mesh_flat.num_ghost_nodes();
  ASSERT_EQ(27, num_nodes);
  
  // there are no ghosts, all faces are considered owned  
  int num_faces = num_owned_faces + source_mesh_flat.num_ghost_faces();
  ASSERT_EQ(36, num_faces);

  // Check coordinates
  // List coordinates of cell 0 - others are equal to this
  // with a shift
  std::vector<Wonton::Point<3>> cell0Coords =
    {{0.0, 0.0, 0.0},  {0.0, 0.0, 0.5},  {0.0, 0.5, 0.0},  {0.0, 0.5, 0.5},
     {0.5, 0.0, 0.0},  {0.5, 0.0, 0.5},  {0.5, 0.5, 0.0},  {0.5, 0.5, 0.5}};
     
  std::vector<int> expOwnedGids = {0, 1, 2, 3, 4, 5, 6, 7};
  
  // Check global IDs
  std::vector<int>& cell_gids = source_mesh_flat.get_global_cell_ids();
  ASSERT_EQ(num_cells, cell_gids.size());
  for (int c=0; c<num_owned_cells; ++c)
    ASSERT_EQ(expOwnedGids[c], cell_gids[c]);

  
  // check coordinates
  for (int c=0; c<num_owned_cells; ++c) {
    std::vector<Wonton::Point<3>> myCoords;
    source_mesh_flat.cell_get_coordinates(c, &myCoords);
    ASSERT_EQ(8, myCoords.size());
    std::sort(myCoords.begin(), myCoords.end());
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
  
  // Check neighbors
  // Each owned cell should have all of the 7 other cells as neighbors
  for (int c=0; c<num_owned_cells; ++c) {
    // Get my 7 neighbors
    std::vector<int> myNeighbors;
    source_mesh_flat.cell_get_node_adj_cells(c, Portage::Entity_type::ALL, &myNeighbors);
    ASSERT_EQ(7, myNeighbors.size());
    // Convert to global IDs
    for (int n=0; n<7; ++n) myNeighbors[n] = cell_gids[myNeighbors[n]];
    // Add my own ID
    myNeighbors.push_back(cell_gids[c]);
    // Now make sure all 8 neighbors are present, in any order
    std::sort(myNeighbors.begin(), myNeighbors.end());
    for (int n=0; n<8; ++n)
      ASSERT_EQ(n, myNeighbors[n]);
  }
  
  // Check field values
  double* ddata1 = nullptr;
  source_state_flat.mesh_get_data(Portage::Entity_kind::CELL, "d1", &ddata1);
  double* ddata2 = nullptr;
  source_state_flat.mesh_get_data(Portage::Entity_kind::CELL, "d2", &ddata2);

  for (int c=0; c<num_owned_cells; ++c) {
    int gid = cell_gids[c];
    int expValue1 = 10. + gid;
    ASSERT_EQ(expValue1, ddata1[c]);
    int expValue2 = 100. + gid*gid;
    ASSERT_EQ(expValue2, ddata2[c]);
  }

}


TEST(MPI_Bounding_Boxes, SimpleTest2D) {

  int commRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &commRank);

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
  Wonton::Jali_Mesh_Wrapper inputMeshWrapper(*source_mesh);
  
  Wonton::Flat_Mesh_Wrapper<> source_mesh_flat;
  source_mesh_flat.initialize(inputMeshWrapper);

  // Create state data for source mesh
  // Add multiple state vector types
  std::vector<int>& gids = source_mesh_flat.get_global_cell_ids();

  // states are merged based on the gid of each cell, so we need to be consistent
  // across ranks. the easiest way to generate fields like this is to make the
  // field a function of gid, that way we are guaranteed to be consistent
  
  double dtest1[16];
  for (int i=0; i<gids.size(); ++i) dtest1[i]=gids[i]+10.;
  double dtest2[16];
  for (int i=0; i<gids.size(); ++i) dtest2[i]=gids[i]*gids[i]+100.;
  
  std::shared_ptr<Jali::State> state(Jali::State::create(source_mesh));
  state->add("d1", source_mesh, Jali::Entity_kind::CELL,
             Jali::Entity_type::ALL, dtest1);
  state->add("d2", source_mesh, Jali::Entity_kind::CELL,
             Jali::Entity_type::ALL, dtest2);

  Wonton::Jali_State_Wrapper wrapper(*state);
  Wonton::Flat_State_Wrapper<Wonton::Flat_Mesh_Wrapper<>> source_state_flat(source_mesh_flat);
  source_state_flat.initialize(wrapper, {"d1", "d2"});
  
  // Target mesh
  std::shared_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 1.0, 1.0, 3, 3);
  Wonton::Jali_Mesh_Wrapper target_mesh_(*target_mesh);
  std::shared_ptr<Jali::State> target_state(Jali::State::create(target_mesh));
  Wonton::Jali_State_Wrapper target_state_(*target_state);

  // Use a bounding box distributor to send the source cells to the target
  // partitions where they are needed
  Portage::MPI_Bounding_Boxes distributor;
  distributor.distribute(source_mesh_flat, source_state_flat, target_mesh_,
                         target_state_);

  // After distributing, the meshes are merged and all entities are considered
  // owned. There is no reason to make a distinction between owned and ghost 
  // entities.
 
  //Reference counts for entities on each rank 
  int ref_owned_cells[4] = {9, 12, 12, 16};
  int ref_owned_nodes[4] = {16, 20, 20, 25};
  int ref_owned_faces[4] = {24, 31, 31, 40};
 
  // all target ranks get all 16 cells
  int num_owned_cells = source_mesh_flat.num_owned_cells();
  ASSERT_EQ(ref_owned_cells[commRank], num_owned_cells);
  
  // all target ranks get all 25 nodes
  int num_owned_nodes = source_mesh_flat.num_owned_nodes();
  ASSERT_EQ(ref_owned_nodes[commRank], num_owned_nodes);

  // all target ranks get all 40 faces
  int num_owned_faces = source_mesh_flat.num_owned_faces();
  ASSERT_EQ(ref_owned_faces[commRank], num_owned_faces);

  // there are no ghosts, all cells are considered owned
  int num_cells = num_owned_cells + source_mesh_flat.num_ghost_cells();
  ASSERT_EQ(ref_owned_cells[commRank], num_cells);

  // there are no ghosts, all nodes are considered owned
  int num_nodes = num_owned_nodes + source_mesh_flat.num_ghost_nodes();
  ASSERT_EQ(ref_owned_nodes[commRank], num_nodes);

  // there are no ghosts, all faces are considered owned
  int num_faces = num_owned_faces + source_mesh_flat.num_ghost_faces();
  ASSERT_EQ(ref_owned_faces[commRank], num_faces);
  
  // Reference List of owned cells on each rank after distribution
  std::vector<int> expOwnedCellGids[4] = 
  {{0,1,2,3,4,6,8,9,12},
   {0,1,2,3,4,5,6,7,8,9,12,13},
   {0,1,2,3,4,6,8,9,10,11,12,14},
   {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}};
  // Check coordinates
  // List coordinates of cell 0 - others are equal to this
  // with a shift
  std::vector<Wonton::Point<2>> cell0Coords =
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
    std::vector<Wonton::Point<2>> myCoords;
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
>>>>>>> 1a42de35621ea75ef9b54d2f223422412b36e781

  // Check global IDs
  std::vector<int>& cell_gids = source_mesh_flat.get_global_cell_ids();
  ASSERT_EQ(num_cells, cell_gids.size());
  for (int c=0; c<num_owned_cells; ++c)
    ASSERT_EQ(expOwnedCellGids[commRank][c], cell_gids[c]);

   //Reference list of vertices
   std::vector<Portage::Point<2>> ref_nodes = 
   {{0.0, 0.0}, {0.25, 0.0}, {0.5, 0.0}, {0.0, 0.25}, {0.25, 0.25},
    {0.5, 0.25}, {0.0, 0.5}, {0.25, 0.5}, {0.5, 0.5}, {0.0, 0.75},
    {0.25, 0.75}, {0.5, 0.75}, {0.0, 1.0}, {0.25, 1.0}, {0.5, 1.0},
    {0.75, 0.0}, {1.0, 0.0}, {0.75, 0.25}, {1.0, 0.25}, {0.75, 0.5},
    {1.0, 0.5}, {0.75, 0.75}, {1.0, 0.75}, {0.75, 1.0}, {1.0, 1.0}};

    std::vector<std::vector<int>> ref_cell_to_nodes = 
    {{0,1,4,3}, {3,4,7,6}, {1,2,5,4}, {4,5,8,7},
     {6,7,10,9}, {9,10,13,12}, {7,8,11,10}, {10,11,14,13},
     {2,15,17,5}, {5,17,19,8}, {15,16,18,17}, {17,18,20,19},
     {8,19,21,11}, {11,21,23,14}, {19,20,22,21}, {21,22,24,23}};  

 
    // Check coordinates 
    for (int c=0; c<num_owned_cells; ++c) {
      std::vector<Portage::Point<2>> mycoords;
      source_mesh_flat.cell_get_coordinates(c, &mycoords);
      ASSERT_EQ(4, mycoords.size());
      int gid = expOwnedCellGids[commRank][c];

      for (int n=0; n<4; ++n) {
        int nid = ref_cell_to_nodes[gid][n];
        ASSERT_EQ(ref_nodes[nid][0], mycoords[n][0]);
        ASSERT_EQ(ref_nodes[nid][1], mycoords[n][1]);
      }
    }
  
  // Check neighbors
  for (int c=0; c<num_owned_cells; ++c) {
    // Get my neighbors
    std::vector<int> myNeighbors;
    source_mesh_flat.cell_get_node_adj_cells(c, Portage::Entity_type::ALL, &myNeighbors);
    int count = myNeighbors.size();
    // Convert to global IDs
    for (int n=0; n<count; ++n) myNeighbors[n] = cell_gids[myNeighbors[n]];
    // Add my own ID
    myNeighbors.push_back(cell_gids[c]);
    count += 1;
    // Convert to cartesian IDs
    for (int n=0; n<count; ++n) {
      int ngid = myNeighbors[n];
      myNeighbors[n] = (ngid & 9) | ((ngid & 4) >> 1) | ((ngid & 2) << 1);
    }
    // Which neighbors do I expect?
    std::vector<int> expNeighbors;
    int gid = cell_gids[c];
    int cid = (gid & 9) | ((gid & 4) >> 1) | ((gid & 2) << 1);
   
    // Loop through all possible neighbors, skipping if on boundary
     if (commRank == 0)
     {
      for (int dx=-1; dx<=+1; ++dx) {
        if (dx == -1 && cid / 4 == 0) continue;
        if (dx == +1 && cid / 4 == 2) continue;
        for (int dy=-1; dy<=+1; ++dy) {
          if (dy == -1 && cid % 4 == 0) continue;
          if (dy == +1 && cid % 4 == 2) continue;
          expNeighbors.push_back(cid + dx * 4 + dy);
        }
      } 
     }
     else if (commRank == 1)
     {
      for (int dx=-1; dx<=+1; ++dx) {
        if (dx == -1 && cid / 4 == 0) continue;
        if (dx == +1 && cid / 4 == 2) continue;
        for (int dy=-1; dy<=+1; ++dy) {
          if (dy == -1 && cid % 4 == 0) continue;
          if (dy == +1 && cid % 4 == 3) continue;
          expNeighbors.push_back(cid + dx * 4 + dy);
        }
      } 
     }
     else if (commRank == 2)
     {
      for (int dx=-1; dx<=+1; ++dx) {
        if (dx == -1 && cid / 4 == 0) continue;
        if (dx == +1 && cid / 4 == 3) continue;
        for (int dy=-1; dy<=+1; ++dy) {
          if (dy == -1 && cid % 4 == 0) continue;
          if (dy == +1 && cid % 4 == 2) continue;
          expNeighbors.push_back(cid + dx * 4 + dy);
        }
      } 
     }
     else if (commRank == 3)
     {
      for (int dx=-1; dx<=+1; ++dx) {
        if (dx == -1 && cid / 4 == 0) continue;
        if (dx == +1 && cid / 4 == 3) continue;
        for (int dy=-1; dy<=+1; ++dy) {
          if (dy == -1 && cid % 4 == 0) continue;
          if (dy == +1 && cid % 4 == 3) continue;
          expNeighbors.push_back(cid + dx * 4 + dy);
        }
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
  source_state_flat.mesh_get_data(Portage::Entity_kind::CELL, "d1", &ddata1);
  double* ddata2 = nullptr;
  source_state_flat.mesh_get_data(Portage::Entity_kind::CELL, "d2", &ddata2);
  
  for (int c=0; c<num_owned_cells; ++c) {
    int gid = cell_gids[c];
    int expValue1 = 10. + gid;
    ASSERT_EQ(expValue1, ddata1[c]);
    int expValue2 = 100. + gid*gid;
    ASSERT_EQ(expValue2, ddata2[c]);
  }
  
}


