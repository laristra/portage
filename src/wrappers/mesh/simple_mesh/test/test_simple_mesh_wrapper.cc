/*---------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/
#include <iostream>
#include <vector>
#include <iterator>

#include "portage/wrappers/mesh/simple_mesh/simple_mesh_wrapper.h"
#include "portage/simple_mesh/simple_mesh.h"
#include "portage/support/portage.h"
#include "portage/support/Point.h"

#include "gtest/gtest.h"
#include "mpi.h"

TEST(Simple_Mesh, OneCell) {
  Portage::Simple_Mesh mesh(0.0, 0.0, 0.0,
                            1.0, 1.0, 1.0,
                            1, 1, 1);
  Portage::Simple_Mesh_Wrapper mesh_wrapper(mesh);

  // Check basic dimensionality
  ASSERT_EQ(mesh_wrapper.space_dimension(), 3);

  // Basic ownership
  ASSERT_EQ(mesh_wrapper.num_owned_cells(), 1);
  ASSERT_EQ(mesh_wrapper.num_owned_nodes(), 8);
  ASSERT_EQ(mesh_wrapper.num_owned_faces(), 6);
  ASSERT_EQ(mesh_wrapper.num_ghost_cells(), 0);
  ASSERT_EQ(mesh_wrapper.num_ghost_nodes(), 0);
  ASSERT_EQ(mesh_wrapper.num_ghost_faces(), 0);

  // Cell information
  ASSERT_EQ(mesh_wrapper.cell_get_type(0),
            Portage::Entity_type::PARALLEL_OWNED);
  ASSERT_EQ(mesh_wrapper.cell_get_element_type(0),
            Portage::Element_type::HEX);

  // Connectivity information
  // Nodes for this cell; the expected are global ids and the ordering is
  // local to the cell
  std::vector<int> expectedNodes = {0, 1, 3, 2, 4, 5, 7, 6};
  std::vector<int> cellnodes;
  mesh_wrapper.cell_get_nodes(0, &cellnodes);
  for (int i(0); i < 8; ++i)
    ASSERT_EQ(expectedNodes[i], cellnodes[i]);
  // Faces and dirs; the expected are global ids
  std::vector<int> expectedFaces = {2, 5, 3, 4, 0, 1};
  std::vector<int> expectedDirs = {-1, 1, 1, -1, -1, 1};
  std::vector<int> cellfaces, cellfaceDirs;
  mesh_wrapper.cell_get_faces_and_dirs(0, &cellfaces, &cellfaceDirs);
  for (int i(0); i < 6; ++i) {
    ASSERT_EQ(expectedFaces[i], cellfaces[i]);
    ASSERT_EQ(expectedDirs[i], cellfaceDirs[i]);
  }
  // Nodes of the faces; the expected are global ids
  std::vector<std::vector<int>> expectedFaceNodes = {{0, 1, 3, 2},
                                                     {4, 5, 7, 6},
                                                     {0, 1, 5, 4},
                                                     {2, 3, 7, 6},
                                                     {0, 2, 6, 4},
                                                     {1, 3, 7, 5}};
  std::vector<int> facenodes;
  for (int f(0); f < 6; ++f) {
    mesh_wrapper.face_get_nodes(f, &facenodes);
    for (int i(0); i < 4; ++i)
      ASSERT_EQ(expectedFaceNodes[f][i], facenodes[i]);
  }

  // Coordinate information
  // Ordering is the GLOBAL ordering
  std::vector<Portage::Point<3>> expectedCoords = {{0.0, 0.0, 0.0},
                                                   {1.0, 0.0, 0.0},
                                                   {0.0, 1.0, 0.0},
                                                   {1.0, 1.0, 0.0},
                                                   {0.0, 0.0, 1.0},
                                                   {1.0, 0.0, 1.0},
                                                   {0.0, 1.0, 1.0},
                                                   {1.0, 1.0, 1.0}};
  Portage::Point<3> nodeCoord;
  for (int i(0); i < 8; ++i) {
    mesh_wrapper.node_get_coordinates(i, &nodeCoord);
    for (int d(0); d < 3; ++d)
      ASSERT_EQ(expectedCoords[i][d], nodeCoord[d]);
  }

  // Check for nodes in cells adjacent to this
  // In a single cell setup, this should only contain the other nodes
  // of the cell.
  std::vector<int> adjnodes;
  for (int i(0); i < 8; ++i) {
    mesh_wrapper.node_get_cell_adj_nodes(i,
					 Portage::Entity_type::PARALLEL_OWNED,
					 &adjnodes);
    ASSERT_EQ(adjnodes.size(), 7);
  }

  // Check for cells adjacent to this cell; adjacency determined from
  // nodes.
  // FOr a single cell setup, this should be empty
  std::vector<int> adjcells;
  mesh_wrapper.cell_get_node_adj_cells(0,
				       Portage::Entity_type::PARALLEL_OWNED,
				       &adjcells);
  ASSERT_EQ(adjcells.size(), 0);
}


TEST(Simple_Mesh, MultiCell) {
  Portage::Simple_Mesh mesh(0.0, 0.0, 0.0,
                            1.0, 1.0, 1.0,
                            2, 2, 2);
  Portage::Simple_Mesh_Wrapper mesh_wrapper(mesh);

  // Check basic dimensionality
  ASSERT_EQ(mesh_wrapper.space_dimension(), 3);

  // Basic ownership
  int ncells = mesh_wrapper.num_owned_cells();
  int nnodes = mesh_wrapper.num_owned_nodes();
  int nfaces = mesh_wrapper.num_owned_faces();

  ASSERT_EQ(ncells, 8);
  ASSERT_EQ(nnodes, 27);
  ASSERT_EQ(nfaces, 36);
  ASSERT_EQ(mesh_wrapper.num_ghost_cells(), 0);
  ASSERT_EQ(mesh_wrapper.num_ghost_nodes(), 0);
  ASSERT_EQ(mesh_wrapper.num_ghost_faces(), 0);

  // Cell information
  for (int i(0); i < ncells; ++i) {
    ASSERT_EQ(mesh_wrapper.cell_get_type(i),
	      Portage::Entity_type::PARALLEL_OWNED);
    ASSERT_EQ(mesh_wrapper.cell_get_element_type(i),
	      Portage::Element_type::HEX);
  }

  // Connectivity information
  // Nodes for each cell; the expected are global ids and the ordering is
  // local to the cell
  std::vector<std::vector<int>> expectedNodes = {
    {0, 1, 4, 3, 9, 10, 13, 12},
    {1, 2, 5, 4, 10, 11, 14, 13},
    {3, 4, 7, 6, 12, 13, 16, 15},
    {4, 5, 8, 7, 13, 14, 17, 16},
    {9, 10, 13, 12, 18, 19, 22, 21},
    {10, 11, 14, 13, 19, 20, 23, 22},
    {12, 13, 16, 15, 21, 22, 25, 24},
    {13, 14, 17, 16, 22, 23, 26, 25}};
  std::vector<int> cellnodes;
  for (int i(0); i < ncells; ++i) {
    mesh_wrapper.cell_get_nodes(i, &cellnodes);
    for (int n(0); n < 8; ++n)
      ASSERT_EQ(expectedNodes[i][n], cellnodes[n]);
  }

  std::vector<std::vector<int>> expectedFaces = {
    {12, 25, 14, 24, 0, 4},
    {13, 26, 15, 25, 1, 5},
    {14, 28, 16, 27, 2, 6},
    {15, 29, 17, 28, 3, 7},
    {18, 31, 20, 30, 4, 8},
    {19, 32, 21, 31, 5, 9},
    {20, 34, 22, 33, 6, 10},
    {21, 35, 23, 34, 7, 11}
  };
  std::vector<int> expectedDirs = {-1, 1, 1, -1, -1, 1};
  std::vector<int> cellfaces, cellfacedirs;
  for (int i(0); i < ncells; ++i) {
    mesh_wrapper.cell_get_faces_and_dirs(i, &cellfaces, &cellfacedirs);
    for (int f(0); f < 6; ++f) {
      ASSERT_EQ(expectedFaces[i][f], cellfaces[f]);
      ASSERT_EQ(expectedDirs[f], cellfacedirs[f]);
    }
  }
  // // Nodes of the faces; the expected are global ids
  // std::vector<std::vector<int>> expectedFaceNodes = {{0, 1, 3, 2},
  //                                                    {4, 5, 7, 6},
  //                                                    {0, 1, 5, 4},
  //                                                    {2, 3, 7, 6},
  //                                                    {0, 2, 6, 4},
  //                                                    {1, 3, 7, 5}};
  // std::vector<int> facenodes;
  // for (int f(0); f < 6; ++f) {
  //   mesh_wrapper.face_get_nodes(f, &facenodes);
  //   for (int i(0); i < 4; ++i)
  //     ASSERT_EQ(expectedFaceNodes[f][i], facenodes[i]);
  // }

  // // Coordinate information
  // // Ordering is the GLOBAL ordering
  // std::vector<Portage::Point<3>> expectedCoords = {{0.0, 0.0, 0.0},
  //                                                  {1.0, 0.0, 0.0},
  //                                                  {0.0, 1.0, 0.0},
  //                                                  {1.0, 1.0, 0.0},
  //                                                  {0.0, 0.0, 1.0},
  //                                                  {1.0, 0.0, 1.0},
  //                                                  {0.0, 1.0, 1.0},
  //                                                  {1.0, 1.0, 1.0}};
  // Portage::Point<3> nodeCoord;
  // for (int i(0); i < 8; ++i) {
  //   mesh_wrapper.node_get_coordinates(i, &nodeCoord);
  //   for (int d(0); d < 3; ++d)
  //     ASSERT_EQ(expectedCoords[i][d], nodeCoord[d]);
  // }

  // // Check for nodes in cells adjacent to this
  // // In a single cell setup, this should only contain the other nodes
  // // of the cell.
  // std::vector<int> adjnodes;
  // for (int i(0); i < 8; ++i) {
  //   mesh_wrapper.node_get_cell_adj_nodes(i,
  // 					 Portage::Entity_type::PARALLEL_OWNED,
  // 					 &adjnodes);
  //   ASSERT_EQ(adjnodes.size(), 7);
  // }

  // // Check for cells adjacent to this cell; adjacency determined from
  // // nodes.
  // // FOr a single cell setup, this should be empty
  // std::vector<int> adjcells;
  // mesh_wrapper.cell_get_node_adj_cells(0,
  // 				       Portage::Entity_type::PARALLEL_OWNED,
  // 				       &adjcells);
  // ASSERT_EQ(adjcells.size(), 0);
}
