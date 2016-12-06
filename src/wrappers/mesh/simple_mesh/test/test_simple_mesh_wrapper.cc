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
  std::vector<int> expectedDirs = {1, 1, -1, -1, -1, 1};
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
  // For a single cell setup, this should be empty
  std::vector<int> adjcells;
  mesh_wrapper.cell_get_node_adj_cells(0,
                                       Portage::Entity_type::PARALLEL_OWNED,
                                       &adjcells);
  ASSERT_EQ(adjcells.size(), 0);
}



TEST(Simple_Mesh, MultiCell) {
  // Create a 2x2x2 mesh
  double xmin(0.0), ymin(0.0), zmin(0.0);
  double xmax(1.0), ymax(1.0), zmax(1.0);
  int nx(2), ny(2), nz(2);
  Portage::Simple_Mesh mesh(xmin, ymin, zmin,
                            xmax, ymax, zmax,
                            nx, ny, nz);
  Portage::Simple_Mesh_Wrapper mesh_wrapper(mesh);

  // Check basic dimensionality
  ASSERT_EQ(mesh_wrapper.space_dimension(), 3);

  // Basic ownership
  int ncells = mesh_wrapper.num_owned_cells();
  int nnodes = mesh_wrapper.num_owned_nodes();
  int nfaces = mesh_wrapper.num_owned_faces();

  ASSERT_EQ(ncells, nx*ny*nz);
  ASSERT_EQ(nnodes, (nx+1)*(ny+1)*(nz+1));
  ASSERT_EQ(nfaces, nx*ny*(nz+1) + nx*(ny+1)*nz + (nx+1)*ny*nz);
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

  // Faces for each cell; the expected are global ids and the ordering is
  // local to the cell
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
  std::vector<int> expectedDirs = {1, 1, -1, -1, -1, 1};
  std::vector<int> cellfaces, cellfacedirs;
  for (int i(0); i < ncells; ++i) {
    mesh_wrapper.cell_get_faces_and_dirs(i, &cellfaces, &cellfacedirs);
    for (int f(0); f < 6; ++f) {
      ASSERT_EQ(expectedFaces[i][f], cellfaces[f]);
      ASSERT_EQ(expectedDirs[f], cellfacedirs[f]);
    }
  }

  // Nodes of the faces; the expected are global ids and the ordering is
  // local to the face
  std::vector<std::vector<int>> expectedFaceNodes = {
    // xy faces
    // z = 0 plane
    {0, 1, 4, 3},
    {1, 2, 5, 4},
    {3, 4, 7, 6},
    {4, 5, 8, 7},
    // z = 0.5 plane
    {9, 10, 13, 12},
    {10, 11, 14, 13},
    {12, 13, 16, 15},
    {13, 14, 17, 16},
    // z = 1.0 plane
    {18, 19, 22, 21},
    {19, 20, 23, 22},
    {21, 22, 25, 24},
    {22, 23, 26, 25},
    // xz faces
    // y = 0 plane bottom
    {0, 1, 10, 9},
    {1, 2, 11, 10},
    // y = 0.5 plane bottom
    {3, 4, 13, 12},
    {4, 5, 14, 13},
    // y = 1.0 plane bottom
    {6, 7, 16, 15},
    {7, 8, 17, 16},
    // y = 0 plane top
    {9, 10, 19, 18},
    {10, 11, 20, 19},
    // y = 0.5 plane top
    {12, 13, 22, 21},
    {13, 14, 23, 22},
    // y = 1.0 plane top
    {15, 16, 25, 24},
    {16, 17, 26, 25},
    // yz faces
    // x = 0 plane bottom front
    {0, 3, 12, 9},
    // x = 0.5 plane bottom front
    {1, 4, 13, 10},
    // x = 1.0 plane bottom front
    {2, 5, 14, 11},
    // x = 0 plane bottom back
    {3, 6, 15, 12},
    // x = 0.5 plane bottom back
    {4, 7, 16, 13},
    // x = 1.0 plane bottom back
    {5, 8, 17, 14},
    // x = 0 plane top front
    {9, 12, 21, 18},
    // x = 0.5 plane top front
    {10, 13, 22, 19},
    // x = 1.0 plane top front
    {11, 14, 23, 20},
    // x = 0 plane top back
    {12, 15, 24, 21},
    // x = 0.5 plane top back
    {13, 16, 25, 22},
    // x = 1.0 plane top back
    {14, 17, 26, 23}
  };
  std::vector<int> faceNodes;
  for (int f(0); f < nfaces; ++f) {
    mesh_wrapper.face_get_nodes(f, &faceNodes);
    for (int i(0); i < 4; ++i)
      ASSERT_EQ(expectedFaceNodes[f][i], faceNodes[i]);
  }

  // Coordinate information
  // Ordering is the global ordering
  Portage::Point<3> expectedCoords, nodeCoords;
  double dx((xmax-xmin)/nx), dy((ymax-ymin)/ny), dz((zmax-zmin)/nz);
  for (int iz(0); iz < nz+1; ++iz)
    for (int iy(0); iy < ny+1; ++iy)
      for (int ix(0); ix < nx+1; ++ix) {
        expectedCoords = {ix*dx + xmin,
                          iy*dy + ymin,
                          iz*dz + zmin};
        int thisNode = ix + iy*(nx+1) + iz*(nx+1)*(ny+1);
        mesh_wrapper.node_get_coordinates(thisNode, &nodeCoords);
        for (int d(0); d < 3; ++d)
          ASSERT_EQ(expectedCoords[d], nodeCoords[d]);
      }

  // Check for nodes in cells adjacent to each node
  // Ordering is cell first, then nodes within a cell
  std::vector<std::vector<int>> adjExpectedNodes = {
    // node 0; just the others in this cell
    {1, 4, 3, 9, 10, 13, 12},
    // node 1; cells 0,1
    {0, 4, 3, 9, 10, 13, 12, 2, 5, 11, 14},
    // node 2; just the others in this cell
    {1, 5, 4, 10, 11, 14, 13},
    // node 3; cells 0,3
    {0, 1, 4, 9, 10, 13, 12, 7, 6, 16, 15},
    // node 4; cells 0,1,2,3
    {0, 1, 3, 9, 10, 13, 12, 2, 5, 11, 14, 7, 6, 16, 15, 8, 17},
    // node 5; cells 1,3
    {1, 2, 4, 10, 11, 14, 13, 8, 7, 17, 16},
    // node 6; just the others in this cell
    {3, 4, 7, 12, 13, 16, 15},
    // node 7; cells 2,3
    {3, 4, 6, 12, 13, 16, 15, 5, 8, 14, 17},
    // node 8; just the others in this cell
    {4, 5, 7, 13, 14, 17, 16},
    // node 9; cells 0,4
    {0, 1, 4, 3, 10, 13, 12, 18, 19, 22, 21},
    // node 10; cells 0,1,4,5
    {0, 1, 4, 3, 9, 13, 12, 2, 5, 11, 14, 18, 19, 22, 21, 20, 23},
    // node 11; cells 1,5
    {1, 2, 5, 4, 10, 14, 13, 19, 20, 23, 22},
    // node 12; cells 0,2,4,6
    {0, 1, 4, 3, 9, 10, 13, 7, 6, 16, 15, 18, 19, 22, 21, 25, 24},
    // node 13; everyone else
    {0, 1, 4, 3, 9, 10, 12, 2, 5, 11, 14, 7, 6, 16, 15, 8, 17,
     18, 19, 22, 21, 20, 23, 25, 24, 26},
    // node 14; cells 1,3,5,7
    {1, 2, 5, 4, 10, 11, 13, 8, 7, 17, 16, 19, 20, 23, 22, 26, 25},
    // node 15; cells 2,6
    {3, 4, 7, 6, 12, 13, 16, 21, 22, 25, 24},
    // node 16; cells 2,3,6,7
    {3, 4, 7, 6, 12, 13, 15, 5, 8, 14, 17, 21, 22, 25, 24, 23, 26},
    // node 17; cells 3,7
    {4, 5, 8, 7, 13, 14, 16, 22, 23, 26, 25},
    // node 18; just the others in this cell
    {9, 10, 13, 12, 19, 22, 21},
    // node 19; cells 4,5
    {9, 10, 13, 12, 18, 22, 21, 11, 14, 20, 23},
    // node 20; just the others in this cell
    {10, 11, 14, 13, 19, 23, 22},
    // node 21; cells 4,6
    {9, 10, 13, 12, 18, 19, 22, 16, 15, 25, 24},
    // node 22; cells 4,5,6,7
    {9, 10, 13, 12, 18, 19, 21, 11, 14, 20, 23, 16, 15, 25, 24, 17, 26},
    // node 23; cells 5,7
    {10, 11, 14, 13, 19, 20, 22, 17, 16, 26, 25},
    // node 24; just the others in this cell
    {12, 13, 16, 15, 21, 22, 25},
    // node 25; cells 6,7
    {12, 13, 16, 15, 21, 22, 24, 14, 17, 23, 26},
    // node 26; just the others in this cell
    {13, 14, 17, 16, 22, 23, 25}
  };
  std::vector<int> adjNodes;
  for (int i(0); i < nnodes; ++i) {
    int adjExpectedNumNodes = adjExpectedNodes[i].size();

    mesh_wrapper.node_get_cell_adj_nodes(i,
                                         Portage::Entity_type::PARALLEL_OWNED,
                                         &adjNodes);
    ASSERT_EQ(adjExpectedNumNodes, adjNodes.size());

    for (int j(0); j < adjExpectedNumNodes; ++j)
      ASSERT_EQ(adjExpectedNodes[i][j], adjNodes[j]);
  }

  // Check for cells adjacent to this cell, determined by nodes
  // For this 2x2x2 case, each cell is attached to all the others
  int expectedNumAdjCells = ncells-1;
  std::vector<int> adjCells;
  for (int i(0); i < ncells; ++i) {
    mesh_wrapper.cell_get_node_adj_cells(i,
                                         Portage::Entity_type::PARALLEL_OWNED,
                                         &adjCells);
    ASSERT_EQ(expectedNumAdjCells, adjCells.size());
  }
}


TEST(Simple_Mesh, AdjCell) {
  // Create a 2x2x3 mesh
  double xmin(0.0), ymin(0.0), zmin(0.0);
  double xmax(1.0), ymax(1.0), zmax(1.0);
  int nx(2), ny(2), nz(3);
  Portage::Simple_Mesh mesh(xmin, ymin, zmin,
                            xmax, ymax, zmax,
                            nx, ny, nz);
  Portage::Simple_Mesh_Wrapper mesh_wrapper(mesh);

  int ncells = mesh_wrapper.num_owned_cells();
  ASSERT_EQ(ncells, nx*ny*nz);

  // Check for cells adjacent to each cell; adjacency determined from nodes
  std::vector<std::vector<int>> expectedAdjCells = {
    // bottom plane
    {1, 2, 3, 4, 5, 6, 7},
    {0, 3, 2, 4, 5, 7, 6},
    {0, 1, 3, 4, 6, 5, 7},
    {0, 1, 2, 4, 5, 6, 7},
    // middle plane
    {0, 1, 5, 2, 3, 6, 7, 8, 9, 10, 11},
    {0, 1, 4, 3, 7, 2, 6, 8, 9, 11, 10},
    {0, 2, 4, 1, 3, 5, 7, 8, 10, 9, 11},
    {0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11},
    // top plane
    {4, 5, 9, 6, 7, 10, 11},
    {4, 5, 8, 7, 11, 6, 10},
    {4, 6, 8, 5, 7, 9, 11},
    {4, 5, 6, 7, 8, 9, 10}
  };
  std::vector<int> adjCells;
  for (int i(0); i < ncells; ++i) {
    int expectedNumAdjCells = expectedAdjCells[i].size();
    mesh_wrapper.cell_get_node_adj_cells(i,
                                         Portage::Entity_type::PARALLEL_OWNED,
                                         &adjCells);
    ASSERT_EQ(expectedNumAdjCells, adjCells.size());
    for (int j(0); j < expectedNumAdjCells; ++j)
      ASSERT_EQ(expectedAdjCells[i][j], adjCells[j]);
  }
}
