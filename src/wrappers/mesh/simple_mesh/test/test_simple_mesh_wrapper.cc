/*---------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/
#include "portage/wrappers/mesh/simple_mesh/simple_mesh_wrapper.h"
#include "portage/simple_mesh/simple_mesh.h"
#include "portage/support/portage.h"
#include "portage/support/Point.h"

#include <iostream>
#include <vector>
#include <iterator>

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
  // Nodes for this cell; the expected are global ids
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
  // Cell node coordinates; this has ordering of cellnodes
  std::vector<Portage::Point<3>> expectedCoords = {{0.0, 0.0, 0.0},
                                                   {1.0, 0.0, 0.0},
                                                   {1.0, 1.0, 0.0},
                                                   {0.0, 1.0, 0.0},
                                                   {0.0, 0.0, 1.0},
                                                   {1.0, 0.0, 1.0},
                                                   {1.0, 1.0, 1.0},
                                                   {0.0, 1.0, 1.0}};
  Portage::Point<3> nodeCoord;
  for (int i(0); i < 8; ++i) {
    mesh_wrapper.node_get_coordinates(i, &nodeCoord);
    for (int d(0); d < 3; ++d)
      ASSERT_EQ(expectedCoords[i][d], nodeCoord[d]);
  }
}
