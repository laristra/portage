



#include <iostream>
#include <vector>

#include "portage/simple_mesh/simple_mesh.h"

#include "portage/support/portage.h"
#include "portage/support/Point.h"

#include "gtest/gtest.h"
#include "mpi.h"

/*!
  @file test_simple_mesh.cc
  @brief Tests for the simple mesh definition in simple_mesh.h
*/
TEST(Simple_Mesh, SingleCell) {
  // Create a single cell mesh
  Portage::Simple_Mesh mesh(0.0, 0.0, 0.0,
                            1.0, 1.0, 1.0,
                            1, 1, 1);
  // Dimensionality
  ASSERT_EQ(mesh.space_dimension(), 3);
  // Cells
  ASSERT_EQ(mesh.num_entities(Portage::Entity_kind::CELL,
                              Portage::Entity_type::PARALLEL_OWNED),
            1);
  ASSERT_EQ(mesh.num_entities(Portage::Entity_kind::CELL,
                              Portage::Entity_type::PARALLEL_GHOST),
            0);
  // Faces
  ASSERT_EQ(mesh.num_entities(Portage::Entity_kind::FACE,
                              Portage::Entity_type::PARALLEL_OWNED),
            6);
  ASSERT_EQ(mesh.num_entities(Portage::Entity_kind::FACE,
                              Portage::Entity_type::PARALLEL_GHOST),
            0);
  // nodes
  ASSERT_EQ(mesh.num_entities(Portage::Entity_kind::NODE,
                              Portage::Entity_type::PARALLEL_OWNED),
            8);
  ASSERT_EQ(mesh.num_entities(Portage::Entity_kind::NODE,
                              Portage::Entity_type::PARALLEL_GHOST),
            0);

  // These are global id's - the global ordering goes x first, then y, then z
  std::vector<int> nodes;
  // For a cell, the ordering is counterclockwise starting at (x0,y0,z0)
  mesh.cell_get_nodes(0, &nodes);
  // So, `nodes` should have the following global id order:
  std::vector<int> expectedNodes = {0, 1, 3, 2, 4, 5, 7, 6};
  for (int i(0); i < 8; ++i)
    ASSERT_EQ(nodes[i], expectedNodes[i]);

  // Coordinate locations of cell node 5: (1.0, 0.0, 1.0)
  Portage::Point<3> pp;
  mesh.node_get_coordinates(nodes[5], &pp);
  ASSERT_EQ(1.0, pp[0]);
  ASSERT_EQ(0.0, pp[1]);
  ASSERT_EQ(1.0, pp[2]);
}

TEST(Simple_Mesh, SmallGrid) {
  // Create a simple 2x2x2 mesh
  Portage::Simple_Mesh mesh(0.0, 0.0, 0.0,
                            1.0, 1.0, 1.0,
                            2, 2, 2);
  // Number of cells
  ASSERT_EQ(mesh.num_entities(Portage::Entity_kind::CELL,
                              Portage::Entity_type::PARALLEL_OWNED),
            8);
  // Number of faces - interior faces only counted once
  ASSERT_EQ(mesh.num_entities(Portage::Entity_kind::FACE,
                              Portage::Entity_type::PARALLEL_OWNED),
            36);
  // Number of nodes
  ASSERT_EQ(mesh.num_entities(Portage::Entity_kind::NODE,
                              Portage::Entity_type::PARALLEL_OWNED),
            27);

  // Check that the global ID of a node in the center - shared by cells -
  // is the same.
  // This is the global ID of the node in the center of the mesh
  int gid_nodeWant = 13;
  std::vector<int> cellnodes;
  // This is the cell-local index which should correspond to the global index
  std::vector<int> cell_local_index = {6, 7, 5, 4, 2, 3, 1, 0};
  for (int i(0); i < 8; ++i) {
    // Get the nodes associated with this cell
    mesh.cell_get_nodes(i, &cellnodes);
    // Make sure that we match the expected global id
    ASSERT_EQ(gid_nodeWant, cellnodes[cell_local_index[i]]);
  }
}
