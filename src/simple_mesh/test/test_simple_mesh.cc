/*
Copyright (c) 2016, Los Alamos National Security, LLC
All rights reserved.

Copyright 2016. Los Alamos National Security, LLC. This software was produced
under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National
Laboratory (LANL), which is operated by Los Alamos National Security, LLC for
the U.S. Department of Energy. The U.S. Government has rights to use,
reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS
NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
derivative works, such modified software should be clearly marked, so as not to
confuse it with the version available from LANL.

Additionally, redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. Neither the name of Los Alamos National Security, LLC, Los Alamos
   National Laboratory, LANL, the U.S. Government, nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL
SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/



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
