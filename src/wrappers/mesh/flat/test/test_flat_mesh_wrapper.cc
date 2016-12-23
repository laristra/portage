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



#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"

#include <algorithm>
#include <iostream>

#include "gtest/gtest.h"
#include "mpi.h"

#include "Mesh.hh"
#include "MeshFactory.hh"

#include "portage/support/Point.h"
#include "portage/wrappers/mesh/flat/flat_mesh_wrapper.h"

using std::abs;

/*!
  @file test_flag_mesh_wrapper.cc
  @brief Unit tests for the flat mesh wrapper class
 */


/*!
  @brief Unit test for basic routines in Flat_Mesh_Wrapper in 3D
 */
TEST(Flat_Mesh_Wrapper, basic_routines_3d) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::EDGE,
                        Jali::Entity_kind::FACE,
                        Jali::Entity_kind::WEDGE,
                        Jali::Entity_kind::CORNER});
  std::shared_ptr<Jali::Mesh> mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2);
  ASSERT_TRUE(mesh != NULL);
  Portage::Jali_Mesh_Wrapper mesh_wrapper(*mesh);
  Portage::Flat_Mesh_Wrapper<> mesh_flat;
  mesh_flat.initialize(mesh_wrapper);

  ASSERT_EQ(3, mesh_flat.space_dimension());

  // Test cells, nodes, and faces of flat mesh
  ASSERT_EQ(8, mesh_flat.num_owned_cells());
  ASSERT_EQ(8, mesh_flat.num_entities(Portage::Entity_kind::CELL));
  ASSERT_EQ(27, mesh_flat.num_owned_nodes());
  ASSERT_EQ(27, mesh_flat.num_entities(Portage::Entity_kind::NODE));
  ASSERT_EQ(36, mesh_flat.num_owned_faces());
  ASSERT_EQ(36, mesh_flat.num_entities(Portage::Entity_kind::FACE));

  // Check coordinates
  // List coordinates of cell 0 - others are equal to this
  // with a shift
  std::vector<Portage::Point<3>> cell0Coords =
    {{0.0, 0.0, 0.0},  {0.5, 0.0, 0.0},  {0.5, 0.5, 0.0},  {0.0, 0.5, 0.0},
     {0.0, 0.0, 0.5},  {0.5, 0.0, 0.5},  {0.5, 0.5, 0.5},  {0.0, 0.5, 0.5}};
  for (int c=0; c<8; ++c) {
    std::vector<Portage::Point<3>> pplist;
    mesh_flat.cell_get_coordinates(c, &pplist);
    ASSERT_EQ(8, pplist.size());
    double dx = (c & 4 ? 0.5 : 0.0);
    double dy = (c & 2 ? 0.5 : 0.0);
    double dz = (c & 1 ? 0.5 : 0.0);
    for (int n=0; n<8; ++n) {
      ASSERT_EQ(cell0Coords[n][0] + dx, pplist[n][0]);
      ASSERT_EQ(cell0Coords[n][1] + dy, pplist[n][1]);
      ASSERT_EQ(cell0Coords[n][2] + dz, pplist[n][2]);
    }
  }

  // Test global ids
  std::vector<int>& gids = mesh_flat.get_global_cell_ids();
  for (int c=0; c<8; ++c)
    ASSERT_EQ(c, gids[c]);

  // Test decompose_cell_into_tets()
  std::vector<std::array<Portage::Point<3>, 4>> tcoords;
  mesh_flat.decompose_cell_into_tets(0, &tcoords, true);
  ASSERT_EQ(tcoords.size(), 24);

  // Test centroids
  Portage::Point<3> centroid;
  mesh_flat.cell_centroid(0, &centroid);
  for (int d=0; d<3; ++d)
    ASSERT_EQ(0.25, centroid[d]);

  // Test cell_volume()
  for (int c=0; c<8; ++c)
    ASSERT_TRUE(abs(mesh_flat.cell_volume(c)-0.125) < 1e-12);

  // Test cell neighbors
  for (int c=0; c<8; ++c) {
    std::vector<int> adjcells;
    mesh_flat.cell_get_node_adj_cells(c, Portage::Entity_type::ALL,
                                      &adjcells);
    ASSERT_EQ(7, adjcells.size());
    std::sort(adjcells.begin(), adjcells.end());
    for (int n=0; n<7; ++n)
      ASSERT_EQ(n + (n >= c ? 1 : 0), adjcells[n]);
  }

}


/*!
  @brief Unit test for basic routines in Flat_Mesh_Wrapper in 2D
 */
TEST(Flat_Mesh_Wrapper, basic_routines_2d) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::EDGE,
                        Jali::Entity_kind::FACE,
                        Jali::Entity_kind::WEDGE,
                        Jali::Entity_kind::CORNER});
  std::shared_ptr<Jali::Mesh> mesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4);
  ASSERT_TRUE(mesh != NULL);
  Portage::Jali_Mesh_Wrapper mesh_wrapper(*mesh);
  Portage::Flat_Mesh_Wrapper<> mesh_flat;
  mesh_flat.initialize(mesh_wrapper);
  ASSERT_EQ(2, mesh_flat.space_dimension());

  // Test cells and nodes of flat mesh
  ASSERT_EQ(16, mesh_flat.num_owned_cells());
  ASSERT_EQ(16, mesh_flat.num_entities(Portage::Entity_kind::CELL));
  ASSERT_EQ(25, mesh_flat.num_owned_nodes());
  ASSERT_EQ(25, mesh_flat.num_entities(Portage::Entity_kind::NODE));

  // Check coordinates
  // List coordinates of cell 0 - others are equal to this
  // with a shift
  std::vector<Portage::Point<2>> cell0Coords =
    {{0.0, 0.0},  {0.25, 0.0},  {0.25, 0.25},  {0.0, 0.25}};
  for (int c=0; c<16; ++c) {
    std::vector<Portage::Point<2>> pplist;
    mesh_flat.cell_get_coordinates(c, &pplist);
    ASSERT_EQ(4, pplist.size());
    double dx = (c / 4) * 0.25;
    double dy = (c % 4) * 0.25;
    for (int n=0; n<4; ++n) {
      ASSERT_EQ(cell0Coords[n][0] + dx, pplist[n][0]);
      ASSERT_EQ(cell0Coords[n][1] + dy, pplist[n][1]);
    }
  }

  // Test global ids
  std::vector<int>& gids = mesh_flat.get_global_cell_ids();
  for (int c=0; c<16; ++c)
    ASSERT_EQ(c, gids[c]);

  // Test centroids
  Portage::Point<2> centroid;
  mesh_flat.cell_centroid(0, &centroid);
  for (int d=0; d<2; ++d)
    ASSERT_EQ(0.125, centroid[d]);

  // Test cell_volume()
  for (int c=0; c<16; ++c)
    ASSERT_TRUE(abs(mesh_flat.cell_volume(c)-0.0625) < 1e-12);

  // Test cell neighbors
  for (int c=0; c<16; ++c) {
    // Get my neighbors
    std::vector<int> adjcells;
    mesh_flat.cell_get_node_adj_cells(c, Portage::Entity_type::ALL,
                                      &adjcells);
    // Add my own ID
    adjcells.push_back(c);
    // Which neighbors do I expect?
    // Loop through all possible neighbors, skipping if on boundary
    std::vector<int> expNeighbors;
    for (int dx=-1; dx<=+1; ++dx) {
      if (dx == -1 && c / 4 == 0) continue;
      if (dx == +1 && c / 4 == 3) continue;
      for (int dy=-1; dy<=+1; ++dy) {
        if (dy == -1 && c % 4 == 0) continue;
        if (dy == +1 && c % 4 == 3) continue;
        expNeighbors.push_back(c + dx * 4 + dy);
      }
    }
    int count = adjcells.size();
    ASSERT_EQ(expNeighbors.size(), count);
    std::sort(adjcells.begin(), adjcells.end());
    for (int n=0; n<count; ++n)
      ASSERT_EQ(expNeighbors[n], adjcells[n]);
  }

}


