/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

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
  @brief Unit test for basic routines in Flat_Mesh_Wrapper
 */
TEST(Flat_Mesh_Wrapper, basic_routines) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::EDGE,
                                  Jali::Entity_kind::FACE,
                                  Jali::Entity_kind::WEDGE,
                                  Jali::Entity_kind::CORNER});
  std::shared_ptr<Jali::Mesh> mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2);
  ASSERT_TRUE(mesh != NULL);
  Portage::Jali_Mesh_Wrapper mesh_wrapper(*mesh);
  Portage::Flat_Mesh_Wrapper<> mesh_flat(8, mesh_wrapper);

  ASSERT_EQ(3, mesh_flat.space_dimension());

  // Test cells of flat mesh
  ASSERT_EQ(8, mesh_flat.num_owned_cells());
  ASSERT_EQ(8, mesh_flat.num_entities(Portage::Entity_kind::CELL));
  std::vector<Portage::Point<3>> pplist;
  mesh_flat.cell_get_coordinates(0, &pplist);
  Portage::Point<3> minpoint = *std::min_element(pplist.begin(), pplist.end());
  Portage::Point<3> maxpoint = *std::max_element(pplist.begin(), pplist.end());
  for (int d=0; d<2; ++d)
  {
    // arithmetic should be exact here, so use normal equality
    ASSERT_EQ(0.0, minpoint[d]);
    ASSERT_EQ(0.5, maxpoint[d]);
  }

  // Test decompose_cell_into_tets()
  std::vector<std::array<Portage::Point<3>, 4>> tcoords;
  mesh_flat.decompose_cell_into_tets(0, &tcoords, true);
  ASSERT_EQ(tcoords.size(), 12);

  // Test centroids
  Portage::Point<3> centroid;
  mesh_flat.cell_centroid(0, &centroid);
  for (int d=0; d<2; ++d)
    ASSERT_EQ(0.25, centroid[d]);

  // Test cell_volume()
  for (int c=0; c<7; ++c)
    ASSERT_TRUE(abs(mesh_flat.cell_volume(c)-0.125) < 1e-12);

  // Test cell neighbors
  std::vector<int> adjcells;
  mesh_flat.cell_get_node_adj_cells(0, Portage::Entity_type::ALL,
                                    &adjcells);
  ASSERT_EQ(7, adjcells.size());
  std::sort(adjcells.begin(), adjcells.end());
  for (int c=0; c<7; ++c)
    ASSERT_EQ(c+1, adjcells[c]);

}
