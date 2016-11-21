/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"

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
  Portage::Flat_Mesh_Wrapper<> mesh_flat;
  mesh_flat.initialize(8, mesh_wrapper);
  std::vector<std::array<Portage::Point<3>, 4>> tcoords;

  // Test decompose_cell_into_tets()
  mesh_flat.decompose_cell_into_tets(0, &tcoords, true);
  EXPECT_EQ(tcoords.size(), 12);

  // Test cell_volume()
  ASSERT_TRUE(abs(mesh_flat.cell_volume(0)-0.125) < 1e-12);
}
