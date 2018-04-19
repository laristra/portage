/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/
/*
 * test_swarm.cc
 *
 *  Created on: Apr 19, 2017
 *      Author: gad
 *      Updated: rao
 */

#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <memory>

#include "portage/wonton/mesh/simple_mesh/simple_mesh_wrapper.h"
#include "portage/wonton/mesh/flat/flat_mesh_wrapper.h"

#include "portage/swarm/swarm.h"
#include "portage/support/Point.h"
#include "portage/support/portage.h"

#include "gtest/gtest.h"

TEST(Swarm, Sanity_Check) {
  Portage::vector<Portage::Point<3>> points(10);

  srand(time(NULL));
  for (int i = 0; i < 10; i++)
    points[i] = Portage::Point<3>(
        (static_cast<double>(rand()) / RAND_MAX),
        (static_cast<double>(rand()) / RAND_MAX),
        (static_cast<double>(rand()) / RAND_MAX));

  auto p_ptr = std::make_shared<Portage::vector<Portage::Point<3>>>(points);

  Portage::Meshfree::Swarm<3> swarm(p_ptr);

  // Did we get back 10 points?
  ASSERT_EQ(10, swarm.num_owned_particles());

  // Check the bounding box and center of the hexahedral "cell"
  // surrounding each point

  for (int i = 0; i < 10; i++) {
    // Are the point coordinates correct?
    auto pt = swarm.get_particle_coordinates(i);
    for (size_t j=0; j < 3; j++) ASSERT_EQ(pt[j], points[i][j]);
  }
}  // TEST


/*!
  @brief Unit test for constructor with Flat_Mesh_Wrapper in 3D using cells
*/
TEST(Swarm, Build_Flat_Mesh_Wrapper_Cell) {
  Portage::Simple_Mesh mesh(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2);
  Wonton::Simple_Mesh_Wrapper mesh_wrapper(mesh);
  Wonton::Flat_Mesh_Wrapper<double> mesh_flat;
  mesh_flat.initialize(mesh_wrapper);

  // create swarm from mesh wrapper cells
  Portage::Meshfree::Swarm<3> swarmc(mesh_flat, Portage::CELL);

  // test size
  ASSERT_EQ(8, swarmc.num_particles());

  // test points
  for (size_t ijk = 0; ijk < 8; ijk++) {
    auto pt = swarmc.get_particle_coordinates(ijk);
    Portage::Point<3> cent;
    mesh_flat.cell_centroid<3>(ijk, &cent);
    for (int i = 0; i < 3; i++) ASSERT_TRUE(pt[i] == cent[i]);
  }
}


/*!
  @brief Unit test for constructor with Flat_Mesh_Wrapper in 3D using cells
*/
TEST(Swarm, Build_Flat_Mesh_Wrapper_Node) {
  Portage::Simple_Mesh mesh(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2);
  Wonton::Simple_Mesh_Wrapper mesh_wrapper(mesh);
  Wonton::Flat_Mesh_Wrapper<double> mesh_flat;
  mesh_flat.initialize(mesh_wrapper);

  // create swarm from mesh wrapper cells
  Portage::Meshfree::Swarm<3> swarmn(mesh_flat, Portage::NODE);

  // test size
  ASSERT_EQ(27, swarmn.num_particles());

  // test points
  for (size_t ijk = 0; ijk < 27; ijk++) {
    auto pt = swarmn.get_particle_coordinates(ijk);
    Portage::Point<3> node;
    mesh_flat.node_get_coordinates(ijk, &node);
    for (int i = 0; i < 3; i++) ASSERT_TRUE(pt[i] == node[i]);
  }
}
