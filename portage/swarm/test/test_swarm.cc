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

#include "gtest/gtest.h"

// portage includes
#include "portage/swarm/swarm.h"
#include "portage/support/portage.h"

// wonton includes
#include "wonton/mesh/simple/simple_mesh.h"
#include "wonton/mesh/simple/simple_mesh_wrapper.h"

TEST(SwarmFactory, Check_1D_random) {
  std::shared_ptr<Portage::Meshfree::Swarm<1>> swarm=Portage::Meshfree::SwarmFactory(-4.,4.,17,0);
  ASSERT_EQ(swarm->num_owned_particles(), 17);
  for (int i=0; i<17; i++) {
    Portage::Point<1> pt = swarm->get_particle_coordinates(i);
    ASSERT_LE(pt[0],  4.);
    ASSERT_GE(pt[0], -4.);
  }
}

TEST(SwarmFactory, Check_2D_random) {
  std::shared_ptr<Portage::Meshfree::Swarm<2>> swarm=Portage::Meshfree::SwarmFactory(-4.,-4.,4.,4.,17*17,0);
  ASSERT_EQ(swarm->num_owned_particles(), 17*17);
  for (int i=0; i<17*17; i++) {
    Portage::Point<2> pt = swarm->get_particle_coordinates(i);
    ASSERT_LE(pt[0],  4.);
    ASSERT_GE(pt[0], -4.);
    ASSERT_LE(pt[1],  4.);
    ASSERT_GE(pt[1], -4.);
  }
}

TEST(SwarmFactory, Check_2D_random_seed) {
  std::shared_ptr<Portage::Meshfree::Swarm<2>> swarm=Portage::Meshfree::SwarmFactory(-4.,-4.,4.,4.,17*17,0,1234);
  std::shared_ptr<Portage::Meshfree::Swarm<2>> swarm2=Portage::Meshfree::SwarmFactory(-4.,-4.,4.,4.,17*17,0,1234);
  ASSERT_EQ(swarm->num_owned_particles(), 17*17);
  for (int i=0; i<17*17; i++) {
    Portage::Point<2> pt = swarm->get_particle_coordinates(i);
    Portage::Point<2> pt2 = swarm->get_particle_coordinates(i);
    Portage::Point<2> pt3 = swarm->get_particle_coordinates(i);
    ASSERT_EQ(pt[0],  pt2[0]);
    ASSERT_EQ(pt[1],  pt2[1]);
  }
}

TEST(SwarmFactory, Check_3D_random) {
  std::shared_ptr<Portage::Meshfree::Swarm<3>> swarm=Portage::Meshfree::SwarmFactory(-4.,-4.,-4.,4.,4.,4.,17*17*17,0);
  ASSERT_EQ(swarm->num_owned_particles(), 17*17*17);
  for (int i=0; i<17*17*17; i++) {
    Portage::Point<3> pt = swarm->get_particle_coordinates(i);
    ASSERT_LE(pt[0],  4.);
    ASSERT_GE(pt[0], -4.);
    ASSERT_LE(pt[1],  4.);
    ASSERT_GE(pt[1], -4.);
    ASSERT_LE(pt[2],  4.);
    ASSERT_GE(pt[2], -4.);
  }
}

TEST(SwarmFactory, Check_1D_regular) {
  std::shared_ptr<Portage::Meshfree::Swarm<1>> swarm=Portage::Meshfree::SwarmFactory(-4.,4.,17,1);
  ASSERT_EQ(swarm->num_owned_particles(), 17);
  for (int i=0; i<17; i++) {
    Portage::Point<1> pt = swarm->get_particle_coordinates(i);
    ASSERT_NEAR(pt[0], -4.+8./16.*i, 1.e-12);
  }
}

TEST(SwarmFactory, Check_2D_regular) {
  std::shared_ptr<Portage::Meshfree::Swarm<2>> swarm=Portage::Meshfree::SwarmFactory(-4.,-4.,4.,4.,17*17,1);
  ASSERT_EQ(swarm->num_owned_particles(), 17*17);
  for (int i=0; i<17; i++) {
    for (int j=0; j<17; j++) {
      Portage::Point<2> pt = swarm->get_particle_coordinates(i*17+j);
      ASSERT_NEAR(pt[0], -4.+8./16.*i, 1.e-12);
      ASSERT_NEAR(pt[1], -4.+8./16.*j, 1.e-12);
    }
  }
}

TEST(SwarmFactory, Check_3D_regular) {
  std::shared_ptr<Portage::Meshfree::Swarm<3>> swarm=Portage::Meshfree::SwarmFactory(-4.,-4.,-4.,4.,4.,4.,17*17*17,1);
  ASSERT_EQ(swarm->num_owned_particles(), 17*17*17);
  for (int i=0; i<17; i++) {
    for (int j=0; j<17; j++) {
      for (int k=0; k<17; k++) {
	Portage::Point<3> pt = swarm->get_particle_coordinates((i*17+j)*17+k);
	ASSERT_NEAR(pt[0], -4.+8./16.*i, 1.e-12);
	ASSERT_NEAR(pt[1], -4.+8./16.*j, 1.e-12);
	ASSERT_NEAR(pt[2], -4.+8./16.*k, 1.e-12);
      }
    }
  }
}

TEST(SwarmFactory, Check_1D_perturbed) {
  std::shared_ptr<Portage::Meshfree::Swarm<1>> swarm=Portage::Meshfree::SwarmFactory(-4.,4.,17,2);
  ASSERT_EQ(swarm->num_owned_particles(), 17);
  for (int i=0; i<17; i++) {
    Portage::Point<1> pt = swarm->get_particle_coordinates(i);
    ASSERT_LE(pt[0],  4.);
    ASSERT_GE(pt[0], -4.);
  }
}

TEST(SwarmFactory, Check_2D_perturbed) {
  std::shared_ptr<Portage::Meshfree::Swarm<2>> swarm=Portage::Meshfree::SwarmFactory(-4.,-4.,4.,4.,17*17,2);
  ASSERT_EQ(swarm->num_owned_particles(), 17*17);
  for (int i=0; i<17; i++) {
    for (int j=0; j<17; j++) {
      Portage::Point<2> pt = swarm->get_particle_coordinates(i*17+j);
      ASSERT_LE(pt[0],  4.);
      ASSERT_GE(pt[0], -4.);
      ASSERT_LE(pt[1],  4.);
      ASSERT_GE(pt[1], -4.);
    }
  }
}

TEST(SwarmFactory, Check_3D_perturbed) {
  std::shared_ptr<Portage::Meshfree::Swarm<3>> swarm=Portage::Meshfree::SwarmFactory(-4.,-4.,-4.,4.,4.,4.,17*17*17,2);
  ASSERT_EQ(swarm->num_owned_particles(), 17*17*17);
  for (int i=0; i<17; i++) {
    for (int j=0; j<17; j++) {
      for (int k=0; k<17; k++) {
	Portage::Point<3> pt = swarm->get_particle_coordinates((i*17+j)*17+k);
	ASSERT_LE(pt[0],  4.);
	ASSERT_GE(pt[0], -4.);
	ASSERT_LE(pt[1],  4.);
	ASSERT_GE(pt[1], -4.);
	ASSERT_LE(pt[2],  4.);
	ASSERT_GE(pt[2], -4.);
      }
    }
  }
}

TEST(Swarm, Sanity_Check_3D) {
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
    Portage::Point<3> p1 = swarm.get_particle_coordinates(i);
    Portage::Point<3> p2 = points[i];
    for (size_t j=0; j < 3; j++) ASSERT_EQ(p1[j], p2[j]);
  }
}  // TEST


/*!
  @brief Unit test for constructor with Simple_Mesh_Wrapper in 3D using cells
*/
TEST(Swarm, Build_Simple_Mesh_Wrapper_Cell) {
  Wonton::Simple_Mesh mesh(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2);
  Wonton::Simple_Mesh_Wrapper mesh_wrapper(mesh);

  // create swarm from mesh wrapper cells
  std::shared_ptr<Portage::Meshfree::Swarm<3>> swarmc_ptr =
    Portage::Meshfree::SwarmFactory<3, Wonton::Simple_Mesh_Wrapper>(
        mesh_wrapper, Portage::Entity_kind::CELL);
  Portage::Meshfree::Swarm<3> &swarmc(*swarmc_ptr);

  // test size
  ASSERT_EQ(8, swarmc.num_particles());

  // test points
  for (size_t ijk = 0; ijk < 8; ijk++) {
    auto pt = swarmc.get_particle_coordinates(ijk);
    Portage::Point<3> cent;
    mesh_wrapper.cell_centroid<3>(ijk, &cent);
    for (int i = 0; i < 3; i++) ASSERT_TRUE(pt[i] == cent[i]);
  }
}


/*!
  @brief Unit test for constructor with Simple_Mesh_Wrapper in 3D using cells
*/
TEST(Swarm, Build_Simple_Mesh_Wrapper_Node) {
  Wonton::Simple_Mesh mesh(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2);
  Wonton::Simple_Mesh_Wrapper mesh_wrapper(mesh);

  // create swarm from mesh wrapper cells
  std::shared_ptr<Portage::Meshfree::Swarm<3>> swarmn_ptr =
    Portage::Meshfree::SwarmFactory<3,Wonton::Simple_Mesh_Wrapper>(mesh_wrapper, Portage::Entity_kind::NODE);
  Portage::Meshfree::Swarm<3> &swarmn(*swarmn_ptr);

  // test size
  ASSERT_EQ(27, swarmn.num_particles());

  // test points
  for (size_t ijk = 0; ijk < 27; ijk++) {
    auto pt = swarmn.get_particle_coordinates(ijk);
    Portage::Point<3> node;
    mesh_wrapper.node_get_coordinates(ijk, &node);
    for (int i = 0; i < 3; i++) ASSERT_TRUE(pt[i] == node[i]);
  }
}
