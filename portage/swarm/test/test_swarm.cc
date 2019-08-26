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
#include "wonton/support/Point.h"


TEST(SwarmFactory, Check_1D_random) {
  std::shared_ptr<Portage::Meshfree::Swarm<1>> swarm=Portage::Meshfree::SwarmFactory(-4.,4.,17,0);
  ASSERT_EQ(swarm->num_owned_particles(), 17);
  for (int i=0; i<17; i++) {
    Wonton::Point<1> pt = swarm->get_particle_coordinates(i);
    ASSERT_LE(pt[0],  4.);
    ASSERT_GE(pt[0], -4.);
  }
}

TEST(SwarmFactory, Check_1D_random_seed) {
  std::shared_ptr<Portage::Meshfree::Swarm<1>> swarm=Portage::Meshfree::SwarmFactory(-4.,4.,17,0,1234);
  std::shared_ptr<Portage::Meshfree::Swarm<1>> swarm2=Portage::Meshfree::SwarmFactory(-4.,4.,17,0,1234);
  ASSERT_EQ(swarm->num_owned_particles(), 17);
  for (int i=0; i<17; i++) {
    Wonton::Point<1> pt = swarm->get_particle_coordinates(i);
    Wonton::Point<1> pt2 = swarm2->get_particle_coordinates(i);
    ASSERT_EQ(pt[0], pt2[0]);
  }
}

TEST(SwarmFactory, Check_2D_random) {
  std::shared_ptr<Portage::Meshfree::Swarm<2>> swarm=Portage::Meshfree::SwarmFactory(-4.,-4.,4.,4.,17*17,0);
  ASSERT_EQ(swarm->num_owned_particles(), 17*17);
  for (int i=0; i<17*17; i++) {
    Wonton::Point<2> pt = swarm->get_particle_coordinates(i);
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
    Wonton::Point<2> pt = swarm->get_particle_coordinates(i);
    Wonton::Point<2> pt2 = swarm2->get_particle_coordinates(i);
    ASSERT_EQ(pt[0],  pt2[0]);
    ASSERT_EQ(pt[1],  pt2[1]);
  }
}

TEST(SwarmFactory, Check_3D_random) {
  std::shared_ptr<Portage::Meshfree::Swarm<3>> swarm=Portage::Meshfree::SwarmFactory(-4.,-4.,-4.,4.,4.,4.,17*17*17,0);
  ASSERT_EQ(swarm->num_owned_particles(), 17*17*17);
  for (int i=0; i<17*17*17; i++) {
    Wonton::Point<3> pt = swarm->get_particle_coordinates(i);
    ASSERT_LE(pt[0],  4.);
    ASSERT_GE(pt[0], -4.);
    ASSERT_LE(pt[1],  4.);
    ASSERT_GE(pt[1], -4.);
    ASSERT_LE(pt[2],  4.);
    ASSERT_GE(pt[2], -4.);
  }
}

TEST(SwarmFactory, Check_3D_random_seed) {
  std::shared_ptr<Portage::Meshfree::Swarm<3>> swarm=Portage::Meshfree::SwarmFactory(-4.,-4.,-4.,4.,4.,4.,17*17*17,0,1234);
  std::shared_ptr<Portage::Meshfree::Swarm<3>> swarm2=Portage::Meshfree::SwarmFactory(-4.,-4.,-4.,4.,4.,4.,17*17*17,0,1234);
  ASSERT_EQ(swarm->num_owned_particles(), 17*17*17);
  for (int i=0; i<17; i++) {
    for (int j=0; j<17; j++) {
      for (int k=0; k<17; k++) {
	Wonton::Point<3> pt = swarm->get_particle_coordinates((i*17+j)*17+k);
	Wonton::Point<3> pt2 = swarm2->get_particle_coordinates((i*17+j)*17+k);
	ASSERT_EQ(pt[0],  pt2[0]);
	ASSERT_EQ(pt[1],  pt2[1]);
	ASSERT_EQ(pt[2],  pt2[2]);
      }
    }
  }
}

TEST(SwarmFactory, Check_1D_regular) {
  std::shared_ptr<Portage::Meshfree::Swarm<1>> swarm=Portage::Meshfree::SwarmFactory(-4.,4.,17,1);
  ASSERT_EQ(swarm->num_owned_particles(), 17);
  for (int i=0; i<17; i++) {
    Wonton::Point<1> pt = swarm->get_particle_coordinates(i);
    ASSERT_NEAR(pt[0], -4.+8./16.*i, 1.e-12);
  }
}

TEST(SwarmFactory, Check_2D_regular) {
  std::shared_ptr<Portage::Meshfree::Swarm<2>> swarm=Portage::Meshfree::SwarmFactory(-4.,-4.,4.,4.,17*17,1);
  ASSERT_EQ(swarm->num_owned_particles(), 17*17);
  for (int i=0; i<17; i++) {
    for (int j=0; j<17; j++) {
      Wonton::Point<2> pt = swarm->get_particle_coordinates(i*17+j);
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
	Wonton::Point<3> pt = swarm->get_particle_coordinates((i*17+j)*17+k);
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
    Wonton::Point<1> pt = swarm->get_particle_coordinates(i);
    ASSERT_LE(pt[0],  4.);
    ASSERT_GE(pt[0], -4.);
  }
}

TEST(SwarmFactory, Check_2D_perturbed) {
  std::shared_ptr<Portage::Meshfree::Swarm<2>> swarm=Portage::Meshfree::SwarmFactory(-4.,-4.,4.,4.,17*17,2);
  ASSERT_EQ(swarm->num_owned_particles(), 17*17);
  for (int i=0; i<17; i++) {
    for (int j=0; j<17; j++) {
      Wonton::Point<2> pt = swarm->get_particle_coordinates(i*17+j);
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
	Wonton::Point<3> pt = swarm->get_particle_coordinates((i*17+j)*17+k);
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

TEST(SwarmFactory, Check_3D_perturbed_seed) {
  std::shared_ptr<Portage::Meshfree::Swarm<3>> swarm=Portage::Meshfree::SwarmFactory(-4.,-4.,-4.,4.,4.,4.,17*17*17,2,1234);
  std::shared_ptr<Portage::Meshfree::Swarm<3>> swarm2=Portage::Meshfree::SwarmFactory(-4.,-4.,-4.,4.,4.,4.,17*17*17,2,1234);
  ASSERT_EQ(swarm->num_owned_particles(), 17*17*17);
  for (int i=0; i<17; i++) {
    for (int j=0; j<17; j++) {
      for (int k=0; k<17; k++) {
	Wonton::Point<3> pt = swarm->get_particle_coordinates((i*17+j)*17+k);
	Wonton::Point<3> pt2 = swarm2->get_particle_coordinates((i*17+j)*17+k);
        ASSERT_EQ(pt[0], pt2[0]);
        ASSERT_EQ(pt[1], pt2[1]);
        ASSERT_EQ(pt[2], pt2[2]);
      }
    }
  }
}

TEST(Swarm, Sanity_Check_3D) {
  Portage::vector<Wonton::Point<3>> points(10);

  srand(time(NULL));
  for (int i = 0; i < 10; i++)
    points[i] = Wonton::Point<3>(
        (static_cast<double>(rand()) / RAND_MAX),
        (static_cast<double>(rand()) / RAND_MAX),
        (static_cast<double>(rand()) / RAND_MAX));

  auto p_ptr = std::make_shared<Portage::vector<Wonton::Point<3>>>(points);

  Portage::Meshfree::Swarm<3> swarm(p_ptr);

  // Did we get back 10 points?
  ASSERT_EQ(10, swarm.num_owned_particles());

  // Check the bounding box and center of the hexahedral "cell"
  // surrounding each point

  for (int i = 0; i < 10; i++) {
    // Are the point coordinates correct?
    Wonton::Point<3> p1 = swarm.get_particle_coordinates(i);
    Wonton::Point<3> p2 = points[i];
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
    Wonton::Point<3> cent;
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
    Wonton::Point<3> node;
    mesh_wrapper.node_get_coordinates(ijk, &node);
    for (int i = 0; i < 3; i++) ASSERT_TRUE(pt[i] == node[i]);
  }
}


/*!
  @brief Unit test for factory with multiple 2D meshes
*/
TEST(Swarm, Multiple_2D) {
  Wonton::Simple_Mesh mesh0(0.0, 0.0, 1.0, 1.0, 4, 4);
  Wonton::Simple_Mesh_Wrapper wrapper0(mesh0);

  Wonton::Simple_Mesh mesh1(1.0, 0.0, 2.0, 1.0, 4, 4);
  Wonton::Simple_Mesh_Wrapper wrapper1(mesh1);

  Wonton::Simple_Mesh mesh2(1.0, 1.0, 2.0, 2.0, 4, 4);
  Wonton::Simple_Mesh_Wrapper wrapper2(mesh2);

  Wonton::Simple_Mesh mesh3(0.0, 1.0, 1.0, 2.0, 4, 4);
  Wonton::Simple_Mesh_Wrapper wrapper3(mesh3);

  std::vector<Wonton::Simple_Mesh_Wrapper*> wrappers={&wrapper0, &wrapper1, &wrapper2, &wrapper3};

  // create swarm from mesh wrapper cells
  std::shared_ptr<Portage::Meshfree::Swarm<2>> swarm_ptr =
    Portage::Meshfree::SwarmFactory<2,Wonton::Simple_Mesh_Wrapper>(wrappers, Portage::Entity_kind::NODE);
  Portage::Meshfree::Swarm<2> &swarm(*swarm_ptr);

  // test size
  ASSERT_EQ(100, swarm.num_particles());

  // test points
  for (size_t ijk = 0; ijk < 100; ijk++) {
    auto pt = swarm.get_particle_coordinates(ijk);
    Wonton::Point<2> node;
    if (ijk < 25) {
      wrappers[0]->node_get_coordinates(ijk   , &node);
      for (int i = 0; i < 2; i++) ASSERT_TRUE(pt[i] == node[i]);
    } else if (ijk < 50) {
      wrappers[1]->node_get_coordinates(ijk-25, &node);
      for (int i = 0; i < 2; i++) ASSERT_TRUE(pt[i] == node[i]);
    } else if (ijk < 75) {
      wrappers[2]->node_get_coordinates(ijk-50, &node);
      for (int i = 0; i < 2; i++) ASSERT_TRUE(pt[i] == node[i]);
    } else {
      wrappers[3]->node_get_coordinates(ijk-75, &node);
      for (int i = 0; i < 2; i++) ASSERT_TRUE(pt[i] == node[i]);
    }
  }
}


/*!
  @brief Unit test for factory with multiple 3D meshes
*/
TEST(Swarm, Multiple_3D) {
  Wonton::Simple_Mesh mesh0(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 4, 4);
  Wonton::Simple_Mesh mesh1(1.0, 0.0, 0.0, 2.0, 1.0, 1.0, 4, 4, 4);
  Wonton::Simple_Mesh mesh2(1.0, 1.0, 0.0, 2.0, 2.0, 1.0, 4, 4, 4);
  Wonton::Simple_Mesh mesh3(0.0, 1.0, 0.0, 1.0, 2.0, 1.0, 4, 4, 4);
  Wonton::Simple_Mesh mesh4(0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 4, 4, 4);
  Wonton::Simple_Mesh mesh5(1.0, 0.0, 1.0, 2.0, 1.0, 2.0, 4, 4, 4);
  Wonton::Simple_Mesh mesh6(1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 4, 4, 4);
  Wonton::Simple_Mesh mesh7(0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 4, 4, 4);
  Wonton::Simple_Mesh_Wrapper wrapper0(mesh0);
  Wonton::Simple_Mesh_Wrapper wrapper1(mesh0);
  Wonton::Simple_Mesh_Wrapper wrapper2(mesh0);
  Wonton::Simple_Mesh_Wrapper wrapper3(mesh0);
  Wonton::Simple_Mesh_Wrapper wrapper4(mesh0);
  Wonton::Simple_Mesh_Wrapper wrapper5(mesh0);
  Wonton::Simple_Mesh_Wrapper wrapper6(mesh0);
  Wonton::Simple_Mesh_Wrapper wrapper7(mesh0);

  std::vector<Wonton::Simple_Mesh_Wrapper*> wrappers={&wrapper0, &wrapper1, &wrapper2, &wrapper3, 
                                                      &wrapper4, &wrapper5, &wrapper6, &wrapper7};

  // create swarm from mesh wrapper cells
  std::shared_ptr<Portage::Meshfree::Swarm<3>> swarm_ptr =
    Portage::Meshfree::SwarmFactory<3,Wonton::Simple_Mesh_Wrapper>(wrappers, Portage::Entity_kind::NODE);
  Portage::Meshfree::Swarm<3> &swarm(*swarm_ptr);

  // test size
  ASSERT_EQ(1000, swarm.num_particles());

  // test points
  for (size_t ijk = 0; ijk < 1000; ijk++) {
    auto pt = swarm.get_particle_coordinates(ijk);
    Wonton::Point<3> node;
    if (ijk < 125) {
      wrappers[0]->node_get_coordinates(ijk    , &node);
      for (int i = 0; i < 3; i++) ASSERT_TRUE(pt[i] == node[i]);
    } else if (ijk < 250) {
      wrappers[1]->node_get_coordinates(ijk-125, &node);
      for (int i = 0; i < 3; i++) ASSERT_TRUE(pt[i] == node[i]);
    } else if (ijk < 375) {
      wrappers[2]->node_get_coordinates(ijk-250, &node);
      for (int i = 0; i < 3; i++) ASSERT_TRUE(pt[i] == node[i]);
    } else if (ijk < 500) {
      wrappers[2]->node_get_coordinates(ijk-375, &node);
      for (int i = 0; i < 3; i++) ASSERT_TRUE(pt[i] == node[i]);
    } else if (ijk < 625) {
      wrappers[2]->node_get_coordinates(ijk-500, &node);
      for (int i = 0; i < 3; i++) ASSERT_TRUE(pt[i] == node[i]);
    } else if (ijk < 750) {
      wrappers[2]->node_get_coordinates(ijk-625, &node);
      for (int i = 0; i < 3; i++) ASSERT_TRUE(pt[i] == node[i]);
    } else if (ijk < 875) {
      wrappers[2]->node_get_coordinates(ijk-750, &node);
      for (int i = 0; i < 3; i++) ASSERT_TRUE(pt[i] == node[i]);
    } else {
      wrappers[3]->node_get_coordinates(ijk-875, &node);
      for (int i = 0; i < 3; i++) ASSERT_TRUE(pt[i] == node[i]);
    }
  }
}
