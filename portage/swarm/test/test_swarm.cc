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

using namespace Portage::Meshfree;

class SwarmTest : public testing::Test {
protected:
  // particles per axis count and bounds, spatial step.
  int const n = 17;
  double const p_min = -4.0;
  double const p_max =  4.0;
  double const step = (p_max - p_min) / (n - 1);
  double const epsilon = 1.E-12;

  // distribution for generated particle fields
  unsigned const seed = 1234;
  unsigned const no_seed = 0;
  int const random = 0;
  int const regular = 1;
  int const perturbed = 2;
};

TEST_F(SwarmTest, Check_1D_random) {

  Swarm<1> swarm(n, random, no_seed, p_min, p_max);

  for (int i = 0; i < n; i++) {
    auto p = swarm.get_particle_coordinates(i);
    ASSERT_LE(p[0], p_max);
    ASSERT_GE(p[0], p_min);
  }
}

TEST_F(SwarmTest, Check_1D_random_seed) {

  Swarm<1> orig(n, random, seed, p_min, p_max);
  Swarm<1> twin(n, random, seed, p_min, p_max);

  for (int i = 0; i < n; i++) {
    auto p = orig.get_particle_coordinates(i);
    auto q = twin.get_particle_coordinates(i);
    ASSERT_DOUBLE_EQ(p[0], q[0]);
  }
}

TEST_F(SwarmTest, Check_2D_random) {

  Swarm<2> swarm(n * n, random, no_seed, p_min, p_max, p_min, p_max);

  for (int i = 0; i < n * n; i++) {
    auto p = swarm.get_particle_coordinates(i);
    ASSERT_LE(p[0], p_max);
    ASSERT_GE(p[0], p_min);
    ASSERT_LE(p[1], p_max);
    ASSERT_GE(p[1], p_min);
  }
}

TEST_F(SwarmTest, Check_2D_random_seed) {

  Swarm<2> orig(n * n, random, seed, p_min, p_max, p_min, p_max);
  Swarm<2> twin(n * n, random, seed, p_min, p_max, p_min, p_max);

  for (int i = 0; i < n * n; i++) {
    auto p = orig.get_particle_coordinates(i);
    auto q = twin.get_particle_coordinates(i);
    ASSERT_DOUBLE_EQ(p[0], q[0]);
    ASSERT_DOUBLE_EQ(p[1], q[1]);
  }
}

TEST_F(SwarmTest, Check_3D_random) {

  Swarm<3> swarm(n * n * n, random, no_seed, p_min, p_max, p_min, p_max, p_max, p_max);

  for (int i = 0; i < n * n * n; i++) {
    auto p = swarm.get_particle_coordinates(i);
    ASSERT_LE(p[0], p_max);
    ASSERT_GE(p[0], p_min);
    ASSERT_LE(p[1], p_max);
    ASSERT_GE(p[1], p_min);
    ASSERT_LE(p[2], p_max);
    ASSERT_GE(p[2], p_min);
  }
}

TEST_F(SwarmTest, Check_3D_random_seed) {

  Swarm<3> orig(n * n * n, random, seed, p_min, p_max, p_min, p_max, p_max, p_max);
  Swarm<3> twin(n * n * n, random, seed, p_min, p_max, p_min, p_max, p_max, p_max);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      for (int k = 0; k < n; k++) {
        int const index = (i * n + j) * n + k;
        auto p = orig.get_particle_coordinates(index);
        auto q = twin.get_particle_coordinates(index);
        ASSERT_DOUBLE_EQ(p[0], q[0]);
        ASSERT_DOUBLE_EQ(p[1], q[1]);
        ASSERT_DOUBLE_EQ(p[2], q[2]);
      }
    }
  }
}

TEST_F(SwarmTest, Check_1D_regular) {

  Swarm<1> swarm(n, regular, no_seed, p_min, p_max);

  for (int i = 0; i < n; i++) {
    auto p = swarm.get_particle_coordinates(i);
    ASSERT_NEAR(p[0], p_min + i * step, epsilon);
  }
}

TEST_F(SwarmTest, Check_2D_regular) {

  Swarm<2> swarm(n * n, regular, no_seed, p_min, p_max, p_min, p_max);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      auto p = swarm.get_particle_coordinates(i * n + j);
      ASSERT_NEAR(p[0], p_min + i * step, epsilon);
      ASSERT_NEAR(p[1], p_min + j * step, epsilon);
    }
  }
}

TEST_F(SwarmTest, Check_3D_regular) {

  Swarm<3> swarm(n * n * n, regular, no_seed, p_min, p_max, p_min, p_max, p_min, p_max);

 // std::cout << "step: "<< step << std::endl;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      for (int k = 0; k < n; k++) {
        auto p = swarm.get_particle_coordinates((i * n + j) * n + k);
        ASSERT_NEAR(p[0], p_min + i * step, epsilon);
        ASSERT_NEAR(p[1], p_min + j * step, epsilon);
        ASSERT_NEAR(p[2], p_min + k * step, epsilon);
      }
    }
  }
}

TEST_F(SwarmTest, Check_1D_perturbed) {

  Swarm<1> swarm(n, perturbed, no_seed, p_min, p_max);

  for (int i = 0; i < n; i++) {
    auto p = swarm.get_particle_coordinates(i);
    ASSERT_LE(p[0],  4.);
    ASSERT_GE(p[0], -4.);
  }
}

TEST_F(SwarmTest, Check_2D_perturbed) {

  Swarm<3> swarm(n * n, perturbed, no_seed, p_min, p_max, p_min, p_max);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      auto p = swarm.get_particle_coordinates(i * n + j);
      ASSERT_LE(p[0], p_max);
      ASSERT_GE(p[0], p_min);
      ASSERT_LE(p[1], p_max);
      ASSERT_GE(p[1], p_min);
    }
  }
}

TEST_F(SwarmTest, Check_3D_perturbed) {

  Swarm<3> swarm(n * n * n, perturbed, no_seed, p_min, p_max, p_min, p_max, p_min, p_max);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      for (int k = 0; k < n; k++) {
        auto p = swarm.get_particle_coordinates((i * n + j) * n + k);
        ASSERT_LE(p[0], p_max);
        ASSERT_GE(p[0], p_min);
        ASSERT_LE(p[1], p_max);
        ASSERT_GE(p[1], p_min);
        ASSERT_LE(p[2], p_max);
        ASSERT_GE(p[2], p_min);
      }
    }
  }
}

TEST_F(SwarmTest, Check_3D_perturbed_seed) {

  Swarm<3> orig(n * n * n, perturbed, seed, p_min, p_max, p_min, p_max, p_min, p_max);
  Swarm<3> twin(n * n * n, perturbed, seed, p_min, p_max, p_min, p_max, p_min, p_max);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      for (int k = 0; k < n; k++) {
        int const index = (i * n + j) * n + k;
        auto p = orig.get_particle_coordinates(index);
        auto q = twin.get_particle_coordinates(index);
        ASSERT_DOUBLE_EQ(p[0], q[0]);
        ASSERT_DOUBLE_EQ(p[1], q[1]);
        ASSERT_DOUBLE_EQ(p[2], q[2]);
      }
    }
  }
}

TEST_F(SwarmTest, Sanity_Check_3D) {

  // generate a random particle field
  std::random_device device;
  std::mt19937 engine { device() };
  std::uniform_real_distribution<double> generator(p_min, p_max);

  Wonton::vector<Wonton::Point<3>> points(n);

  for (auto&& current : points) {
    // point coordinates are not always initialized in order
    // so enforce random number picking sequence
    double const noise[] = { generator(engine),
                             generator(engine),
                             generator(engine) };

    current = Wonton::Point<3>(noise[0], noise[1], noise[2]);
  }

  // create particle field from point list
  Swarm<3> swarm(points);

  ASSERT_EQ(swarm.num_owned_particles(), n);

  for (int i = 0; i < n; i++) {
    Wonton::Point<3> p = swarm.get_particle_coordinates(i);
    Wonton::Point<3> q = points[i];
    for (int j = 0; j < 3; j++)
      ASSERT_DOUBLE_EQ(p[j], q[j]);
  }
}

TEST_F(SwarmTest, Build_Simple_Mesh_Wrapper_Cell) {

  // generate a cartesian grid using simple mesh
  Wonton::Simple_Mesh grid(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2);
  Wonton::Simple_Mesh_Wrapper mesh(grid);

  // create particle field from generated mesh
  Swarm<3> swarm(mesh, Wonton::CELL);

  ASSERT_EQ(swarm.num_particles(), 8);

  for (int i = 0; i < 8; ++i) {
    auto p = swarm.get_particle_coordinates(i);
    Wonton::Point<3> centroid;
    mesh.cell_centroid(i, &centroid);
    for (int j = 0; j < 3; j++)
      ASSERT_DOUBLE_EQ(p[j], centroid[j]);
  }
}

TEST_F(SwarmTest, Build_Simple_Mesh_Wrapper_Node) {

  // generate a cartesian grid using simple mesh
  Wonton::Simple_Mesh grid(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2);
  Wonton::Simple_Mesh_Wrapper mesh(grid);

  // create particle field from generated mesh
  Swarm<3> swarm(mesh, Wonton::NODE);

  ASSERT_EQ(swarm.num_particles(), 27);

  for (int i = 0; i < 27; ++i) {
    auto p = swarm.get_particle_coordinates(i);
    Wonton::Point<3> q;
    mesh.node_get_coordinates(i, &q);
    for (int j = 0; j < 3; j++)
      ASSERT_DOUBLE_EQ(p[j], q[j]);
  }
}

TEST_F(SwarmTest, Multiple_2D) {

  using Mesh = Wonton::Simple_Mesh;
  using Wrapper = Wonton::Simple_Mesh_Wrapper;

  // generate a set of meshes
  Mesh mesh0(0.0, 0.0, 1.0, 1.0, 4, 4);
  Mesh mesh1(1.0, 0.0, 2.0, 1.0, 4, 4);
  Mesh mesh2(1.0, 1.0, 2.0, 2.0, 4, 4);
  Mesh mesh3(0.0, 1.0, 1.0, 2.0, 4, 4);

  Wrapper wrapper0(mesh0);
  Wrapper wrapper1(mesh1);
  Wrapper wrapper2(mesh2);
  Wrapper wrapper3(mesh3);

  std::vector<Wrapper*> wrappers = { &wrapper0, &wrapper1, &wrapper2, &wrapper3 };

  // create particle field from list of meshes
  Swarm<2> swarm(wrappers, Wonton::NODE);

  ASSERT_EQ(swarm.num_particles(), 100);

  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 25; ++j) {
      auto p = swarm.get_particle_coordinates(i * 25 + j);
      Wonton::Point<2> q;
      wrappers[i]->node_get_coordinates(j, &q);
      for (int k = 0; k < 2; ++k)
        ASSERT_DOUBLE_EQ(p[k], q[k]);
    }
  }
}


TEST_F(SwarmTest, Multiple_3D) {

  using Mesh = Wonton::Simple_Mesh;
  using Wrapper = Wonton::Simple_Mesh_Wrapper;

  Mesh mesh0(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 4, 4);
  Mesh mesh1(1.0, 0.0, 0.0, 2.0, 1.0, 1.0, 4, 4, 4);
  Mesh mesh2(1.0, 1.0, 0.0, 2.0, 2.0, 1.0, 4, 4, 4);
  Mesh mesh3(0.0, 1.0, 0.0, 1.0, 2.0, 1.0, 4, 4, 4);
  Mesh mesh4(0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 4, 4, 4);
  Mesh mesh5(1.0, 0.0, 1.0, 2.0, 1.0, 2.0, 4, 4, 4);
  Mesh mesh6(1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 4, 4, 4);
  Mesh mesh7(0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 4, 4, 4);

  // create a set of meshes
  Wrapper wrapper0(mesh0);
  Wrapper wrapper1(mesh1);
  Wrapper wrapper2(mesh2);
  Wrapper wrapper3(mesh3);
  Wrapper wrapper4(mesh4);
  Wrapper wrapper5(mesh5);
  Wrapper wrapper6(mesh6);
  Wrapper wrapper7(mesh7);

  std::vector<Wrapper*> wrappers = { &wrapper0, &wrapper1, &wrapper2, &wrapper3,
                                     &wrapper4, &wrapper5, &wrapper6, &wrapper7 };

  // create particle field from list of meshes
  Swarm<3> swarm(wrappers, Wonton::NODE);

  ASSERT_EQ(swarm.num_particles(), 1000);

  for (int i = 0; i < 8; ++i) {
    for (int j = 0; j < 125; ++j) {
      auto p = swarm.get_particle_coordinates(i * 125 + j);
      Wonton::Point<3> q;
      wrappers[i]->node_get_coordinates(j, &q);
      for (int k = 0; k < 3; ++k)
        ASSERT_DOUBLE_EQ(p[k], q[k]);
    }
  }
}
