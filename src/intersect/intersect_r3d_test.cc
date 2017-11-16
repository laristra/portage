/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/
#include "intersect_r3d.h"
#include "gtest/gtest.h"
#include "MeshFactory.hh"
#include "portage/wonton/mesh/jali/jali_mesh_wrapper.h"
#include <math.h>

TEST(intersectR3D, simple1) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::EDGE, Jali::Entity_kind::FACE,
                        Jali::Entity_kind::WEDGE, Jali::Entity_kind::CORNER});
  std::shared_ptr<Jali::Mesh> sm = mf(0, 0, 0, 2, 2, 2, 1, 1, 1);
  std::shared_ptr<Jali::Mesh> tm = mf(1, 1, 1, 2, 2, 2, 1, 1, 1);
  const Portage::Jali_Mesh_Wrapper s(*sm);
  const Portage::Jali_Mesh_Wrapper t(*tm);

  const double eps = 1e-12;
  const Portage::IntersectR3D<Portage::Jali_Mesh_Wrapper> isect{s, t};
  const std::vector<std::vector<double>> moments = isect(0, 0);
  for (int i = 0; i < moments.size(); i++) {
    for (int j = 0; j < moments[i].size(); j++) {
      std::cout << "i, j, m " << i << ", " << j << ", " << moments[i][j]
                << std::endl;
    }
  }

  ASSERT_TRUE(moments.size() == 1);
  ASSERT_TRUE(moments[0].size() == 4);

  ASSERT_TRUE(std::abs(moments[0][0] - 1) < eps);
  ASSERT_TRUE(std::abs(moments[0][1] / moments[0][0] - 1.5) < eps);
  ASSERT_TRUE(std::abs(moments[0][2] / moments[0][0] - 1.5) < eps);
  ASSERT_TRUE(std::abs(moments[0][3] / moments[0][0] - 1.5) < eps);
}

TEST(intersectR3D, simple2) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::EDGE, Jali::Entity_kind::FACE,
                        Jali::Entity_kind::WEDGE, Jali::Entity_kind::CORNER});
  std::shared_ptr<Jali::Mesh> sm = mf(0, 0, 0, 2, 2, 2, 1, 1, 1);
  std::shared_ptr<Jali::Mesh> tm = mf(0, 0, 0, 2, 2, 2, 1, 1, 1);
  const Portage::Jali_Mesh_Wrapper s(*sm);
  const Portage::Jali_Mesh_Wrapper t(*tm);

  const double eps = 1e-12;
  const Portage::IntersectR3D<Portage::Jali_Mesh_Wrapper> isect{s, t};
  const std::vector<std::vector<double>> moments = isect(0, 0);
  for (int i = 0; i < moments.size(); i++) {
    for (int j = 0; j < moments[i].size(); j++) {
      std::cout << "i, j, m " << i << ", " << j << ", " << moments[i][j]
                << std::endl;
    }
  }

  ASSERT_TRUE(moments.size() == 1);
  ASSERT_TRUE(moments[0].size() == 4);

  ASSERT_TRUE(std::abs(moments[0][0] - 8) < eps);
  ASSERT_TRUE(std::abs(moments[0][1] / moments[0][0] - 1) < eps);
  ASSERT_TRUE(std::abs(moments[0][2] / moments[0][0] - 1) < eps);
  ASSERT_TRUE(std::abs(moments[0][3] / moments[0][0] - 1) < eps);
}

TEST(intersectR3D, simple3) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::EDGE, Jali::Entity_kind::FACE,
                        Jali::Entity_kind::WEDGE, Jali::Entity_kind::CORNER});
  std::shared_ptr<Jali::Mesh> sm = mf(0, 0, 0, 3, 3, 3, 1, 1, 1);
  std::shared_ptr<Jali::Mesh> tm = mf(1, 1, 1, 2, 2, 2, 1, 1, 1);
  const Portage::Jali_Mesh_Wrapper s(*sm);
  const Portage::Jali_Mesh_Wrapper t(*tm);

  const double eps = 1e-12;
  const Portage::IntersectR3D<Portage::Jali_Mesh_Wrapper> isect{s, t};
  const std::vector<std::vector<double>> moments = isect(0, 0);
  for (int i = 0; i < moments.size(); i++) {
    for (int j = 0; j < moments[i].size(); j++) {
      std::cout << "i, j, m " << i << ", " << j << ", " << moments[i][j]
                << std::endl;
    }
  }

  ASSERT_TRUE(moments.size() == 1);
  ASSERT_TRUE(moments[0].size() == 4);

  ASSERT_TRUE(std::abs(moments[0][0] - 1) < eps);
  ASSERT_TRUE(std::abs(moments[0][1] / moments[0][0] - 1.5) < eps);
  ASSERT_TRUE(std::abs(moments[0][2] / moments[0][0] - 1.5) < eps);
  ASSERT_TRUE(std::abs(moments[0][3] / moments[0][0] - 1.5) < eps);
}

TEST(intersectR3D, simple4) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::EDGE, Jali::Entity_kind::FACE,
                        Jali::Entity_kind::WEDGE, Jali::Entity_kind::CORNER});
  std::shared_ptr<Jali::Mesh> sm = mf(0, 0, 0, 2, 2, 2, 1, 1, 1);
  std::shared_ptr<Jali::Mesh> tm = mf(1, 1, 1, 3, 3, 3, 1, 1, 1);
  const Portage::Jali_Mesh_Wrapper s(*sm);
  const Portage::Jali_Mesh_Wrapper t(*tm);

  const double eps = 1e-12;
  const Portage::IntersectR3D<Portage::Jali_Mesh_Wrapper> isect{s, t};
  const std::vector<std::vector<double>> moments = isect(0, 0);
  for (int i = 0; i < moments.size(); i++) {
    for (int j = 0; j < moments[i].size(); j++) {
      std::cout << "i, j, m " << i << ", " << j << ", " << moments[i][j]
                << std::endl;
    }
  }

  ASSERT_TRUE(moments.size() == 1);
  ASSERT_TRUE(moments[0].size() == 4);

  ASSERT_TRUE(std::abs(moments[0][0] - 1) < eps);
  ASSERT_TRUE(std::abs(moments[0][1] / moments[0][0] - 1.5) < eps);
  ASSERT_TRUE(std::abs(moments[0][2] / moments[0][0] - 1.5) < eps);
  ASSERT_TRUE(std::abs(moments[0][3] / moments[0][0] - 1.5) < eps);
}

TEST(intersectR3D, simple5) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::EDGE, Jali::Entity_kind::FACE,
                        Jali::Entity_kind::WEDGE, Jali::Entity_kind::CORNER});
  std::shared_ptr<Jali::Mesh> sm = mf(0, 0, 0, 10, 10, 10, 1, 1, 1);
  std::shared_ptr<Jali::Mesh> tm = mf(-5, -5, -5, 5, 5, 5, 1, 1, 1);
  const Portage::Jali_Mesh_Wrapper s(*sm);
  const Portage::Jali_Mesh_Wrapper t(*tm);

  const double eps = 1e-12;
  const Portage::IntersectR3D<Portage::Jali_Mesh_Wrapper> isect{s, t};
  const std::vector<std::vector<double>> moments = isect(0, 0);
  for (int i = 0; i < moments.size(); i++) {
    for (int j = 0; j < moments[i].size(); j++) {
      std::cout << "i, j, m " << i << ", " << j << ", " << moments[i][j]
                << std::endl;
    }
  }

  ASSERT_TRUE(moments.size() == 1);
  ASSERT_TRUE(moments[0].size() == 4);

  ASSERT_TRUE(std::abs(moments[0][0] - 125) < eps);
  ASSERT_TRUE(std::abs(moments[0][1] / moments[0][0] - 2.5) < eps);
  ASSERT_TRUE(std::abs(moments[0][2] / moments[0][0] - 2.5) < eps);
  ASSERT_TRUE(std::abs(moments[0][3] / moments[0][0] - 2.5) < eps);
}

TEST(intersectR3D, simple6) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::EDGE, Jali::Entity_kind::FACE,
                        Jali::Entity_kind::WEDGE, Jali::Entity_kind::CORNER});
  std::shared_ptr<Jali::Mesh> sm = mf(0, 0, 0, 10, 10, 10, 5, 5, 5);
  std::shared_ptr<Jali::Mesh> tm = mf(0, 0, 0, 10, 10, 10, 2, 2, 2);
  const Portage::Jali_Mesh_Wrapper s(*sm);
  const Portage::Jali_Mesh_Wrapper t(*tm);

  const double eps = 1e-12;
  const Portage::IntersectR3D<Portage::Jali_Mesh_Wrapper> isect{s, t};
  std::vector<std::vector<double>> moments = isect(0, 0);
  for (int i = 0; i < moments.size(); i++) {
    for (int j = 0; j < moments[i].size(); j++) {
      std::cout << "i, j, m " << i << ", " << j << ", " << moments[i][j]
                << std::endl;
    }
  }

  ASSERT_TRUE(moments.size() == 1);
  ASSERT_TRUE(moments[0].size() == 4);

  ASSERT_TRUE(std::abs(moments[0][0] - 8) < eps);
  ASSERT_TRUE(std::abs(moments[0][1] / moments[0][0] - 1) < eps);
  ASSERT_TRUE(std::abs(moments[0][2] / moments[0][0] - 1) < eps);
  ASSERT_TRUE(std::abs(moments[0][3] / moments[0][0] - 1) < eps);

  moments = isect(1, 0);
  for (int i = 0; i < moments.size(); i++) {
    for (int j = 0; j < moments[i].size(); j++) {
      std::cout << "i, j, m " << i << ", " << j << ", " << moments[i][j]
                << std::endl;
    }
  }

  ASSERT_TRUE(moments.size() == 1);
  ASSERT_TRUE(moments[0].size() == 4);

  ASSERT_TRUE(std::abs(moments[0][0] - 8) < eps);
  ASSERT_TRUE(std::abs(moments[0][1] / moments[0][0] - 1) < eps);
  ASSERT_TRUE(std::abs(moments[0][2] / moments[0][0] - 1) < eps);
  ASSERT_TRUE(std::abs(moments[0][3] / moments[0][0] - 3) < eps);

  moments = isect(2, 0);
  for (int i = 0; i < moments.size(); i++) {
    for (int j = 0; j < moments[i].size(); j++) {
      std::cout << "i, j, m " << i << ", " << j << ", " << moments[i][j]
                << std::endl;
    }
  }

  ASSERT_TRUE(moments.size() == 1);
  ASSERT_TRUE(moments[0].size() == 4);

  ASSERT_TRUE(std::abs(moments[0][0] - 4) < eps);
  ASSERT_TRUE(std::abs(moments[0][1] / moments[0][0] - 1) < eps);
  ASSERT_TRUE(std::abs(moments[0][2] / moments[0][0] - 1) < eps);
  ASSERT_TRUE(std::abs(moments[0][3] / moments[0][0] - 4.5) < eps);
}

TEST(intersectR3D, simple7) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::EDGE, Jali::Entity_kind::FACE,
                        Jali::Entity_kind::WEDGE, Jali::Entity_kind::CORNER});
  std::shared_ptr<Jali::Mesh> sm = mf(-2, -2, -2, 0, 0, 0, 1, 1, 1);
  std::shared_ptr<Jali::Mesh> tm = mf(-1, -1, -1, 0, 0, 0, 1, 1, 1);
  const Portage::Jali_Mesh_Wrapper s(*sm);
  const Portage::Jali_Mesh_Wrapper t(*tm);

  const double eps = 1e-12;
  const Portage::IntersectR3D<Portage::Jali_Mesh_Wrapper> isect{s, t};
  const std::vector<std::vector<double>> moments = isect(0, 0);
  for (int i = 0; i < moments.size(); i++) {
    for (int j = 0; j < moments[i].size(); j++) {
      std::cout << "i, j, m " << i << ", " << j << ", " << moments[i][j]
                << std::endl;
    }
  }

  ASSERT_TRUE(moments.size() == 1);
  ASSERT_TRUE(moments[0].size() == 4);

  ASSERT_TRUE(std::abs(moments[0][0] - 1) < eps);
  ASSERT_TRUE(std::abs(moments[0][1] / moments[0][0] - (-0.5)) < eps);
  ASSERT_TRUE(std::abs(moments[0][2] / moments[0][0] - (-0.5)) < eps);
  ASSERT_TRUE(std::abs(moments[0][3] / moments[0][0] - (-0.5)) < eps);
}

TEST(intersectR3D, simple8) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::EDGE, Jali::Entity_kind::FACE,
                        Jali::Entity_kind::WEDGE, Jali::Entity_kind::CORNER});
  std::shared_ptr<Jali::Mesh> sm = mf(-4, -4, -4, 0, 0, 0, 1, 1, 1);
  std::shared_ptr<Jali::Mesh> tm = mf(-3, -3, -3, 0, 0, 0, 1, 1, 1);
  const Portage::Jali_Mesh_Wrapper s(*sm);
  const Portage::Jali_Mesh_Wrapper t(*tm);

  const double eps = 1e-12;
  const Portage::IntersectR3D<Portage::Jali_Mesh_Wrapper> isect{s, t};
  const std::vector<std::vector<double>> moments = isect(0, 0);
  for (int i = 0; i < moments.size(); i++) {
    for (int j = 0; j < moments[i].size(); j++) {
      std::cout << "i, j, m " << i << ", " << j << ", " << moments[i][j]
                << std::endl;
    }
  }

  ASSERT_TRUE(moments.size() == 1);
  ASSERT_TRUE(moments[0].size() == 4);

  ASSERT_TRUE(std::abs(moments[0][0] - 27) < eps);
  ASSERT_TRUE(std::abs(moments[0][1] / moments[0][0] - (-1.5)) < eps);
  ASSERT_TRUE(std::abs(moments[0][2] / moments[0][0] - (-1.5)) < eps);
  ASSERT_TRUE(std::abs(moments[0][3] / moments[0][0] - (-1.5)) < eps);
}

TEST(intersectR3D, simple9) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::EDGE, Jali::Entity_kind::FACE,
                        Jali::Entity_kind::WEDGE, Jali::Entity_kind::CORNER});
  std::shared_ptr<Jali::Mesh> sm = mf(-4, -3, -2, 0, 1, 2, 1, 1, 1);
  std::shared_ptr<Jali::Mesh> tm = mf(-3, -2, -1, 0, 1, 2, 1, 1, 1);
  const Portage::Jali_Mesh_Wrapper s(*sm);
  const Portage::Jali_Mesh_Wrapper t(*tm);

  const double eps = 1e-12;
  const Portage::IntersectR3D<Portage::Jali_Mesh_Wrapper> isect{s, t};
  const std::vector<std::vector<double>> moments = isect(0, 0);
  for (int i = 0; i < moments.size(); i++) {
    for (int j = 0; j < moments[i].size(); j++) {
      std::cout << "i, j, m " << i << ", " << j << ", " << moments[i][j]
                << std::endl;
    }
  }

  ASSERT_TRUE(moments.size() == 1);
  ASSERT_TRUE(moments[0].size() == 4);

  ASSERT_TRUE(std::abs(moments[0][0] - 27) < eps);
  ASSERT_TRUE(std::abs(moments[0][1] / moments[0][0] - (-1.5)) < eps);
  ASSERT_TRUE(std::abs(moments[0][2] / moments[0][0] - (-0.5)) < eps);
  ASSERT_TRUE(std::abs(moments[0][3] / moments[0][0] - (+0.5)) < eps);
}

// test intersecting a cube, read in from a file,
// against itself
TEST(intersectR3D, cube_self_intersection) {

  // create the mesh factory
  Jali::MeshFactory mesh_factory(MPI_COMM_WORLD);

  // tell the mesh factory what kind of entitites to include in the
  // description
  mesh_factory.included_entities(
      {Jali::Entity_kind::EDGE, Jali::Entity_kind::FACE,
       Jali::Entity_kind::WEDGE, Jali::Entity_kind::CORNER});

  // read the source exodus file
  std::shared_ptr<Jali::Mesh> sm = mesh_factory("../test_data/1x3d.exo");
  ASSERT_TRUE(sm != NULL);

  // read the target exodus file
  std::shared_ptr<Jali::Mesh> tm = mesh_factory("../test_data/1x3d.exo");
  ASSERT_TRUE(tm != NULL);

  const Portage::Jali_Mesh_Wrapper s(*sm);
  const Portage::Jali_Mesh_Wrapper t(*tm);

  const double eps = 1e-12;
  const Portage::IntersectR3D<Portage::Jali_Mesh_Wrapper> isect{s, t};
  const std::vector<std::vector<double>> moments = isect(0, 0);
  for (int i = 0; i < moments.size(); i++) {
    for (int j = 0; j < moments[i].size(); j++) {
      std::cout << "i, j, m " << i << ", " << j << ", " << moments[i][j]
                << std::endl;
    }
  }

  ASSERT_TRUE(moments.size() == 1);
  ASSERT_TRUE(moments[0].size() == 4);

  ASSERT_TRUE(std::abs(moments[0][0] - 1) < eps);
  ASSERT_TRUE(std::abs(moments[0][1] / moments[0][0] - (.5)) < eps);
  ASSERT_TRUE(std::abs(moments[0][2] / moments[0][0] - (.5)) < eps);
  ASSERT_TRUE(std::abs(moments[0][3] / moments[0][0] - (.5)) < eps);
}

// test intersecting a cube, read in from a file,
// against another cube intersecting half it's volume
TEST(intersectR3D, cube_0_5) {

  // create the mesh factory
  Jali::MeshFactory mesh_factory(MPI_COMM_WORLD);

  // tell the mesh factory what kind of entitites to include in the
  // description
  mesh_factory.included_entities(
      {Jali::Entity_kind::EDGE, Jali::Entity_kind::FACE,
       Jali::Entity_kind::WEDGE, Jali::Entity_kind::CORNER});

  // read the source exodus file
  std::shared_ptr<Jali::Mesh> sm = mesh_factory("../test_data/cube_0.exo");
  ASSERT_TRUE(sm != NULL);

  // read the target exodus file
  std::shared_ptr<Jali::Mesh> tm = mesh_factory("../test_data/cube_0.5.exo");
  ASSERT_TRUE(tm != NULL);

  const Portage::Jali_Mesh_Wrapper s(*sm);
  const Portage::Jali_Mesh_Wrapper t(*tm);

  const double eps = 1e-12;
  const Portage::IntersectR3D<Portage::Jali_Mesh_Wrapper> isect{s, t};
  const std::vector<std::vector<double>> moments = isect(0, 0);
  for (int i = 0; i < moments.size(); i++) {
    for (int j = 0; j < moments[i].size(); j++) {
      std::cout << "i, j, m " << i << ", " << j << ", " << moments[i][j]
                << std::endl;
    }
  }

  ASSERT_TRUE(moments.size() == 1);
  ASSERT_TRUE(moments[0].size() == 4);

  ASSERT_TRUE(std::abs(moments[0][0] - .5) < eps);
  ASSERT_TRUE(std::abs(moments[0][1] / moments[0][0] - (.75)) < eps);
  ASSERT_TRUE(std::abs(moments[0][2] / moments[0][0] - (.5)) < eps);
  ASSERT_TRUE(std::abs(moments[0][3] / moments[0][0] - (.5)) < eps);
}

// in this test, the cubes don't intersect at all
TEST(intersectR3D, cube_no_intersect) {

  // create the mesh factory
  Jali::MeshFactory mesh_factory(MPI_COMM_WORLD);

  // tell the mesh factory what kind of entitites to include in the description
  mesh_factory.included_entities(
      {Jali::Entity_kind::EDGE, Jali::Entity_kind::FACE,
       Jali::Entity_kind::WEDGE, Jali::Entity_kind::CORNER});

  // read the source exodus file
  std::shared_ptr<Jali::Mesh> sm = mesh_factory("../test_data/cube_0.exo");
  ASSERT_TRUE(sm != NULL);

  // read the target exodus file
  std::shared_ptr<Jali::Mesh> tm = mesh_factory("../test_data/cube_2.exo");
  ASSERT_TRUE(tm != NULL);

  const Portage::Jali_Mesh_Wrapper s(*sm);
  const Portage::Jali_Mesh_Wrapper t(*tm);

  const double eps = 1e-12;
  const Portage::IntersectR3D<Portage::Jali_Mesh_Wrapper> isect{s, t};
  const std::vector<std::vector<double>> moments = isect(0, 0);
  for (int i = 0; i < moments.size(); i++) {
    for (int j = 0; j < moments[i].size(); j++) {
      std::cout << "i, j, m " << i << ", " << j << ", " << moments[i][j]
                << std::endl;
    }
  }

  ASSERT_TRUE(moments.size() == 1);
  ASSERT_TRUE(moments[0].size() == 4);

  ASSERT_TRUE(std::abs(moments[0][0] - 0) < eps);
  ASSERT_TRUE(isnan(moments[0][1] / moments[0][0]));
  ASSERT_TRUE(isnan(moments[0][2] / moments[0][0]));
  ASSERT_TRUE(isnan(moments[0][3] / moments[0][0]));
}

// in this test, the cubes share a face but have zero intersection volume
TEST(intersectR3D, cube_0_1) {

  // create the mesh factory
  Jali::MeshFactory mesh_factory(MPI_COMM_WORLD);

  // tell the mesh factory what kind of entitites to include in the
  // description
  mesh_factory.included_entities(
      {Jali::Entity_kind::EDGE, Jali::Entity_kind::FACE,
       Jali::Entity_kind::WEDGE, Jali::Entity_kind::CORNER});

  // read the source exodus file
  std::shared_ptr<Jali::Mesh> sm = mesh_factory("../test_data/cube_0.exo");
  ASSERT_TRUE(sm != NULL);

  // read the target exodus file
  std::shared_ptr<Jali::Mesh> tm = mesh_factory("../test_data/cube_1.exo");
  ASSERT_TRUE(tm != NULL);

  const Portage::Jali_Mesh_Wrapper s(*sm);
  const Portage::Jali_Mesh_Wrapper t(*tm);

  const double eps = 1e-12;
  const Portage::IntersectR3D<Portage::Jali_Mesh_Wrapper> isect{s, t};
  const std::vector<std::vector<double>> moments = isect(0, 0);
  for (int i = 0; i < moments.size(); i++) {
    for (int j = 0; j < moments[i].size(); j++) {
      std::cout << "i, j, m " << i << ", " << j << ", " << moments[i][j]
                << std::endl;
    }
  }

  ASSERT_TRUE(moments.size() == 1);
  ASSERT_TRUE(moments[0].size() == 4);

  ASSERT_TRUE(std::abs(moments[0][0] - 0) < eps);
  ASSERT_TRUE(std::abs(moments[0][1] - 0) < eps);
  ASSERT_TRUE(std::abs(moments[0][2] - 0) < eps);
  ASSERT_TRUE(std::abs(moments[0][3] - 0) < eps);
}
