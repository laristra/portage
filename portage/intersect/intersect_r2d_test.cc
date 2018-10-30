/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include "gtest/gtest.h"

// portage includes
#include "portage/intersect/intersect_r2d.h"
#include "portage/support/portage.h"

// wonton includes
#include "wonton/mesh/simple/simple_mesh.h"
#include "wonton/mesh/simple/simple_mesh_wrapper.h"
#include "wonton/state/simple/simple_state.h"
#include "wonton/state/simple/simple_state_wrapper.h"

/*!
 * @brief Intersect two cells on two single cell meshes to compute moments.
 * Intersect two cells contained in the mesh:
 * (0, 0) (2, 0) (2, 2) (0,2) with (1,1) (2,1) (2,2) (1,2)
 * Results should be an area of 1 and a centroid of 1.5, 1.5.
 */
TEST(intersectR2D, simple1) {
  auto sourcemesh = std::make_shared<Wonton::Simple_Mesh>(0, 0, 2, 2, 1, 1);
  auto targetmesh = std::make_shared<Wonton::Simple_Mesh>(1, 1, 2, 2, 1, 1);
  auto sourcestate = std::make_shared<Wonton::Simple_State>(sourcemesh);
  const Wonton::Simple_Mesh_Wrapper sm(*sourcemesh);
  const Wonton::Simple_Mesh_Wrapper tm(*targetmesh);
  const Wonton::Simple_State_Wrapper ss(*sourcestate);

  Portage::IntersectR2D<Portage::Entity_kind::CELL, Wonton::Simple_Mesh_Wrapper,
                        Wonton::Simple_State_Wrapper,
                        Wonton::Simple_Mesh_Wrapper>
      isect{sm, ss, tm};

  std::vector<int> srccells({0});

  std::vector<Portage::Weights_t> srcwts = isect(0, srccells);
  ASSERT_EQ(1, srcwts.size());
  int srcent = srcwts[0].entityID;
  std::vector<double> moments = srcwts[0].weights;
  for (int j = 0; j < moments.size(); j++)
    std::cout << "i, j, m " << srcent << ", " << j << ", " << moments[j]
              << std::endl;

  double eps = 1.0e-12;
  ASSERT_NEAR(moments[0], 1, eps);
  ASSERT_NEAR(moments[1], 1.5, eps);
  ASSERT_NEAR(moments[2], 1.5, eps);
}
