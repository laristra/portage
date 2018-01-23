/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include "intersect_r2d.h"
#include "gtest/gtest.h"
#include "portage/support/portage.h"
#include "portage/wonton/mesh/simple_mesh/simple_mesh_wrapper.h"

/*!
 * @brief Intersect two cells on two single cell meshes to compute moments.
 * Intersect two cells contained in the mesh:
 * (0, 0) (2, 0) (2, 2) (0,2) with (1,1) (2,1) (2,2) (1,2)
 * Results should be an area of 1 and a centroid of 1.5, 1.5.
 */
TEST(intersectR2D, simple1) {
  Portage::Simple_Mesh sm{0, 0, 2, 2, 1, 1};
  Portage::Simple_Mesh tm{1, 1, 2, 2, 1, 1};
  const Wonton::Simple_Mesh_Wrapper s(sm);
  const Wonton::Simple_Mesh_Wrapper t(tm);

  Portage::IntersectR2D<Portage::Entity_kind::CELL, Wonton::Simple_Mesh_Wrapper,
                        Wonton::Simple_Mesh_Wrapper>
      isect{s, t};

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
