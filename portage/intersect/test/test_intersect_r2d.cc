/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include "gtest/gtest.h"

// wonton includes
#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"
#include "wonton/mesh/simple/simple_mesh.h"
#include "wonton/mesh/simple/simple_mesh_wrapper.h"
#include "wonton/state/simple/simple_state.h"
#include "wonton/state/simple/simple_state_wrapper.h"

// portage includes
#include "portage/support/portage.h"
// templated struct for invoking 2D/3D intersect uniformly (see last test)
// include the include files intersect_r2d.h and intersect_r3d.h
#include "portage/intersect/intersect_rNd.h"

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

  Portage::NumericTolerances_t num_tols = Portage::DEFAULT_NUMERIC_TOLERANCES<2>;

  Portage::IntersectRND<2>::Intersect<Portage::Entity_kind::CELL,
                                      Wonton::Simple_Mesh_Wrapper,
                                      Wonton::Simple_State_Wrapper,
                                      Wonton::Simple_Mesh_Wrapper>
      isect{sm, ss, tm, num_tols};

  std::vector<int> srccells({0});

  std::vector<Portage::Weights_t> srcwts = isect(0, srccells);
  ASSERT_EQ(unsigned(1), srcwts.size());
  int srcent = srcwts[0].entityID;
  std::vector<double> moments = srcwts[0].weights;
  int const num_moments = moments.size();
  for (int j = 0; j < num_moments; j++)
    std::cout << "i, j, m " << srcent << ", " << j << ", " << moments[j] << std::endl;

  double const eps = 1.E-12;
  ASSERT_NEAR(moments[0], 1.0, eps);
  ASSERT_NEAR(moments[1], 1.5, eps);
  ASSERT_NEAR(moments[2], 1.5, eps);
}
