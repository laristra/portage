/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include "gtest/gtest.h"

// wonton includes
#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"
#include "wonton/support/CoordinateSystem.h"
#include "wonton/mesh/simple/simple_mesh.h"
#include "wonton/mesh/simple/simple_mesh_wrapper.h"
#include "wonton/state/simple/simple_state.h"
#include "wonton/state/simple/simple_state_wrapper.h"

// portage includes
#include "portage/support/portage.h"
#include "portage/intersect/intersect_r2d.h"

/*!
 * @brief Intersect two cells on two single cell meshes to compute moments.
 * Intersect two cells contained in the mesh:
 * (0, 0) (2, 0) (2, 2) (0,2) with (1,1) (2,1) (2,2) (1,2)
 * XY: results should be an area of 1 and a centroid of 1.5, 1.5.
 * RZ: results should be an area of 1.5 and a centroid of 14/9, 1.5.
 */

class intersectR2D : public ::testing::TestWithParam<Wonton::CoordSysType> {};

TEST_P(intersectR2D, simple1) {
  auto sys = GetParam();

  auto sourcemesh = std::make_shared<Wonton::Simple_Mesh>(0, 0, 2, 2, 1, 1);
  auto targetmesh = std::make_shared<Wonton::Simple_Mesh>(1, 1, 2, 2, 1, 1);

  const Wonton::Simple_Mesh_Wrapper sm(*sourcemesh, true, true, true, sys);
  const Wonton::Simple_Mesh_Wrapper tm(*targetmesh, true, true, true, sys);

  auto sourcestate = std::make_shared<Wonton::Simple_State>(sourcemesh);
  const Wonton::Simple_State_Wrapper ss(*sourcestate);

  Portage::NumericTolerances_t num_tols = Portage::DEFAULT_NUMERIC_TOLERANCES<2>;

  Portage::IntersectR2D<Portage::Entity_kind::CELL,
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
    std::cout << "i, j, m: " << srcent << ", " << j << ", " << moments[j] << std::endl;

  double area;
  Wonton::Point<2> centroid;

  if (sys == Wonton::CoordSysType::Cartesian) {
    area = 1.0;
    centroid = Wonton::Point<2>{1.5, 1.5};
  } else {
    area = 2 * M_PI * 1.5;
    centroid = Wonton::Point<2>{14./9, 1.5};
  }

  double const eps = 1.E-12;
  ASSERT_NEAR(moments[0], area, eps);
  ASSERT_NEAR(moments[1], area * centroid[0], eps);
  ASSERT_NEAR(moments[2], area * centroid[1], eps);
}


INSTANTIATE_TEST_CASE_P(
  intersectR2DAll,
  intersectR2D,
  ::testing::Values(Wonton::CoordSysType::Cartesian, Wonton::CoordSysType::CylindricalAxisymmetric));

