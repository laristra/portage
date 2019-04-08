/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <vector>

#include "gtest/gtest.h"

#include "wonton/mesh/simple/simple_mesh.h"
#include "wonton/mesh/simple/simple_mesh_wrapper.h"
#include "portage/support/weight.h"
#include "portage/support/faceted_setup.h"

const unsigned int NCELLS=5;

TEST(Faceted_Setup, Simple2D) {
  Wonton::Simple_Mesh mesh(0., 0., 1., .1, NCELLS, NCELLS); // high aspect ratio
  Wonton::Simple_Mesh_Wrapper wrapper(mesh);
  std::vector<std::vector<std::vector<double>>> smoothing
    (NCELLS*NCELLS,std::vector<std::vector<double>>(4,std::vector<double>(3)));
  Portage::Meshfree::Weight::faceted_setup_cell<2,
    Wonton::Simple_Mesh_Wrapper>(wrapper, smoothing, 1.0);
  double dx1=1./NCELLS, dx0=.1/NCELLS;
  for (int i=0; i<NCELLS; i++) {
    ASSERT_EQ(smoothing[i].size(), 4);
    ASSERT_NEAR(smoothing[i][0][0],  0.0, 1.e-12);
    ASSERT_NEAR(smoothing[i][0][1], -1.0, 1.e-12);
    ASSERT_NEAR(smoothing[i][0][2],  dx0, 1.e-12);
    ASSERT_NEAR(smoothing[i][1][0],  1.0, 1.e-12);
    ASSERT_NEAR(smoothing[i][1][1],  0.0, 1.e-12);
    ASSERT_NEAR(smoothing[i][1][2],  dx1, 1.e-12); 
    ASSERT_NEAR(smoothing[i][2][0],  0.0, 1.e-12);
    ASSERT_NEAR(smoothing[i][2][1],  1.0, 1.e-12);
    ASSERT_NEAR(smoothing[i][2][2],  dx0, 1.e-12); 
    ASSERT_NEAR(smoothing[i][3][0], -1.0, 1.e-12);
    ASSERT_NEAR(smoothing[i][3][1],  0.0, 1.e-12);
    ASSERT_NEAR(smoothing[i][3][2],  dx1, 1.e-12); 
  }
}
