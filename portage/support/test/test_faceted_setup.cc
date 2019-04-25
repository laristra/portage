/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <vector>
#include <algorithm>

#include "gtest/gtest.h"

#include "wonton/mesh/simple/simple_mesh.h"
#include "wonton/mesh/simple/simple_mesh_wrapper.h"
#include "wonton/support/Matrix.h"
#include "portage/support/weight.h"
#include "portage/support/faceted_setup.h"

const unsigned int NCELLS=5;

TEST(Faceted_Setup, Simple2D) {
  Wonton::Simple_Mesh mesh(0., 0., 1., .1, NCELLS, NCELLS); // high aspect ratio
  Wonton::Simple_Mesh_Wrapper wrapper(mesh);
  std::vector<std::vector<std::vector<double>>> smoothing
    (NCELLS*NCELLS,std::vector<std::vector<double>>(4,std::vector<double>(3)));
  std::vector<Wonton::Point<2>> extents;
  double factor = 1.5;

  Portage::Meshfree::Weight::faceted_setup_cell
    <2,Wonton::Simple_Mesh_Wrapper>(wrapper, smoothing, extents, factor);

  double dx0=1.*factor/NCELLS, dx1=.1*factor/NCELLS;
  ASSERT_EQ(smoothing.size(), NCELLS*NCELLS);
  ASSERT_EQ(extents.size(), NCELLS*NCELLS);
  for (int i=0; i<NCELLS; i++) {
    ASSERT_EQ(smoothing[i].size(), 4);
    ASSERT_NEAR(smoothing[i][0][0],  0.0, 1.e-12);
    ASSERT_NEAR(smoothing[i][0][1], -1.0, 1.e-12);
    ASSERT_NEAR(smoothing[i][0][2],  dx1, 1.e-12);
    ASSERT_NEAR(smoothing[i][1][0],  1.0, 1.e-12);
    ASSERT_NEAR(smoothing[i][1][1],  0.0, 1.e-12);
    ASSERT_NEAR(smoothing[i][1][2],  dx0, 1.e-12); 
    ASSERT_NEAR(smoothing[i][2][0],  0.0, 1.e-12);
    ASSERT_NEAR(smoothing[i][2][1],  1.0, 1.e-12);
    ASSERT_NEAR(smoothing[i][2][2],  dx1, 1.e-12); 
    ASSERT_NEAR(smoothing[i][3][0], -1.0, 1.e-12);
    ASSERT_NEAR(smoothing[i][3][1],  0.0, 1.e-12);
    ASSERT_NEAR(smoothing[i][3][2],  dx0, 1.e-12); 

    ASSERT_NEAR(extents[i][0], 2*dx0, 1.e-12);
    ASSERT_NEAR(extents[i][1], 2*dx1, 1.e-12);
  }
}

TEST(Faceted_Setup, Simple3D) {
  Wonton::Simple_Mesh mesh(0., 0., 0., 1., .1, .01, NCELLS, NCELLS, NCELLS); // high aspect ratio
  Wonton::Simple_Mesh_Wrapper wrapper(mesh);
  std::vector<std::vector<std::vector<double>>> smoothing
    (NCELLS*NCELLS*NCELLS,std::vector<std::vector<double>>(6,std::vector<double>(4)));
  std::vector<Wonton::Point<3>> extents;
  double factor = 1.5;

  Portage::Meshfree::Weight::faceted_setup_cell
    <3,Wonton::Simple_Mesh_Wrapper>(wrapper, smoothing, extents, factor);

  double dx0=1.*factor/NCELLS, dx1=.1*factor/NCELLS, dx2=.01*factor/NCELLS;
  ASSERT_EQ(smoothing.size(), NCELLS*NCELLS*NCELLS);
  ASSERT_EQ(extents.size(), NCELLS*NCELLS*NCELLS);
  for (int i=0; i<NCELLS; i++) {
    ASSERT_EQ(smoothing[i].size(), 6);
    ASSERT_NEAR(smoothing[i][0][0],  0.0, 1.e-12);
    ASSERT_NEAR(smoothing[i][0][1], -1.0, 1.e-12);
    ASSERT_NEAR(smoothing[i][0][2],  0.0, 1.e-12);
    ASSERT_NEAR(smoothing[i][0][3],  dx1, 1.e-12);
    ASSERT_NEAR(smoothing[i][1][0],  1.0, 1.e-12);
    ASSERT_NEAR(smoothing[i][1][1],  0.0, 1.e-12);
    ASSERT_NEAR(smoothing[i][1][2],  0.0, 1.e-12);
    ASSERT_NEAR(smoothing[i][1][3],  dx0, 1.e-12); 
    ASSERT_NEAR(smoothing[i][2][0],  0.0, 1.e-12);
    ASSERT_NEAR(smoothing[i][2][1],  1.0, 1.e-12);
    ASSERT_NEAR(smoothing[i][2][2],  0.0, 1.e-12);
    ASSERT_NEAR(smoothing[i][2][3],  dx1, 1.e-12); 
    ASSERT_NEAR(smoothing[i][3][0], -1.0, 1.e-12);
    ASSERT_NEAR(smoothing[i][3][1],  0.0, 1.e-12);
    ASSERT_NEAR(smoothing[i][3][2],  0.0, 1.e-12);
    ASSERT_NEAR(smoothing[i][3][3],  dx0, 1.e-12); 
    ASSERT_NEAR(smoothing[i][4][0],  0.0, 1.e-12);
    ASSERT_NEAR(smoothing[i][4][1],  0.0, 1.e-12);
    ASSERT_NEAR(smoothing[i][4][2], -1.0, 1.e-12);
    ASSERT_NEAR(smoothing[i][4][3],  dx2, 1.e-12); 
    ASSERT_NEAR(smoothing[i][5][0],  0.0, 1.e-12);
    ASSERT_NEAR(smoothing[i][5][1],  0.0, 1.e-12);
    ASSERT_NEAR(smoothing[i][5][2],  1.0, 1.e-12);
    ASSERT_NEAR(smoothing[i][5][3],  dx2, 1.e-12); 

    ASSERT_NEAR(extents[i][0], 2*dx0, 1.e-12);
    ASSERT_NEAR(extents[i][1], 2*dx1, 1.e-12);
    ASSERT_NEAR(extents[i][2], 2*dx2, 1.e-12);
  }
}

TEST(Faceted_Setup, Simple2D_Tilted) {
  Wonton::Simple_Mesh mesh(0., 0., 1., .1, NCELLS, NCELLS); // high aspect ratio
  // rotate Pi/6 radians and translate by {10,20}.
  // by Mathematica
  const double a=sqrt(3.)/2, b=1./2.;
  std::vector<std::vector<double>> affinev(2,std::vector<double>(3));
  affinev[0] = {a,-b, 10.};
  affinev[1] = {b, a, 20.};
  Wonton::Matrix affine(affinev);
  mesh.transform<2>(affine);
  // set smoothing
  Wonton::Simple_Mesh_Wrapper wrapper(mesh);
  std::vector<std::vector<std::vector<double>>> smoothing
    (NCELLS*NCELLS,std::vector<std::vector<double>>(4,std::vector<double>(3)));
  std::vector<Wonton::Point<2>> extents;
  double factor = 1.5;

  Portage::Meshfree::Weight::faceted_setup_cell
    <2,Wonton::Simple_Mesh_Wrapper>(wrapper, smoothing, extents, factor);

  double dx0=1.*factor/NCELLS, dx1=.1*factor/NCELLS;
  const std::vector<std::vector<double>> xpts=
    {{0, a*dx0, a*dx0 - b*dx1, -(b*dx1)}, {0, b*dx0, b*dx0 + a*dx1, a*dx1}};
  const std::vector<double> xmin={*std::min_element(xpts[0].begin(),xpts[0].end()), 
                                  *std::min_element(xpts[1].begin(),xpts[1].end())};
  const std::vector<double> xmax={*std::max_element(xpts[0].begin(),xpts[0].end()), 
                                  *std::max_element(xpts[1].begin(),xpts[1].end())};
  ASSERT_EQ(smoothing.size(), NCELLS*NCELLS);
  ASSERT_EQ(extents.size(), NCELLS*NCELLS);
  for (int i=0; i<NCELLS; i++) {
    ASSERT_EQ(smoothing[i].size(), 4);
    ASSERT_NEAR(smoothing[i][0][0],  b,  1.e-12);
    ASSERT_NEAR(smoothing[i][0][1],  -a,  1.e-12);
    ASSERT_NEAR(smoothing[i][0][2],  dx1,   1.e-12);
    ASSERT_NEAR(smoothing[i][1][0],  a,  1.e-12);
    ASSERT_NEAR(smoothing[i][1][1],  b,  1.e-12);
    ASSERT_NEAR(smoothing[i][1][2],  dx0,   1.e-12); 
    ASSERT_NEAR(smoothing[i][2][0],  -b,  1.e-12);
    ASSERT_NEAR(smoothing[i][2][1],  a,  1.e-12);
    ASSERT_NEAR(smoothing[i][2][2],  dx1,   1.e-12); 
    ASSERT_NEAR(smoothing[i][3][0],  -a,  1.e-12);
    ASSERT_NEAR(smoothing[i][3][1],  -b,  1.e-12);
    ASSERT_NEAR(smoothing[i][3][2],  dx0,   1.e-12); 

    ASSERT_NEAR(extents[i][0], 2*(xmax[0]-xmin[0]), 1.e-12);
    ASSERT_NEAR(extents[i][1], 2*(xmax[1]-xmin[1]), 1.e-12);
  }
}

TEST(Faceted_Setup, Simple3D_Tilted) {
  Wonton::Simple_Mesh mesh(0., 0., 0., 1., .1, .01, NCELLS, NCELLS, NCELLS); // high aspect ratio
  // rotate Pi/2 radians about the axis {1,1,1} and translate by {10,20,30}.
  // by Mathematica
  const double a=1./sqrt(3.), b=1./3.;
  std::vector<std::vector<double>> affinev(3,std::vector<double>(4));
  affinev[0] = {b,b-a,b+a, 10.};
  affinev[1] = {b+a,b,b-a, 20.};
  affinev[2] = {b-a,b+a,b, 30.};
  Wonton::Matrix affine(affinev);
  mesh.transform<3>(affine);
  // set smoothing
  Wonton::Simple_Mesh_Wrapper wrapper(mesh);
  std::vector<std::vector<std::vector<double>>> smoothing
    (NCELLS*NCELLS*NCELLS,std::vector<std::vector<double>>(6,std::vector<double>(4)));
  std::vector<Wonton::Point<3>> extents;
  double factor = 1.5;

  Portage::Meshfree::Weight::faceted_setup_cell
    <3,Wonton::Simple_Mesh_Wrapper>(wrapper, smoothing, extents, factor);

  double dx0=1.*factor/NCELLS, dx1=.1*factor/NCELLS, dx2=.01*factor/NCELLS;
  const std::vector<std::vector<double>> xpts=
{{0, b*dx0, b*dx0 + (-a + b)*dx1, (-a + b)*dx1, (a + b)*dx2, b*dx0 + (a + b)*dx2, 
  b*dx0 + (-a + b)*dx1 + (a + b)*dx2, (-a + b)*dx1 + (a + b)*dx2}, 
 {0, (a + b)*dx0, (a + b)*dx0 + b*dx1, b*dx1, (-a + b)*dx2, (a + b)*dx0 + (-a + b)*dx2, 
  (a + b)*dx0 + b*dx1 + (-a + b)*dx2, b*dx1 + (-a + b)*dx2}, 
 {0, (-a + b)*dx0, (-a + b)*dx0 + (a + b)*dx1, (a + b)*dx1, b*dx2, (-a + b)*dx0 + b*dx2, 
     (-a + b)*dx0 + (a + b)*dx1 + b*dx2, (a + b)*dx1 + b*dx2}};
  const std::vector<double> xmin={*std::min_element(xpts[0].begin(),xpts[0].end()), 
                                  *std::min_element(xpts[1].begin(),xpts[1].end()), 
                                  *std::min_element(xpts[2].begin(),xpts[2].end())};
  const std::vector<double> xmax={*std::max_element(xpts[0].begin(),xpts[0].end()), 
                                  *std::max_element(xpts[1].begin(),xpts[1].end()), 
                                  *std::max_element(xpts[2].begin(),xpts[2].end())};
  ASSERT_EQ(smoothing.size(), NCELLS*NCELLS*NCELLS);
  ASSERT_EQ(extents.size(), NCELLS*NCELLS*NCELLS);
  for (int i=0; i<NCELLS; i++) {
    ASSERT_EQ(smoothing[i].size(), 6);
    ASSERT_NEAR(smoothing[i][0][0],  a-b, 1.e-11);
    ASSERT_NEAR(smoothing[i][0][1],  -b, 1.e-11);
    ASSERT_NEAR(smoothing[i][0][2],  -a-b, 1.e-11);
    ASSERT_NEAR(smoothing[i][0][3],  dx1, 1.e-11);
    ASSERT_NEAR(smoothing[i][1][0],  b, 1.e-11);
    ASSERT_NEAR(smoothing[i][1][1],  a+b, 1.e-11);
    ASSERT_NEAR(smoothing[i][1][2],  -a+b, 1.e-11);
    ASSERT_NEAR(smoothing[i][1][3],  dx0, 1.e-11); 
    ASSERT_NEAR(smoothing[i][2][0],  -a+b, 1.e-11);
    ASSERT_NEAR(smoothing[i][2][1],  b, 1.e-11);
    ASSERT_NEAR(smoothing[i][2][2],  a+b, 1.e-11);
    ASSERT_NEAR(smoothing[i][2][3],  dx1, 1.e-11); 
    ASSERT_NEAR(smoothing[i][3][0],  -b, 1.e-11);
    ASSERT_NEAR(smoothing[i][3][1],  -a-b, 1.e-11);
    ASSERT_NEAR(smoothing[i][3][2],  a-b, 1.e-11);
    ASSERT_NEAR(smoothing[i][3][3],  dx0, 1.e-11); 
    ASSERT_NEAR(smoothing[i][4][0],  -a-b, 1.e-11);
    ASSERT_NEAR(smoothing[i][4][1],  a-b, 1.e-11);
    ASSERT_NEAR(smoothing[i][4][2],  -b, 1.e-11);
    ASSERT_NEAR(smoothing[i][4][3],  dx2, 1.e-11); 
    ASSERT_NEAR(smoothing[i][5][0],  a+b, 1.e-11);
    ASSERT_NEAR(smoothing[i][5][1],  -a+b, 1.e-11);
    ASSERT_NEAR(smoothing[i][5][2],  b, 1.e-11);
    ASSERT_NEAR(smoothing[i][5][3],  dx2, 1.e-11); 

    ASSERT_NEAR(extents[i][0], 2*(xmax[0]-xmin[0]), 1.e-12);
    ASSERT_NEAR(extents[i][1], 2*(xmax[1]-xmin[1]), 1.e-12);
    ASSERT_NEAR(extents[i][2], 2*(xmax[2]-xmin[2]), 1.e-12);
  }
}
