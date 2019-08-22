/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <vector>
#include <algorithm>
#include <memory>
#include <cassert>

#include "gtest/gtest.h"

#include "wonton/mesh/simple/simple_mesh.h"
#include "wonton/mesh/simple/simple_mesh_wrapper.h"
#include "wonton/state/simple/simple_state.h"
#include "wonton/state/simple/simple_state_wrapper.h"
#include "wonton/support/Matrix.h"
#include "portage/support/portage.h"
#include "portage/support/weight.h"
#include "portage/support/faceted_setup.h"

const unsigned int NCELLS=5;

// Checks the face normals and distances of a 2D axis-aligned brick mesh. 
TEST(Faceted_Setup, Simple2D) {
  Wonton::Simple_Mesh mesh(0., 0., 1., .1, NCELLS, NCELLS); // high aspect ratio
  Wonton::Simple_Mesh_Wrapper wrapper(mesh);
  Portage::vector<std::vector<std::vector<double>>> smoothing
    (NCELLS*NCELLS,std::vector<std::vector<double>>(4,std::vector<double>(3)));
  Portage::vector<Wonton::Point<2>> extents;
  double factor = 1.5, bfactor = 1.5;

  Portage::Meshfree::Weight::faceted_setup_cell
    <2,Wonton::Simple_Mesh_Wrapper>(wrapper, smoothing, extents, factor, bfactor);

  double dx0=1.*factor/NCELLS, dx1=.1*factor/NCELLS;
  ASSERT_EQ(smoothing.size(), NCELLS*NCELLS);
  ASSERT_EQ(extents.size(), NCELLS*NCELLS);
  for (int i=0; i<NCELLS; i++) {
    std::vector<std::vector<double>> h = smoothing[i];
    ASSERT_EQ(h.size(), 4);
    ASSERT_NEAR(h[0][0],  0.0, 1.e-12);
    ASSERT_NEAR(h[0][1], -1.0, 1.e-12);
    ASSERT_NEAR(h[0][2],  dx1, 1.e-12);
    ASSERT_NEAR(h[1][0],  1.0, 1.e-12);
    ASSERT_NEAR(h[1][1],  0.0, 1.e-12);
    ASSERT_NEAR(h[1][2],  dx0, 1.e-12); 
    ASSERT_NEAR(h[2][0],  0.0, 1.e-12);
    ASSERT_NEAR(h[2][1],  1.0, 1.e-12);
    ASSERT_NEAR(h[2][2],  dx1, 1.e-12); 
    ASSERT_NEAR(h[3][0], -1.0, 1.e-12);
    ASSERT_NEAR(h[3][1],  0.0, 1.e-12);
    ASSERT_NEAR(h[3][2],  dx0, 1.e-12); 

    Wonton::Point<2> dx=extents[i];
    ASSERT_NEAR(dx[0], 2*dx0, 1.e-12);
    ASSERT_NEAR(dx[1], 2*dx1, 1.e-12);
  }
}

// Checks the face normals and distances of a 3D axis-aligned brick mesh. 
TEST(Faceted_Setup, Simple3D) {
  Wonton::Simple_Mesh mesh(0., 0., 0., 1., .1, .01, NCELLS, NCELLS, NCELLS); // high aspect ratio
  Wonton::Simple_Mesh_Wrapper wrapper(mesh);
  Portage::vector<std::vector<std::vector<double>>> smoothing
    (NCELLS*NCELLS*NCELLS,std::vector<std::vector<double>>(6,std::vector<double>(4)));
  Portage::vector<Wonton::Point<3>> extents;
  double factor = 1.5, bfactor = 1.5;

  Portage::Meshfree::Weight::faceted_setup_cell
    <3,Wonton::Simple_Mesh_Wrapper>(wrapper, smoothing, extents, factor, bfactor);

  double dx0=1.*factor/NCELLS, dx1=.1*factor/NCELLS, dx2=.01*factor/NCELLS;
  ASSERT_EQ(smoothing.size(), NCELLS*NCELLS*NCELLS);
  ASSERT_EQ(extents.size(), NCELLS*NCELLS*NCELLS);
  for (int i=0; i<NCELLS; i++) {
    std::vector<std::vector<double>> h = smoothing[i];
    ASSERT_EQ(h.size(), 6);
    ASSERT_NEAR(h[0][0],  0.0, 1.e-12);
    ASSERT_NEAR(h[0][1], -1.0, 1.e-12);
    ASSERT_NEAR(h[0][2],  0.0, 1.e-12);
    ASSERT_NEAR(h[0][3],  dx1, 1.e-12);
    ASSERT_NEAR(h[1][0],  1.0, 1.e-12);
    ASSERT_NEAR(h[1][1],  0.0, 1.e-12);
    ASSERT_NEAR(h[1][2],  0.0, 1.e-12);
    ASSERT_NEAR(h[1][3],  dx0, 1.e-12); 
    ASSERT_NEAR(h[2][0],  0.0, 1.e-12);
    ASSERT_NEAR(h[2][1],  1.0, 1.e-12);
    ASSERT_NEAR(h[2][2],  0.0, 1.e-12);
    ASSERT_NEAR(h[2][3],  dx1, 1.e-12); 
    ASSERT_NEAR(h[3][0], -1.0, 1.e-12);
    ASSERT_NEAR(h[3][1],  0.0, 1.e-12);
    ASSERT_NEAR(h[3][2],  0.0, 1.e-12);
    ASSERT_NEAR(h[3][3],  dx0, 1.e-12); 
    ASSERT_NEAR(h[4][0],  0.0, 1.e-12);
    ASSERT_NEAR(h[4][1],  0.0, 1.e-12);
    ASSERT_NEAR(h[4][2], -1.0, 1.e-12);
    ASSERT_NEAR(h[4][3],  dx2, 1.e-12); 
    ASSERT_NEAR(h[5][0],  0.0, 1.e-12);
    ASSERT_NEAR(h[5][1],  0.0, 1.e-12);
    ASSERT_NEAR(h[5][2],  1.0, 1.e-12);
    ASSERT_NEAR(h[5][3],  dx2, 1.e-12); 

    Wonton::Point<3> dx=extents[i];
    ASSERT_NEAR(dx[0], 2*dx0, 1.e-12);
    ASSERT_NEAR(dx[1], 2*dx1, 1.e-12);
    ASSERT_NEAR(dx[2], 2*dx2, 1.e-12);
  }
}

// Checks the face normals and distances of a 2D brick mesh tilted 30 degrees 
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
  Portage::vector<std::vector<std::vector<double>>> smoothing
    (NCELLS*NCELLS,std::vector<std::vector<double>>(4,std::vector<double>(3)));
  Portage::vector<Wonton::Point<2>> extents;
  double factor = 1.5, bfactor = 1.5;

  Portage::Meshfree::Weight::faceted_setup_cell
    <2,Wonton::Simple_Mesh_Wrapper>(wrapper, smoothing, extents, factor, bfactor);

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
    std::vector<std::vector<double>> h = smoothing[i];
    ASSERT_EQ(h.size(), 4);
    ASSERT_NEAR(h[0][0],  b,  1.e-12);
    ASSERT_NEAR(h[0][1],  -a,  1.e-12);
    ASSERT_NEAR(h[0][2],  dx1,   1.e-12);
    ASSERT_NEAR(h[1][0],  a,  1.e-12);
    ASSERT_NEAR(h[1][1],  b,  1.e-12);
    ASSERT_NEAR(h[1][2],  dx0,   1.e-12); 
    ASSERT_NEAR(h[2][0],  -b,  1.e-12);
    ASSERT_NEAR(h[2][1],  a,  1.e-12);
    ASSERT_NEAR(h[2][2],  dx1,   1.e-12); 
    ASSERT_NEAR(h[3][0],  -a,  1.e-12);
    ASSERT_NEAR(h[3][1],  -b,  1.e-12);
    ASSERT_NEAR(h[3][2],  dx0,   1.e-12); 

    Wonton::Point<2> dx=extents[i];
    ASSERT_NEAR(dx[0], 2*(xmax[0]-xmin[0]), 1.e-12);
    ASSERT_NEAR(dx[1], 2*(xmax[1]-xmin[1]), 1.e-12);
  }
}

// Checks the face normals and distances of a 3D brick mesh tilted 30 degrees 
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
  Portage::vector<std::vector<std::vector<double>>> smoothing
    (NCELLS*NCELLS*NCELLS,std::vector<std::vector<double>>(6,std::vector<double>(4)));
  Portage::vector<Wonton::Point<3>> extents;
  double factor = 1.5, bfactor = 1.5;

  Portage::Meshfree::Weight::faceted_setup_cell
    <3,Wonton::Simple_Mesh_Wrapper>(wrapper, smoothing, extents, factor, bfactor);

  double dx0=1.*factor/NCELLS, dx1=.1*factor/NCELLS, dx2=.01*factor/NCELLS;
  const std::vector<std::vector<double>> xpts= // Mathematica output
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
    std::vector<std::vector<double>> h = smoothing[i];
    ASSERT_EQ(h.size(), 6);
    ASSERT_NEAR(h[0][0],  a-b, 1.e-11);
    ASSERT_NEAR(h[0][1],  -b, 1.e-11);
    ASSERT_NEAR(h[0][2],  -a-b, 1.e-11);
    ASSERT_NEAR(h[0][3],  dx1, 1.e-11);
    ASSERT_NEAR(h[1][0],  b, 1.e-11);
    ASSERT_NEAR(h[1][1],  a+b, 1.e-11);
    ASSERT_NEAR(h[1][2],  -a+b, 1.e-11);
    ASSERT_NEAR(h[1][3],  dx0, 1.e-11); 
    ASSERT_NEAR(h[2][0],  -a+b, 1.e-11);
    ASSERT_NEAR(h[2][1],  b, 1.e-11);
    ASSERT_NEAR(h[2][2],  a+b, 1.e-11);
    ASSERT_NEAR(h[2][3],  dx1, 1.e-11); 
    ASSERT_NEAR(h[3][0],  -b, 1.e-11);
    ASSERT_NEAR(h[3][1],  -a-b, 1.e-11);
    ASSERT_NEAR(h[3][2],  a-b, 1.e-11);
    ASSERT_NEAR(h[3][3],  dx0, 1.e-11); 
    ASSERT_NEAR(h[4][0],  -a-b, 1.e-11);
    ASSERT_NEAR(h[4][1],  a-b, 1.e-11);
    ASSERT_NEAR(h[4][2],  -b, 1.e-11);
    ASSERT_NEAR(h[4][3],  dx2, 1.e-11); 
    ASSERT_NEAR(h[5][0],  a+b, 1.e-11);
    ASSERT_NEAR(h[5][1],  -a+b, 1.e-11);
    ASSERT_NEAR(h[5][2],  b, 1.e-11);
    ASSERT_NEAR(h[5][3],  dx2, 1.e-11); 

    Wonton::Point<3> dx=extents[i];
    ASSERT_NEAR(dx[0], 2*(xmax[0]-xmin[0]), 1.e-12);
    ASSERT_NEAR(dx[1], 2*(xmax[1]-xmin[1]), 1.e-12);
    ASSERT_NEAR(dx[2], 2*(xmax[2]-xmin[2]), 1.e-12);
  }
}



// Checks the face normals and distances of a 2D axis-aligned brick mesh. 
TEST(Faceted_Setup, Simple2DFace) {
  Wonton::Simple_Mesh mesh(0., 0., 1., .1, NCELLS, NCELLS); // high aspect ratio
  Wonton::Simple_Mesh_Wrapper wrapper(mesh);
  Portage::vector<std::vector<std::vector<double>>> smoothing
    (NCELLS*NCELLS,std::vector<std::vector<double>>(4,std::vector<double>(3)));
  Portage::vector<Wonton::Point<2>> extents;
  double factor = 1.5, bfactor=0.5;

  Portage::Meshfree::Weight::faceted_setup_cell
    <2,Wonton::Simple_Mesh_Wrapper>(wrapper, smoothing, extents, factor, bfactor);

  double dx0=1.*factor/NCELLS, dx1=.1*factor/NCELLS;
  ASSERT_EQ(smoothing.size(), NCELLS*NCELLS);
  ASSERT_EQ(extents.size(), NCELLS*NCELLS);
  for (int i=0; i<NCELLS; i++) {
    std::vector<std::vector<double>> h = smoothing[i];
    ASSERT_EQ(h.size(), 4);
    if (not wrapper.on_exterior_boundary(Wonton::CELL, i)) {
      ASSERT_NEAR(h[0][0],  0.0, 1.e-12);
      ASSERT_NEAR(h[0][1], -1.0, 1.e-12);
      ASSERT_NEAR(h[0][2],  dx1, 1.e-12);
      ASSERT_NEAR(h[1][0],  1.0, 1.e-12);
      ASSERT_NEAR(h[1][1],  0.0, 1.e-12);
      ASSERT_NEAR(h[1][2],  dx0, 1.e-12); 
      ASSERT_NEAR(h[2][0],  0.0, 1.e-12);
      ASSERT_NEAR(h[2][1],  1.0, 1.e-12);
      ASSERT_NEAR(h[2][2],  dx1, 1.e-12); 
      ASSERT_NEAR(h[3][0], -1.0, 1.e-12);
      ASSERT_NEAR(h[3][1],  0.0, 1.e-12);
      ASSERT_NEAR(h[3][2],  dx0, 1.e-12); 

      Wonton::Point<2> dx=extents[i];
      ASSERT_NEAR(dx[0], 2*dx0, 1.e-12);
      ASSERT_NEAR(dx[1], 2*dx1, 1.e-12);
    } else {
      std::vector<int> faces, dirs;
      wrapper.cell_get_faces_and_dirs(i, &faces, &dirs);
      if (wrapper.on_exterior_boundary(Wonton::FACE, faces[0])) {
        ASSERT_NEAR(h[0][2],  dx1*bfactor/factor, 1.e-12);
      }
      if (wrapper.on_exterior_boundary(Wonton::FACE, faces[1])) {
        ASSERT_NEAR(h[1][2],  dx0*bfactor/factor, 1.e-12);
      }
      if (wrapper.on_exterior_boundary(Wonton::FACE, faces[2])) {
        ASSERT_NEAR(h[2][2],  dx1*bfactor/factor, 1.e-12);
      }
      if (wrapper.on_exterior_boundary(Wonton::FACE, faces[2])) {
        ASSERT_NEAR(h[3][2],  dx0*bfactor/factor, 1.e-12);
      }
    }
  }
}



// Checks the face normals and distances of a 2D axis-aligned brick mesh with internal boundaries.
TEST(Faceted_Setup, InternalBoundary) {
  const size_t NCELLS = 4;
  std::shared_ptr<Wonton::Simple_Mesh> mesh_ptr = 
    std::make_shared<Wonton::Simple_Mesh>(-1., -1., 1., 1., NCELLS, NCELLS);
  Wonton::Simple_Mesh &mesh = *mesh_ptr;
  Wonton::Simple_Mesh_Wrapper mwrapper(mesh);
  Portage::vector<std::vector<std::vector<double>>> smoothing;
  Portage::vector<Wonton::Point<2>> extents;
  double factor = 1.25, bfactor=0.5, dx=0.5;

  size_t ncells2d = mesh.num_entities(Wonton::CELL, Wonton::PARALLEL_OWNED);
  assert(ncells2d == NCELLS*NCELLS);
  std::vector<double> values(ncells2d, 0.0);
  for (int i=0; i<ncells2d; i++) {
    Wonton::Point<2> pnt;
    mwrapper.cell_centroid(i, &pnt);
    if      (pnt[0]<0. and pnt[1]<0.) values[i] = 1.;
    else if (pnt[0]>0. and pnt[1]>0.) values[i] = 1.;
  }
  double *valptr = values.data();

  Wonton::Simple_State state(mesh_ptr); 
  std::vector<double> &added = state.add("indicate", Wonton::CELL, valptr);
  Wonton::Simple_State_Wrapper swrapper(state);

  Portage::Meshfree::Weight::faceted_setup_cell
    (mwrapper, swrapper, "indicate", 0.25, 
     smoothing, extents, 
     factor, bfactor);

  for (int i=0; i<ncells2d; i++) {
    Wonton::Point<2> pnt;
    mwrapper.cell_centroid(i, &pnt);
    std::vector<std::vector<double>> h = smoothing[i];
    if (pnt[0] == -0.25 and pnt[1] == -0.75) {
      ASSERT_EQ(h[0][2],  factor*dx);
      ASSERT_EQ(h[1][2], bfactor*dx);
      ASSERT_EQ(h[2][2],  factor*dx);
      ASSERT_EQ(h[3][2],  factor*dx);
    } else if (pnt[0] == -0.25 and pnt[1] == -0.25) {
      ASSERT_EQ(h[0][2],  factor*dx);
      ASSERT_EQ(h[1][2], bfactor*dx);
      ASSERT_EQ(h[2][2], bfactor*dx);
      ASSERT_EQ(h[3][2],  factor*dx);
    } else if (pnt[0] == -0.25 and pnt[1] ==  0.25) {
      ASSERT_EQ(h[0][2], bfactor*dx);
      ASSERT_EQ(h[1][2], bfactor*dx);
      ASSERT_EQ(h[2][2],  factor*dx);
      ASSERT_EQ(h[3][2],  factor*dx);
    } else if (pnt[0] == -0.25 and pnt[1] ==  0.75) {
      ASSERT_EQ(h[0][2],  factor*dx);
      ASSERT_EQ(h[1][2], bfactor*dx);
      ASSERT_EQ(h[2][2],  factor*dx);
      ASSERT_EQ(h[3][2],  factor*dx);

    } else if (pnt[0] ==  0.25 and pnt[1] == -0.75) {
      ASSERT_EQ(h[0][2],  factor*dx);
      ASSERT_EQ(h[1][2],  factor*dx);
      ASSERT_EQ(h[2][2],  factor*dx);
      ASSERT_EQ(h[3][2], bfactor*dx);
    } else if (pnt[0] ==  0.25 and pnt[1] == -0.25) {
      ASSERT_EQ(h[0][2],  factor*dx);
      ASSERT_EQ(h[1][2],  factor*dx);
      ASSERT_EQ(h[2][2], bfactor*dx);
      ASSERT_EQ(h[3][2], bfactor*dx);
    } else if (pnt[0] ==  0.25 and pnt[1] ==  0.25) {
      ASSERT_EQ(h[0][2], bfactor*dx);
      ASSERT_EQ(h[1][2],  factor*dx);
      ASSERT_EQ(h[2][2],  factor*dx);
      ASSERT_EQ(h[3][2], bfactor*dx);
    } else if (pnt[0] ==  0.25 and pnt[1] ==  0.75) {
      ASSERT_EQ(h[0][2],  factor*dx);
      ASSERT_EQ(h[1][2],  factor*dx);
      ASSERT_EQ(h[2][2],  factor*dx);
      ASSERT_EQ(h[3][2], bfactor*dx);

    } else if (pnt[0] ==  0.75 and pnt[1] == -0.25) {
      ASSERT_EQ(h[0][2],  factor*dx);
      ASSERT_EQ(h[1][2],  factor*dx);
      ASSERT_EQ(h[2][2], bfactor*dx);
      ASSERT_EQ(h[3][2],  factor*dx);

    } else if (pnt[0] ==  0.75 and pnt[1] ==  0.25) {
      ASSERT_EQ(h[0][2], bfactor*dx);
      ASSERT_EQ(h[1][2],  factor*dx);
      ASSERT_EQ(h[2][2],  factor*dx);
      ASSERT_EQ(h[3][2],  factor*dx);

    } else if (pnt[0] == -0.75 and pnt[1] == -0.25) {
      ASSERT_EQ(h[0][2],  factor*dx);
      ASSERT_EQ(h[1][2],  factor*dx);
      ASSERT_EQ(h[2][2], bfactor*dx);
      ASSERT_EQ(h[3][2],  factor*dx);

    } else if (pnt[0] == -0.75 and pnt[1] ==  0.25) {
      ASSERT_EQ(h[0][2], bfactor*dx);
      ASSERT_EQ(h[1][2],  factor*dx);
      ASSERT_EQ(h[2][2],  factor*dx);
      ASSERT_EQ(h[3][2],  factor*dx);

    } else {
      ASSERT_EQ(h[0][2],  factor*dx);
      ASSERT_EQ(h[1][2],  factor*dx);
      ASSERT_EQ(h[2][2],  factor*dx);
      ASSERT_EQ(h[3][2],  factor*dx);
    }
  }
}
    
