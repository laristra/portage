/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/


#include <iostream>

#include "gtest/gtest.h"

// wonton includes
#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"
#include "wonton/support/Vector.h"
#include "wonton/mesh/simple/simple_mesh.h"
#include "wonton/mesh/simple/simple_mesh_wrapper.h"
#include "wonton/state/simple/simple_state.h"
#include "wonton/state/simple/simple_state_wrapper.h"

// portage includes
#include "portage/interpolate/quadfit.h"
#include "portage/support/portage.h"


/// Test quadfit computation for cell centered fields

TEST(Quadfit, Fields_Cell_Ctr) {

  // Create a 4 cell mesh
  std::shared_ptr<Wonton::Simple_Mesh> mesh1 =
      std::make_shared<Wonton::Simple_Mesh>(0.0, 0.0, 1.0, 1.0, 4, 4);
  ASSERT_TRUE(mesh1 != nullptr);

  // create the wrapper
  Wonton::Simple_Mesh_Wrapper meshWrapper(*mesh1);

  // Create a state object
  Wonton::Simple_State mystate(mesh1);

  // Create a state Wrapper
  Wonton::Simple_State_Wrapper stateWrapper(mystate);

  const int nc1 = meshWrapper.num_owned_cells();

  // Define three state vectors, one with constant value, the next
  // with a linear function that is x+2y, the final x*x+y*y

  std::vector<double> data1(nc1, 1.25);

  // add the data vector to the state
  mystate.add("cellvars1", Portage::Entity_kind::CELL, &(data1[0]));

  // create the second vector
  std::vector<double> data2(nc1);

  // set the data (x+2*y)
  for (int c = 0; c < nc1; c++) {
    Wonton::Point<2> ccen;
    meshWrapper.cell_centroid(c, &ccen);
    data2[c] = ccen[0] + 2 * ccen[1];
  }

  // add the second data vector to the state
  mystate.add("cellvars2", Portage::Entity_kind::CELL, &(data2[0]));

  // create the third vector
  std::vector<double> data3(nc1);

  // set the data (x*x+y+y)
  for (int c = 0; c < nc1; c++) {
    Wonton::Point<2> ccen;
    meshWrapper.cell_centroid(c, &ccen);
    data3[c] = ccen[0]*ccen[0] + ccen[1]*ccen[1];
  }

  // add the second data vector to the state
  mystate.add("cellvars3", Portage::Entity_kind::CELL, &(data3[0]));

  // Create Quadfit objects

  Portage::Limited_Quadfit<2, Portage::Entity_kind::CELL,
                           Wonton::Simple_Mesh_Wrapper,
                           Wonton::Simple_State_Wrapper>
      qfitcalc1(meshWrapper, stateWrapper, "cellvars1", 
                Portage::NOLIMITER, Portage::BND_NOLIMITER);
  Portage::Limited_Quadfit<2, Portage::Entity_kind::CELL,
                           Wonton::Simple_Mesh_Wrapper,
                           Wonton::Simple_State_Wrapper>
      qfitcalc2(meshWrapper, stateWrapper, "cellvars2", 
                Portage::NOLIMITER, Portage::BND_NOLIMITER);
  Portage::Limited_Quadfit<2, Portage::Entity_kind::CELL,
                           Wonton::Simple_Mesh_Wrapper,
                           Wonton::Simple_State_Wrapper>
      qfitcalc3(meshWrapper, stateWrapper, "cellvars3", 
                Portage::NOLIMITER, Portage::BND_NOLIMITER);
  Portage::Limited_Quadfit<2, Portage::Entity_kind::CELL,
                           Wonton::Simple_Mesh_Wrapper,
                           Wonton::Simple_State_Wrapper>
      qfitcalc4(meshWrapper, stateWrapper, "cellvars1",
                Portage::BARTH_JESPERSEN, Portage::BND_NOLIMITER);
  Portage::Limited_Quadfit<2, Portage::Entity_kind::CELL,
                           Wonton::Simple_Mesh_Wrapper,
                           Wonton::Simple_State_Wrapper>
      qfitcalc5(meshWrapper, stateWrapper, "cellvars2",
                Portage::BARTH_JESPERSEN, Portage::BND_NOLIMITER);
  Portage::Limited_Quadfit<2, Portage::Entity_kind::CELL,
                           Wonton::Simple_Mesh_Wrapper,
                           Wonton::Simple_State_Wrapper>
      qfitcalc6(meshWrapper, stateWrapper, "cellvars3",
                Portage::BARTH_JESPERSEN, Portage::BND_NOLIMITER);

  // Compute the quadfit for each of these fields

  Wonton::Vector<5> qfit;

  // Verify the quadfit values
  // For field 1 (constant), it is is 0,0
  // For field 2 (x+2y), it is (1,2)
  // For field 3 x^2+y^2 it is (2x, 2y, 1, 0, 1)

  for (int c = 0; c < nc1; ++c) {
    // Create boundary_cell logical variable

    bool boundary_cell = false;
    std::vector<int> cfaces, cellfaceDirs;
    meshWrapper.cell_get_faces_and_dirs(c, &cfaces, &cellfaceDirs);
    for (auto f : cfaces) {
      std::vector<int> fcells;
      meshWrapper.face_get_cells(f, Portage::Entity_type::ALL, &fcells);
      if (fcells.size() == 1) {
        boundary_cell = true;
        break;
      }
    }

    // unlimited quadfit of constant function

    qfit = qfitcalc1(c);
    ASSERT_NEAR(0.0, qfit[0], 1.0e-10);
    ASSERT_NEAR(0.0, qfit[1], 1.0e-10);
    ASSERT_NEAR(0.0, qfit[2], 1.0e-10);
    ASSERT_NEAR(0.0, qfit[3], 1.0e-10);
    ASSERT_NEAR(0.0, qfit[4], 1.0e-10);

    // unlimited quadfit of linear function

    qfit = qfitcalc2(c);
    ASSERT_NEAR(1.0, qfit[0], 1.0e-10);
    ASSERT_NEAR(2.0, qfit[1], 1.0e-10);
    ASSERT_NEAR(0.0, qfit[2], 1.0e-10);
    ASSERT_NEAR(0.0, qfit[3], 1.0e-10);
    ASSERT_NEAR(0.0, qfit[4], 1.0e-10);

    // unlimited quadfit of quadratic function
    //
    // Avoid boundary cells.  LimitedQuadfit won't crash at boundaries, however it
    // won't give the correct answer in the ASSERT statements for the quadratic
    // function test.

    if (!boundary_cell) {
      qfit = qfitcalc3(c);
      Wonton::Point<2> ccen;
      meshWrapper.cell_centroid(c, &ccen);
      double cx = ccen[0]; // x
      double cy = ccen[1]; // y
      ASSERT_NEAR(2.0*cx, qfit[0], 1.0e-10); // partial of f wrt x
      ASSERT_NEAR(2.0*cy, qfit[1], 1.0e-10); // partial of f wrt y
      ASSERT_NEAR(1.0, qfit[2], 1.0e-10); // (2nd partial of f wrt x)/2
      ASSERT_NEAR(0.0, qfit[3], 1.0e-10); // (2nd partials of f wrt x and y)/2
      ASSERT_NEAR(1.0, qfit[4], 1.0e-10); // (2nd partial of f wrt y)/2
    }

    // limited quadfit of constant function

    qfit = qfitcalc4(c);
    ASSERT_NEAR(0.0, qfit[0], 1.0e-10);
    ASSERT_NEAR(0.0, qfit[1], 1.0e-10);
    ASSERT_NEAR(0.0, qfit[2], 1.0e-10);
    ASSERT_NEAR(0.0, qfit[3], 1.0e-10);
    ASSERT_NEAR(0.0, qfit[4], 1.0e-10);

    // limited quadfit of linear function
    if (!boundary_cell) {
      qfit = qfitcalc5(c);
      ASSERT_NEAR(1.0, qfit[0], 1.0e-10);
      ASSERT_NEAR(2.0, qfit[1], 1.0e-10);
      ASSERT_NEAR(0.0, qfit[2], 1.0e-10);
      ASSERT_NEAR(0.0, qfit[3], 1.0e-10);
      ASSERT_NEAR(0.0, qfit[4], 1.0e-10);
    }
    // limited quadfit of quadratic function
    //
    // Avoid boundary cells.  LimitedQuadfit won't crash at boundaries, but it
    // won't give the correct answer in the ASSERT statements for the quadratic
    // function test.

    if (!boundary_cell) {
      qfit = qfitcalc6(c);
      Wonton::Point<2> ccen;
      meshWrapper.cell_centroid(c, &ccen);
      double cx = ccen[0]; // x
      double cy = ccen[1]; // y
      ASSERT_NEAR(2.0*cx, qfit[0], 1.0e-10); // partial of f wrt x
      ASSERT_NEAR(2.0*cy, qfit[1], 1.0e-10); // partial of f wrt y
      ASSERT_NEAR(1.0, qfit[2], 1.0e-10); // (2nd partial of f wrt x)/2
      ASSERT_NEAR(0.0, qfit[3], 1.0e-10); // (2nd partials of f wrt x and y)/2
      ASSERT_NEAR(1.0, qfit[4], 1.0e-10); // (2nd partial of f wrt y)/2
    }

  }
}

/// Test quadfit computation with node centered fields

TEST(Quadfit, Fields_Node_Ctr) {

  // Create a 4 cell mesh
  std::shared_ptr<Wonton::Simple_Mesh> mesh1 =
      std::make_shared<Wonton::Simple_Mesh>(0.0, 0.0, 1.0, 1.0, 4, 4);
  ASSERT_TRUE(mesh1 != nullptr);

  // create the wrapper
  Wonton::Simple_Mesh_Wrapper meshWrapper(*mesh1);

  // Create a state object
  Wonton::Simple_State mystate(mesh1);

  // Create a state Wrapper
  Wonton::Simple_State_Wrapper stateWrapper(mystate);

  const int nn1 = meshWrapper.num_owned_nodes();

  // Define three state vectors, one with constant value, the next
  // with a linear function that is x+2y, the final x*x+y*y

  std::vector<double> data1(nn1, 1.5);

  // add the data vector to the state
  mystate.add("nodevars1", Portage::Entity_kind::NODE, &(data1[0]));

  std::vector<double> data2(nn1);

  for (int n = 0; n < nn1; ++n) {
    Wonton::Point<2> nodexy;
    meshWrapper.node_get_coordinates(n, &nodexy);
    data2[n] = 3 * nodexy[0] + nodexy[1];
  }

  // add the data vector to the state
  mystate.add("nodevars2", Portage::Entity_kind::NODE, &(data2[0]));

  std::vector<double> data3(nn1);

  for (int n = 0; n < nn1; ++n) {
    Wonton::Point<2> nodexy;
    meshWrapper.node_get_coordinates(n, &nodexy);
    data3[n] = nodexy[0]*nodexy[0] + nodexy[1]*nodexy[1];
  }

  // add the data vector to the state
  mystate.add("nodevars3", Portage::Entity_kind::NODE, &(data3[0]));

  // Create Quadfit calculater objects

  Portage::Limited_Quadfit<2, Portage::Entity_kind::NODE,
                           Wonton::Simple_Mesh_Wrapper,
                           Wonton::Simple_State_Wrapper>
      qfitcalc1(meshWrapper, stateWrapper, "nodevars1", 
                Portage::NOLIMITER, Portage::BND_NOLIMITER);
  Portage::Limited_Quadfit<2, Portage::Entity_kind::NODE,
                           Wonton::Simple_Mesh_Wrapper,
                           Wonton::Simple_State_Wrapper>
      qfitcalc2(meshWrapper, stateWrapper, "nodevars2", 
                Portage::NOLIMITER, Portage::BND_NOLIMITER);
  Portage::Limited_Quadfit<2, Portage::Entity_kind::NODE,
                           Wonton::Simple_Mesh_Wrapper,
                           Wonton::Simple_State_Wrapper>
      qfitcalc3(meshWrapper, stateWrapper, "nodevars3", 
                Portage::NOLIMITER, Portage::BND_NOLIMITER);
  Portage::Limited_Quadfit<2, Portage::Entity_kind::NODE,
                           Wonton::Simple_Mesh_Wrapper,
                           Wonton::Simple_State_Wrapper>
      qfitcalc4(meshWrapper, stateWrapper, "nodevars1",
                Portage::BARTH_JESPERSEN, Portage::BND_NOLIMITER);
  Portage::Limited_Quadfit<2, Portage::Entity_kind::NODE,
                           Wonton::Simple_Mesh_Wrapper,
                           Wonton::Simple_State_Wrapper>
      qfitcalc5(meshWrapper, stateWrapper, "nodevars2",
                Portage::BARTH_JESPERSEN, Portage::BND_NOLIMITER);
  Portage::Limited_Quadfit<2, Portage::Entity_kind::NODE,
                           Wonton::Simple_Mesh_Wrapper,
                           Wonton::Simple_State_Wrapper>
      qfitcalc6(meshWrapper, stateWrapper, "nodevars3",
                Portage::BARTH_JESPERSEN, Portage::BND_NOLIMITER);

  // Make sure we retrieved the correct quadfit value for each node
  // For field 1, it is a constant
  // For field 2, it is a linear function
  // For field 3, it is a quadratic function

  Wonton::Vector<5> qfit;

  for (int n = 0; n < nn1; ++n) {
    bool boundary_node = false;
    std::vector<int> nodecells;
    meshWrapper.node_get_cells(n, Portage::Entity_type::ALL, &nodecells);
    for (auto nc : nodecells) {
      std::vector<int> cfaces, cellfaceDirs;
      meshWrapper.cell_get_faces_and_dirs(nc, &cfaces, &cellfaceDirs);
      for (auto f : cfaces) {
        std::vector<int> fcells;
        meshWrapper.face_get_cells(f, Portage::Entity_type::ALL, &fcells);
        if (fcells.size() == 1) {
          boundary_node = true;
          break;
        }
      }
      if (boundary_node)
        break;
    }

    // if (!boundary_node) {
    // unlimited quadfit of constant function
    qfit = qfitcalc1(n);
    ASSERT_NEAR(0.0, qfit[0], 1.0e-10);
    ASSERT_NEAR(0.0, qfit[1], 1.0e-10);
    ASSERT_NEAR(0.0, qfit[2], 1.0e-10);
    ASSERT_NEAR(0.0, qfit[3], 1.0e-10);
    ASSERT_NEAR(0.0, qfit[4], 1.0e-10);

    // unlimited quadfit of linear function

    qfit = qfitcalc2(n);
    ASSERT_NEAR(3.0, qfit[0], 1.0e-10);
    ASSERT_NEAR(1.0, qfit[1], 1.0e-10);
    ASSERT_NEAR(0.0, qfit[2], 1.0e-10);
    ASSERT_NEAR(0.0, qfit[3], 1.0e-10);
    ASSERT_NEAR(0.0, qfit[4], 1.0e-10);

    // ulimited quadfit of quadratic function
    //
    // Avoid boundary nodes.  LimitedQuadfit won't crash at boundaries, but it
    // won't give the correct answer in the ASSERT statements for the quadratic
    // function test.

    if (!boundary_node) {
      qfit = qfitcalc3(n);
      Wonton::Point<2> nodexy;
      meshWrapper.node_get_coordinates(n, &nodexy);
      double nx = nodexy[0]; // x
      double ny = nodexy[1]; // y
      ASSERT_NEAR(2.0*nx, qfit[0], 1.0e-10); // partial of f wrt x
      ASSERT_NEAR(2.0*ny, qfit[1], 1.0e-10); // partial of f wrt y
      ASSERT_NEAR(1.0, qfit[2], 1.0e-10); // (2nd partial of f wrt x)/2
      ASSERT_NEAR(0.0, qfit[3], 1.0e-10); // (2nd partial of f wrt x and y)/2
      ASSERT_NEAR(1.0, qfit[4], 1.0e-10); // (2nd partial of f wrt y)/2
    }

    // limited quadfit of constant function

    qfit = qfitcalc4(n);
    ASSERT_NEAR(0.0, qfit[0], 1.0e-10);
    ASSERT_NEAR(0.0, qfit[1], 1.0e-10);
    ASSERT_NEAR(0.0, qfit[2], 1.0e-10);
    ASSERT_NEAR(0.0, qfit[3], 1.0e-10);
    ASSERT_NEAR(0.0, qfit[4], 1.0e-10);

    // unlimited quadfit of linear function
    if (!boundary_node) {
      qfit = qfitcalc5(n);
      ASSERT_NEAR(3.0, qfit[0], 1.0e-10);
      ASSERT_NEAR(1.0, qfit[1], 1.0e-10);
      ASSERT_NEAR(0.0, qfit[2], 1.0e-10);
      ASSERT_NEAR(0.0, qfit[3], 1.0e-10);
      ASSERT_NEAR(0.0, qfit[4], 1.0e-10);
    }

    // limited quadfit of quadratic function

    // Avoid boundary nodes.  LimitedQuadfit won't crash at boundaries, but it
    // won't give the correct answer in the ASSERT statements for the quadratic
    // function test.

    if (!boundary_node) {
      qfit = qfitcalc6(n);
      Wonton::Point<2> nodexy;
      meshWrapper.node_get_coordinates(n, &nodexy);
      double nx = nodexy[0]; // x
      double ny = nodexy[1]; // y
      ASSERT_NEAR(2.0*nx, qfit[0], 1.0e-10); // partial of f wrt x
      ASSERT_NEAR(2.0*ny, qfit[1], 1.0e-10); // partial of f wrt y
      ASSERT_NEAR(1.0, qfit[2], 1.0e-10); // (2nd partial of f wrt x)/2
      ASSERT_NEAR(0.0, qfit[3], 1.0e-10); // (2nd partials of f wrt x and y)/2
      ASSERT_NEAR(1.0, qfit[4], 1.0e-10); // (2nd partial of f wrt y)/2
    }
  }
}
