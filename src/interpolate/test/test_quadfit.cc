/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include "portage/interpolate/quadfit.h"

#include <iostream>

#include "gtest/gtest.h"
#ifdef ENABLE_MPI
#include "mpi.h"
#endif

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliState.h"
#include "JaliStateVector.h"

#include "portage/support/portage.h"
#include "portage/support/Vector.h"
#include "portage/wonton/mesh/jali/jali_mesh_wrapper.h"
#include "portage/wonton/state/jali/jali_state_wrapper.h"
#include "portage/wonton/mesh/AuxMeshTopology.h"

/// Test quadfit computation for cell centered fields

TEST(Quadfit, Fields_Cell_Ctr) {

  // Make a 4x4 mesh

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  if (Jali::framework_available(Jali::MSTK))
    mf.framework(Jali::MSTK);
  mf.included_entities({Jali::Entity_kind::EDGE,
                        Jali::Entity_kind::FACE,
                        Jali::Entity_kind::WEDGE,
                        Jali::Entity_kind::CORNER});

  std::shared_ptr<Jali::Mesh> mesh1 = mf(0.0, 0.0, 1.0, 1.0, 4, 4);
  ASSERT_TRUE(mesh1 != nullptr);

  // Create a state object and add the first two vectors to it

  Jali::State mystate(mesh1);

  // Define three state vectors, one with constant value, one with
  // a linear function that is x+2y, and one with quadratic, x^2+y^2.

  int nc1 = mesh1->num_entities(Jali::Entity_kind::CELL,
                                Jali::Entity_type::PARALLEL_OWNED);

  // constant function
  std::vector<double> data1(nc1, 1.25);
  Jali::StateVector<double> myvec1("cellvars1", mesh1,
                                   Jali::Entity_kind::CELL,
                                   Jali::Entity_type::PARALLEL_OWNED,
                                   &(data1[0]));
  Jali::StateVector<double> &addvec1 = mystate.add(myvec1);

  // linear function
  std::vector<double> data2(nc1);
  for (int c = 0; c < nc1; c++) {
    JaliGeometry::Point ccen = mesh1->cell_centroid(c);
    data2[c] = ccen[0]+2*ccen[1];
  }

  Jali::StateVector<double> myvec2("cellvars2", mesh1,
                                   Jali::Entity_kind::CELL,
                                   Jali::Entity_type::PARALLEL_OWNED,
                                   &(data2[0]));
  Jali::StateVector<double> &addvec2 = mystate.add(myvec2);

  // quadratic function
  std::vector<double> data3(nc1);
  for (int c = 0; c < nc1; c++) {
    JaliGeometry::Point ccen = mesh1->cell_centroid(c);
    data3[c] = ccen[0]*ccen[0]+ccen[1]*ccen[1];
  }

  Jali::StateVector<double> myvec3("cellvars3", mesh1,
                                   Jali::Entity_kind::CELL,
                                   Jali::Entity_type::PARALLEL_OWNED,
                                   &(data3[0]));
  Jali::StateVector<double> &addvec3 = mystate.add(myvec3);

  Wonton::Jali_Mesh_Wrapper meshwrapper(*mesh1);
  Wonton::Jali_State_Wrapper statewrapper(mystate);

  // Create Quadfit objects

  Portage::Limited_Quadfit<2, Portage::CELL,
                           Wonton::Jali_Mesh_Wrapper,
                           Wonton::Jali_State_Wrapper>
      qfitcalc1(meshwrapper, statewrapper, "cellvars1", Portage::NOLIMITER);
  Portage::Limited_Quadfit<2, Portage::CELL,
                           Wonton::Jali_Mesh_Wrapper,
                           Wonton::Jali_State_Wrapper>
      qfitcalc2(meshwrapper, statewrapper, "cellvars2", Portage::NOLIMITER);
  Portage::Limited_Quadfit<2, Portage::CELL,
                           Wonton::Jali_Mesh_Wrapper,
                           Wonton::Jali_State_Wrapper>
      qfitcalc3(meshwrapper, statewrapper, "cellvars3", Portage::NOLIMITER);
  Portage::Limited_Quadfit<2, Portage::CELL,
                           Wonton::Jali_Mesh_Wrapper,
                           Wonton::Jali_State_Wrapper>
      qfitcalc4(meshwrapper, statewrapper, "cellvars1",
                Portage::BARTH_JESPERSEN);
  Portage::Limited_Quadfit<2, Portage::CELL,
                           Wonton::Jali_Mesh_Wrapper,
                           Wonton::Jali_State_Wrapper>
      qfitcalc5(meshwrapper, statewrapper, "cellvars2",
                Portage::BARTH_JESPERSEN);
  Portage::Limited_Quadfit<2, Portage::CELL,
                           Wonton::Jali_Mesh_Wrapper,
                           Wonton::Jali_State_Wrapper>
      qfitcalc6(meshwrapper, statewrapper, "cellvars3",
                Portage::BARTH_JESPERSEN);


  // Compute the quadfit for each of these fields

  Portage::Vector<5> qfit;

  // Verify the quadfit values
  // For field 1 (constant), it is is 0,0
  // For field 2 (x+2y), it is (1,2)
  // For field 3 x^2+y^2 it is (2x, 2y, 1, 0, 1)

  for (int c = 0; c < nc1; ++c) {
    // Create boundary_cell logical variable

    bool boundary_cell = false;
    std::vector<int> cfaces;
    mesh1->cell_get_faces(c, &cfaces);
    for (auto f : cfaces) {
      std::vector<int> fcells;
      mesh1->face_get_cells(f, Jali::Entity_type::ALL, &fcells);
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
      JaliGeometry::Point ccen = mesh1->cell_centroid(c);
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
      JaliGeometry::Point ccen = mesh1->cell_centroid(c);
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

  // Make a 3x3 mesh

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  if (Jali::framework_available(Jali::MSTK))
    mf.framework(Jali::MSTK);
  mf.included_entities({Jali::Entity_kind::EDGE,
                        Jali::Entity_kind::FACE,
                        Jali::Entity_kind::WEDGE,
                        Jali::Entity_kind::CORNER});

  std::shared_ptr<Jali::Mesh> mesh1 = mf(0.0, 0.0, 1.0, 1.0, 4, 4);
  ASSERT_TRUE(mesh1 != nullptr);

  // Create a state object and add the first two vectors to it

  Jali::State mystate(mesh1);

  // Define three state vectors, one with constant value, the other
  // with a linear function

  int nn1 = mesh1->num_entities(Jali::Entity_kind::NODE,
                                Jali::Entity_type::PARALLEL_OWNED);

  // constant function
  std::vector<double> data1(nn1, 1.5);

  Jali::StateVector<double> myvec1("nodevars1", mesh1,
                                   Jali::Entity_kind::NODE,
                                   Jali::Entity_type::PARALLEL_OWNED,
                                   &(data1[0]));
  Jali::StateVector<double> &addvec1 = mystate.add(myvec1);

  // linear function
  std::vector<double> data2(nn1);
  for (int n = 0; n < nn1; ++n) {
    JaliGeometry::Point nodexy;
    mesh1->node_get_coordinates(n, &nodexy);
    data2[n] = 3*nodexy[0]+nodexy[1];
  }
  Jali::StateVector<double> myvec2("nodevars2", mesh1,
                                   Jali::Entity_kind::NODE,
                                   Jali::Entity_type::PARALLEL_OWNED,
                                   &(data2[0]));
  Jali::StateVector<double> &addvec2 = mystate.add(myvec2);

  // quadratic function
  std::vector<double> data3(nn1);
  for (int n = 0; n < nn1; ++n) {
    JaliGeometry::Point nodexy;
    mesh1->node_get_coordinates(n, &nodexy);
    data3[n] = nodexy[0]*nodexy[0] + nodexy[1]*nodexy[1];
    //printf ("data3[%d] = %4.3f\n", n, data3[n]);
  }

  Jali::StateVector<double> myvec3("nodevars3", mesh1,
                                   Jali::Entity_kind::NODE,
                                   Jali::Entity_type::PARALLEL_OWNED,
                                   &(data3[0]));
  Jali::StateVector<double> &addvec3 = mystate.add(myvec3);

  // Create Quadfit calculater objects

  Wonton::Jali_Mesh_Wrapper meshwrapper(*mesh1);
  Wonton::Jali_State_Wrapper statewrapper(mystate);

  Portage::Limited_Quadfit<2, Portage::NODE,
                           Wonton::Jali_Mesh_Wrapper,
                           Wonton::Jali_State_Wrapper>
      qfitcalc1(meshwrapper, statewrapper, "nodevars1", Portage::NOLIMITER);
  Portage::Limited_Quadfit<2, Portage::NODE,
                           Wonton::Jali_Mesh_Wrapper,
                           Wonton::Jali_State_Wrapper>
      qfitcalc2(meshwrapper, statewrapper, "nodevars2", Portage::NOLIMITER);
  Portage::Limited_Quadfit<2, Portage::NODE,
                           Wonton::Jali_Mesh_Wrapper,
                           Wonton::Jali_State_Wrapper>
      qfitcalc3(meshwrapper, statewrapper, "nodevars3", Portage::NOLIMITER);
  Portage::Limited_Quadfit<2, Portage::NODE,
                           Wonton::Jali_Mesh_Wrapper,
                           Wonton::Jali_State_Wrapper>
      qfitcalc4(meshwrapper, statewrapper, "nodevars1",
                Portage::BARTH_JESPERSEN);
  Portage::Limited_Quadfit<2, Portage::NODE,
                           Wonton::Jali_Mesh_Wrapper,
                           Wonton::Jali_State_Wrapper>
      qfitcalc5(meshwrapper, statewrapper, "nodevars2",
                Portage::BARTH_JESPERSEN);
  Portage::Limited_Quadfit<2, Portage::NODE,
                           Wonton::Jali_Mesh_Wrapper,
                           Wonton::Jali_State_Wrapper>
      qfitcalc6(meshwrapper, statewrapper, "nodevars3",
                Portage::BARTH_JESPERSEN);


  // Make sure we retrieved the correct quadfit value for each node
  // For field 1, it is a constant
  // For field 2, it is a linear function
  // For field 3, it is a quadratic function

  Portage::Vector<5> qfit;

  for (int n = 0; n < nn1; ++n) {
    bool boundary_node = false;
    std::vector<int> nodecells;
    mesh1->node_get_cells(n, Jali::Entity_type::ALL, &nodecells);
    for (auto nc : nodecells) {
      std::vector<int> cfaces;
      mesh1->cell_get_faces(nc, &cfaces);
      for (auto f : cfaces) {
        std::vector<int> fcells;
        mesh1->face_get_cells(f, Jali::Entity_type::ALL, &fcells);
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
      JaliGeometry::Point nodexy;
      mesh1->node_get_coordinates(n, &nodexy);
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
      JaliGeometry::Point nodexy;
      mesh1->node_get_coordinates(n, &nodexy);
      double nx = nodexy[0]; // x
      double ny = nodexy[1]; // y
      ASSERT_NEAR(2.0*nx, qfit[0], 1.0e-10); // partial of f wrt x
      ASSERT_NEAR(2.0*ny, qfit[1], 1.0e-10); // partial of f wrt y
      ASSERT_NEAR(1.0, qfit[2], 1.0e-10); // (2nd partial of f wrt x)/2
      ASSERT_NEAR(0.0, qfit[3], 1.0e-10); // (2nd partials of f wrt x and y)/2
      ASSERT_NEAR(1.0, qfit[4], 1.0e-10); // (2nd partial of f wrt y)/2
    }
    //}
  }
}
