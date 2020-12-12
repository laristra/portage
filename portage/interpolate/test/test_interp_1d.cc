/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/


#include <cmath>
#include <iostream>
#include <memory>
#include <portage/search/search_kdtree.h>
#include <portage/intersect/intersect_r2d.h>
#include <portage/driver/mmdriver.h>

#include "gtest/gtest.h"

// wonton includes
#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"
#include "wonton/mesh/simple/simple_mesh.h"
#include "wonton/mesh/simple/simple_mesh_wrapper.h"
#include "wonton/state/simple/simple_state.h"
#include "wonton/state/simple/simple_state_wrapper.h"

// portage includes
#include "portage/interpolate/interpolate_2nd_order.h"
#include "portage/intersect/simple_intersect_for_tests.h"
#include "portage/support/portage.h"
#include "portage/driver/coredriver.h"

double TOL = 1e-12;

// Unlimited 2nd-order interpolation of polynomial cell-centered fields
class Order2Test : public ::testing::TestWithParam<int> {
 protected:
  int itest;
};


// Tests 1-2 Cartesian coordinates for constant and linear functions
// Tests 3-4 Cylindrical coordinates for constant and linear functions

TEST_P(Order2Test, SimpleMesh1D) {
  int itest = GetParam();

  // Create simple meshes
  std::shared_ptr<Wonton::Simple_Mesh> srcmesh = std::make_shared<Wonton::Simple_Mesh>(0.0, 1.0, 4);
  std::shared_ptr<Wonton::Simple_Mesh> tgtmesh = std::make_shared<Wonton::Simple_Mesh>(0.0, 1.0, 8);

  // Create mesh wrappers
  auto CoordSys = Wonton::CoordSysType::Cartesian;
  if (itest == 3 || itest == 4) CoordSys = Wonton::CoordSysType::CylindricalRadial;
  if (itest == 5 || itest == 6) CoordSys = Wonton::CoordSysType::SphericalRadial;

  Wonton::Simple_Mesh_Wrapper srcmesh_wrapper(*srcmesh, true, true, true, CoordSys);
  Wonton::Simple_Mesh_Wrapper tgtmesh_wrapper(*tgtmesh, true, true, true, CoordSys);

  // count cells
  int ncells_src = srcmesh_wrapper.num_owned_cells();
  int ncells_tgt = tgtmesh_wrapper.num_owned_cells();

  // create a state object
  Wonton::Simple_State srcstate(srcmesh);
  Wonton::Simple_State tgtstate(tgtmesh);

  // Define a state vector with constant value and add it to the source state
  double mass0(0.0);
  Wonton::Point<1> xc;
  std::vector<double> data(ncells_src);

  for (int c = 0; c < ncells_src; ++c) {
    double vol = srcmesh_wrapper.cell_volume(c);
    if (itest == 1) {
      data[c] = 1.25;
    } else if (itest == 2 || itest == 3) {
      srcmesh_wrapper.cell_centroid(c, &xc);
      data[c] = xc[0];
    } else if (itest == 4 || itest == 5) {
      Wonton::Point<1> a, b;
      std::vector<int> nodes;
      srcmesh_wrapper.cell_get_nodes(c, &nodes);
      srcmesh_wrapper.node_get_coordinates(nodes[0], &a);
      srcmesh_wrapper.node_get_coordinates(nodes[1], &b);
      data[c] = fabs(std::pow(b[0], 3) - std::pow(a[0], 3)) / (3 * vol);
    } else if (itest == 6) {
      Wonton::Point<1> a, b;
      std::vector<int> nodes;
      srcmesh_wrapper.cell_get_nodes(c, &nodes);
      srcmesh_wrapper.node_get_coordinates(nodes[0], &a);
      srcmesh_wrapper.node_get_coordinates(nodes[1], &b);
      data[c] = fabs(std::pow(b[0], 4) - std::pow(a[0], 4)) / (4 * vol);
    }
    mass0 += data[c] * vol;
  }
  srcstate.add("cellvars", Wonton::Entity_kind::CELL, &(data[0]));

  // create state wrapper
  Wonton::Simple_State_Wrapper srcstate_wrapper(srcstate);
  Wonton::Simple_State_Wrapper tgtstate_wrapper(tgtstate);

  // get Wonton::Points
  std::vector<std::vector<Wonton::Point<1>>> src_cell_coords(ncells_src);
  std::vector<std::vector<Wonton::Point<1>>> tgt_cell_coords(ncells_tgt);

  for (int c = 0; c < ncells_src; ++c)
    srcmesh_wrapper.cell_get_coordinates(c, &(src_cell_coords[c]));
  for (int c = 0; c < ncells_tgt; ++c) 
    tgtmesh_wrapper.cell_get_coordinates(c, &(tgt_cell_coords[c]));

  // Interpolate from source to target mesh using the independent 
  // of cell intersections in simple_intersect_for_tests.h
  std::vector<double> outvals(ncells_tgt);
  std::vector<std::vector<Portage::Weights_t>> sources_and_weights(ncells_tgt);

  for (int c = 0; c < ncells_tgt; ++c) {
    std::vector<int> xcells;
    std::vector<std::vector<double>> xwts;

    // Compute the moments
    if (itest == 1 || itest == 2) {
      BOX_INTERSECT::intersection_moments<1>(
          tgt_cell_coords[c], src_cell_coords, &xcells, &xwts);
    } else if (itest == 3 || itest == 4) {
      BOX_INTERSECT::intersection_moments<1, Wonton::CylindricalRadialCoordinates>(
          tgt_cell_coords[c], src_cell_coords, &xcells, &xwts);
    } else if (itest == 5 || itest == 6) {
      BOX_INTERSECT::intersection_moments<1, Wonton::SphericalRadialCoordinates>(
          tgt_cell_coords[c], src_cell_coords, &xcells, &xwts);
    }

    // Pack the results into a vector of true Portage::Weights_t
    int num_intersect_cells = xcells.size();
    std::vector<Portage::Weights_t> wtsvec(num_intersect_cells);

    for (int i = 0; i < num_intersect_cells; ++i) {
      wtsvec[i].entityID = xcells[i];
      wtsvec[i].weights = xwts[i];
    }

    // Put the weights in final form
    sources_and_weights[c] = wtsvec;
  }

  // use default tolerances
  Portage::NumericTolerances_t num_tols = Portage::DEFAULT_NUMERIC_TOLERANCES<1>;

  // compute gradient field to pass to the interpolator
  using Driver = Portage::CoreDriver<1, Wonton::Entity_kind::CELL,
                                     Wonton::Simple_Mesh_Wrapper,
                                     Wonton::Simple_State_Wrapper>;

  Driver driver(srcmesh_wrapper, srcstate_wrapper,
                tgtmesh_wrapper, tgtstate_wrapper);

  auto gradients = driver.compute_source_gradient("cellvars");

  // Create Interpolation object
  Portage::Interpolate_2ndOrder<1, Wonton::Entity_kind::CELL,
                                Wonton::Simple_Mesh_Wrapper,
                                Wonton::Simple_Mesh_Wrapper,
                                Wonton::Simple_State_Wrapper,
                                Wonton::Simple_State_Wrapper,
                                double>
      interpolator(srcmesh_wrapper, tgtmesh_wrapper, srcstate_wrapper,
                   num_tols);

  interpolator.set_interpolation_variable("cellvars", &gradients);

  Wonton::transform(tgtmesh_wrapper.begin(Wonton::Entity_kind::CELL),
                    tgtmesh_wrapper.end(Wonton::Entity_kind::CELL),
                    sources_and_weights.begin(),
                    outvals.begin(), interpolator);

  // Make sure we retrieved the correct value for each cell on the target
  double stdval(0.0), mass1(0.0);
  for (int c = 0; c < ncells_tgt; ++c) {
    double vol = tgtmesh_wrapper.cell_volume(c);
    if (itest == 1) {
      stdval = 1.25;
    } else if (itest == 2 || itest == 3) {
      tgtmesh_wrapper.cell_centroid(c, &xc);
      stdval = xc[0];
    } else if (itest == 4 || itest == 5) {
      Wonton::Point<1> a, b;
      std::vector<int> nodes;
      tgtmesh_wrapper.cell_get_nodes(c, &nodes);
      tgtmesh_wrapper.node_get_coordinates(nodes[0], &a);
      tgtmesh_wrapper.node_get_coordinates(nodes[1], &b);
      stdval = fabs(std::pow(b[0], 3) - std::pow(a[0], 3)) / (3 * vol);
    } else if (itest == 6) {
      Wonton::Point<1> a, b;
      std::vector<int> nodes;
      tgtmesh_wrapper.cell_get_nodes(c, &nodes);
      tgtmesh_wrapper.node_get_coordinates(nodes[0], &a);
      tgtmesh_wrapper.node_get_coordinates(nodes[1], &b);
      stdval = fabs(std::pow(b[0], 4) - std::pow(a[0], 4)) / (4 * vol);
    } 
    // there is no clear definition of linearity preservation in curvilinear coordinates
    // instaead, we use a mass conservation test.
    if (itest < 4) {
      ASSERT_NEAR(stdval, outvals[c], TOL);
    }
    mass1 += stdval * vol;
  }

  ASSERT_NEAR(mass0, mass1, TOL);
}

INSTANTIATE_TEST_CASE_P(
  Order2TestAll,
  Order2Test,
  ::testing::Values(1,2,3,4,5,6));
