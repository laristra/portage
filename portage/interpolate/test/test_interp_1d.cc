/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/


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


TEST_P(Order2Test, SimpleMesh1D) {
  int itest = GetParam();

  // Create simple meshes
  std::shared_ptr<Wonton::Simple_Mesh> srcmesh = std::make_shared<Wonton::Simple_Mesh>(0.0, 1.0, 4);
  std::shared_ptr<Wonton::Simple_Mesh> tgtmesh = std::make_shared<Wonton::Simple_Mesh>(0.0, 1.0, 7);

  // Create mesh wrappers
  Wonton::Simple_Mesh_Wrapper srcmesh_wrapper(*srcmesh);
  Wonton::Simple_Mesh_Wrapper tgtmesh_wrapper(*tgtmesh);

  // count cells
  int ncells_src = srcmesh_wrapper.num_owned_cells();
  int ncells_tgt = tgtmesh_wrapper.num_owned_cells();

  // create a state object
  Wonton::Simple_State srcstate(srcmesh);
  Wonton::Simple_State tgtstate(tgtmesh);

  // Define a state vector with constant value and add it to the source state
  Wonton::Point<1> xc;
  std::vector<double> data(ncells_src);

  for (int c = 0; c < ncells_src; ++c) {
    if (itest == 1) {
      data[c] = 1.25;
    } else if (itest == 2) {
      srcmesh_wrapper.cell_centroid(c, &xc);
      data[c] = xc[0];
    }
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
    BOX_INTERSECT::intersection_moments<1>(tgt_cell_coords[c],
                                           src_cell_coords,
                                           &xcells, &xwts);

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
  double stdval;
  for (int c = 0; c < ncells_tgt; ++c) {
    if (itest == 1) {
      stdval = 1.25;
    } else if (itest == 2) {
      tgtmesh_wrapper.cell_centroid(c, &xc);
      stdval = xc[0];
    } 
    ASSERT_NEAR(stdval, outvals[c], TOL);
  }
}

INSTANTIATE_TEST_CASE_P(
  Order2TestAll,
  Order2Test,
  ::testing::Values(2));
