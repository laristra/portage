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
#include "wonton/mesh/simple/simple_mesh.h"
#include "wonton/mesh/simple/simple_mesh_wrapper.h"
#include "wonton/state/simple/simple_state.h"
#include "wonton/state/simple/simple_state_wrapper.h"

// portage includes
#include "portage/interpolate/interpolate_1st_order.h"
#include "portage/intersect/simple_intersect_for_tests.h"
#include "portage/support/portage.h"

double TOL = 1e-12;

class Order1Test : public ::testing::TestWithParam<int> {
 protected:
  int itest;
};

/// First order interpolation of constant cell-centered field in 2D

TEST_P(Order1Test, Cell_Constant_2D) {
  int itest = GetParam();
  int nx;
  Wonton::Vector<3> coefs;
  if (itest == 1) {
    nx = 5;
    coefs = { 1.25, 0.0, 0.0 }; 
  } else if (itest == 2) {
    nx = 2;
    coefs = { 0.0, 1.0, 1.0 }; 
  }

  // Create simple meshes
  std::shared_ptr<Wonton::Simple_Mesh> source_mesh =
    std::make_shared<Wonton::Simple_Mesh>(0.0, 0.0, 1.0, 1.0, 4, 4);
  std::shared_ptr<Wonton::Simple_Mesh> target_mesh =
    std::make_shared<Wonton::Simple_Mesh>(0.0, 0.0, 1.0, 1.0, nx, nx);

  // Create mesh wrappers
  Wonton::Simple_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Wonton::Simple_Mesh_Wrapper targetMeshWrapper(*target_mesh);

  const int ncells_source = sourceMeshWrapper.num_owned_cells();
  const int ncells_target = targetMeshWrapper.num_owned_cells();

  // Create a state object
  Wonton::Simple_State source_state(source_mesh);

  // Define a state vector with constant value and add it to the source state
  std::vector<double> data(ncells_source);
  for (int c = 0; c < ncells_source; ++c) {
    Wonton::Point<2> xc;
    sourceMeshWrapper.cell_centroid(c, &xc);
    data[c] = coefs[0] + xc[0] * coefs[1] + xc[1] * coefs[2];
  }
  source_state.add("cellvars", Wonton::Entity_kind::CELL, &(data[0]));

  // Create state wrapper
  Wonton::Simple_State_Wrapper sourceStateWrapper(source_state);

  // Gather the cell coordinates as Portage Points for source and target meshes
  // for intersection. The outer vector is the cells, the inner vector is the
  // points of the vertices of that cell.
  std::vector<std::vector<Wonton::Point<2>>> source_cell_coords(ncells_source);
  std::vector<std::vector<Wonton::Point<2>>> target_cell_coords(ncells_target);

  for (int c = 0; c < ncells_source; ++c)
    sourceMeshWrapper.cell_get_coordinates(c, &(source_cell_coords[c]));
  for (int c = 0; c < ncells_target; ++c)
    targetMeshWrapper.cell_get_coordinates(c, &(target_cell_coords[c]));

  // Interpolate from source to target mesh using the independent calculation
  // in simple_intersect_for_tests.h
  std::vector<double> outvals(ncells_target);
  std::vector<std::vector<Portage::Weights_t>> sources_and_weights(ncells_target);

  // Loop over target cells
  for (int c = 0; c < ncells_target; ++c) {
    std::vector<int> xcells;
    std::vector<std::vector<double>> xwts;

    // Compute the moments
    // xcells is the source cell indices that intersect
    // xwts is the moments vector for each cell that intersects
    BOX_INTERSECT::intersection_moments<2>(target_cell_coords[c],
                                           source_cell_coords,
                                           &xcells, &xwts);

    // Pack the results into a vector of true Portage::Weights_t
    int const num_intersect_cells = xcells.size();
    std::vector<Portage::Weights_t> wtsvec(num_intersect_cells);
    for (int i = 0; i < num_intersect_cells; ++i) {
      wtsvec[i].entityID = xcells[i];
      wtsvec[i].weights = xwts[i];
    }

    // Put the weights in final form
    sources_and_weights[c] = wtsvec;
  }

  // create interpolation object
  Portage::NumericTolerances_t num_tols = Portage::DEFAULT_NUMERIC_TOLERANCES<2>;
 
  Portage::Interpolate_1stOrder<2, Wonton::Entity_kind::CELL,
                                Wonton::Simple_Mesh_Wrapper,
                                Wonton::Simple_Mesh_Wrapper,
                                Wonton::Simple_State_Wrapper,
                                Wonton::Simple_State_Wrapper,
                                double,
                                Portage::DummyInterfaceReconstructor,
                                void,
                                void,
                                Wonton::CylindricalAxisymmetricCoordinates>
      interpolator(sourceMeshWrapper, targetMeshWrapper, sourceStateWrapper,
                   num_tols);

  interpolator.set_interpolation_variable("cellvars");

  Wonton::transform(targetMeshWrapper.begin(Wonton::Entity_kind::CELL),
                     targetMeshWrapper.end(Wonton::Entity_kind::CELL),
                     sources_and_weights.begin(),
                     outvals.begin(), interpolator);

  // Make sure we retrieved the correct value for each cell on the target
  double stdval;
  for (int c = 0; c < ncells_target; ++c) {
    Wonton::Point<2> xc;
    targetMeshWrapper.cell_centroid(c, &xc);
    stdval = coefs[0] + xc[0] * coefs[1] + xc[1] * coefs[2];
    ASSERT_NEAR(stdval, outvals[c], TOL);
  }
}


INSTANTIATE_TEST_CASE_P(
  Order1TestAll,
  Order1Test,
  ::testing::Values(1, 2));

