/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/


#include <iostream>

#include "gtest/gtest.h"

// portage includes
#include "portage/interpolate/interpolate_1st_order.h"
#include "portage/interpolate/interpolate_2nd_order.h"
#include "portage/intersect/simple_intersect_for_tests.h"
#include "portage/support/portage.h"

// wonton includes
#include "wonton/mesh/simple/simple_mesh.h"
#include "wonton/mesh/simple/simple_mesh_wrapper.h"
#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"
#include "wonton/support/Point.h"

// Jali includes
#include "JaliState.h"
#include "MeshFactory.hh"

double TOL = 1e-12;

/// Second order interpolation of cell-centered vector field with no
/// limiter in 2D

TEST(Interpolate_Vector_2nd, Cell_Ctr_Const_NOLIMITER_2D) {

  // Create simple meshes
  auto source_mesh = Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 1.0, 1.0, 4, 4);
  auto target_mesh = Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 1.0, 1.0, 5, 5);

  // Create mesh wrappers
  Wonton::Jali_Mesh_Wrapper source_mesh_wrapper(*source_mesh);
  Wonton::Jali_Mesh_Wrapper target_mesh_wrapper(*target_mesh);

  // count cells
  const int ncells_source = source_mesh_wrapper.num_owned_cells();
  const int ncells_target = target_mesh_wrapper.num_owned_cells();

  // Create a state object
  std::shared_ptr<Jali::State> source_state = Jali::State::create(source_mesh);

  // Define a state vector with constant value and add it to the source state
  std::vector<Wonton::Point<2>> data(ncells_source, Wonton::Point<2>(1.0, 2.0));
  source_state->add("cellvec", source_mesh, Jali::Entity_kind::CELL,
                    Jali::Entity_type::ALL, &(data[0]));

  // Create state wrapper
  Wonton::Jali_State_Wrapper source_state_wrapper(*source_state);

  // Gather the cell coordinates as Portage Points for source and target meshes
  // for intersection. The outer vector is the cells, the inner vector is the
  // points of the vertices of that cell.
  std::vector<std::vector<Wonton::Point<2>>> source_cell_coords(ncells_source);
  std::vector<std::vector<Wonton::Point<2>>> target_cell_coords(ncells_target);

  // Actually get the Wonton::Points
  for (int c = 0; c < ncells_source; ++c)
    source_mesh_wrapper.cell_get_coordinates(c, &(source_cell_coords[c]));
  for (int c = 0; c < ncells_target; ++c)
    target_mesh_wrapper.cell_get_coordinates(c, &(target_cell_coords[c]));

  // Interpolate from source to target mesh using the independent calculation
  // in simple_intersect_for_tests.h
  std::vector<Wonton::Point<2>> outvals(ncells_target);
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
    std::vector<Portage::Weights_t> wtsvec(xcells.size());
    for (int i = 0; i < xcells.size(); ++i) {
      wtsvec[i].entityID = xcells[i];
      wtsvec[i].weights = xwts[i];
    }

    // Put the weights in final form
    sources_and_weights[c] = wtsvec;
  }

  // Now do it the Portage way

  // use default tolerances
  Portage::NumericTolerances_t num_tols;
  num_tols.use_default();

  // Create Interpolation object
  // Portage::Interpolate_1stOrder<2, Wonton::Entity_kind::CELL,
  Portage::Interpolate_2ndOrder<2, Wonton::Entity_kind::CELL,
                                Wonton::Jali_Mesh_Wrapper,
                                Wonton::Jali_Mesh_Wrapper,
                                Wonton::Jali_State_Wrapper>
      interpolator(source_mesh_wrapper, target_mesh_wrapper, source_state_wrapper,
                   num_tols);

  interpolator.set_interpolation_variable("cellvec");

  Portage::transform(target_mesh_wrapper.begin(Wonton::Entity_kind::CELL),
                     target_mesh_wrapper.end(Wonton::Entity_kind::CELL),
                     sources_and_weights.begin(),
                     outvals.begin(), interpolator);

  // Make sure we retrieved the correct value for each cell on the target
  // const double stdval = data[0];
  // for (int c = 0; c < ncells_target; ++c)
  //   ASSERT_NEAR(stdval, outvals[c],TOL);
}

