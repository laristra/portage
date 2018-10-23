/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/


#include <iostream>

#include "gtest/gtest.h"

// portage includes
#include "portage/interpolate/interpolate_1st_order.h"
#include "portage/interpolate/test/simple_intersect_for_tests.h"
#include "portage/support/portage.h"

// wonton includes
#include "wonton/mesh/simple/simple_mesh.h"
#include "wonton/mesh/simple/simple_mesh_wrapper.h"
#include "wonton/state/simple/simple_state.h"
#include "wonton/state/simple/simple_state_wrapper.h"

double TOL = 1e-12;

/// First order interpolation of constant cell-centered field in 2D

TEST(Interpolate_1st_Order, Cell_Ctr_Const_2D) {

  // Create simple meshes

  std::shared_ptr<Wonton::Simple_Mesh> source_mesh =
    std::make_shared<Wonton::Simple_Mesh>(0.0, 0.0, 1.0, 1.0, 4, 4);
  std::shared_ptr<Wonton::Simple_Mesh> target_mesh =
    std::make_shared<Wonton::Simple_Mesh>(0.0, 0.0, 1.0, 1.0, 5, 5);

  // Create mesh wrappers

  Wonton::Simple_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Wonton::Simple_Mesh_Wrapper targetMeshWrapper(*target_mesh);

  // count cells

  const int ncells_source =
      sourceMeshWrapper.num_owned_cells();
  const int ncells_target =
      targetMeshWrapper.num_owned_cells();

  // Create a state object

  Wonton::Simple_State source_state(source_mesh);

  // Define a state vector with constant value and add it to the source state

  std::vector<double> data(ncells_source, 1.25);
  source_state.add("cellvars", Portage::Entity_kind::CELL, &(data[0]));

  // Create state wrapper

  Wonton::Simple_State_Wrapper sourceStateWrapper(source_state);

  // Gather the cell coordinates as Portage Points for source and target meshes
  // for intersection. The outer vector is the cells, the inner vector is the
  // points of the vertices of that cell.

  std::vector<std::vector<Portage::Point<2>>>
      source_cell_coords(ncells_source);
  std::vector<std::vector<Portage::Point<2>>>
      target_cell_coords(ncells_target);

  // Actually get the Portage::Points

  for (int c = 0; c < ncells_source; ++c)
    sourceMeshWrapper.cell_get_coordinates(c, &(source_cell_coords[c]));
  for (int c = 0; c < ncells_target; ++c)
    targetMeshWrapper.cell_get_coordinates(c, &(target_cell_coords[c]));

  // Interpolate from source to target mesh using the independent calculation
  // in simple_intersect_for_tests.h

  std::vector<double> outvals(ncells_target);
  std::vector<std::vector<Portage::Weights_t>>
      sources_and_weights(ncells_target);

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

  // Create Interpolation object

  Portage::Interpolate_1stOrder<2, Portage::Entity_kind::CELL,
                                Wonton::Simple_Mesh_Wrapper,
                                Wonton::Simple_Mesh_Wrapper,
                                Wonton::Simple_State_Wrapper>
      interpolater(sourceMeshWrapper, targetMeshWrapper, sourceStateWrapper);

  interpolater.set_interpolation_variable("cellvars");


  Portage::transform(targetMeshWrapper.begin(Portage::Entity_kind::CELL),
                     targetMeshWrapper.end(Portage::Entity_kind::CELL),
                     sources_and_weights.begin(),
                     outvals.begin(), interpolater);

  // Make sure we retrieved the correct value for each cell on the target
  const double stdval = data[0];
  for (int c = 0; c < ncells_target; ++c) ASSERT_NEAR(stdval, outvals[c], TOL);
}


/// First order interpolation of linear cell-centered field in 2D

TEST(Interpolate_1st_Order, Cell_Ctr_Lin_2D) {

  // Create simple meshes

  std::shared_ptr<Wonton::Simple_Mesh> source_mesh =
    std::make_shared<Wonton::Simple_Mesh>(0.0, 0.0, 1.0, 1.0, 4, 4);
  std::shared_ptr<Wonton::Simple_Mesh> target_mesh =
    std::make_shared<Wonton::Simple_Mesh>(0.0, 0.0, 1.0, 1.0, 2, 2);

  // Create mesh wrappers

  Wonton::Simple_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Wonton::Simple_Mesh_Wrapper targetMeshWrapper(*target_mesh);

  // count cells

  const int ncells_source =
      sourceMeshWrapper.num_owned_cells();
  const int ncells_target =
      targetMeshWrapper.num_owned_cells();

  // Create a state object

  Wonton::Simple_State source_state(source_mesh);

// Define a state vector with linear value and add it to the source state

  std::vector<double> data(ncells_source);
  for (int c = 0; c < ncells_source; ++c) {
    Portage::Point<2> cen;
    sourceMeshWrapper.cell_centroid(c, &cen);
    data[c] = cen[0]+cen[1];
  }
  source_state.add("cellvars", Portage::Entity_kind::CELL, &(data[0]));

  // Create state wrapper

  Wonton::Simple_State_Wrapper sourceStateWrapper(source_state);

  // Gather the cell coordinates as Portage Points for source and target meshes
  // for intersection. The outer vector is the cells, the inner vector is the
  // points of the vertices of that cell.

  std::vector<std::vector<Portage::Point<2>>>
      source_cell_coords(ncells_source);
  std::vector<std::vector<Portage::Point<2>>>
      target_cell_coords(ncells_target);

  // Actually get the Portage::Points

  for (int c = 0; c < ncells_source; ++c)
    sourceMeshWrapper.cell_get_coordinates(c, &(source_cell_coords[c]));
  for (int c = 0; c < ncells_target; ++c)
    targetMeshWrapper.cell_get_coordinates(c, &(target_cell_coords[c]));

  // Interpolate from source to target mesh using the independent calculation
  // in simple_intersect_for_tests.h

  std::vector<double> outvals(ncells_target);
  std::vector<std::vector<Portage::Weights_t>>
      sources_and_weights(ncells_target);

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

  // Create Interpolation object

  Portage::Interpolate_1stOrder<2, Portage::Entity_kind::CELL,
                                Wonton::Simple_Mesh_Wrapper,
                                Wonton::Simple_Mesh_Wrapper,
                                Wonton::Simple_State_Wrapper>
      interpolater(sourceMeshWrapper, targetMeshWrapper, sourceStateWrapper);

  interpolater.set_interpolation_variable("cellvars");


  Portage::transform(targetMeshWrapper.begin(Portage::Entity_kind::CELL),
                     targetMeshWrapper.end(Portage::Entity_kind::CELL),
                     sources_and_weights.begin(),
                     outvals.begin(), interpolater);

  // Make sure we retrieved the correct value for each cell on the target
  // NOTE: EVEN THOUGH 1ST ORDER INTERPOLATION ALGORITHM DOES NOT IN
  // GENERAL PRESERVE A LINEAR FIELD THE SPECIAL STRUCTURE OF THE SOURCE
  // AND TARGET MESHES ENSURES THAT THE LINEAR FIELD IS INTERPOLATED CORRECTLY
  // IN THIS TEST

  std::vector<double> stdvals(ncells_target);
  for (int c = 0; c < ncells_target; ++c) {
    Portage::Point<2> cen;
    targetMeshWrapper.cell_centroid(c, &cen);
    stdvals[c] = cen[0]+cen[1];
  }

  for (int c = 0; c < ncells_target; ++c)
    ASSERT_NEAR(stdvals[c], outvals[c], TOL);
}

/// First order interpolation of constant node-centered field in 2D

TEST(Interpolate_1st_Order, Node_Ctr_Const_2D) {

  // Create simple meshes

  std::shared_ptr<Wonton::Simple_Mesh> source_mesh =
    std::make_shared<Wonton::Simple_Mesh>(0.0, 0.0, 1.0, 1.0, 4, 4);
  std::shared_ptr<Wonton::Simple_Mesh> target_mesh =
    std::make_shared<Wonton::Simple_Mesh>(0.0, 0.0, 1.0, 1.0, 5, 5);

  // Create mesh wrappers

  Wonton::Simple_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Wonton::Simple_Mesh_Wrapper targetMeshWrapper(*target_mesh);

  // count nodes

  const int nnodes_source =
      sourceMeshWrapper.num_owned_nodes();
  const int nnodes_target =
      targetMeshWrapper.num_owned_nodes();

  // Create a state object

  Wonton::Simple_State source_state(source_mesh);

  // Define a state vector with constant value and add it to the source state

  std::vector<double> data(nnodes_source, 1.25);
  source_state.add("nodevars", Portage::Entity_kind::NODE, &(data[0]));

  // Create state wrapper

  Wonton::Simple_State_Wrapper sourceStateWrapper(source_state);

  // Gather the cell coordinates as Portage Points for source and target meshes
  // for intersection. The outer vector is the cells, the inner vector is the
  // points of the vertices of that cell.

  std::vector<std::vector<Portage::Point<2>>>
      source_dualcell_coords(nnodes_source);
  std::vector<std::vector<Portage::Point<2>>>
      target_dualcell_coords(nnodes_target);

  // Actually get the Portage::Points for the dual cells

  for (int n = 0; n < nnodes_source; ++n)
	  sourceMeshWrapper.dual_cell_get_coordinates(n, &source_dualcell_coords[n]);
  for (int n = 0; n < nnodes_target; ++n)
	  targetMeshWrapper.dual_cell_get_coordinates(n, &target_dualcell_coords[n]);


  // Interpolate from source to target mesh using the independent calculation
  // in simple_intersect_for_tests.h

  std::vector<double> outvals(nnodes_target);
  std::vector<std::vector<Portage::Weights_t>>
      sources_and_weights(nnodes_target);

  // Loop over target nodes

  for (int c = 0; c < nnodes_target; ++c) {

    std::vector<int> xcells;
    std::vector<std::vector<double>> xwts;

		// Compute the moments
		// xcells is the source cell indices that intersect
		// xwts is the moments vector for each cell that intersects

    BOX_INTERSECT::intersection_moments<2>(target_dualcell_coords[c],
                                        source_dualcell_coords,
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

  // Create Interpolation object

  Portage::Interpolate_1stOrder<2, Portage::Entity_kind::NODE,
                                Wonton::Simple_Mesh_Wrapper,
                                Wonton::Simple_Mesh_Wrapper,
                                Wonton::Simple_State_Wrapper>
      interpolater(sourceMeshWrapper, targetMeshWrapper, sourceStateWrapper);

  interpolater.set_interpolation_variable("nodevars");


  Portage::transform(targetMeshWrapper.begin(Portage::Entity_kind::NODE),
                     targetMeshWrapper.end(Portage::Entity_kind::NODE),
                     sources_and_weights.begin(),
                     outvals.begin(), interpolater);

  // Make sure we retrieved the correct value for each cell on the target
  const double stdval = data[0];
  for (int c = 0; c < nnodes_target; ++c) ASSERT_NEAR(stdval, outvals[c],TOL);


}



/// First order interpolation of constant cell-centered field in 3D

TEST(Interpolate_1st_Order, Cell_Ctr_Const_3D) {
 // Create simple meshes

  std::shared_ptr<Wonton::Simple_Mesh> source_mesh =
    std::make_shared<Wonton::Simple_Mesh>(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 4, 4);
  std::shared_ptr<Wonton::Simple_Mesh> target_mesh =
    std::make_shared<Wonton::Simple_Mesh>(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 5, 5, 5);

  // Create mesh wrappers

  Wonton::Simple_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Wonton::Simple_Mesh_Wrapper targetMeshWrapper(*target_mesh);

  // count cells

  const int ncells_source =
      sourceMeshWrapper.num_owned_cells();
  const int ncells_target =
      targetMeshWrapper.num_owned_cells();

  // Create a state object

  Wonton::Simple_State source_state(source_mesh);

  // Define a state vector with constant value and add it to the source state

  std::vector<double> data(ncells_source, 1.25);
  source_state.add("cellvars", Portage::Entity_kind::CELL, &(data[0]));

  // Create state wrapper

  Wonton::Simple_State_Wrapper sourceStateWrapper(source_state);

  // Gather the cell coordinates as Portage Points for source and target meshes
  // for intersection. The outer vector is the cells, the inner vector is the
  // points of the vertices of that cell.

  std::vector<std::vector<Portage::Point<3>>>
      source_cell_coords(ncells_source);
  std::vector<std::vector<Portage::Point<3>>>
      target_cell_coords(ncells_target);

  // Actually get the Portage::Points

  for (int c = 0; c < ncells_source; ++c)
    sourceMeshWrapper.cell_get_coordinates(c, &(source_cell_coords[c]));
  for (int c = 0; c < ncells_target; ++c)
    targetMeshWrapper.cell_get_coordinates(c, &(target_cell_coords[c]));

  // Interpolate from source to target mesh using the independent calculation
  // in simple_intersect_for_tests.h

  std::vector<double> outvals(ncells_target);
  std::vector<std::vector<Portage::Weights_t>>
      sources_and_weights(ncells_target);

  // Loop over target cells

  for (int c = 0; c < ncells_target; ++c) {

    std::vector<int> xcells;
    std::vector<std::vector<double>> xwts;

		// Compute the moments
		// xcells is the source cell indices that intersect
		// xwts is the moments vector for each cell that intersects

    BOX_INTERSECT::intersection_moments<3>(target_cell_coords[c],
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

  // Create Interpolation object

  Portage::Interpolate_1stOrder<3, Portage::Entity_kind::CELL,
                                Wonton::Simple_Mesh_Wrapper,
                                Wonton::Simple_Mesh_Wrapper,
                                Wonton::Simple_State_Wrapper>
      interpolater(sourceMeshWrapper, targetMeshWrapper, sourceStateWrapper);

  interpolater.set_interpolation_variable("cellvars");


  Portage::transform(targetMeshWrapper.begin(Portage::Entity_kind::CELL),
                     targetMeshWrapper.end(Portage::Entity_kind::CELL),
                     sources_and_weights.begin(),
                     outvals.begin(), interpolater);

  // Make sure we retrieved the correct value for each cell on the target
  const double stdval = data[0];
  for (int c = 0; c < ncells_target; ++c)
    ASSERT_NEAR(stdval, outvals[c], TOL);
}

/// First order interpolation of linear cell-centered field in 3D

TEST(Interpolate_1st_Order, Cell_Ctr_Lin_3D) {
  // Create simple meshes

  std::shared_ptr<Wonton::Simple_Mesh> source_mesh =
    std::make_shared<Wonton::Simple_Mesh>(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 4, 4);
  std::shared_ptr<Wonton::Simple_Mesh> target_mesh =
    std::make_shared<Wonton::Simple_Mesh>(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2);

  // Create mesh wrappers

  Wonton::Simple_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Wonton::Simple_Mesh_Wrapper targetMeshWrapper(*target_mesh);

  // count cells

  const int ncells_source =
      sourceMeshWrapper.num_owned_cells();
  const int ncells_target =
      targetMeshWrapper.num_owned_cells();

  // Create a state object

  Wonton::Simple_State source_state(source_mesh);

// Define a state vector with constant value and add it to the source state

  std::vector<double> data(ncells_source);
  for (int c = 0; c < ncells_source; ++c) {
    Portage::Point<3> cen;
    sourceMeshWrapper.cell_centroid(c, &cen);
    data[c] = cen[0]+cen[1];
  }
  source_state.add("cellvars", Portage::Entity_kind::CELL, &(data[0]));

  // Create state wrapper

  Wonton::Simple_State_Wrapper sourceStateWrapper(source_state);

  // Gather the cell coordinates as Portage Points for source and target meshes
  // for intersection. The outer vector is the cells, the inner vector is the
  // points of the vertices of that cell.

  std::vector<std::vector<Portage::Point<3>>>
      source_cell_coords(ncells_source);
  std::vector<std::vector<Portage::Point<3>>>
      target_cell_coords(ncells_target);

  // Actually get the Portage::Points

  for (int c = 0; c < ncells_source; ++c)
    sourceMeshWrapper.cell_get_coordinates(c, &(source_cell_coords[c]));
  for (int c = 0; c < ncells_target; ++c)
    targetMeshWrapper.cell_get_coordinates(c, &(target_cell_coords[c]));

  // Interpolate from source to target mesh using the independent calculation
  // in simple_intersect_for_tests.h

  std::vector<double> outvals(ncells_target);
  std::vector<std::vector<Portage::Weights_t>>
      sources_and_weights(ncells_target);

  // Loop over target cells

  for (int c = 0; c < ncells_target; ++c) {

    std::vector<int> xcells;
    std::vector<std::vector<double>> xwts;

		// Compute the moments
		// xcells is the source cell indices that intersect
		// xwts is the moments vector for each cell that intersects

    BOX_INTERSECT::intersection_moments<3>(target_cell_coords[c],
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

  // Create Interpolation object

  Portage::Interpolate_1stOrder<3, Portage::Entity_kind::CELL,
                                Wonton::Simple_Mesh_Wrapper,
                                Wonton::Simple_Mesh_Wrapper,
                                Wonton::Simple_State_Wrapper>
      interpolater(sourceMeshWrapper, targetMeshWrapper, sourceStateWrapper);

  interpolater.set_interpolation_variable("cellvars");


  Portage::transform(targetMeshWrapper.begin(Portage::Entity_kind::CELL),
                     targetMeshWrapper.end(Portage::Entity_kind::CELL),
                     sources_and_weights.begin(),
                     outvals.begin(), interpolater);

  // Make sure we retrieved the correct value for each cell on the target
  // NOTE: EVEN THOUGH 1ST ORDER INTERPOLATION ALGORITHM DOES NOT IN
  // GENERAL PRESERVE A LINEAR FIELD THE SPECIAL STRUCTURE OF THE SOURCE
  // AND TARGET MESHES ENSURES THAT THE LINEAR FIELD IS INTERPOLATED CORRECTLY
  // IN THIS TEST

  std::vector<double> stdvals(ncells_target);
  for (int c = 0; c < ncells_target; ++c) {
    Portage::Point<3> cen;
    targetMeshWrapper.cell_centroid(c, &cen);
    stdvals[c] = cen[0]+cen[1];
  }

  for (int c = 0; c < ncells_target; ++c)
    ASSERT_NEAR(stdvals[c], outvals[c], TOL);
}


/// First order interpolation of constant node-centered field in 3D

TEST(Interpolate_1st_Order, Node_Ctr_Const_3D) {

  // Create simple meshes

  std::shared_ptr<Wonton::Simple_Mesh> source_mesh =
    std::make_shared<Wonton::Simple_Mesh>(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 4, 4);
  std::shared_ptr<Wonton::Simple_Mesh> target_mesh =
    std::make_shared<Wonton::Simple_Mesh>(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 5, 5, 5);

  // Create mesh wrappers

  Wonton::Simple_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Wonton::Simple_Mesh_Wrapper targetMeshWrapper(*target_mesh);

  // count nodes

  const int nnodes_source =
      sourceMeshWrapper.num_owned_nodes();
  const int nnodes_target =
      targetMeshWrapper.num_owned_nodes();

  // Create a state object

  Wonton::Simple_State source_state(source_mesh);

  // Define a state vector with constant value and add it to the source state

  std::vector<double> data(nnodes_source, 1.25);
  source_state.add("nodevars", Portage::Entity_kind::NODE, &(data[0]));

  // Create state wrapper

  Wonton::Simple_State_Wrapper sourceStateWrapper(source_state);

  // Gather the cell coordinates as Portage Points for source and target meshes
  // for intersection. The outer vector is the cells, the inner vector is the
  // points of the vertices of that cell.

  std::vector<std::vector<Portage::Point<3>>>
      source_dualcell_coords(nnodes_source);
  std::vector<std::vector<Portage::Point<3>>>
      target_dualcell_coords(nnodes_target);

  // Actually get the Portage::Points for the dual cells

  for (int n = 0; n < nnodes_source; ++n)
	  sourceMeshWrapper.dual_cell_get_coordinates(n, &source_dualcell_coords[n]);
  for (int n = 0; n < nnodes_target; ++n)
	  targetMeshWrapper.dual_cell_get_coordinates(n, &target_dualcell_coords[n]);

  // Interpolate from source to target mesh using the independent calculation
  // in simple_intersect_for_tests.h

  std::vector<double> outvals(nnodes_target);
  std::vector<std::vector<Portage::Weights_t>>
      sources_and_weights(nnodes_target);

  // Loop over target nodes

  for (int c = 0; c < nnodes_target; ++c) {

    std::vector<int> xcells;
    std::vector<std::vector<double>> xwts;

		// Compute the moments
		// xcells is the source cell indices that intersect
		// xwts is the moments vector for each cell that intersects

    BOX_INTERSECT::intersection_moments<3>(target_dualcell_coords[c],
                                        source_dualcell_coords,
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

  // Create Interpolation object

  Portage::Interpolate_1stOrder<3, Portage::Entity_kind::NODE,
                                Wonton::Simple_Mesh_Wrapper,
                                Wonton::Simple_Mesh_Wrapper,
                                Wonton::Simple_State_Wrapper>
      interpolater(sourceMeshWrapper, targetMeshWrapper, sourceStateWrapper);

  interpolater.set_interpolation_variable("nodevars");


  Portage::transform(targetMeshWrapper.begin(Portage::Entity_kind::NODE),
                     targetMeshWrapper.end(Portage::Entity_kind::NODE),
                     sources_and_weights.begin(),
                     outvals.begin(), interpolater);

  // Make sure we retrieved the correct value for each cell on the target
  const double stdval = data[0];
  for (int c = 0; c < nnodes_target; ++c) ASSERT_NEAR(stdval, outvals[c],TOL);
}
