/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/


#include <iostream>

#include "gtest/gtest.h"

// Jali includes
#include "Mesh.hh"
#include "JaliState.h"
#include "MeshFactory.hh"


// wonton includes
#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"
#include "wonton/support/Vector.h"
#include "wonton/support/Matrix.h"
#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"

// portage includes
#include "portage/interpolate/interpolate_1st_order.h"
#include "portage/intersect/simple_intersect_for_tests.h"
#include "portage/intersect/dummy_interface_reconstructor.h"
#include "portage/support/portage.h"

// generic structure for invoking 1st & 2nd order interpolate (see last test)
// includes the include files for 1st and 2nd order interpolator
#include "portage/interpolate/interpolate_nth_order.h"

double TOL = 1e-12;

/// First order interpolation of constant cell-centered field in 2D

TEST(Interpolate_1st_Order_Vec, Cell_Ctr_Const_2D) {

  // Create meshes

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  std::shared_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4);
  std::shared_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 1.0, 1.0, 2, 2);
  std::shared_ptr<Jali::State> source_state = Jali::State::create(source_mesh);
  
  // Create mesh wrappers

  Wonton::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Wonton::Jali_Mesh_Wrapper targetMeshWrapper(*target_mesh);

  // count cells

  const int ncells_source =
      sourceMeshWrapper.num_owned_cells();
  const int ncells_target =
      targetMeshWrapper.num_owned_cells();

  // Define a state vector with constant value and add it to the source state

  std::vector<Wonton::Vector<2>> data(ncells_source, {1.25, 1.25});
  source_state->add("cellvars", source_mesh, Jali::Entity_kind::CELL,
                    Jali::Entity_type::ALL, &(data[0]));

  
  // Create state wrappers
  Wonton::Jali_State_Wrapper sourceStateWrapper(*source_state);
  
  // Gather the cell coordinates as Portage Points for source and target meshes
  // for intersection. The outer vector is the cells, the inner vector is the
  // points of the vertices of that cell.

  std::vector<std::vector<Wonton::Point<2>>>
      source_cell_coords(ncells_source);
  std::vector<std::vector<Wonton::Point<2>>>
      target_cell_coords(ncells_target);

  // Actually get the Wonton::Points

  for (int c = 0; c < ncells_source; ++c)
    sourceMeshWrapper.cell_get_coordinates(c, &(source_cell_coords[c]));
  for (int c = 0; c < ncells_target; ++c)
    targetMeshWrapper.cell_get_coordinates(c, &(target_cell_coords[c]));

  // Interpolate from source to target mesh using the independent calculation
  // in simple_intersect_for_tests.h

  std::vector<Wonton::Vector<2>> outvals(ncells_target);
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

  // use default tolerances
  Portage::NumericTolerances_t num_tols = Portage::DEFAULT_NUMERIC_TOLERANCES<2>;

  // Create Interpolation object

  Portage::Interpolate_1stOrder<2, Wonton::Entity_kind::CELL,
                                Wonton::Jali_Mesh_Wrapper,
                                Wonton::Jali_Mesh_Wrapper,
                                Wonton::Jali_State_Wrapper,
                                Wonton::Jali_State_Wrapper,
                                Wonton::Vector<2>
                                >
      interpolator(sourceMeshWrapper, targetMeshWrapper, sourceStateWrapper,
                   num_tols);

  interpolator.set_interpolation_variable("cellvars");


  Wonton::transform(targetMeshWrapper.begin(Wonton::Entity_kind::CELL),
                     targetMeshWrapper.end(Wonton::Entity_kind::CELL),
                     sources_and_weights.begin(),
                     outvals.begin(), interpolator);

  // Make sure we retrieved the correct value for each cell on the target
  const Wonton::Vector<2> stdval = data[0];
  for (int c = 0; c < ncells_target; ++c)
    for (int i = 0; i < 2; i++)
      ASSERT_NEAR(stdval[i], outvals[c][i], TOL);
}


/// First order interpolation of linear cell-centered field in 2D

TEST(Interpolate_1st_Order_Vec, Cell_Ctr_Lin_2D) {

  // Create simple meshes

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  std::shared_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4);
  std::shared_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 1.0, 1.0, 2, 2);
  std::shared_ptr<Jali::State> source_state = Jali::State::create(source_mesh);
  
  // Create mesh wrappers

  Wonton::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Wonton::Jali_Mesh_Wrapper targetMeshWrapper(*target_mesh);

  // count cells

  const int ncells_source =
      sourceMeshWrapper.num_owned_cells();
  const int ncells_target =
      targetMeshWrapper.num_owned_cells();

  // Define a state vector with linear value and add it to the source state

  std::vector<Wonton::Vector<2>> data;
  for (int c = 0; c < ncells_source; ++c) {
    Wonton::Point<2> cen;
    sourceMeshWrapper.cell_centroid(c, &cen);
    data.emplace_back(cen[0],cen[1]);
  }
  source_state->add("cellvars", source_mesh, Jali::Entity_kind::CELL,
                    Jali::Entity_type::ALL, &(data[0]));

  // Create state wrappers
  Wonton::Jali_State_Wrapper sourceStateWrapper(*source_state);

  // Gather the cell coordinates as Portage Points for source and target meshes
  // for intersection. The outer vector is the cells, the inner vector is the
  // points of the vertices of that cell.

  std::vector<std::vector<Wonton::Point<2>>>
      source_cell_coords(ncells_source);
  std::vector<std::vector<Wonton::Point<2>>>
      target_cell_coords(ncells_target);

  // Actually get the Wonton::Points

  for (int c = 0; c < ncells_source; ++c)
    sourceMeshWrapper.cell_get_coordinates(c, &(source_cell_coords[c]));
  for (int c = 0; c < ncells_target; ++c)
    targetMeshWrapper.cell_get_coordinates(c, &(target_cell_coords[c]));

  // Interpolate from source to target mesh using the independent calculation
  // in simple_intersect_for_tests.h

  std::vector<Wonton::Vector<2>> outvals(ncells_target);
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

  // use default tolerances
  Portage::NumericTolerances_t num_tols = Portage::DEFAULT_NUMERIC_TOLERANCES<2>;

  // Create Interpolation object

  Portage::Interpolate_1stOrder<2, Wonton::Entity_kind::CELL,
                                Wonton::Jali_Mesh_Wrapper,
                                Wonton::Jali_Mesh_Wrapper,
                                Wonton::Jali_State_Wrapper,
                                Wonton::Jali_State_Wrapper,
                                Wonton::Vector<2>>
      interpolator(sourceMeshWrapper, targetMeshWrapper, sourceStateWrapper,
                   num_tols);

  interpolator.set_interpolation_variable("cellvars");

  outvals[0] = interpolator(0,sources_and_weights[0]);
  
  std::transform(targetMeshWrapper.begin(Wonton::Entity_kind::CELL),
                     targetMeshWrapper.end(Wonton::Entity_kind::CELL),
                     sources_and_weights.begin(),
                     outvals.begin(), interpolator);

  // Make sure we retrieved the correct value for each cell on the target
  // NOTE: EVEN THOUGH 1ST ORDER INTERPOLATION ALGORITHM DOES NOT IN
  // GENERAL PRESERVE A LINEAR FIELD THE SPECIAL STRUCTURE OF THE SOURCE
  // AND TARGET MESHES ENSURES THAT THE LINEAR FIELD IS INTERPOLATED CORRECTLY
  // IN THIS TEST

  std::vector<Wonton::Vector<2>> stdvals;
  for (int c = 0; c < ncells_target; ++c) {
    Wonton::Point<2> cen;
    targetMeshWrapper.cell_centroid(c, &cen);
    stdvals.emplace_back(cen[0],cen[1]);
  }

  for (int c = 0; c < ncells_target; ++c)
    for (int i = 0; i < 2; i++)
      ASSERT_NEAR(stdvals[c][i], outvals[c][i], TOL);
}



/// First order interpolation of constant node-centered vector field in 2D

TEST(Interpolate_1st_Order_Vec, Node_Ctr_Const_2D) {

  // Create meshes

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  std::shared_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4);
  std::shared_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 1.0, 1.0, 5, 5);
  std::shared_ptr<Jali::State> source_state = Jali::State::create(source_mesh);

  // Create mesh wrappers

  Wonton::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Wonton::Jali_Mesh_Wrapper targetMeshWrapper(*target_mesh);

  // count nodes

  const int nnodes_source =
      sourceMeshWrapper.num_owned_nodes();
  const int nnodes_target =
      targetMeshWrapper.num_owned_nodes();

  // Define a state vector with constant value and add it to the source state

  std::vector<Wonton::Vector<2>> data(nnodes_source, {1.25, 1.25});
  source_state->add("nodevars", source_mesh, Jali::Entity_kind::NODE,
                    Jali::Entity_type::ALL, &(data[0]));

  // Create state wrapper

  Wonton::Jali_State_Wrapper sourceStateWrapper(*source_state);

  // Gather the cell coordinates as Portage Points for source and target meshes
  // for intersection. The outer vector is the cells, the inner vector is the
  // points of the vertices of that cell.

  std::vector<std::vector<Wonton::Point<2>>>
      source_dualcell_coords(nnodes_source);
  std::vector<std::vector<Wonton::Point<2>>>
      target_dualcell_coords(nnodes_target);

  // Actually get the Wonton::Points for the dual cells

  for (int n = 0; n < nnodes_source; ++n)
    sourceMeshWrapper.dual_cell_get_coordinates(n, &source_dualcell_coords[n]);
  for (int n = 0; n < nnodes_target; ++n)
    targetMeshWrapper.dual_cell_get_coordinates(n, &target_dualcell_coords[n]);


  // Interpolate from source to target mesh using the independent calculation
  // in simple_intersect_for_tests.h

  std::vector<Wonton::Vector<2>> outvals(nnodes_target);
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

  // use default tolerances
  Portage::NumericTolerances_t num_tols = Portage::DEFAULT_NUMERIC_TOLERANCES<2>;

  // Create Interpolation object

  Portage::Interpolate_1stOrder<2, Wonton::Entity_kind::NODE,
                                Wonton::Jali_Mesh_Wrapper,
                                Wonton::Jali_Mesh_Wrapper,
                                Wonton::Jali_State_Wrapper,
                                Wonton::Jali_State_Wrapper,
                                Wonton::Vector<2>>
      interpolator(sourceMeshWrapper, targetMeshWrapper, sourceStateWrapper,
                   num_tols);

  interpolator.set_interpolation_variable("nodevars");


  Wonton::transform(targetMeshWrapper.begin(Wonton::Entity_kind::NODE),
                     targetMeshWrapper.end(Wonton::Entity_kind::NODE),
                     sources_and_weights.begin(),
                     outvals.begin(), interpolator);

  // Make sure we retrieved the correct value for each cell on the target
  const Wonton::Vector<2> stdval = data[0];
  for (int c = 0; c < nnodes_target; ++c)
    for (int i = 0; i < 2; i++)
      ASSERT_NEAR(stdval[i], outvals[c][i], TOL);
}



// First order interpolation of linear node-centered "Matrix/Tensor" field in 3D
//
// This will not work until (1) Matrix class is templated on
// dimensions (2) can be initialized using T(0.0) operation (2) Jali
// state vectors disable operator<< for vectors whose elements don't
// support them (or Matrix class supports operator<< or we switch to
// the new simple state manager)
//
// But it is being left as an example for review

// TEST(Interpolate_1st_Order_mat, Node_Ctr_Const_3D) {

//   Jali::MeshFactory mf(MPI_COMM_WORLD);
//   std::shared_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 4, 4);
//   std::shared_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 5, 5, 5);
//   std::shared_ptr<Jali::State> source_state = Jali::State::create(source_mesh);
  
//   // Create mesh wrappers

//   Wonton::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
//   Wonton::Jali_Mesh_Wrapper targetMeshWrapper(*target_mesh);

//   // count nodes

//   const int nnodes_source =
//       sourceMeshWrapper.num_owned_nodes();
//   const int nnodes_target =
//       targetMeshWrapper.num_owned_nodes();

//   // Define a state vector with a constant value rectangular Matrix
//   // and add it to the source state

//   std::vector<Wonton::Matrix> data(nnodes_source, Wonton::Matrix(3, 2, 42));
//   source_state->add("nodevars", source_mesh, Jali::Entity_kind::NODE,
//                     Jali::Entity_type::ALL, &(data[0]));

//   // Create state wrapper

//   Wonton::Jali_State_Wrapper sourceStateWrapper(*source_state);

//   // Gather the cell coordinates as Portage Points for source and target meshes
//   // for intersection. The outer vector is the cells, the inner vector is the
//   // points of the vertices of that cell.

//   std::vector<std::vector<Wonton::Point<3>>>
//       source_dualcell_coords(nnodes_source);
//   std::vector<std::vector<Wonton::Point<3>>>
//       target_dualcell_coords(nnodes_target);

//   // Actually get the Wonton::Points for the dual cells

//   for (int n = 0; n < nnodes_source; ++n)
//     sourceMeshWrapper.dual_cell_get_coordinates(n, &source_dualcell_coords[n]);
//   for (int n = 0; n < nnodes_target; ++n)
//     targetMeshWrapper.dual_cell_get_coordinates(n, &target_dualcell_coords[n]);

//   // Interpolate from source to target mesh using the independent calculation
//   // in simple_intersect_for_tests.h

//   std::vector<Wonton::Matrix> outvals(nnodes_target); 
//   std::vector<std::vector<Portage::Weights_t>>
//       sources_and_weights(nnodes_target);

//   // Loop over target nodes

//   for (int c = 0; c < nnodes_target; ++c) {

//     std::vector<int> xcells;
//     std::vector<std::vector<double>> xwts;

//     // Compute the moments
//     // xcells is the source cell indices that intersect
//     // xwts is the moments vector for each cell that intersects

//     BOX_INTERSECT::intersection_moments<3>(target_dualcell_coords[c],
//                                         source_dualcell_coords,
//                                         &xcells, &xwts);

//     // Pack the results into a vector of true Portage::Weights_t

//     std::vector<Portage::Weights_t> wtsvec(xcells.size());
//     for (int i = 0; i < xcells.size(); ++i) {
//       wtsvec[i].entityID = xcells[i];
//       wtsvec[i].weights = xwts[i];
//     }

//     // Put the weights in final form

//     sources_and_weights[c] = wtsvec;
//   }

//   // Now do it the Portage way

//   // use default tolerances
//   Portage::NumericTolerances_t num_tols = Portage::DEFAULT_NUMERIC_TOLERANCES<3>;

//   // Create Interpolation object

//   Portage::Interpolate_1stOrder<3, Wonton::Entity_kind::NODE,
//                                 Wonton::Jali_Mesh_Wrapper,
//                                 Wonton::Jali_Mesh_Wrapper,
//                                 Wonton::Jali_State_Wrapper,
//                                 Wonton::Jali_State_Wrapper,
//                                 Wonton::Matrix>
//       interpolator(sourceMeshWrapper, targetMeshWrapper, sourceStateWrapper,
//                    num_tols);

//   interpolator.set_interpolation_variable("nodevars");


//   std::transform(targetMeshWrapper.begin(Wonton::Entity_kind::NODE),
//                      targetMeshWrapper.end(Wonton::Entity_kind::NODE),
//                      sources_and_weights.begin(),
//                      outvals.begin(), interpolator);

//   // Make sure we retrieved the correct value for each cell on the target
//   const Wonton::Matrix stdval = data[0];
//   for (int c = 0; c < nnodes_target; ++c)
//     for (int i = 0; i < 2; i++)
//       for (int j = 0; j < 3; j++)  {
//         double exp_mat_val = stdval[i][j];
//         double mat_val = outvals[c][i][j];
//         ASSERT_NEAR(exp_mat_val, mat_val, TOL);
//       }
// }
