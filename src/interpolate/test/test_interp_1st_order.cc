/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include "portage/driver/driver.h"
#include "portage/interpolate/interpolate_1st_order.h"

#include <iostream>
#include <memory>

#include "gtest/gtest.h"
#include "mpi.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "FrameworkTraits.hh"
#include "MeshFramework.hh"
#include "JaliState.h"
#include "JaliStateVector.h"

#include "portage/support/portage.h"
#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"
#include "portage/wrappers/state/jali/jali_state_wrapper.h"
#include "portage/interpolate/test/simple_intersect_for_tests.h"


double TOL = 1e-12;


/// First order interpolation of constant cell-centered field in 2D

TEST(Interpolate_1st_Order, Cell_Ctr_Const_2D) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);

  std::shared_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4);
  std::shared_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 1.0, 1.0, 5, 5);

  const int ncells_source =
      source_mesh->num_entities(Jali::Entity_kind::CELL,
                                Jali::Entity_type::PARALLEL_OWNED);
  const int ncells_target =
      target_mesh->num_entities(Jali::Entity_kind::CELL,
                                Jali::Entity_type::PARALLEL_OWNED);

  // Create a state object

  Jali::State source_state(source_mesh);

  // Define a state vectors with constant value

  std::vector<double> data(ncells_source, 1.25);
  Jali::StateVector<double> myvec("cellvars", source_mesh,
                                  Jali::Entity_kind::CELL,
                                  Jali::Entity_type::PARALLEL_OWNED,
                                  &(data[0]));
  source_state.add(myvec);

  Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Portage::Jali_Mesh_Wrapper targetMeshWrapper(*target_mesh);
  Portage::Jali_State_Wrapper sourceStateWrapper(source_state);

  // Gather the cell coordinates for source and target meshes for intersection

  std::vector<std::vector<JaliGeometry::Point>>
      source_cell_coords(ncells_source);
  std::vector<std::vector<JaliGeometry::Point>>
      target_cell_coords(ncells_target);

  for (int c = 0; c < ncells_source; ++c)
    source_mesh->cell_get_coordinates(c, &(source_cell_coords[c]));
  for (int c = 0; c < ncells_target; ++c)
    target_mesh->cell_get_coordinates(c, &(target_cell_coords[c]));

  // Interpolate from source to target mesh

  std::vector<double> outvals(ncells_target);
  std::vector<std::vector<Portage::Weights_t>>
      sources_and_weights(ncells_target);
  for (int c = 0; c < ncells_target; ++c) {
    std::vector<int> xcells;
    std::vector<std::vector<double>> xwts;

    BOX_INTERSECT::intersection_moments(target_cell_coords[c],
                                        source_cell_coords,
                                        &xcells, &xwts);

    std::vector<Portage::Weights_t> wtsvec(xcells.size());
    for (int i = 0; i < xcells.size(); ++i) {
      wtsvec[i].entityID = xcells[i];
      wtsvec[i].weights = xwts[i];
    }
    sources_and_weights[c] = wtsvec;
  }

  // Create Interpolation object

  Portage::Interpolate_1stOrder<Portage::Jali_Mesh_Wrapper,
                                Portage::Jali_Mesh_Wrapper,
                                Portage::Jali_State_Wrapper,
                                Portage::CELL, 2>
      interpolater(sourceMeshWrapper, targetMeshWrapper, sourceStateWrapper);

  interpolater.set_interpolation_variable("cellvars");

  Jali::Entity_ID_List const& targetcells =
      target_mesh->cells<Jali::Entity_type::ALL>();
  Portage::transform(targetcells.begin(), targetcells.end(),
                     sources_and_weights.begin(),
                     outvals.begin(), interpolater);

  // Make sure we retrieved the correct value for each cell on the target
  const double stdval = data[0];
  for (int c = 0; c < ncells_target; ++c)
    ASSERT_DOUBLE_EQ(stdval, outvals[c]);
}


/// First order interpolation of linear cell-centered field in 2D

TEST(Interpolate_1st_Order, Cell_Ctr_Lin_2D) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);

  std::shared_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4);
  std::shared_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 1.0, 1.0, 2, 2);

  const int ncells_source =
      source_mesh->num_entities(Jali::Entity_kind::CELL,
                                Jali::Entity_type::PARALLEL_OWNED);
  const int ncells_target =
      target_mesh->num_entities(Jali::Entity_kind::CELL,
                                Jali::Entity_type::PARALLEL_OWNED);

  // Create a state object

  Jali::State source_state(source_mesh);

  // Define a state vector representing a linear function

  std::vector<double> data(ncells_source);
  for (int c = 0; c < ncells_source; ++c) {
    JaliGeometry::Point cen;
    source_mesh->cell_centroid(c, &cen);
    data[c] = cen[0]+cen[1];
  }
  Jali::StateVector<double> myvec("cellvars", source_mesh,
                                  Jali::Entity_kind::CELL,
                                  Jali::Entity_type::PARALLEL_OWNED,
                                  &(data[0]));
  source_state.add(myvec);

  // Create Interpolation objects

  Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Portage::Jali_Mesh_Wrapper targetMeshWrapper(*target_mesh);
  Portage::Jali_State_Wrapper sourceStateWrapper(source_state);

  // Gather the cell coordinates for source and target meshes for intersection

  std::vector<std::vector<JaliGeometry::Point>>
      source_cell_coords(ncells_source);
  std::vector<std::vector<JaliGeometry::Point>>
      target_cell_coords(ncells_target);

  for (int c = 0; c < ncells_source; ++c)
    source_mesh->cell_get_coordinates(c, &(source_cell_coords[c]));
  for (int c = 0; c < ncells_target; ++c)
    target_mesh->cell_get_coordinates(c, &(target_cell_coords[c]));

  // Interpolate from source to target mesh

  std::vector<double> outvals(ncells_target);
  std::vector<std::vector<Portage::Weights_t>>
      sources_and_weights(ncells_target);

  for (int c = 0; c < ncells_target; ++c) {
    std::vector<int> xcells;
    std::vector<std::vector<double>> xwts;

    BOX_INTERSECT::intersection_moments(target_cell_coords[c],
                                        source_cell_coords,
                                        &xcells, &xwts);

    std::vector<Portage::Weights_t> wtsvec(xcells.size());
    for (int i = 0; i < xcells.size(); ++i) {
      wtsvec[i].entityID = xcells[i];
      wtsvec[i].weights = xwts[i];
    }
    sources_and_weights[c] = wtsvec;
  }

  Portage::Interpolate_1stOrder<Portage::Jali_Mesh_Wrapper,
                                Portage::Jali_Mesh_Wrapper,
                                Portage::Jali_State_Wrapper,
                                Portage::CELL, 2>
      interpolater(sourceMeshWrapper, targetMeshWrapper, sourceStateWrapper);

  interpolater.set_interpolation_variable("cellvars");

  Jali::Entity_ID_List const& targetcells =
      target_mesh->cells<Jali::Entity_type::ALL>();
  Portage::transform(targetcells.begin(), targetcells.end(),
                     sources_and_weights.begin(),
                     outvals.begin(), interpolater);

  // Make sure we retrieved the correct value for each cell on the target
  // NOTE: EVEN THOUGH 1ST ORDER INTERPOLATION ALGORITHM DOES NOT IN
  // GENERAL PRESERVE A LINEAR FIELD THE SPECIAL STRUCTURE OF THE SOURCE
  // AND TARGET MESHES ENSURES THAT THE LINEAR FIELD IS INTERPOLATED CORRECTLY
  // IN THIS TEST

  std::vector<double> stdvals(ncells_target);
  for (int c = 0; c < ncells_target; ++c) {
    JaliGeometry::Point cen;
    target_mesh->cell_centroid(c, &cen);
    stdvals[c] = cen[0]+cen[1];
  }

  for (int c = 0; c < ncells_target; ++c)
    ASSERT_DOUBLE_EQ(stdvals[c], outvals[c]);
}


/// First order interpolation of constant node-centered field in 2D

TEST(Interpolate_1st_Order, Node_Ctr_Const_2D) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);
  mf.included_entities({Jali::Entity_kind::EDGE,
                        Jali::Entity_kind::FACE,
                        Jali::Entity_kind::WEDGE,
                        Jali::Entity_kind::CORNER});

  std::shared_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4);
  std::shared_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 1.0, 1.0, 5, 5);

  const int nnodes_source =
      source_mesh->num_entities(Jali::Entity_kind::NODE,
                                Jali::Entity_type::PARALLEL_OWNED);
  const int nnodes_target =
      target_mesh->num_entities(Jali::Entity_kind::NODE,
                                Jali::Entity_type::PARALLEL_OWNED);


  // Create a state object and add the first two vectors to it

  Jali::State source_state(source_mesh);

  // Define two state vectors, one with constant value, the other
  // with a linear function

  std::vector<double> data(nnodes_source, 1.5);
  Jali::StateVector<double> myvec("nodevars", source_mesh,
                                  Jali::Entity_kind::NODE,
                                  Jali::Entity_type::PARALLEL_OWNED,
                                  &(data[0]));
  source_state.add(myvec);

  // Create Interpolation objects

  Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Portage::Jali_Mesh_Wrapper targetMeshWrapper(*target_mesh);
  Portage::Jali_State_Wrapper sourceStateWrapper(source_state);

  Portage::MeshWrapperDual<Portage::Jali_Mesh_Wrapper>
      sourceDualWrapper(sourceMeshWrapper);
  Portage::MeshWrapperDual<Portage::Jali_Mesh_Wrapper>
      targetDualWrapper(targetMeshWrapper);

  // Gather the dual cell coordinates for source and target meshes for
  // intersection

  std::vector<std::vector<JaliGeometry::Point>>
      source_dualcell_coords(nnodes_source);
  std::vector<std::vector<JaliGeometry::Point>>
      target_dualcell_coords(nnodes_target);

  // Because the meshes are rectangular we can get away with examining
  // the coordinates of the corners instead of the wedges

  // Also, because we will use only the bounding box of each cell to
  // do the search and intersection we can get away with adding all
  // the coordinates of the corners to list including duplicates

  std::vector<std::vector<Portage::Weights_t>>
      sources_and_weights(nnodes_target);

  for (int n = 0; n < nnodes_source; ++n) {
    std::vector<JaliGeometry::Point> dualcoords;
    std::vector<int> corners;
    source_mesh->node_get_corners(n, Jali::Entity_type::ALL, &corners);

    for (auto cn : corners) {
      std::vector<JaliGeometry::Point> cncoords;
      source_mesh->corner_get_coordinates(cn, &cncoords);
      for (auto coord : cncoords)
        source_dualcell_coords[n].push_back(coord);
    }
  }

  for (int n = 0; n < nnodes_target; ++n) {
    std::vector<JaliGeometry::Point> dualcoords;
    std::vector<int> corners;
    target_mesh->node_get_corners(n, Jali::Entity_type::ALL, &corners);

    for (auto cn : corners) {
      std::vector<JaliGeometry::Point> cncoords;
      target_mesh->corner_get_coordinates(cn, &cncoords);
      for (auto coord : cncoords)
        target_dualcell_coords[n].push_back(coord);
    }
  }

  for (int n = 0; n < nnodes_target; ++n) {
    std::vector<int> xcells;
    std::vector<std::vector<double>> xwts;

    BOX_INTERSECT::intersection_moments(target_dualcell_coords[n],
                                        source_dualcell_coords,
                                        &xcells, &xwts);


    std::vector<Portage::Weights_t> wtsvec(xcells.size());
    for (int i = 0; i < xcells.size(); ++i) {
      wtsvec[i].entityID = xcells[i];
      wtsvec[i].weights = xwts[i];
    }
    sources_and_weights[n] = wtsvec;
  }

  // Interpolate from source to target mesh

  std::vector<double> outvals(nnodes_target);

  Portage::Interpolate_1stOrder<
    Portage::MeshWrapperDual<Portage::Jali_Mesh_Wrapper>,
    Portage::MeshWrapperDual<Portage::Jali_Mesh_Wrapper>,
    Portage::Jali_State_Wrapper,
    Portage::NODE, 2>
      interpolater(sourceDualWrapper, targetDualWrapper, sourceStateWrapper);

  interpolater.set_interpolation_variable("nodevars");

  Jali::Entity_ID_List const& targetnodes =
      target_mesh->nodes<Jali::Entity_type::ALL>();

  Portage::transform(targetnodes.begin(), targetnodes.end(),
                     sources_and_weights.begin(),
                     outvals.begin(), interpolater);

  // Make sure we retrieved the correct value for each cell on the target

  const double stdval = data[0];
  for (int n = 0; n < nnodes_target; ++n)
    ASSERT_NEAR(stdval, outvals[n], TOL);
}


/// First order interpolation of constant cell-centered field in 3D

TEST(Interpolate_1st_Order, Cell_Ctr_Const_3D) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);

  std::shared_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                               4, 4, 4);
  std::shared_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                               5, 5, 5);

  const int ncells_source =
      source_mesh->num_entities(Jali::Entity_kind::CELL,
                                Jali::Entity_type::PARALLEL_OWNED);
  const int ncells_target =
      target_mesh->num_entities(Jali::Entity_kind::CELL,
                                Jali::Entity_type::PARALLEL_OWNED);

  // Create a state object

  Jali::State source_state(source_mesh);

  // Define a state vectors with constant value

  std::vector<double> data(ncells_source, 1.25);
  Jali::StateVector<double> myvec("cellvars", source_mesh,
                                  Jali::Entity_kind::CELL,
                                  Jali::Entity_type::PARALLEL_OWNED,
                                  &(data[0]));
  source_state.add(myvec);

  Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Portage::Jali_Mesh_Wrapper targetMeshWrapper(*target_mesh);
  Portage::Jali_State_Wrapper sourceStateWrapper(source_state);

  // Gather the cell coordinates for source and target meshes for intersection

  std::vector<std::vector<JaliGeometry::Point>>
      source_cell_coords(ncells_source);
  std::vector<std::vector<JaliGeometry::Point>>
      target_cell_coords(ncells_target);

  for (int c = 0; c < ncells_source; ++c)
    source_mesh->cell_get_coordinates(c, &(source_cell_coords[c]));
  for (int c = 0; c < ncells_target; ++c)
    target_mesh->cell_get_coordinates(c, &(target_cell_coords[c]));

  // Interpolate from source to target mesh

  Jali::Entity_ID_List const& targetcells =
      target_mesh->cells<Jali::Entity_type::ALL>();

  std::vector<double> outvals(ncells_target);
  std::vector<std::vector<Portage::Weights_t>>
      sources_and_weights(ncells_target);

  for (auto const& c : targetcells) {
    std::vector<int> xcells;
    std::vector<std::vector<double>> xwts;

    BOX_INTERSECT::intersection_moments(target_cell_coords[c],
                                        source_cell_coords,
                                        &xcells, &xwts);


    std::vector<Portage::Weights_t> wtsvec(xcells.size());
    for (int i = 0; i < xcells.size(); ++i) {
      wtsvec[i].entityID = xcells[i];
      wtsvec[i].weights = xwts[i];
    }
    sources_and_weights[c] = wtsvec;
  }

  // Create Interpolation object

  Portage::Interpolate_1stOrder<Portage::Jali_Mesh_Wrapper,
                                Portage::Jali_Mesh_Wrapper,
                                Portage::Jali_State_Wrapper,
                                Portage::CELL, 3>
      interpolater(sourceMeshWrapper, targetMeshWrapper, sourceStateWrapper);

  interpolater.set_interpolation_variable("cellvars");
  
  Portage::transform(targetcells.begin(), targetcells.end(),
                     sources_and_weights.begin(),
                     outvals.begin(), interpolater);
  
  // Make sure we retrieved the correct value for each cell on the target
  const double stdval = data[0];
  for (int c = 0; c < ncells_target; ++c)
    ASSERT_DOUBLE_EQ(stdval, outvals[c]);
}


/// First order interpolation of linear cell-centered field in 3D

TEST(Interpolate_1st_Order, Cell_Ctr_Lin_3D) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);

  std::shared_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                               4, 4, 4);
  std::shared_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                               2, 2, 2);

  const int ncells_source =
      source_mesh->num_entities(Jali::Entity_kind::CELL,
                                Jali::Entity_type::PARALLEL_OWNED);
  const int ncells_target =
      target_mesh->num_entities(Jali::Entity_kind::CELL,
                                Jali::Entity_type::PARALLEL_OWNED);

  // Create a state object

  Jali::State source_state(source_mesh);

  // Define a state vector representing a linear function

  std::vector<double> data(ncells_source);
  for (int c = 0; c < ncells_source; ++c) {
    JaliGeometry::Point cen;
    source_mesh->cell_centroid(c, &cen);
    data[c] = cen[0]+cen[1];
  }
  Jali::StateVector<double> myvec("cellvars", source_mesh,
                                  Jali::Entity_kind::CELL,
                                  Jali::Entity_type::PARALLEL_OWNED,
                                  &(data[0]));
  source_state.add(myvec);

  Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Portage::Jali_Mesh_Wrapper targetMeshWrapper(*target_mesh);
  Portage::Jali_State_Wrapper sourceStateWrapper(source_state);

  // Gather the cell coordinates for source and target meshes for intersection

  std::vector<std::vector<JaliGeometry::Point>>
      source_cell_coords(ncells_source);
  std::vector<std::vector<JaliGeometry::Point>>
      target_cell_coords(ncells_target);

  for (int c = 0; c < ncells_source; ++c)
    source_mesh->cell_get_coordinates(c, &(source_cell_coords[c]));
  for (int c = 0; c < ncells_target; ++c)
    target_mesh->cell_get_coordinates(c, &(target_cell_coords[c]));

  // Interpolate from source to target mesh

  Jali::Entity_ID_List const& targetcells =
      target_mesh->cells<Jali::Entity_type::ALL>();
  std::vector<double> outvals(ncells_target);
  std::vector<std::vector<Portage::Weights_t>>
      sources_and_weights(ncells_target);

  for (auto const& c : targetcells) {
    std::vector<int> xcells;
    std::vector<std::vector<double>> xwts;

    BOX_INTERSECT::intersection_moments(target_cell_coords[c],
                                        source_cell_coords,
                                        &xcells, &xwts);


    std::vector<Portage::Weights_t> wtsvec(xcells.size());
    for (int i = 0; i < xcells.size(); ++i) {
      wtsvec[i].entityID = xcells[i];
      wtsvec[i].weights = xwts[i];
    }
    sources_and_weights[c] = wtsvec;
  }

  // Create Interpolation objects

  Portage::Interpolate_1stOrder<Portage::Jali_Mesh_Wrapper,
                                Portage::Jali_Mesh_Wrapper,
                                Portage::Jali_State_Wrapper,
                                Portage::CELL, 3>
      interpolater(sourceMeshWrapper, targetMeshWrapper, sourceStateWrapper);

  interpolater.set_interpolation_variable("cellvars");

  Portage::transform(targetcells.begin(), targetcells.end(),
                     sources_and_weights.begin(),
                     outvals.begin(), interpolater);

  // Make sure we retrieved the correct value for each cell on the target
  // NOTE: EVEN THOUGH 1ST ORDER INTERPOLATION ALGORITHM DOES NOT IN
  // GENERAL PRESERVE A LINEAR FIELD THE SPECIAL STRUCTURE OF THE SOURCE
  // AND TARGET MESHES ENSURES THAT THE LINEAR FIELD IS INTERPOLATING CORRECTLY
  // IN THIS TEST

  std::vector<double> stdvals(ncells_target);
  for (int c = 0; c < ncells_target; ++c) {
    JaliGeometry::Point cen;
    target_mesh->cell_centroid(c, &cen);
    stdvals[c] = cen[0]+cen[1];
  }

  for (int c = 0; c < ncells_target; ++c)
    ASSERT_DOUBLE_EQ(stdvals[c], outvals[c]);
}


/// First order interpolation of constant node-centered field in 3D

TEST(Interpolate_1st_Order, Node_Ctr_Const_3D) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);
  mf.included_entities({Jali::Entity_kind::EDGE,
                        Jali::Entity_kind::FACE,
                        Jali::Entity_kind::WEDGE,
                        Jali::Entity_kind::CORNER});

  std::shared_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                               4, 4, 4);
  std::shared_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                               5, 5, 5);

  const int nnodes_source =
      source_mesh->num_entities(Jali::Entity_kind::NODE,
                                Jali::Entity_type::PARALLEL_OWNED);
  const int nnodes_target =
      target_mesh->num_entities(Jali::Entity_kind::NODE,
                                Jali::Entity_type::PARALLEL_OWNED);

  // Create a state object and add the first two vectors to it

  Jali::State source_state(source_mesh);

  // Define two state vectors, one with constant value, the other
  // with a linear function

  std::vector<double> data(nnodes_source, 1.5);
  Jali::StateVector<double> myvec("nodevars", source_mesh,
                                  Jali::Entity_kind::NODE,
                                  Jali::Entity_type::PARALLEL_OWNED,
                                  &(data[0]));
  source_state.add(myvec);

  Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Portage::Jali_Mesh_Wrapper targetMeshWrapper(*target_mesh);
  Portage::Jali_State_Wrapper sourceStateWrapper(source_state);


  Portage::MeshWrapperDual<Portage::Jali_Mesh_Wrapper>
      sourceDualWrapper(sourceMeshWrapper);
  Portage::MeshWrapperDual<Portage::Jali_Mesh_Wrapper>
      targetDualWrapper(targetMeshWrapper);

  // Interpolate from source to target mesh

  Jali::Entity_ID_List const& targetnodes =
      target_mesh->nodes<Jali::Entity_type::ALL>();

  std::vector<double> outvals(nnodes_target);
  std::vector<std::vector<Portage::Weights_t>>
      sources_and_weights(nnodes_target);

  // Gather the dual cell coordinates for source and target meshes for
  // intersection

  std::vector<std::vector<JaliGeometry::Point>>
      source_dualcell_coords(nnodes_source);
  std::vector<std::vector<JaliGeometry::Point>>
      target_dualcell_coords(nnodes_target);

  // Because the meshes are rectangular we can get away with examining
  // the coordinates of the corners instead of the wedges

  // Also, because we will use only the bounding box of each cell to
  // do the search and intersection we can get away with adding all
  // the coordinates of the corners to list including duplicates

  for (int n = 0; n < nnodes_source; ++n) {
    std::vector<JaliGeometry::Point> dualcoords;
    std::vector<int> corners;
    source_mesh->node_get_corners(n, Jali::Entity_type::ALL, &corners);

    for (auto cn : corners) {
      std::vector<JaliGeometry::Point> cncoords;
      source_mesh->corner_get_coordinates(cn, &cncoords);
      for (auto coord : cncoords)
        source_dualcell_coords[n].push_back(coord);
    }
  }

  for (int n = 0; n < nnodes_target; ++n) {
    std::vector<JaliGeometry::Point> dualcoords;
    std::vector<int> corners;
    target_mesh->node_get_corners(n, Jali::Entity_type::ALL, &corners);

    for (auto cn : corners) {
      std::vector<JaliGeometry::Point> cncoords;
      target_mesh->corner_get_coordinates(cn, &cncoords);
      for (auto coord : cncoords)
        target_dualcell_coords[n].push_back(coord);
    }
  }

  for (int n = 0; n < nnodes_target; ++n) {
    std::vector<int> xcells;
    std::vector<std::vector<double>> xwts;

    BOX_INTERSECT::intersection_moments(target_dualcell_coords[n],
                                        source_dualcell_coords,
                                        &xcells, &xwts);


    std::vector<Portage::Weights_t> wtsvec(xcells.size());
    for (int i = 0; i < xcells.size(); ++i) {
      wtsvec[i].entityID = xcells[i];
      wtsvec[i].weights = xwts[i];
    }
    sources_and_weights[n] = wtsvec;
  }

  // Create Interpolation object

  Portage::Interpolate_1stOrder<
    Portage::MeshWrapperDual<Portage::Jali_Mesh_Wrapper>,
    Portage::MeshWrapperDual<Portage::Jali_Mesh_Wrapper>,
    Portage::Jali_State_Wrapper,
    Portage::NODE, 3>
      interpolater(sourceDualWrapper, targetDualWrapper, sourceStateWrapper);

  interpolater.set_interpolation_variable("nodevars");

  Portage::transform(targetnodes.begin(), targetnodes.end(),
                     sources_and_weights.begin(),
                     outvals.begin(), interpolater);

  // Make sure we retrieved the correct value for each cell on the target

  const double stdval = data[0];
  for (int n = 0; n < nnodes_target; ++n)
    ASSERT_NEAR(stdval, outvals[n], TOL);
}

