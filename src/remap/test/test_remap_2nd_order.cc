/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include <iostream>

#include "portage/support/portage.h"
#include "portage/remap/remap_2nd_order.h"
#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"
#include "portage/wrappers/state/jali/jali_state_wrapper.h"
#include "portage/driver/driver.h"


#include "gtest/gtest.h"
#include "mpi.h"

// Jali includes
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "FrameworkTraits.hh"
#include "MeshFramework.hh"
#include "JaliState.h"
#include "JaliStateVector.h"
#include "Point.hh"

// Local include
#include "simple_intersect_for_tests.h"


/// Second order remap of constant cell-centered field with no limiter in 2D

TEST(Remap_2nd_Order, Cell_Ctr_Const_No_Limiter_2D) {

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);

  std::unique_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4);
  std::unique_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 1.0, 1.0, 5, 5);
  Jali::State source_state(source_mesh.get());
  Jali::State target_state(target_mesh.get());

  // Define two state vectors, one with constant value and the other
  // with a linear function that is x+y

  int ncells_source = source_mesh->num_entities(Jali::CELL, Jali::OWNED);
  std::vector<double> data(ncells_source, 1.25);
  Jali::StateVector<double> myvec("cellvars", Jali::CELL, source_mesh.get(),
                                  &(data[0]));
  Jali::StateVector<double> &addvec = source_state.add(myvec);

  int ncells_target = target_mesh->num_entities(Jali::CELL, Jali::OWNED);


  // Create Remap objects
  Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Portage::Jali_State_Wrapper sourceStateWrapper(source_state);

  Portage::Remap_2ndOrder<Portage::Jali_Mesh_Wrapper,
                          Portage::Jali_State_Wrapper,
                          Portage::CELL>
      remapper(sourceMeshWrapper, sourceStateWrapper, "cellvars",
               Portage::NOLIMITER);


  // Gather the cell coordinates for source and target meshes for intersection

  std::vector<std::vector<JaliGeometry::Point>>
      source_cell_coords(ncells_source);
  std::vector<std::vector<JaliGeometry::Point>>
      target_cell_coords(ncells_target);

  for (int c = 0; c < ncells_source; c++)
    source_mesh->cell_get_coordinates(c, &(source_cell_coords[c]));
  for (int c = 0; c < ncells_target; c++)
    target_mesh->cell_get_coordinates(c, &(target_cell_coords[c]));

  // Remap from source to target mesh

  std::vector<double> outvals(ncells_target);  // field values on target mesh
  for (int c = 0; c < ncells_target; ++c) {
    std::vector<int> xcells;
    std::vector<std::vector<double>> xweights;

    BOX_INTERSECT::intersection_moments(target_cell_coords[c],
                                        source_cell_coords,
                                        &xcells, &xweights);

    std::pair< std::vector<int> const &,
               std::vector<std::vector<double>> const &>
        cells_and_weights(xcells, xweights);

    outvals[c] = remapper(cells_and_weights);
  }

  // Make sure we retrieved the correct value for each cell on the target

  double stdval = data[0];
  for (int c = 0; c < ncells_target; ++c)
    ASSERT_NEAR(stdval, outvals[c], 1.0e-10);
}


/// Second order remap of linear cell-centered field with no limiting in 2D

TEST(Remap_2nd_Order, Cell_Ctr_Lin_No_Limiter_2D) {

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);
  
  std::unique_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4);
  std::unique_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 1.0, 1.0, 5, 5);
  Jali::State source_state(source_mesh.get());
  Jali::State target_state(target_mesh.get());
  
  // Define a state vectors, with a linear function that is x+y

  int ncells_source = source_mesh->num_entities(Jali::CELL, Jali::OWNED);
  std::vector<double> data(ncells_source);
  for (int c = 0; c < ncells_source; c++) {
    JaliGeometry::Point ccen = source_mesh->cell_centroid(c);
    data[c] = ccen[0]+ccen[1];
  }


  Jali::StateVector<double> myvec("cellvars", Jali::CELL, source_mesh.get(),
                                  &(data[0]));
  Jali::StateVector<double> &addvec = source_state.add(myvec);

  // Create Remap objects

  Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Portage::Jali_State_Wrapper sourceStateWrapper(source_state);

  Portage::Remap_2ndOrder<Portage::Jali_Mesh_Wrapper,
                          Portage::Jali_State_Wrapper,
                          Portage::CELL>
      remapper(sourceMeshWrapper, sourceStateWrapper, "cellvars",
               Portage::NOLIMITER);

  // Gather the cell coordinates for source and target meshes for intersection

  int ncells_target = target_mesh->num_entities(Jali::CELL, Jali::OWNED);
  std::vector<std::vector<JaliGeometry::Point>>
      source_cell_coords(ncells_source);
  std::vector<std::vector<JaliGeometry::Point>>
      target_cell_coords(ncells_target);

  for (int c = 0; c < ncells_source; c++)
    source_mesh->cell_get_coordinates(c, &(source_cell_coords[c]));
  for (int c = 0; c < ncells_target; c++)
    target_mesh->cell_get_coordinates(c, &(target_cell_coords[c]));


  // Remap from source to target mesh
  std::vector<double> outvals(ncells_target);

  for (int c = 0; c < ncells_target; ++c) {
    std::vector<int> xcells;
    std::vector<std::vector<double>> xweights;

    BOX_INTERSECT::intersection_moments(target_cell_coords[c],
                                        source_cell_coords,
                                        &xcells, &xweights);

    std::pair< std::vector<int> const &,
               std::vector<std::vector<double>> const &>
        cells_and_weights(xcells, xweights);

    outvals[c] = remapper(cells_and_weights);
  }

  // Make sure we retrieved the correct value for each cell on the target

  std::vector<double> stdvals(ncells_target);
  for (int c = 0; c < ncells_target; ++c) {
    JaliGeometry::Point ccen = target_mesh->cell_centroid(c);
    stdvals[c] = ccen[0]+ccen[1];
  }
  for (int c = 0; c < ncells_target; ++c)
    ASSERT_NEAR(stdvals[c], outvals[c], 1.0e-10);
}


/*!  @brief Second order remap of linear cell-centered field with
  Barth-Jespersen limiting in 2D */

TEST(Remap_2nd_Order, Cell_Ctr_Lin_BJ_Limiter_2D) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);

  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);

  std::unique_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4);
  std::unique_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 1.0, 1.0, 5, 5);
  Jali::State source_state(source_mesh.get());
  Jali::State target_state(target_mesh.get());

  // Define a state vectors, with a linear function that is x+y

  int ncells_source = source_mesh->num_entities(Jali::CELL, Jali::OWNED);
  std::vector<double> data(ncells_source);
  double minval = 1e+10, maxval = -1e+10;
  for (int c = 0; c < ncells_source; c++) {
    JaliGeometry::Point ccen = source_mesh->cell_centroid(c);
    if (ccen[0] < 0.5)
      data[c] = ccen[0]+ccen[1];
    else
      data[c] = 100*ccen[0];
    if (data[c] < minval) minval = data[c];
    if (data[c] > maxval) maxval = data[c];
  }

  Jali::StateVector<double> myvec("cellvars", Jali::CELL, source_mesh.get(),
                                  &(data[0]));
  Jali::StateVector<double> &addvec = source_state.add(myvec);

  // Create Remap objects - one with no limiter and one with limiter
  Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Portage::Jali_State_Wrapper sourceStateWrapper(source_state);


  Portage::Remap_2ndOrder<Portage::Jali_Mesh_Wrapper,
                          Portage::Jali_State_Wrapper,
                          Portage::CELL>
      remapper1(sourceMeshWrapper, sourceStateWrapper, "cellvars",
               Portage::NOLIMITER);
  Portage::Remap_2ndOrder<Portage::Jali_Mesh_Wrapper,
                          Portage::Jali_State_Wrapper,
                          Portage::CELL>
      remapper2(sourceMeshWrapper, sourceStateWrapper, "cellvars",
               Portage::BARTH_JESPERSEN);


  // Gather the cell coordinates for the source and target meshes for
  // intersection

  int ncells_target = target_mesh->num_entities(Jali::CELL, Jali::OWNED);
  std::vector<std::vector<JaliGeometry::Point>>
      source_cell_coords(ncells_source);
  std::vector<std::vector<JaliGeometry::Point>>
      target_cell_coords(ncells_target);

  for (int c = 0; c < ncells_source; c++)
    source_mesh->cell_get_coordinates(c, &(source_cell_coords[c]));
  for (int c = 0; c < ncells_target; c++)
    target_mesh->cell_get_coordinates(c, &(target_cell_coords[c]));


  std::vector<double> outvals1(ncells_target);
  std::vector<double> outvals2(ncells_target);

  // Remap from source to target mesh
  for (int c = 0; c < ncells_target; ++c) {
    std::vector<int> xcells;
    std::vector<std::vector<double>> xweights;

    BOX_INTERSECT::intersection_moments(target_cell_coords[c],
                                        source_cell_coords,
                                        &xcells, &xweights);

    std::pair< std::vector<int> const &,
               std::vector<std::vector<double>> const &>
        cells_and_weights(xcells, xweights);

    outvals1[c] = remapper1(cells_and_weights);
    outvals2[c] = remapper2(cells_and_weights);
  }

  // Check if we violated the bounds on at least one cell in the
  // unlimited remap and if we respected the bounds in all cells in
  // the limited case

  bool outofbounds_unlimited = false;
  for (int c = 0; c < ncells_target; ++c) {
    if (outvals1[c] < minval  || outvals1[c] > maxval) {
      outofbounds_unlimited = true;
      break;
    }
  }

  bool inbounds_limited = true;
  for (int c = 0; c < ncells_target; ++c) {
    if (outvals2[c] < minval  || outvals2[c] > maxval) {
      inbounds_limited = false;
      break;
    }
  }

  EXPECT_TRUE(outofbounds_unlimited && inbounds_limited);
}


/// Second order remap of constant node-centered field with no limiting in 2D

TEST(Remap_2nd_Order, Node_Ctr_Const_No_Limiter) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);

  std::unique_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4, NULL,
                               true, true, true, true);
  std::unique_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 1.0, 1.0, 5, 5, NULL,
                               true, true, true, true);

  int nnodes_source = source_mesh->num_entities(Jali::NODE, Jali::OWNED);
  int nnodes_target = target_mesh->num_entities(Jali::NODE, Jali::OWNED);

  // Create a state object and add the first two vectors to it

  Jali::State source_state(source_mesh.get());


  // Define two state vectors, one with constant value, the other
  // with a linear function

  std::vector<double> data(nnodes_source, 1.5);
  Jali::StateVector<double> myvec("nodevars", Jali::NODE, source_mesh.get(),
                                  &(data[0]));
  Jali::StateVector<double> &addvec = source_state.add(myvec);

  // Create Remap objects
  Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Portage::Jali_State_Wrapper sourceStateWrapper(source_state);

  Portage::Remap_2ndOrder<Portage::Jali_Mesh_Wrapper,
                          Portage::Jali_State_Wrapper,
                          Portage::NODE>
      remapper(sourceMeshWrapper, sourceStateWrapper, "nodevars",
               Portage::NOLIMITER);

  // Remap from source to target mesh

  std::vector<double> outvals(nnodes_target);


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

  for (int n = 0; n < nnodes_source; n++) {
    std::vector<JaliGeometry::Point> dualcoords;
    std::vector<int> corners;
    source_mesh->node_get_corners(n, Jali::ALL, &corners);

    for (auto cn : corners) {
      std::vector<JaliGeometry::Point> cncoords;
      source_mesh->corner_get_coordinates(cn, &cncoords);
      for (auto coord : cncoords)
        source_dualcell_coords[n].push_back(coord);
    }
  }

  for (int n = 0; n < nnodes_target; n++) {
    std::vector<JaliGeometry::Point> dualcoords;
    std::vector<int> corners;
    target_mesh->node_get_corners(n, Jali::ALL, &corners);

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
    
    std::pair< std::vector<int> const &,
               std::vector<std::vector<double>> const &>
        nodes_and_weights(xcells, xwts);

    outvals[n] = remapper(nodes_and_weights);
  }

  // Make sure we retrieved the correct value for each cell on the target

  double stdval = data[0];
  for (int n = 0; n < nnodes_target; ++n)
    ASSERT_DOUBLE_EQ(stdval, outvals[n]);
}


/// Second order remap of linear node-centered field with no limiting in 2D

TEST(Remap_2nd_Order, Node_Ctr_Lin_No_Limiter) {

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);

  std::unique_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4, NULL,
                                               true, true, true, true);
  std::unique_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 1.0, 1.0, 5, 5, NULL,
                                               true, true, true, true);

  int nnodes_source = source_mesh->num_entities(Jali::NODE, Jali::OWNED);
  int nnodes_target = target_mesh->num_entities(Jali::NODE, Jali::OWNED);

  // Create a state object and add the first two vectors to it

  Jali::State source_state(source_mesh.get());

  // Define a state vectors representing a linear function

  std::vector<double> data(nnodes_source);
  for (int n = 0; n < nnodes_source; n++) {
    JaliGeometry::Point coord;
    source_mesh->node_get_coordinates(n, &coord);
    data[n] = coord[0]+coord[1];
  }
  Jali::StateVector<double> myvec("nodevars", Jali::NODE, source_mesh.get(),
                                  &(data[0]));
  Jali::StateVector<double> &addvec = source_state.add(myvec);


  // Create a mesh wrapper

  Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Portage::Jali_Mesh_Wrapper targetMeshWrapper(*target_mesh);

  // Create Remap objects

  Portage::Remap_2ndOrder<Portage::Jali_Mesh_Wrapper,
                          Portage::Jali_State_Wrapper,
                          Portage::NODE>
      remapper(sourceMeshWrapper, source_state, "nodevars",
               Portage::NOLIMITER);


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

  for (int n = 0; n < nnodes_source; n++) {
    std::vector<JaliGeometry::Point> dualcoords;
    std::vector<int> corners;
    source_mesh->node_get_corners(n, Jali::ALL, &corners);

    for (auto cn : corners) {
      std::vector<JaliGeometry::Point> cncoords;
      source_mesh->corner_get_coordinates(cn, &cncoords);
      for (auto coord : cncoords)
        source_dualcell_coords[n].push_back(coord);
    }
  }

  for (int n = 0; n < nnodes_target; n++) {
    std::vector<JaliGeometry::Point> dualcoords;
    std::vector<int> corners;
    target_mesh->node_get_corners(n, Jali::ALL, &corners);

    for (auto cn : corners) {
      std::vector<JaliGeometry::Point> cncoords;
      target_mesh->corner_get_coordinates(cn, &cncoords);
      for (auto coord : cncoords)
        target_dualcell_coords[n].push_back(coord);
    }
  }


  // Remap from source to target mesh

  std::vector<double> outvals(nnodes_target);

  for (int n = 0; n < nnodes_target; ++n) {
    std::vector<int> xcells;
    std::vector<std::vector<double>> xwts;

    BOX_INTERSECT::intersection_moments(target_dualcell_coords[n],
                                        source_dualcell_coords,
                                        &xcells, &xwts);

    std::pair< std::vector<int> const &,
               std::vector<std::vector<double>> const &>
        nodes_and_weights(xcells, xwts);

    outvals[n] = remapper(nodes_and_weights);
  }

  // Make sure we retrieved the correct value for each node on the
  // target Second order remapping should preserve a linear field
  // exactly but node-centered remapping has some quirks - the field
  // does not get preserved exactly at the boundary because the source
  // values for boundary dual cells are not specified inside the dual
  // cells but at one of their vertices or edges. So, check only
  // interior nodes

  std::vector<double> stdvals(nnodes_target);
  for (int n = 0; n < nnodes_target; ++n) {
    JaliGeometry::Point coord;
    target_mesh->node_get_coordinates(n, &coord);
    if (fabs(coord[0]) < 1e-16 || fabs(1-coord[0]) < 1e-16 ||
        fabs(coord[1]) < 1e-16 || fabs(1-coord[1]) < 1.e-16)
      continue;
    stdvals[n] = coord[0]+coord[1];
    EXPECT_DOUBLE_EQ(stdvals[n], outvals[n]);
  }
}


/// Second order remap of constant cell-centered field with no limiting in 3D

TEST(Remap_2nd_Order, Cell_Ctr_Const_No_Limiter_3D) {

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);

  std::unique_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                               4, 4, 4);
  std::unique_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                               5, 5, 5);
  Jali::State source_state(source_mesh.get());
  Jali::State target_state(target_mesh.get());

  // Define two state vectors, one with constant value and the other
  // with a linear function that is x+y

  int ncells_source = source_mesh->num_entities(Jali::CELL, Jali::OWNED);
  std::vector<double> data(ncells_source, 1.25);
  Jali::StateVector<double> myvec("cellvars", Jali::CELL, source_mesh.get(),
                                  &(data[0]));
  Jali::StateVector<double> &addvec = source_state.add(myvec);

  int ncells_target = target_mesh->num_entities(Jali::CELL, Jali::OWNED);


  // Create Remap objects
  Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Portage::Jali_State_Wrapper sourceStateWrapper(source_state);

  Portage::Remap_2ndOrder<Portage::Jali_Mesh_Wrapper,
                          Portage::Jali_State_Wrapper,
                          Portage::CELL>
      remapper(sourceMeshWrapper, sourceStateWrapper, "cellvars",
               Portage::NOLIMITER);


  // Gather the cell coordinates for source and target meshes for intersection

  std::vector<std::vector<JaliGeometry::Point>>
      source_cell_coords(ncells_source);
  std::vector<std::vector<JaliGeometry::Point>>
      target_cell_coords(ncells_target);

  for (int c = 0; c < ncells_source; c++)
    source_mesh->cell_get_coordinates(c, &(source_cell_coords[c]));
  for (int c = 0; c < ncells_target; c++)
    target_mesh->cell_get_coordinates(c, &(target_cell_coords[c]));

  // Remap from source to target mesh

  std::vector<double> outvals(ncells_target);  // field values on target mesh
  for (int c = 0; c < ncells_target; ++c) {
    std::vector<int> xcells;
    std::vector<std::vector<double>> xweights;

    BOX_INTERSECT::intersection_moments(target_cell_coords[c],
                                        source_cell_coords,
                                        &xcells, &xweights);

    std::pair< std::vector<int> const &,
               std::vector<std::vector<double>> const &>
        cells_and_weights(xcells, xweights);

    outvals[c] = remapper(cells_and_weights);
  }

  // Make sure we retrieved the correct value for each cell on the target

  double stdval = data[0];
  for (int c = 0; c < ncells_target; ++c)
    ASSERT_NEAR(stdval, outvals[c], 1.0e-10);
}


/// Second order remap of linear cell-centered field with no limiting in 3D

TEST(Remap_2nd_Order, Cell_Ctr_Lin_No_Limiter_3D) {

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);

  std::unique_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                               4, 4, 4);
  std::unique_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                               5, 5, 5);
  Jali::State source_state(source_mesh.get());
  Jali::State target_state(target_mesh.get());

  // Define a state vectors, with a linear function that is x+y

  int ncells_source = source_mesh->num_entities(Jali::CELL, Jali::OWNED);
  std::vector<double> data(ncells_source);
  for (int c = 0; c < ncells_source; c++) {
    JaliGeometry::Point ccen = source_mesh->cell_centroid(c);
    data[c] = ccen[0]+ccen[1]+ccen[2];
  }


  Jali::StateVector<double> myvec("cellvars", Jali::CELL, source_mesh.get(),
                                  &(data[0]));
  Jali::StateVector<double> &addvec = source_state.add(myvec);

  // Create Remap objects

  Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Portage::Jali_State_Wrapper sourceStateWrapper(source_state);

  Portage::Remap_2ndOrder<Portage::Jali_Mesh_Wrapper,
                          Portage::Jali_State_Wrapper,
                          Portage::CELL>
      remapper(sourceMeshWrapper, sourceStateWrapper, "cellvars",
               Portage::NOLIMITER);

  // Gather the cell coordinates for source and target meshes for intersection

  int ncells_target = target_mesh->num_entities(Jali::CELL, Jali::OWNED);
  std::vector<std::vector<JaliGeometry::Point>>
      source_cell_coords(ncells_source);
  std::vector<std::vector<JaliGeometry::Point>>
      target_cell_coords(ncells_target);

  for (int c = 0; c < ncells_source; c++)
    source_mesh->cell_get_coordinates(c, &(source_cell_coords[c]));
  for (int c = 0; c < ncells_target; c++)
    target_mesh->cell_get_coordinates(c, &(target_cell_coords[c]));


  // Remap from source to target mesh
  std::vector<double> outvals(ncells_target);

  for (int c = 0; c < ncells_target; ++c) {
    std::vector<int> xcells;
    std::vector<std::vector<double>> xweights;

    BOX_INTERSECT::intersection_moments(target_cell_coords[c],
                                        source_cell_coords,
                                        &xcells, &xweights);

    std::pair<std::vector<int> const &,
              std::vector<std::vector<double>> const &>
        cells_and_weights(xcells, xweights);

    outvals[c] = remapper(cells_and_weights);
  }

  // Make sure we retrieved the correct value for each cell on the target

  std::vector<double> stdvals(ncells_target);
  for (int c = 0; c < ncells_target; ++c) {
    JaliGeometry::Point ccen = target_mesh->cell_centroid(c);
    stdvals[c] = ccen[0]+ccen[1]+ccen[2];
  }
  for (int c = 0; c < ncells_target; ++c)
    ASSERT_NEAR(stdvals[c], outvals[c], 1.0e-10);
}


/// Second order remap of discontinuous cell-centered field with
/// Barth-Jespersen limiting in 3D

TEST(Remap_2nd_Order, Cell_Ctr_BJ_Limiter_3D) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);

  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);

  std::unique_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                               4, 4, 4);
  std::unique_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                               5, 5, 5);
  Jali::State source_state(source_mesh.get());
  Jali::State target_state(target_mesh.get());

  // Define a state vectors, with a linear function that is x+y

  int ncells_source = source_mesh->num_entities(Jali::CELL, Jali::OWNED);
  std::vector<double> data(ncells_source);
  double minval = 1e+10, maxval = -1e+10;
  for (int c = 0; c < ncells_source; c++) {
    JaliGeometry::Point ccen = source_mesh->cell_centroid(c);
    if (ccen[0] < 0.5)
      data[c] = ccen[0]+ccen[1]+ccen[2];
    else
      data[c] = 100*ccen[0];
    if (data[c] < minval) minval = data[c];
    if (data[c] > maxval) maxval = data[c];
  }

  Jali::StateVector<double> myvec("cellvars", Jali::CELL, source_mesh.get(),
                                  &(data[0]));
  Jali::StateVector<double> &addvec = source_state.add(myvec);

  // Create Remap objects - one with no limiter and one with limiter
  Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Portage::Jali_State_Wrapper sourceStateWrapper(source_state);


  Portage::Remap_2ndOrder<Portage::Jali_Mesh_Wrapper,
                          Portage::Jali_State_Wrapper,
                          Portage::CELL>
      remapper1(sourceMeshWrapper, sourceStateWrapper, "cellvars",
               Portage::NOLIMITER);
  Portage::Remap_2ndOrder<Portage::Jali_Mesh_Wrapper,
                          Portage::Jali_State_Wrapper,
                          Portage::CELL>
      remapper2(sourceMeshWrapper, sourceStateWrapper, "cellvars",
               Portage::BARTH_JESPERSEN);


  // Gather the cell coordinates for the source and target meshes for
  // intersection

  int ncells_target = target_mesh->num_entities(Jali::CELL, Jali::OWNED);
  std::vector<std::vector<JaliGeometry::Point>>
      source_cell_coords(ncells_source);
  std::vector<std::vector<JaliGeometry::Point>>
      target_cell_coords(ncells_target);

  for (int c = 0; c < ncells_source; c++)
    source_mesh->cell_get_coordinates(c, &(source_cell_coords[c]));
  for (int c = 0; c < ncells_target; c++)
    target_mesh->cell_get_coordinates(c, &(target_cell_coords[c]));


  std::vector<double> outvals1(ncells_target);
  std::vector<double> outvals2(ncells_target);

  // Remap from source to target mesh
  for (int c = 0; c < ncells_target; ++c) {
    std::vector<int> xcells;
    std::vector<std::vector<double>> xweights;

    BOX_INTERSECT::intersection_moments(target_cell_coords[c],
                                        source_cell_coords,
                                        &xcells, &xweights);

    std::pair< std::vector<int> const &,
               std::vector<std::vector<double>> const &>
        cells_and_weights(xcells, xweights);

    outvals1[c] = remapper1(cells_and_weights);
    outvals2[c] = remapper2(cells_and_weights);
  }

  // Check if we violated the bounds on at least one cell in the
  // unlimited remap and if we respected the bounds in all cells in
  // the limited case

  bool outofbounds_unlimited = false;
  for (int c = 0; c < ncells_target; ++c) {
    if (outvals1[c] < minval  || outvals1[c] > maxval) {
      outofbounds_unlimited = true;
      break;
    }
  }

  bool inbounds_limited = true;
  for (int c = 0; c < ncells_target; ++c) {
    if (outvals2[c] < minval && minval-outvals2[c] > 1.0e-10) {
      inbounds_limited = false;
      break;
    }
    else if (outvals2[c] > maxval && outvals2[c]-maxval > 1.0e-10) {
      inbounds_limited = false;
      break;
    }
  }

  EXPECT_TRUE(outofbounds_unlimited);
  EXPECT_TRUE(inbounds_limited);
}


/// Second order remap of constant node-centered field with no limiting in 3D

TEST(Remap_2nd_Order, Node_Ctr_Const_No_Limiter_3D) {

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);

  std::unique_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                               4, 4, 4, NULL,
                                               true, true, true, true);
  std::unique_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                               5, 5, 5, NULL,
                                               true, true, true, true);

  int nnodes_source = source_mesh->num_entities(Jali::NODE, Jali::OWNED);
  int nnodes_target = target_mesh->num_entities(Jali::NODE, Jali::OWNED);

  // Create a state object and add the first two vectors to it

  Jali::State source_state(source_mesh.get());


  // Define a state vectors representing a linear function

  double nodeval = 3.0;
  std::vector<double> data(nnodes_source, nodeval);
  Jali::StateVector<double> myvec("nodevars", Jali::NODE, source_mesh.get(),
                                  &(data[0]));
  Jali::StateVector<double> &addvec = source_state.add(myvec);


  // Create a mesh wrapper

  Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Portage::Jali_Mesh_Wrapper targetMeshWrapper(*target_mesh);

  // Create Remap objects

  Portage::Remap_2ndOrder<Portage::Jali_Mesh_Wrapper,
                          Portage::Jali_State_Wrapper,
                          Portage::NODE>
      remapper(sourceMeshWrapper, source_state, "nodevars", Portage::NOLIMITER);


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

  for (int n = 0; n < nnodes_source; n++) {
    std::vector<JaliGeometry::Point> dualcoords;
    std::vector<int> corners;
    source_mesh->node_get_corners(n, Jali::ALL, &corners);

    for (auto cn : corners) {
      std::vector<JaliGeometry::Point> cncoords;
      source_mesh->corner_get_coordinates(cn, &cncoords);
      for (auto coord : cncoords)
        source_dualcell_coords[n].push_back(coord);
    }
  }

  for (int n = 0; n < nnodes_target; n++) {
    std::vector<JaliGeometry::Point> dualcoords;
    std::vector<int> corners;
    target_mesh->node_get_corners(n, Jali::ALL, &corners);

    for (auto cn : corners) {
      std::vector<JaliGeometry::Point> cncoords;
      target_mesh->corner_get_coordinates(cn, &cncoords);
      for (auto coord : cncoords)
        target_dualcell_coords[n].push_back(coord);
    }
  }


  // Remap from source to target mesh

  std::vector<double> outvals(nnodes_target);

  for (int n = 0; n < nnodes_target; ++n) {
    std::vector<int> xcells;
    std::vector<std::vector<double>> xwts;

    BOX_INTERSECT::intersection_moments(target_dualcell_coords[n],
                                        source_dualcell_coords,
                                        &xcells, &xwts);

    std::pair< std::vector<int> const &,
               std::vector<std::vector<double>> const &>
        nodes_and_weights(xcells, xwts);

    outvals[n] = remapper(nodes_and_weights);
  }

  // Make sure we retrieved the correct value for each node on the
  // target Second order remapping should preserve a linear field
  // exactly but node-centered remapping has some quirks - the field
  // does not get preserved exactly at the boundary because the source
  // values for boundary dual cells are not specified inside the dual
  // cells but at one of their vertices or edges. So, check only
  // interior nodes

  for (int n = 0; n < nnodes_target; ++n) {
    JaliGeometry::Point coord;
    target_mesh->node_get_coordinates(n, &coord);
    //    if (fabs(coord[0]) < 1e-16 || fabs(1-coord[0]) < 1e-16 ||
    //        fabs(coord[1]) < 1e-16 || fabs(1-coord[1]) < 1.e-16 ||
    //        fabs(coord[2]) < 1e-16 || fabs(1-coord[2]) < 1.e-16)
    //      continue;
    EXPECT_DOUBLE_EQ(nodeval, outvals[n]);
  }
}



TEST(Remap_2nd_Order, Node_Ctr_Lin_No_Limiter_3D) {

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);

  std::unique_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                               4, 4, 4, NULL,
                                               true, true, true, true);
  std::unique_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                               5, 5, 5, NULL,
                                               true, true, true, true);

  int nnodes_source = source_mesh->num_entities(Jali::NODE, Jali::OWNED);
  int nnodes_target = target_mesh->num_entities(Jali::NODE, Jali::OWNED);

  // Create a state object and add the first two vectors to it

  Jali::State source_state(source_mesh.get());


  // Define a state vectors representing a linear function

  std::vector<double> data(nnodes_source);
  for (int n = 0; n < nnodes_source; n++) {
    JaliGeometry::Point coord;
    source_mesh->node_get_coordinates(n, &coord);
    data[n] = coord[0]+coord[1]+coord[2];
  }
  Jali::StateVector<double> myvec("nodevars", Jali::NODE, source_mesh.get(),
                                  &(data[0]));
  Jali::StateVector<double> &addvec = source_state.add(myvec);


  // Create a mesh wrapper

  Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Portage::Jali_Mesh_Wrapper targetMeshWrapper(*target_mesh);

  // Create Remap objects

  Portage::Remap_2ndOrder<Portage::Jali_Mesh_Wrapper,
                          Portage::Jali_State_Wrapper,
                          Portage::NODE>
      remapper(sourceMeshWrapper, source_state, "nodevars",
               Portage::NOLIMITER);


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

  for (int n = 0; n < nnodes_source; n++) {
    std::vector<JaliGeometry::Point> dualcoords;
    std::vector<int> corners;
    source_mesh->node_get_corners(n, Jali::ALL, &corners);

    for (auto cn : corners) {
      std::vector<JaliGeometry::Point> cncoords;
      source_mesh->corner_get_coordinates(cn, &cncoords);
      for (auto coord : cncoords)
        source_dualcell_coords[n].push_back(coord);
    }
  }

  for (int n = 0; n < nnodes_target; n++) {
    std::vector<JaliGeometry::Point> dualcoords;
    std::vector<int> corners;
    target_mesh->node_get_corners(n, Jali::ALL, &corners);

    for (auto cn : corners) {
      std::vector<JaliGeometry::Point> cncoords;
      target_mesh->corner_get_coordinates(cn, &cncoords);
      for (auto coord : cncoords)
        target_dualcell_coords[n].push_back(coord);
    }
  }


  // Remap from source to target mesh

  std::vector<double> outvals(nnodes_target);

  for (int n = 0; n < nnodes_target; ++n) {
    std::vector<int> xcells;
    std::vector<std::vector<double>> xwts;

    BOX_INTERSECT::intersection_moments(target_dualcell_coords[n],
                                        source_dualcell_coords,
                                        &xcells, &xwts);

    std::pair< std::vector<int> const &,
               std::vector<std::vector<double>> const &>
        nodes_and_weights(xcells, xwts);

    outvals[n] = remapper(nodes_and_weights);
  }

  // Make sure we retrieved the correct value for each node on the
  // target Second order remapping should preserve a linear field
  // exactly but node-centered remapping has some quirks - the field
  // does not get preserved exactly at the boundary because the source
  // values for boundary dual cells are not specified inside the dual
  // cells but at one of their vertices or edges. So, check only
  // interior nodes

  std::vector<double> stdvals(nnodes_target);
  for (int n = 0; n < nnodes_target; ++n) {
    JaliGeometry::Point coord;
    target_mesh->node_get_coordinates(n, &coord);
    if (fabs(coord[0]) < 1e-16 || fabs(1-coord[0]) < 1e-16 ||
        fabs(coord[1]) < 1e-16 || fabs(1-coord[1]) < 1.e-16 ||
        fabs(coord[2]) < 1e-16 || fabs(1-coord[2]) < 1.e-16)
      continue;
    stdvals[n] = coord[0]+coord[1]+coord[2];
    EXPECT_DOUBLE_EQ(stdvals[n], outvals[n]);
  }
}


/// Second order remap of discontinuous node-centered field with
/// Barth-Jespersen limiting in 3D

TEST(Remap_2nd_Order, Node_Ctr_BJ_Limiter_3D) {

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);

  std::unique_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4, NULL,
                               true, true, true, true);
  std::unique_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 1.0, 1.0, 5, 5, NULL,
                               true, true, true, true);

  int nnodes_source = source_mesh->num_entities(Jali::NODE, Jali::OWNED);
  int nnodes_target = target_mesh->num_entities(Jali::NODE, Jali::OWNED);

  // Create a state object and add the first two vectors to it

  Jali::State source_state(source_mesh.get());


  // Define a state vector representing two piecewise linear functions,
  // where the right half has 100 times the value it would
  // have in the left linear function

  std::vector<double> data(nnodes_source);
  double minval = 1e+10, maxval = -1e+10;
  for (int n = 0; n < nnodes_source; n++) {
    JaliGeometry::Point coord;
    source_mesh->node_get_coordinates(n, &coord);
    if (coord[0] >= 0.5)
      data[n] = 100*(coord[0]+coord[1]+coord[2]);
    else
      data[n] = coord[0]+coord[1]+coord[2];
    if (data[n] < minval) minval = data[n];
    if (data[n] > maxval) maxval = data[n];
  }
  Jali::StateVector<double> myvec("nodevars", Jali::NODE, source_mesh.get(),
                                  &(data[0]));
  Jali::StateVector<double> &addvec = source_state.add(myvec);


  // Create a mesh wrapper

  Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Portage::Jali_Mesh_Wrapper targetMeshWrapper(*target_mesh);

  // Create Remap objects - one with no limiter and one with limiter

  Portage::Remap_2ndOrder<Portage::Jali_Mesh_Wrapper,
                          Portage::Jali_State_Wrapper,
                          Portage::NODE>
      remapper1(sourceMeshWrapper, source_state, "nodevars",
                Portage::NOLIMITER);
  Portage::Remap_2ndOrder<Portage::Jali_Mesh_Wrapper,
                          Portage::Jali_State_Wrapper,
                          Portage::NODE>
      remapper2(sourceMeshWrapper, source_state, "nodevars",
                Portage::BARTH_JESPERSEN);


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

  for (int n = 0; n < nnodes_source; n++) {
    std::vector<JaliGeometry::Point> dualcoords;
    std::vector<int> corners;
    source_mesh->node_get_corners(n, Jali::ALL, &corners);

    for (auto cn : corners) {
      std::vector<JaliGeometry::Point> cncoords;
      source_mesh->corner_get_coordinates(cn, &cncoords);
      for (auto coord : cncoords)
        source_dualcell_coords[n].push_back(coord);
    }
  }

  for (int n = 0; n < nnodes_target; n++) {
    std::vector<JaliGeometry::Point> dualcoords;
    std::vector<int> corners;
    target_mesh->node_get_corners(n, Jali::ALL, &corners);

    for (auto cn : corners) {
      std::vector<JaliGeometry::Point> cncoords;
      target_mesh->corner_get_coordinates(cn, &cncoords);
      for (auto coord : cncoords)
        target_dualcell_coords[n].push_back(coord);
    }
  }


  // Remap from source to target mesh

  std::vector<double> outvals1(nnodes_target);
  std::vector<double> outvals2(nnodes_target);

  for (int n = 0; n < nnodes_target; ++n) {
    std::vector<int> xcells;
    std::vector<std::vector<double>> xwts;

    BOX_INTERSECT::intersection_moments(target_dualcell_coords[n],
                                        source_dualcell_coords,
                                        &xcells, &xwts);

    std::pair< std::vector<int> const &,
               std::vector<std::vector<double>> const &>
        nodes_and_weights(xcells, xwts);

    outvals1[n] = remapper1(nodes_and_weights);
    outvals2[n] = remapper2(nodes_and_weights);
  }

  // Check if we violated the bounds on at least one node in the
  // unlimited remap and if we respected the bounds on all nodes in
  // the limited case

  bool outofbounds_unlimited = false;
  for (int n = 0; n < nnodes_target; ++n) {
    if (outvals1[n] < minval  || outvals1[n] > maxval) {
      outofbounds_unlimited = true;
      break;
    }
  }

  bool inbounds_limited = true;
  for (int n = 0; n < nnodes_target; ++n) {
    if (outvals2[n] < minval  || outvals2[n] > maxval) {
      inbounds_limited = false;
      break;
    }
  }

  EXPECT_TRUE(outofbounds_unlimited);
  EXPECT_TRUE(inbounds_limited);
}
