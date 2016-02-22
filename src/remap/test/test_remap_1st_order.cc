/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/
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
#include "portage/remap/remap_1st_order.h"
#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"
#include "portage/wrappers/state/jali/jali_state_wrapper.h"
#include "portage/remap/test/simple_intersect_for_tests.h"


/// First order remap of constant cell-centered field in 2D

TEST(Remap_1st_Order, Cell_Ctr_Const_2D) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);

  std::unique_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4);
  std::unique_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 1.0, 1.0, 5, 5);

  const int ncells_source = source_mesh->num_entities(Jali::CELL, Jali::OWNED);
  const int ncells_target = target_mesh->num_entities(Jali::CELL, Jali::OWNED);

  // Create a state object

  Jali::State source_state(source_mesh.get());

  // Define a state vectors with constant value

  std::vector<double> data(ncells_source, 1.25);
  Jali::StateVector<double> myvec("cellvars", Jali::CELL,
                                  source_mesh.get(), &(data[0]));
  source_state.add(myvec);

  // Create Remap object

  Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Portage::Jali_State_Wrapper sourceStateWrapper(source_state);

  Portage::Remap_1stOrder<Portage::Jali_Mesh_Wrapper,
                          Portage::Jali_State_Wrapper,
                          Portage::CELL>
      remapper(sourceMeshWrapper, sourceStateWrapper, "cellvars");

  // Gather the cell coordinates for source and target meshes for intersection

  std::vector<std::vector<JaliGeometry::Point>>
      source_cell_coords(ncells_source);
  std::vector<std::vector<JaliGeometry::Point>>
      target_cell_coords(ncells_target);

  for (int c = 0; c < ncells_source; ++c)
    source_mesh->cell_get_coordinates(c, &(source_cell_coords[c]));
  for (int c = 0; c < ncells_target; ++c)
    target_mesh->cell_get_coordinates(c, &(target_cell_coords[c]));

  // Remap from source to target mesh

  std::vector<double> outvals(ncells_target);

  for (int c = 0; c < ncells_target; ++c) {
    std::vector<int> xcells;
    std::vector<std::vector<double>> xwts;

    BOX_INTERSECT::intersection_moments(target_cell_coords[c],
                                        source_cell_coords,
                                        &xcells, &xwts);

    std::pair< std::vector<int> const &,
               std::vector< std::vector<double> > const & >
        cells_and_weights(xcells, xwts);

    outvals[c] = remapper(cells_and_weights);
  }

  // Make sure we retrieved the correct value for each cell on the target
  const double stdval = data[0];
  for (int c = 0; c < ncells_target; ++c)
    ASSERT_DOUBLE_EQ(stdval, outvals[c]);
}


/// First order remap of linear cell-centered field in 2D

TEST(Remap_1st_Order, Cell_Ctr_Lin_2D) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);

  std::unique_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4);
  std::unique_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 1.0, 1.0, 2, 2);

  const int ncells_source = source_mesh->num_entities(Jali::CELL, Jali::OWNED);
  const int ncells_target = target_mesh->num_entities(Jali::CELL, Jali::OWNED);

  // Create a state object

  Jali::State source_state(source_mesh.get());

  // Define a state vector representing a linear function

  std::vector<double> data(ncells_source);
  for (int c = 0; c < ncells_source; ++c) {
    JaliGeometry::Point cen;
    source_mesh->cell_centroid(c, &cen);
    data[c] = cen[0]+cen[1];
  }
  Jali::StateVector<double> myvec("cellvars", Jali::CELL,
                                  source_mesh.get(), &(data[0]));
  source_state.add(myvec);

  // Create Remap objects

  Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Portage::Jali_State_Wrapper sourceStateWrapper(source_state);

  Portage::Remap_1stOrder<Portage::Jali_Mesh_Wrapper,
                          Portage::Jali_State_Wrapper,
                          Portage::CELL>
      remapper(sourceMeshWrapper, sourceStateWrapper, "cellvars");

  // Gather the cell coordinates for source and target meshes for intersection

  std::vector<std::vector<JaliGeometry::Point>>
      source_cell_coords(ncells_source);
  std::vector<std::vector<JaliGeometry::Point>>
      target_cell_coords(ncells_target);

  for (int c = 0; c < ncells_source; ++c)
    source_mesh->cell_get_coordinates(c, &(source_cell_coords[c]));
  for (int c = 0; c < ncells_target; ++c)
    target_mesh->cell_get_coordinates(c, &(target_cell_coords[c]));

  // Remap from source to target mesh

  std::vector<double> outvals(ncells_target);

  for (int c = 0; c < ncells_target; ++c) {
    std::vector<int> xcells;
    std::vector<std::vector<double>> xwts;

    BOX_INTERSECT::intersection_moments(target_cell_coords[c],
                                        source_cell_coords,
                                        &xcells, &xwts);

    std::pair< std::vector<int> const &,
               std::vector< std::vector<double> > const & >
        cells_and_weights(xcells, xwts);

    outvals[c] = remapper(cells_and_weights);
  }

  // Make sure we retrieved the correct value for each cell on the target
  // NOTE: EVEN THOUGH 1ST ORDER REMAPPING ALGORITHM DOES NOT IN
  // GENERAL PRESERVE A LINEAR FIELD THE SPECIAL STRUCTURE OF THE SOURCE
  // AND TARGET MESHES ENSURES THAT THE LINEAR FIELD IS REMAPPED CORRECTLY
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


/// First order remap of constant node-centered field in 2D

TEST(Remap_1st_Order, Node_Ctr_Const_2D) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);

  std::unique_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4, NULL,
                                               true, true, true, true);
  std::unique_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 1.0, 1.0, 5, 5, NULL,
                                               true, true, true, true);

  const int nnodes_source = source_mesh->num_entities(Jali::NODE, Jali::OWNED);
  const int nnodes_target = target_mesh->num_entities(Jali::NODE, Jali::OWNED);


  // Create a state object and add the first two vectors to it

  Jali::State source_state(source_mesh.get());

  // Define two state vectors, one with constant value, the other
  // with a linear function

  std::vector<double> data(nnodes_source, 1.5);
  Jali::StateVector<double> myvec("nodevars", Jali::NODE,
                                  source_mesh.get(), &(data[0]));
  source_state.add(myvec);

  // Create Remap objects

  Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Portage::Jali_State_Wrapper sourceStateWrapper(source_state);

  Portage::Remap_1stOrder<Portage::Jali_Mesh_Wrapper,
                          Portage::Jali_State_Wrapper,
                          Portage::NODE>
      remapper(sourceMeshWrapper, sourceStateWrapper, "nodevars");

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

  for (int n = 0; n < nnodes_source; ++n) {
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

  for (int n = 0; n < nnodes_target; ++n) {
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
               std::vector< std::vector<double> > const & >
        nodes_and_weights(xcells, xwts);

    outvals[n] = remapper(nodes_and_weights);
  }

  // Make sure we retrieved the correct value for each cell on the target

  const double stdval = data[0];
  for (int n = 0; n < nnodes_target; ++n)
    ASSERT_DOUBLE_EQ(stdval, outvals[n]);
}


/// First order remap of constant cell-centered field in 3D

TEST(Remap_1st_Order, Cell_Ctr_Const_3D) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);

  std::unique_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                               4, 4, 4);
  std::unique_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                               5, 5, 5);

  const int ncells_source = source_mesh->num_entities(Jali::CELL, Jali::OWNED);
  const int ncells_target = target_mesh->num_entities(Jali::CELL, Jali::OWNED);

  // Create a state object

  Jali::State source_state(source_mesh.get());

  // Define a state vectors with constant value

  std::vector<double> data(ncells_source, 1.25);
  Jali::StateVector<double> myvec("cellvars", Jali::CELL,
                                  source_mesh.get(), &(data[0]));
  source_state.add(myvec);

  // Create Remap object

  Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Portage::Jali_State_Wrapper sourceStateWrapper(source_state);

  Portage::Remap_1stOrder<Portage::Jali_Mesh_Wrapper,
                          Portage::Jali_State_Wrapper,
                          Portage::CELL>
      remapper(sourceMeshWrapper, sourceStateWrapper, "cellvars");

  // Gather the cell coordinates for source and target meshes for intersection

  std::vector<std::vector<JaliGeometry::Point>>
      source_cell_coords(ncells_source);
  std::vector<std::vector<JaliGeometry::Point>>
      target_cell_coords(ncells_target);

  for (int c = 0; c < ncells_source; ++c)
    source_mesh->cell_get_coordinates(c, &(source_cell_coords[c]));
  for (int c = 0; c < ncells_target; ++c)
    target_mesh->cell_get_coordinates(c, &(target_cell_coords[c]));

  // Remap from source to target mesh

  std::vector<double> outvals(ncells_target);

  for (int c = 0; c < ncells_target; ++c) {
    std::vector<int> xcells;
    std::vector<std::vector<double>> xwts;

    BOX_INTERSECT::intersection_moments(target_cell_coords[c],
                                        source_cell_coords,
                                        &xcells, &xwts);

    std::pair< std::vector<int> const &,
               std::vector< std::vector<double> > const & >
        cells_and_weights(xcells, xwts);

    outvals[c] = remapper(cells_and_weights);
}

  // Make sure we retrieved the correct value for each cell on the target
  const double stdval = data[0];
  for (int c = 0; c < ncells_target; ++c)
    ASSERT_DOUBLE_EQ(stdval, outvals[c]);
}


/// First order remap of linear cell-centered field in 3D

TEST(Remap_1st_Order, Cell_Ctr_Lin_3D) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);

  std::unique_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                               4, 4, 4);
  std::unique_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                               2, 2, 2);

  const int ncells_source = source_mesh->num_entities(Jali::CELL, Jali::OWNED);
  const int ncells_target = target_mesh->num_entities(Jali::CELL, Jali::OWNED);

  // Create a state object

  Jali::State source_state(source_mesh.get());

  // Define a state vector representing a linear function

  std::vector<double> data(ncells_source);
  for (int c = 0; c < ncells_source; ++c) {
    JaliGeometry::Point cen;
    source_mesh->cell_centroid(c, &cen);
    data[c] = cen[0]+cen[1];
  }
  Jali::StateVector<double> myvec("cellvars", Jali::CELL,
                                  source_mesh.get(), &(data[0]));
  source_state.add(myvec);

  // Create Remap objects

  Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Portage::Jali_State_Wrapper sourceStateWrapper(source_state);

  Portage::Remap_1stOrder<Portage::Jali_Mesh_Wrapper,
                          Portage::Jali_State_Wrapper,
                          Portage::CELL>
      remapper(sourceMeshWrapper, sourceStateWrapper, "cellvars");

  // Gather the cell coordinates for source and target meshes for intersection

  std::vector<std::vector<JaliGeometry::Point>>
      source_cell_coords(ncells_source);
  std::vector<std::vector<JaliGeometry::Point>>
      target_cell_coords(ncells_target);

  for (int c = 0; c < ncells_source; ++c)
    source_mesh->cell_get_coordinates(c, &(source_cell_coords[c]));
  for (int c = 0; c < ncells_target; ++c)
    target_mesh->cell_get_coordinates(c, &(target_cell_coords[c]));

  // Remap from source to target mesh

  std::vector<double> outvals(ncells_target);

  for (int c = 0; c < ncells_target; ++c) {
    std::vector<int> xcells;
    std::vector<std::vector<double>> xwts;

    BOX_INTERSECT::intersection_moments(target_cell_coords[c],
                                        source_cell_coords,
                                        &xcells, &xwts);

    std::pair< std::vector<int> const &,
               std::vector< std::vector<double> > const & >
        cells_and_weights(xcells, xwts);

    outvals[c] = remapper(cells_and_weights);
  }

  // Make sure we retrieved the correct value for each cell on the target
  // NOTE: EVEN THOUGH 1ST ORDER REMAPPING ALGORITHM DOES NOT IN
  // GENERAL PRESERVE A LINEAR FIELD THE SPECIAL STRUCTURE OF THE SOURCE
  // AND TARGET MESHES ENSURES THAT THE LINEAR FIELD IS REMAPPED CORRECTLY
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


/// First order remap of constant node-centered field in 3D

TEST(Remap_1st_Order, Node_Ctr_Const_3D) {
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

  const int nnodes_source = source_mesh->num_entities(Jali::NODE, Jali::OWNED);
  const int nnodes_target = target_mesh->num_entities(Jali::NODE, Jali::OWNED);

  // Create a state object and add the first two vectors to it

  Jali::State source_state(source_mesh.get());

  // Define two state vectors, one with constant value, the other
  // with a linear function

  std::vector<double> data(nnodes_source, 1.5);
  Jali::StateVector<double> myvec("nodevars", Jali::NODE,
                                  source_mesh.get(), &(data[0]));
  source_state.add(myvec);

  // Create Remap objects

  Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Portage::Jali_State_Wrapper sourceStateWrapper(source_state);

  Portage::Remap_1stOrder<Portage::Jali_Mesh_Wrapper,
                          Portage::Jali_State_Wrapper,
                          Portage::NODE>
      remapper(sourceMeshWrapper, sourceStateWrapper, "nodevars");

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

  for (int n = 0; n < nnodes_source; ++n) {
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

  for (int n = 0; n < nnodes_target; ++n) {
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
               std::vector< std::vector<double> > const & >
        nodes_and_weights(xcells, xwts);

    outvals[n] = remapper(nodes_and_weights);
  }

  // Make sure we retrieved the correct value for each cell on the target

  const double stdval = data[0];
  for (int n = 0; n < nnodes_target; ++n)
    ASSERT_DOUBLE_EQ(stdval, outvals[n]);
}

