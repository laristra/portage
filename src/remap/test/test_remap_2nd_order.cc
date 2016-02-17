/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/
#include <iostream>
#include <memory>

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
#include "portage/remap/test/simple_intersect_for_tests.h"

// Remap of constant, cell-centered field with no limiting - 2D

TEST(Remap_2nd_Order, Cell_Ctr_Const_No_Limiter_2D) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);

  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);

  auto source_mesh = std::unique_ptr<Jali::Mesh>(mf(0.0, 0.0, 1.0, 1.0, 4, 4));
  auto target_mesh = std::unique_ptr<Jali::Mesh>(mf(0.0, 0.0, 1.0, 1.0, 5, 5));
  Jali::State source_state(source_mesh.get());
  Jali::State target_state(target_mesh.get());

  // Define two state vectors, one with constant value and the other
  // with a linear function that is x+y

  const int ncells_source = source_mesh->num_entities(Jali::CELL, Jali::OWNED);
  std::vector<double> data(ncells_source, 1.25);
  Jali::StateVector<double> myvec("cellvars", Jali::CELL,
                                  source_mesh.get(), &(data[0]));
  source_state.add(myvec);

  const int ncells_target = target_mesh->num_entities(Jali::CELL, Jali::OWNED);

  // Create Remap objects
  Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Portage::Jali_State_Wrapper sourceStateWrapper(source_state);

  Portage::Remap_2ndOrder<Portage::Jali_Mesh_Wrapper,
                          Portage::Jali_State_Wrapper,
                          Portage::Entity_kind>
      remapper(sourceMeshWrapper, sourceStateWrapper, Portage::CELL, "cellvars",
               Portage::NOLIMITER);

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

  std::vector<double> outvals(ncells_target);  // field values on target mesh
  for (int c = 0; c < ncells_target; ++c) {
    std::vector<int> xcells;
    std::vector<std::vector<double>> xweights;

    intersection_moments(target_cell_coords[c], source_cell_coords, &xcells,
                          &xweights);

    std::pair< std::vector<int> const &,
               std::vector< std::vector<double> > const & >
        cells_and_weights(xcells, xweights);

    outvals[c] = remapper(cells_and_weights);
  }

  // Make sure we retrieved the correct value for each cell on the target

  const double stdval = data[0];
  for (int c = 0; c < ncells_target; ++c)
    ASSERT_NEAR(stdval, outvals[c], 1.0e-10);
}


// Remap of linear, cell-centered field with no limiting - 2D

TEST(Remap_2nd_Order, Cell_Ctr_Lin_No_Limiter_2D) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);

  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);

  auto source_mesh = std::unique_ptr<Jali::Mesh>(mf(0.0, 0.0, 1.0, 1.0, 4, 4));
  auto target_mesh = std::unique_ptr<Jali::Mesh>(mf(0.0, 0.0, 1.0, 1.0, 5, 5));
  Jali::State source_state(source_mesh.get());
  Jali::State target_state(target_mesh.get());

  // Define a state vectors, with a linear function that is x+y

  const int ncells_source = source_mesh->num_entities(Jali::CELL, Jali::OWNED);
  std::vector<double> data(ncells_source);
  for (int c = 0; c < ncells_source; ++c) {
    JaliGeometry::Point ccen = source_mesh->cell_centroid(c);
    data[c] = ccen[0]+ccen[1];
  }

  Jali::StateVector<double> myvec("cellvars", Jali::CELL,
                                  source_mesh.get(), &(data[0]));
  source_state.add(myvec);

  // Create Remap objects

  Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Portage::Jali_State_Wrapper sourceStateWrapper(source_state);

  Portage::Remap_2ndOrder<Portage::Jali_Mesh_Wrapper,
                          Portage::Jali_State_Wrapper,
                          Portage::Entity_kind>
      remapper(sourceMeshWrapper, sourceStateWrapper, Portage::CELL, "cellvars",
               Portage::NOLIMITER);

  // Gather the cell coordinates for source and target meshes for intersection

  const int ncells_target = target_mesh->num_entities(Jali::CELL, Jali::OWNED);
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
    std::vector<std::vector<double>> xweights;

    intersection_moments(target_cell_coords[c], source_cell_coords,
                         &xcells, &xweights);

    // std::cerr << "Target Cell " << c << ":" << std::endl;
    // for (int i = 0; i < xcells.size(); i++)
    //   std::cerr << "  Source Cell " << xcells[i] <<
    //       "  Xsection Vol " << xweights[i][0] << " Xsection Centroid " <<
    //       xweights[i][1]/xweights[i][0] << "," <<
    //       xweights[i][2]/xweights[i][0] <<
    //       std::endl;

    std::pair< std::vector<int> const &,
               std::vector< std::vector<double> > const & >
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


// Remap of constant, cell-centered field with Barth-Jespersen limiting - 2D

TEST(Remap_2nd_Order, Cell_Ctr_Lin_BJ_Limiter_2D) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);

  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);

  auto source_mesh = std::unique_ptr<Jali::Mesh>(mf(0.0, 0.0, 1.0, 1.0, 4, 4));
  auto target_mesh = std::unique_ptr<Jali::Mesh>(mf(0.0, 0.0, 1.0, 1.0, 5, 5));
  Jali::State source_state(source_mesh.get());
  Jali::State target_state(target_mesh.get());

  // Define a state vectors, with a linear function that is x+y

  int ncells_source = source_mesh->num_entities(Jali::CELL, Jali::OWNED);
  std::vector<double> data(ncells_source);
  const double minval = 1e+10, maxval = -1e+10;
  for (int c = 0; c < ncells_source; ++c) {
    JaliGeometry::Point ccen = source_mesh->cell_centroid(c);
    if (ccen[0] < 0.5)
      data[c] = ccen[0]+ccen[1];
    else
      data[c] = 100*ccen[0];
    if (data[c] < minval) minval = data[c];
    if (data[c] > maxval) maxval = data[c];
  }

  Jali::StateVector<double> myvec("cellvars", Jali::CELL,
                                  source_mes.get()h, &(data[0]));
  source_state.add(myvec);

  // Create Remap objects - one with no limiter and one with limiter
  Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Portage::Jali_State_Wrapper sourceStateWrapper(source_state);

  Portage::Remap_2ndOrder<Portage::Jali_Mesh_Wrapper,
                          Portage::Jali_State_Wrapper,
                          Portage::Entity_kind>
      remapper1(sourceMeshWrapper, sourceStateWrapper, Portage::CELL,
                "cellvars", Portage::NOLIMITER);
  Portage::Remap_2ndOrder<Portage::Jali_Mesh_Wrapper,
                          Portage::Jali_State_Wrapper,
                          Portage::Entity_kind>
      remapper2(sourceMeshWrapper, sourceStateWrapper, Portage::CELL,
                "cellvars", Portage::BARTH_JESPERSEN);

  // Gather the cell coordinates for the source and target meshes for
  // intersection

  const int ncells_target = target_mesh->num_entities(Jali::CELL, Jali::OWNED);
  std::vector<std::vector<JaliGeometry::Point>>
      source_cell_coords(ncells_source);
  std::vector<std::vector<JaliGeometry::Point>>
      target_cell_coords(ncells_target);

  for (int c = 0; c < ncells_source; ++c)
    source_mesh->cell_get_coordinates(c, &(source_cell_coords[c]));
  for (int c = 0; c < ncells_target; ++c)
    target_mesh->cell_get_coordinates(c, &(target_cell_coords[c]));

  std::vector<double> outvals1(ncells_target);
  std::vector<double> outvals2(ncells_target);

  // Remap from source to target mesh
  for (int c = 0; c < ncells_target; ++c) {
    std::vector<int> xcells;
    std::vector<std::vector<double>> xweights;

    intersection_moments(target_cell_coords[c], source_cell_coords,
                         &xcells, &xweights);

    // std::cerr << "Target Cell " << c << ":" << std::endl;
    // for (int i = 0; i < xcells.size(); i++)
    //   std::cerr << "  Source Cell " << xcells[i] <<
    //       "  Xsection Vol " << xweights[i][0] << " Xsection Centroid " <<
    //       xweights[i][1]/xweights[i][0] << "," <<
    //       xweights[i][2]/xweights[i][0] <<
    //       std::endl;

    std::pair< std::vector<int> const &,
               std::vector< std::vector<double> > const & >
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


// **** WILL NOT WORK UNTIL WE ELIMINATE USE OF MeshDualWrapper AND
// FIX THE GRADIENT AND REMAP CODE TO WORK EXPLICITLY WITH NODE
// CENTERED QUANTITIES *****


// TEST(Remap_2nd_Order, Node_Ctr_Const_No_Limiter) {

//   Jali::MeshFactory mf(MPI_COMM_WORLD);
//   Jali::FrameworkPreference pref;
//   pref.push_back(Jali::MSTK);
//   if (Jali::framework_available(Jali::MSTK))
//     mf.preference(pref);

//   Jali::Mesh *source_mesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4,
//                                NULL, true, true, true, true);
//   Jali::Mesh *target_mesh = mf(0.0, 0.0, 1.0, 1.0, 5, 5,
//                                NULL, true, true, true, true);

//   int nnodes_source = source_mesh->num_entities(Jali::NODE, Jali::OWNED);
//   int nnodes_target = target_mesh->num_entities(Jali::NODE, Jali::OWNED);

//   // Create a state object and add the first two vectors to it

//   Jali::State source_state(source_mesh);


//   // Define two state vectors, one with constant value, the other
//   // with a linear function

//   std::vector<double> data(nnodes_source, 1.5);
//   Jali::StateVector<double> myvec("nodevars", Jali::NODE, source_mesh,
//                                   &(data[0]));
//   Jali::StateVector<double> &addvec = source_state.add(myvec);

//   // Create Remap objects
//   Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
//   Portage::Jali_State_Wrapper sourceStateWrapper(source_state);

//   Portage::Remap_2ndOrder<Portage::Jali_Mesh_Wrapper,
//                           Portage::Jali_State_Wrapper,
//                           Portage::Entity_kind>
//       remapper(sourceMeshWrapper, sourceStateWrapper,
//                Portage::NODE, "nodevars",
//                Portage::NOLIMITER);

//   // Remap from source to target mesh

//   std::vector<double> outvals(nnodes_target);


//   // Gather the dual cell coordinates for source and target meshes for
//      intersection

//   std::vector<std::vector<JaliGeometry::Point>>
//                 source_dualcell_coords(nnodes_source);
//   std::vector<std::vector<JaliGeometry::Point>>
//                 target_dualcell_coords(nnodes_target);

//   // Because the meshes are rectangular we can get away with examining
//   // the coordinates of the corners instead of the wedges

//   // Also, because we will use only the bounding box of each cell to
//   // do the search and intersection we can get away with adding all
//   // the coordinates of the corners to list including duplicates

//   for (int n = 0; n < nnodes_source; n++) {
//     std::vector<JaliGeometry::Point> dualcoords;
//     std::vector<int> corners;
//     source_mesh->node_get_corners(n, Jali::ALL,&corners);

//     for (auto cn : corners) {
//       std::vector<JaliGeometry::Point> cncoords;
//       source_mesh->corner_get_coordinates(cn,&cncoords);
//       for (auto coord : cncoords)
//         source_dualcell_coords[n].push_back(coord);
//     }
//   }

//   for (int n = 0; n < nnodes_target; n++) {
//     std::vector<JaliGeometry::Point> dualcoords;
//     std::vector<int> corners;
//     target_mesh->node_get_corners(n, Jali::ALL,&corners);

//     for (auto cn : corners) {
//       std::vector<JaliGeometry::Point> cncoords;
//       target_mesh->corner_get_coordinates(cn,&cncoords);
//       for (auto coord : cncoords)
//         target_dualcell_coords[n].push_back(coord);
//     }
//   }


//   for (int n = 0; n < nnodes_target; ++n) {
//     std::vector<int> xcells;
//     std::vector<std::vector<double>> xwts;

//     intersection_moments(target_dualcell_coords[n], source_dualcell_coords,
//                          &xcells, &xwts);

//     std::pair< std::vector<int> const &,
//                std::vector< std::vector<double> > const & >
//         nodes_and_weights(xcells, xwts);

//     outvals[n] = remapper(nodes_and_weights);
//   }

//   // Make sure we retrieved the correct value for each cell on the target

//   double stdval = data[0];
//   for (int n = 0; n < nnodes_target; ++n)
//     ASSERT_DOUBLE_EQ(stdval, outvals[n]);

// }


// TEST(Remap_2nd_Order, Node_Ctr_Lin) {

//   Jali::MeshFactory mf(MPI_COMM_WORLD);
//   Jali::FrameworkPreference pref;
//   pref.push_back(Jali::MSTK);
//   if (Jali::framework_available(Jali::MSTK))
//     mf.preference(pref);

//   Jali::Mesh *source_mesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4,
//                                NULL, true, true, true, true);
//   Jali::Mesh *target_mesh = mf(0.0, 0.0, 1.0, 1.0, 5, 5,
//                                NULL, true, true, true, true);

//   int nnodes_source = source_mesh->num_entities(Jali::NODE, Jali::OWNED);
//   int nnodes_target = target_mesh->num_entities(Jali::NODE, Jali::OWNED);

//   // Create a state object and add the first two vectors to it

//   Jali::State source_state(source_mesh);


//   // Define two state vectors, one with constant value, the other
//   // with a linear function

//   std::vector<double> data(nnodes_source);
//   for (int n = 0; n < nnodes_source; n++) {
//     JaliGeometry::Point coord;
//     source_mesh->node_get_coordinates(n,&coord);
//     data[n] = coord[0]+coord[1];
//   }
//   Jali::StateVector<double> myvec("nodevars", Jali::NODE,
//                                   source_mesh,&(data[0]));
//   Jali::StateVector<double> &addvec = source_state.add(myvec);


//   // Create a Dual mesh wrapper

//   Portage::MeshWrapper sourceMeshWrapper(source_mesh);
//   Portage::MeshWrapper targetMeshWrapper(target_mesh);

//   // Create Remap objects

//   Portage::Remap_2ndOrder<Portage::Jali_Mesh_Wrapper,
//                           Portage::Jali_State_Wrapper,
//                           Portage::Entity_kind>
//       remapper(sourceMeshWrapper, source_state, Portage::NODE,"nodevars",
//                Portage::NOLIMITER);


//   // Gather the dual cell coordinates for source and target meshes for
//   // intersection

//   std::vector<std::vector<JaliGeometry::Point>>
//                      source_dualcell_coords(nnodes_source);
//   std::vector<std::vector<JaliGeometry::Point>>
//                      target_dualcell_coords(nnodes_target);

//   // Because the meshes are rectangular we can get away with examining
//   // the coordinates of the corners instead of the wedges

//   // Also, because we will use only the bounding box of each cell to
//   // do the search and intersection we can get away with adding all
//   // the coordinates of the corners to list including duplicates

//   for (int n = 0; n < nnodes_source; n++) {
//     std::vector<JaliGeometry::Point> dualcoords;
//     std::vector<int> corners;
//     source_mesh->node_get_corners(n, Jali::ALL,&corners);

//     for (auto cn : corners) {
//       std::vector<JaliGeometry::Point> cncoords;
//       source_mesh->corner_get_coordinates(cn,&cncoords);
//       for (auto coord : cncoords)
//         source_dualcell_coords[n].push_back(coord);
//     }
//   }

//   for (int n = 0; n < nnodes_target; n++) {
//     std::vector<JaliGeometry::Point> dualcoords;
//     std::vector<int> corners;
//     target_mesh->node_get_corners(n, Jali::ALL,&corners);

//     for (auto cn : corners) {
//       std::vector<JaliGeometry::Point> cncoords;
//       target_mesh->corner_get_coordinates(cn,&cncoords);
//       for (auto coord : cncoords)
//         target_dualcell_coords[n].push_back(coord);
//     }
//   }


//   // Remap from source to target mesh

//   std::vector<double> outvals(nnodes_target);

//   for (int n = 0; n < nnodes_target; ++n) {
//     std::vector<int> xcells;
//     std::vector<std::vector<double>> xwts;

//     intersection_moments(target_dualcell_coords[n], source_dualcell_coords,
//                          &xcells, &xwts);

//     std::pair< std::vector<int> const &,
//                std::vector< std::vector<double> > const & >
//         nodes_and_weights(xcells, xwts);

//     outvals[n] = remapper(nodes_and_weights);
//   }

//   // Make sure we retrieved the correct value for each cell on the target
//   // For field 1, it is a constant
//   // For field 2, it is a linear function
//   // NOTE: Even though 1st order remapping algorithm does not in
//   // general preserve a linear field, the special structure of the
//   // source and target mesh ensures that the linear field is remapped
//   // correctly in this test

//   std::vector<double> stdvals(nnodes_target);
//   for (int n = 0; n < nnodes_target; ++n) {
//     JaliGeometry::Point coord;
//     target_mesh->node_get_coordinates(n,&coord);
//     stdvals[n] = coord[0]+coord[1];
//   }
//   for (int n = 0; n < nnodes_target; ++n)
//     EXPECT_DOUBLE_EQ(stdvals[n], outvals[n]);

// }


// Remap of constant, cell-centered field with no limiting - 3D

TEST(Remap_2nd_Order, Cell_Ctr_Const_No_Limiter_3D) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);

  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);

  auto source_mesh = std::unique_ptr<Jali::Mesh>(mf(0.0, 0.0, 0.0,
                                                    1.0, 1.0, 1.0,
                                                    4, 4, 4));
  auto target_mesh = std::unique_ptr<Jali::Mesh>(mf(0.0, 0.0, 0.0,
                                                    1.0, 1.0, 1.0,
                                                    5, 5, 5));
  Jali::State source_state(source_mesh.get());
  Jali::State target_state(target_mesh.get());

  // Define two state vectors, one with constant value and the other
  // with a linear function that is x+y

  const int ncells_source = source_mesh->num_entities(Jali::CELL, Jali::OWNED);
  std::vector<double> data(ncells_source, 1.25);
  Jali::StateVector<double> myvec("cellvars", Jali::CELL,
                                  source_mesh.get(), &(data[0]));
  source_state.add(myvec);

  const int ncells_target = target_mesh->num_entities(Jali::CELL, Jali::OWNED);

  // Create Remap objects
  Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Portage::Jali_State_Wrapper sourceStateWrapper(source_state);

  Portage::Remap_2ndOrder<Portage::Jali_Mesh_Wrapper,
                          Portage::Jali_State_Wrapper,
                          Portage::Entity_kind>
      remapper(sourceMeshWrapper, sourceStateWrapper, Portage::CELL, "cellvars",
               Portage::NOLIMITER);

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

  std::vector<double> outvals(ncells_target);  // field values on target mesh
  for (int c = 0; c < ncells_target; ++c) {
    std::vector<int> xcells;
    std::vector<std::vector<double>> xweights;

    intersection_moments(target_cell_coords[c], source_cell_coords,
                         &xcells, &xweights);

    std::pair< std::vector<int> const &,
               std::vector< std::vector<double> > const & >
        cells_and_weights(xcells, xweights);

    outvals[c] = remapper(cells_and_weights);
  }

  // Make sure we retrieved the correct value for each cell on the target

  const double stdval = data[0];
  for (int c = 0; c < ncells_target; ++c)
    ASSERT_NEAR(stdval, outvals[c], 1.0e-10);
}


// Remap of linear, cell-centered field with no limiting - 3D

TEST(Remap_2nd_Order, Cell_Ctr_Lin_No_Limiter_3D) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);

  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);

  auto source_mesh = std::unique_ptr<Jali::Mesh>(mf(0.0, 0.0, 0.0,
                                                    1.0, 1.0, 1.0,
                                                    4, 4, 4));
  auto target_mesh = std::unique_ptr<Jali::Mesh>(mf(0.0, 0.0, 0.0,
                                                    1.0, 1.0, 1.0,
                                                    5, 5, 5));
  Jali::State source_state(source_mesh.get());
  Jali::State target_state(target_mesh.get());

  // Define a state vectors, with a linear function that is x+y

  const int ncells_source = source_mesh->num_entities(Jali::CELL, Jali::OWNED);
  std::vector<double> data(ncells_source);
  for (int c = 0; c < ncells_source; ++c) {
    JaliGeometry::Point ccen = source_mesh->cell_centroid(c);
    data[c] = ccen[0]+ccen[1]+ccen[2];
  }

  Jali::StateVector<double> myvec("cellvars", Jali::CELL,
                                  source_mesh.get(), &(data[0]));
  source_state.add(myvec);

  // Create Remap objects

  Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Portage::Jali_State_Wrapper sourceStateWrapper(source_state);

  Portage::Remap_2ndOrder<Portage::Jali_Mesh_Wrapper,
                          Portage::Jali_State_Wrapper,
                          Portage::Entity_kind>
      remapper(sourceMeshWrapper, sourceStateWrapper, Portage::CELL,
               "cellvars", Portage::NOLIMITER);

  // Gather the cell coordinates for source and target meshes for intersection

  const int ncells_target = target_mesh->num_entities(Jali::CELL, Jali::OWNED);
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
    std::vector<std::vector<double>> xweights;

    intersection_moments(target_cell_coords[c], source_cell_coords,
                         &xcells, &xweights);

    // std::cerr << "Target Cell " << c << ":" << std::endl;
    // for (int i = 0; i < xcells.size(); i++)
    //   std::cerr << "  Source Cell " << xcells[i] <<
    //       "  Xsection Vol " << xweights[i][0] << " Xsection Centroid " <<
    //       xweights[i][1]/xweights[i][0] << "," <<
    //       xweights[i][2]/xweights[i][0] <<
    //       std::endl;

    std::pair< std::vector<int> const &,
               std::vector< std::vector<double> > const & >
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


// Remap of linear, cell-centered field with Barth-Jespersen limiting - 3D

TEST(Remap_2nd_Order, Cell_Ctr_Lin_BJ_Limiter_3D) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);

  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);

  auto source_mesh = std::unique_ptr<Jali::Mesh>(mf(0.0, 0.0, 0.0,
                                                    1.0, 1.0, 1.0,
                                                    4, 4, 4));
  auto target_mesh = std::unique_ptr<Jali::Mesh>(mf(0.0, 0.0, 0.0,
                                                    1.0, 1.0, 1.0,
                                                    5, 5, 5));
  Jali::State source_state(source_mesh.get());
  Jali::State target_state(target_mesh.get());

  // Define a state vectors, with a linear function that is x+y

  const int ncells_source = source_mesh->num_entities(Jali::CELL, Jali::OWNED);
  std::vector<double> data(ncells_source);
  const double minval = 1e+10, maxval = -1e+10;
  for (int c = 0; c < ncells_source; ++c) {
    JaliGeometry::Point ccen = source_mesh->cell_centroid(c);
    if (ccen[0] < 0.5)
      data[c] = ccen[0]+ccen[1]+ccen[2];
    else
      data[c] = 100*ccen[0];
    if (data[c] < minval) minval = data[c];
    if (data[c] > maxval) maxval = data[c];
  }

  Jali::StateVector<double> myvec("cellvars", Jali::CELL,
                                  source_mesh.get(), &(data[0]));
  source_state.add(myvec);

  // Create Remap objects - one with no limiter and one with limiter
  Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Portage::Jali_State_Wrapper sourceStateWrapper(source_state);

  Portage::Remap_2ndOrder<Portage::Jali_Mesh_Wrapper,
                          Portage::Jali_State_Wrapper,
                          Portage::Entity_kind>
      remapper1(sourceMeshWrapper, sourceStateWrapper, Portage::CELL,
                "cellvars", Portage::NOLIMITER);
  Portage::Remap_2ndOrder<Portage::Jali_Mesh_Wrapper,
                          Portage::Jali_State_Wrapper,
                          Portage::Entity_kind>
      remapper2(sourceMeshWrapper, sourceStateWrapper, Portage::CELL,
                "cellvars", Portage::BARTH_JESPERSEN);

  // Gather the cell coordinates for the source and target meshes for
  // intersection

  const int ncells_target = target_mesh->num_entities(Jali::CELL, Jali::OWNED);
  std::vector<std::vector<JaliGeometry::Point>>
      source_cell_coords(ncells_source);
  std::vector<std::vector<JaliGeometry::Point>>
      target_cell_coords(ncells_target);

  for (int c = 0; c < ncells_source; ++c)
    source_mesh->cell_get_coordinates(c, &(source_cell_coords[c]));
  for (int c = 0; c < ncells_target; ++c)
    target_mesh->cell_get_coordinates(c, &(target_cell_coords[c]));

  std::vector<double> outvals1(ncells_target);
  std::vector<double> outvals2(ncells_target);

  // Remap from source to target mesh
  for (int c = 0; c < ncells_target; ++c) {
    std::vector<int> xcells;
    std::vector<std::vector<double>> xweights;

    intersection_moments(target_cell_coords[c], source_cell_coords,
                         &xcells, &xweights);

    // std::cerr << "Target Cell " << c << ":" << std::endl;
    // for (int i = 0; i < xcells.size(); i++)
    //   std::cerr << "  Source Cell " << xcells[i] <<
    //       "  Xsection Vol " << xweights[i][0] << " Xsection Centroid " <<
    //       xweights[i][1]/xweights[i][0] << "," <<
    //       xweights[i][2]/xweights[i][0] <<
    //       std::endl;

    std::pair< std::vector<int> const &,
               std::vector< std::vector<double> > const & >
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
    } else if (outvals2[c] > maxval && outvals2[c]-maxval > 1.0e-10) {
      inbounds_limited = false;
      break;
    }
  }

  EXPECT_TRUE(outofbounds_unlimited && inbounds_limited);
}

