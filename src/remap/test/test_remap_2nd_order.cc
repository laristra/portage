/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include "portage/support/portage.h"
#include "portage/remap/remap_2nd_order.h"
#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"
#include "portage/wrappers/state/jali/jali_state_wrapper.h"

#include <iostream>

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

//! \todo LOOK AT TEST MORE CLOSELY, IT WAS PASSING EVEN WITH ERRONEOUS FORMULAS

TEST(Remap_2nd_Order, Fields_Cell_Ctr) {

  // Make a 4x4 source mesh and a 2x2 target mesh - so each cell of
  // the target mesh contains four cells of the source mesh

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);

  Jali::Mesh *mesh1 = mf(0.0,0.0,1.0,1.0,4,4);
  ASSERT_TRUE(mesh1 != NULL);

  Jali::Mesh *mesh2 = mf(0.0,0.0,1.0,1.0,2,2);
  ASSERT_TRUE(mesh2 != NULL);


  // Create a state object and add the first two vectors to it

  Jali::State mystate(mesh1);

  // Define two state vectors, one with constant value and the other
  // with a linear function that is x+y 

  int nc1 = mesh1->num_entities(Jali::CELL,Jali::OWNED);
  std::vector<double> data1(nc1);
  for (int c = 0; c < nc1; ++c)
    data1[c] = 1.25;
  Jali::StateVector<double> myvec1("cellvars1",Jali::CELL,mesh1,&(data1[0]));
  Jali::StateVector<double> &addvec1 = mystate.add(myvec1);

  
  std::vector<double> data2(nc1);
  for (int c = 0; c < nc1; c++) {
    JaliGeometry::Point ccen = mesh1->cell_centroid(c);
    data2[c] = ccen[0]+ccen[1];
  }

  Jali::StateVector<double> myvec2("cellvars2",Jali::CELL,mesh1,&(data2[0]));
  Jali::StateVector<double> &addvec2 = mystate.add(myvec2);

  // Create Remap objects

  Portage::Remap_2ndOrder<Portage::Jali_Mesh_Wrapper,Portage::Jali_State_Wrapper,Portage::Entity_kind> 
      remapper1(*mesh1,mystate,Portage::CELL,"cellvars1",Portage::NOLIMITER);
  Portage::Remap_2ndOrder<Portage::Jali_Mesh_Wrapper,Portage::Jali_State_Wrapper,Portage::Entity_kind> 
      remapper2(*mesh1,mystate,Portage::CELL,"cellvars2",Portage::NOLIMITER);

  // Remap from source to target mesh

  double outvals1[4];  // field values on target mesh
  double outvals2[4];

    // Since we know the structure of the two meshes, we can
    // enumerate which source cells intersect a given target cell and
    // what their intersection areas (weights) are

  std::vector<std::vector<int>> all_source_cells =
      {{0,1,4,5},   {2,3,6,7},
       {8,9,12,13}, {10,11,14,15}};
  std::vector<std::vector<std::vector<double>>> all_weights =
      {{{0.25,0.03125,0.03125}, {0.25,0.03125,0.09375}, 
        {0.25,0.09375,0.03125}, {0.25,0.09375,0.09375}},
       {{0.25,0.03125,0.15625}, {0.25,0.03125,0.21875}, 
        {0.25,0.09375,0.15625}, {0.25,0.09375,0.21875}},
       {{0.25,0.15625,0.03125}, {0.25,0.15625,0.09375}, 
        {0.25,0.21875,0.03125}, {0.25,0.21875,0.09375}},
       {{0.25,0.15625,0.15625}, {0.25,0.15625,0.21875}, 
        {0.25,0.21875,0.15625}, {0.25,0.21875,0.21875}}};

  for (int i = 0; i < 4; ++i) {
    std::pair< std::vector<int> const &, 
               std::vector< std::vector<double> > const & > 
        cells_and_weights(all_source_cells[i],all_weights[i]);
    outvals1[i] = remapper1(cells_and_weights);
    outvals2[i] = remapper2(cells_and_weights);
  }

  // Make sure we retrieved the correct value for each cell on the target
  // For field 1, it is a constant
  // For field 2, it is a linear function

  double stdval1 = data1[0];
  int nc2 = mesh2->num_entities(Jali::CELL,Jali::OWNED);
  std::vector<double> stdvals2(nc2);
  for (int c = 0; c < nc2; ++c) {
    JaliGeometry::Point ccen = mesh2->cell_centroid(c);
    stdvals2[c] = ccen[0]+ccen[1];
  }
  for (int i = 0; i < 4; ++i) {
    ASSERT_NEAR(stdval1, outvals1[i], 1.0e-10);
    ASSERT_NEAR(stdvals2[i], outvals2[i], 1.0e-10);
  }

}


// TEST(Remap_2nd_Order, Fields_Node_Ctr) {

//   // Make a 3x3 source mesh and a 1x1 target mesh - so each node of
//   // the target mesh corresponds to four nodes of the source mesh

//   Jali::MeshFactory mf(MPI_COMM_WORLD);

//   Jali::FrameworkPreference pref;
//   pref.push_back(Jali::MSTK);
//   if (Jali::framework_available(Jali::MSTK))
//     mf.preference(pref);

//   Jali::Mesh *mesh1 = mf(0.0,0.0,1.0,1.0,3,3);
//   ASSERT_TRUE(mesh1 != NULL);

//   Jali::Mesh *mesh2 = mf(0.0,0.0,1.0,1.0,1,1);
//   ASSERT_TRUE(mesh2 != NULL);

//   // Create a state object and add the first two vectors to it

//   Jali::State mystate(mesh1);

//   // Define two state vectors, one with constant value, the other
//   // with a linear function

//   int nn1 = mesh1->num_entities(Jali::CELL,Jali::OWNED);

//   std::vector<double> data1(nn1);
//   for (int n = 0; n < nn1; ++n) data1[n] = 1.5;

//   Jali::StateVector<double> myvec1("nodevars1",Jali::NODE,mesh1,&(data1[0]));
//   Jali::StateVector<double> &addvec1 = mystate.add(myvec1);

//   std::vector<double> data2(nn1);
//   for (int n = 0; n < nn1; ++n) {
//     JaliGeometry::Point nodexy;
//     mesh1->node_get_coordinates(n,&nodexy);
//     data2[n] = nodexy[0]+nodexy[1];
//   }
//   Jali::StateVector<double> myvec2("nodevars2",Jali::NODE,mesh1,&(data2[0]));
//   Jali::StateVector<double> &addvec2 = mystate.add(myvec2);

//   // Create Remap objects

//   Portage::Remap_2ndOrder<Portage::Jali_Mesh_Wrapper,Portage::Jali_State_Wrapper,Portage::Entity_kind> 
//       remapper1(*mesh1,mystate,Portage::NODE,"nodevars1",Portage::NOLIMITER);
//   Portage::Remap_2ndOrder<Portage::Jali_Mesh_Wrapper,Portage::Jali_State_Wrapper,Portage::Entity_kind> 
//       remapper2(*mesh1,mystate,Portage::NODE,"nodevars2",Portage::NOLIMITER);

//   // Remap from source to target mesh

//   double outvals1[4];  // field values on target mesh
//   double outvals2[4];

//   // Since we know the structure of the two meshes, we can
//   // enumerate which source nodes contribute to a given target node
//   // and what their intersection areas (weights) are

//   std::vector<std::vector<int>> all_source_nodes =
//       {{0,1,4,5},   {2,3,6,7},
//        {8,9,12,13}, {10,11,14,15}};
//   std::vector<std::vector<std::vector<double>>> all_weights =
//       {{{1./36, (1./36)*(1./12),  (1./36)*(1./12)}, 
//         {2./36, (2./36)*(1./3),   (2./36)*(1./12)}, 
//         {2./36, (2./36)*(1./12),  (2./36)*(1./3)}, 
//         {4./36, (4./36)*(1./3),   (4./36)*(1./3)}},
//        {{2./36, (2./36)*(2./3),   (2./36)*(1./12)}, 
//         {1./36, (1./36)*(11./12), (1./36)*(1./12)}, 
//         {4./36, (4./36)*(2./3),   (4./36)*(1./3)}, 
//         {2./36, (1./36)*(11./12), (1./36)*(1./3)}},
//        {{2./36, (2./36)*(1./12),  (2./36)*(2./3)},
//         {4./36, (4./36)*(1./3),   (4./36)*(2./3)}, 
//         {1./36, (1./36)*(1./12),  (1./36)*(11./12)}, 
//         {2./36, (2./36)*(1./3),   (2./36)*(11./12)}},
//        {{4./36, (4./36)*(2./3),   (4./36)*(2./3)}, 
//         {2./36, (2./36)*(11./12), (2./36)*(2./3)}, 
//         {2./36, (2./36)*(2./3),   (2./36)*(11./12)}, 
//         {1./36, (1./36)*(11./12), (1./36)*(11./12)}}};

//   int nn2 = mesh2->num_entities(Jali::NODE,Jali::OWNED);
//   for (int i = 0; i < nn2; ++i) {
//     std::pair< std::vector<int> const &, 
//                std::vector< std::vector<double> > const & > 
//         nodes_and_weights(all_source_nodes[i],all_weights[i]);
//     outvals1[i] = remapper1(nodes_and_weights);
//     outvals2[i] = remapper2(nodes_and_weights);
//   }

//   // Make sure we retrieved the correct value for each cell on the target
//   // For field 1, it is a constant
//   // For field 2, it is a linear function

//   double stdval1 = data1[0];
//   std::vector<double> stdvals2(nn2);
//   for (int n = 0; n < nn2; ++n) {
//     JaliGeometry::Point nodexy;
//     mesh2->node_get_coordinates(n,&nodexy);
//     stdvals2[n] = nodexy[0]+nodexy[1];
//   }
//   for (int i = 0; i < 4; ++i) {
//     ASSERT_NEAR(stdval1, outvals1[i], 1.0e-10);
//     ASSERT_NEAR(stdvals2[i], outvals2[i], 1.0e-10);
//   }

// }

  
