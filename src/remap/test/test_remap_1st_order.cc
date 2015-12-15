/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include "portage/support/portage.h"
#include "portage/remap/remap_1st_order.h"
#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"
#include "portage/wrappers/state/jali/jali_state_wrapper.h"

#include <iostream>

#include "gtest/gtest.h"
#include "mpi.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "FrameworkTraits.hh"
#include "MeshFramework.hh"
#include "JaliState.h"
#include "JaliStateVector.h"


TEST(Remap_1st_Order, Fields_Cell_Ctr) {

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

  // Define two state vectors, one with constant value, the other
  // with a linear function

  std::vector<double> data1 = {1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,
                               1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25}; 
  Jali::StateVector<double> myvec1("cellvars1",Jali::CELL,mesh1,&(data1[0]));
  Jali::StateVector<double> &addvec1 = mystate.add(myvec1);

  std::vector<double> data2 = {0.,1.,2.,3., 1.,2.,3.,4.,
                               2.,3.,4.,5., 3.,4.,5.,6.,}; 
  Jali::StateVector<double> myvec2("cellvars2",Jali::CELL,mesh1,&(data2[0]));
  Jali::StateVector<double> &addvec2 = mystate.add(myvec2);

  // Create Remap objects

  Portage::Remap_1stOrder<Portage::Jali_Mesh_Wrapper,Portage::Jali_State_Wrapper,Portage::Entity_kind> 
      remapper1(*mesh1,mystate,Portage::CELL,"cellvars1");
  Portage::Remap_1stOrder<Portage::Jali_Mesh_Wrapper,Portage::Jali_State_Wrapper,Portage::Entity_kind> 
      remapper2(*mesh1,mystate,Portage::CELL,"cellvars2");

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
      {{{0.25}, {0.25}, {0.25}, {0.25}},
       {{0.25}, {0.25}, {0.25}, {0.25}},
       {{0.25}, {0.25}, {0.25}, {0.25}},
       {{0.25}, {0.25}, {0.25}, {0.25}}};

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
  // NOTE: Even though 1st order remapping algorithm does not in
  // general preserve a linear field, the special structure of the
  // source and target mesh ensures that the linear field is remapped
  // correctly in this test

  double stdval1 = data1[0];
  double stdvals2[4] = {1., 3., 3., 5.};
  for (int i = 0; i < 4; ++i) {
    ASSERT_EQ(stdval1, outvals1[i]);
    ASSERT_EQ(stdvals2[i], outvals2[i]);
  }

}


TEST(Remap_1st_Order, Fields_Node_Ctr) {

  // Make a 3x3 source mesh and a 1x1 target mesh - so each node of
  // the target mesh corresponds to four nodes of the source mesh

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);

  Jali::Mesh *mesh1 = mf(0.0,0.0,1.0,1.0,3,3);
  ASSERT_TRUE(mesh1 != NULL);

  Jali::Mesh *mesh2 = mf(0.0,0.0,1.0,1.0,1,1);
  ASSERT_TRUE(mesh2 != NULL);

  // Create a state object and add the first two vectors to it

  Jali::State mystate(mesh1);

  // Define two state vectors, one with constant value, the other
  // with a linear function

  std::vector<double> data1 = {1.5,1.5,1.5,1.5, 1.5,1.5,1.5,1.5,
                               1.5,1.5,1.5,1.5, 1.5,1.5,1.5,1.5}; 
  Jali::StateVector<double> myvec1("nodevars1",Jali::NODE,mesh1,&(data1[0]));
  Jali::StateVector<double> &addvec1 = mystate.add(myvec1);

  std::vector<double> data2 = {0.,1.,2.,3., 1.,2.,3.,4.,
                               2.,3.,4.,5., 3.,4.,5.,6.,}; 
  Jali::StateVector<double> myvec2("nodevars2",Jali::NODE,mesh1,&(data2[0]));
  Jali::StateVector<double> &addvec2 = mystate.add(myvec2);

  // Create Remap objects

  Portage::Remap_1stOrder<Portage::Jali_Mesh_Wrapper,Portage::Jali_State_Wrapper,Portage::Entity_kind> 
      remapper1(*mesh1,mystate,Portage::NODE,"nodevars1");
  Portage::Remap_1stOrder<Portage::Jali_Mesh_Wrapper,Portage::Jali_State_Wrapper,Portage::Entity_kind> 
      remapper2(*mesh1,mystate,Portage::NODE,"nodevars2");

  // Remap from source to target mesh

  double outvals1[4];  // field values on target mesh
  double outvals2[4];

  // Since we know the structure of the two meshes, we can
  // enumerate which source nodes contribute to a given target node
  // and what their intersection areas (weights) are

  std::vector<std::vector<int>> all_source_nodes =
      {{0,1,4,5},   {2,3,6,7},
       {8,9,12,13}, {10,11,14,15}};
  std::vector<std::vector<std::vector<double>>> all_weights =
      {{{1./36.}, {2./36.}, {2./36.}, {4./36.}},
       {{2./36.}, {1./36.}, {4./36.}, {2./36.}},
       {{2./36.}, {4./36.}, {1./36.}, {2./36.}},
       {{4./36.}, {2./36.}, {2./36.}, {1./36.}}};

  for (int i = 0; i < 4; ++i) {
    std::pair< std::vector<int> const &, 
               std::vector< std::vector<double> > const & > 
        nodes_and_weights(all_source_nodes[i],all_weights[i]);
    outvals1[i] = remapper1(nodes_and_weights);
    outvals2[i] = remapper2(nodes_and_weights);
  }

  // Make sure we retrieved the correct value for each cell on the target
  // For field 1, it is a constant
  // For field 2, it is a linear function
  // NOTE: Even though 1st order remapping algorithm does not in
  // general preserve a linear field, the special structure of the
  // source and target mesh ensures that the linear field is remapped
  // correctly in this test

  double stdval1 = data1[0];
  double stdvals2[4] = {(4./3.), 3., 3., (14./3.)};
  for (int i = 0; i < 4; ++i) {
    ASSERT_DOUBLE_EQ(stdval1, outvals1[i]);
    ASSERT_DOUBLE_EQ(stdvals2[i], outvals2[i]);
  }

}

  
