/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include "portage/remap/remap_1st_order.h"
#include "portage/state/state.h"
#include "portage/state/state_vector.h"

#include <iostream>

#include "gtest/gtest.h"
#include "mpi.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "FrameworkTraits.hh"
#include "MeshFramework.hh"


TEST(Remap_1st_Order,Constant_Field_Test1) {

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

  Portage::State mystate(mesh1);

  // Define state vector, "density", on source mesh, with the same value
  // on all the cells.

  std::vector<double> data1 = {1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,
                               1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25}; 
  Portage::StateVector myvec1("cellvars",Jali::CELL,mesh1,&(data1[0]));
  Portage::StateVector &addvec1 = mystate.add(myvec1);

  // Create a Remap object

  Portage::Remap_1stOrder remapper(*mesh1,mystate);

  // Remap from source to target mesh

  double outvals[4] = {0.0,0.0,0.0,0.0};  // field values on target mesh

  int all_source_cells[4][4] = {{0,1,4,5},{2,3,6,7},
                                {8,9,12,13},{10,11,14,15}};
  double all_weights[4][4] = {{0.25,0.25,0.25,0.25},{0.25,0.25,0.25,0.25},
                              {0.25,0.25,0.25,0.25},{0.25,0.25,0.25,0.25}};

  for (int i = 0; i < 4; ++i) {

    // Since we are know the structure of the two meshes, we can
    // enumerate which source cells intersect a given target cell and
    // what their intersection areas (weights) are

    Jali::Entity_ID_List source_cells(4);
    std::vector<double> weights(4);
    for (int j = 0; j < 4; ++j) {
      source_cells[j] = all_source_cells[i][j];
      weights[j] = all_weights[i][j];
    }

    outvals[i] = remapper(i,"cellvars",source_cells,weights);
  }

  // Make sure we retrieved a constant value for each cell on the target

  for (int i = 0; i < 4; ++i) {
    ASSERT_EQ(data1[0],outvals[i]);
  }

}

  

  

