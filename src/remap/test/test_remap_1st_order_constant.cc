/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include "portage/remap/remap_1st_order.h"
#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"
#include "portage/wrappers/state/jali/jali_state.h"
#include "portage/wrappers/state/jali/jali_state_vector.h"
#include "portage/wrappers/state/jali/jali_state_wrapper.h"

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

  Jali::State mystate(mesh1);

  // Define state vector, "density", on source mesh, with the same value
  // on all the cells.

  std::string varname("cellvars");
  std::vector<double> data1 = {1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,
                               1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25}; 
  Jali::StateVector myvec1("cellvars",Jali::CELL,mesh1,&(data1[0]));
  Jali::StateVector &addvec1 = mystate.add(myvec1);

  // Create a Remap object

  Portage::Remap_1stOrder<Jali_Mesh_Wrapper,Jali_State_Wrapper,Jali::Entity_kind> 
      remapper(*mesh1,mystate,Jali::CELL,"cellvars");

  // Remap from source to target mesh

  double outvals[4] = {0.0,0.0,0.0,0.0};  // field values on target mesh

  int all_source_cells[4][4] = {{0,1,4,5},{2,3,6,7},
                                {8,9,12,13},{10,11,14,15}};
  double all_weights[4][4] = {{0.25,0.25,0.25,0.25},{0.25,0.25,0.25,0.25},
                              {0.25,0.25,0.25,0.25},{0.25,0.25,0.25,0.25}};

  for (int i = 0; i < 4; ++i) {

    // Since we know the structure of the two meshes, we can
    // enumerate which source cells intersect a given target cell and
    // what their intersection areas (weights) are

    Jali::Entity_ID_List source_cells(4);
    std::vector< std::vector<double> > weight_vectors(4);
    for (int j = 0; j < 4; ++j) {
      std::vector<double> wtvec;
      wtvec.emplace_back(all_weights[i][j]);
      weight_vectors[j] = wtvec;
      source_cells[j] = all_source_cells[i][j];
    }

    std::pair< std::vector<int> const &, 
               std::vector< std::vector<double> > const & > 
        cells_and_weights(source_cells,weight_vectors);

    outvals[i] = remapper(cells_and_weights);
  }

  // Make sure we retrieved a constant value for each cell on the target

  for (int i = 0; i < 4; ++i) {
    ASSERT_EQ(data1[0],outvals[i]);
  }

}

  

  
