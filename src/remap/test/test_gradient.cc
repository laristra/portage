/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include "portage/support/portage.h"
#include "portage/remap/gradient.h"
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


TEST(Gradient, Fields_Cell_Ctr) {

  // Make a 4x4 mesh 

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);

  Jali::Mesh *mesh1 = mf(0.0,0.0,1.0,1.0,4,4);
  ASSERT_TRUE(mesh1 != NULL);

  // Create a state object and add the first two vectors to it

  Jali::State mystate(mesh1);

  // Define three state vectors, one with constant value and the other
  // with a linear function that is x+2y

  int nc1 = mesh1->num_entities(Jali::CELL,Jali::OWNED);
  std::vector<double> data1(nc1);
  for (int c = 0; c < nc1; ++c)
    data1[c] = 1.25;
  Jali::StateVector<double> myvec1("cellvars1",Jali::CELL,mesh1,&(data1[0]));
  Jali::StateVector<double> &addvec1 = mystate.add(myvec1);

  
  std::vector<double> data2(nc1);
  for (int c = 0; c < nc1; c++) {
    JaliGeometry::Point ccen = mesh1->cell_centroid(c);
    data2[c] = ccen[0]+2*ccen[1];
  }

  Jali::StateVector<double> myvec2("cellvars2",Jali::CELL,mesh1,&(data2[0]));
  Jali::StateVector<double> &addvec2 = mystate.add(myvec2);


  // Create Remap objects

  Portage::Limited_Gradient<Portage::Jali_Mesh_Wrapper,Portage::Jali_State_Wrapper,Portage::Entity_kind> 
      gradcalc1(*mesh1,mystate,Portage::CELL,"cellvars1",Portage::NOLIMITER);
  Portage::Limited_Gradient<Portage::Jali_Mesh_Wrapper,Portage::Jali_State_Wrapper,Portage::Entity_kind> 
      gradcalc2(*mesh1,mystate,Portage::CELL,"cellvars2",Portage::NOLIMITER);


  // Compute the gradient for each of these fields

  std::vector<double> grad1(2);
  std::vector<double> grad2(2);

  // Verify the gradient values
  // For field 1 (constant), it is is 0,0
  // For field 2 (x+2y), it is (1,2)

  for (int c = 0; c < nc1; ++c) {
    grad1 = gradcalc1(c);
    ASSERT_NEAR(0.0,grad1[0],1.0e-10);
    ASSERT_NEAR(0.0,grad1[1],1.0e-10);

    grad2 = gradcalc2(c);
    ASSERT_NEAR(1.0,grad2[0],1.0e-10);
    ASSERT_NEAR(2.0,grad2[1],1.0e-10);
  }

}


TEST(Gradient, Fields_Node_Ctr) {

  // Make a 3x3 mesh 

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);

  Jali::Mesh *mesh1 = mf(0.0,0.0,1.0,1.0,3,3);
  ASSERT_TRUE(mesh1 != NULL);

  // Create a state object and add the first two vectors to it

  Jali::State mystate(mesh1);

  // Define three state vectors, one with constant value, the other
  // with a linear function

  int nn1 = mesh1->num_entities(Jali::NODE,Jali::OWNED);

  std::vector<double> data1(nn1);
  for (int n = 0; n < nn1; ++n) data1[n] = 1.5;

  Jali::StateVector<double> myvec1("nodevars1",Jali::NODE,mesh1,&(data1[0]));
  Jali::StateVector<double> &addvec1 = mystate.add(myvec1);

  std::vector<double> data2(nn1);
  for (int n = 0; n < nn1; ++n) {
    JaliGeometry::Point nodexy;
    mesh1->node_get_coordinates(n,&nodexy);
    data2[n] = 3*nodexy[0]+nodexy[1];
  }
  Jali::StateVector<double> myvec2("nodevars2",Jali::NODE,mesh1,&(data2[0]));
  Jali::StateVector<double> &addvec2 = mystate.add(myvec2);

  // Create Gradient calculater objects

  Portage::Limited_Gradient<Portage::Jali_Mesh_Wrapper,Portage::Jali_State_Wrapper,Portage::Entity_kind> 
      gradcalc1(*mesh1,mystate,Portage::NODE,"nodevars1",Portage::NOLIMITER);
  Portage::Limited_Gradient<Portage::Jali_Mesh_Wrapper,Portage::Jali_State_Wrapper,Portage::Entity_kind> 
      gradcalc2(*mesh1,mystate,Portage::NODE,"nodevars2",Portage::NOLIMITER);


  // Make sure we retrieved the correct gradient value for each node
  // For field 1, it is a constant
  // For field 2, it is a linear function

  std::vector<double> grad1(2);
  std::vector<double> grad2(2);

  for (int c = 0; c < nn1; ++c) {
    grad1 = gradcalc1(c);
    ASSERT_NEAR(0.0,grad1[0],1.0e-10);
    ASSERT_NEAR(0.0,grad1[1],1.0e-10);

    grad2 = gradcalc2(c);
    ASSERT_NEAR(3.0,grad2[0],1.0e-10);
    ASSERT_NEAR(1.0,grad2[1],1.0e-10);
  }

}

  
