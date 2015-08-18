/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include "portage/state/state.h"
#include "portage/state/state_vector.h"

#include <iostream>

#include "gtest/gtest.h"
#include "mpi.h"

#include "Mesh.hh"
#include "MeshFactory.hh"

TEST(State, DefineState) {

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  Jali::Mesh *mesh1 = mf(0.0,0.0,1.0,1.0,2,2);

  ASSERT_TRUE(mesh1 != NULL);

  // Define two state vectors

  std::vector<double> data1 = {1.0,3.0,2.5,4.5}; 
  Portage::StateVector myvec1("cellvars",Jali::CELL,mesh1,&(data1[0]));

  std::vector<double> data2 = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0};
  Portage::StateVector myvec2("nodevars",Jali::NODE,mesh1,&(data2[0]));

  // Define another mesh and another statevector on that mesh

  Jali::Mesh *mesh2 = mf(0.0,0.0,1.0,1.0,3,3);
  
  std::vector<double> data3 = {1.0,3.0,2.5,4.5,1.0,2.0}; 
  Portage::StateVector myvec3("cellvars2",Jali::CELL,mesh2,&(data3[0]));
    

  // Create a state object and add the first two vectors to it

  Portage::State mystate(mesh1);

  int add_status;
  Portage::StateVector addvec = mystate.add(myvec1);
  ASSERT_EQ(addvec.size(),myvec1.size());
  for (int i = 0; i < addvec.size(); ++i)
    ASSERT_EQ(addvec[i],myvec1[i]);

  addvec = mystate.add("nodevars",Jali::NODE,&(data2[0]));
  ASSERT_EQ(addvec.size(),myvec2.size());
  for (int i = 0; i < addvec.size(); ++i)
    ASSERT_EQ(addvec[i],myvec2[i]);


  // Try to add the third vector (defined on a different mesh) to it - it 
  // should copy the data but be assigned to mesh1 instead of mesh2

  addvec = mystate.add(myvec3);
  ASSERT_NE(addvec.mesh(),myvec3.mesh());


  // Now retrieve the state vectors from the state object in different ways

  Portage::State::const_iterator itc;
  
  // Make sure we can retrieve the object by name

  itc = mystate.find("cellvars",Jali::CELL);
  ASSERT_NE(mystate.end(),itc);

  // Make sure the object we retrieved is identical to the one we put in

  Portage::StateVector myvec1_copy = *itc;
  
  ASSERT_EQ(myvec1.size(),myvec1_copy.size());
  for (int i = 0; i < myvec1.size(); ++i)
    ASSERT_EQ(myvec1[i],myvec1_copy[i]);
  
  // Make sure the code fails if we ask for the right name but wrong entity type

  itc = mystate.find("cellvars",Jali::FACE);
  ASSERT_EQ(mystate.end(),itc);


  // Try to retrieve a different vector by name

  itc = mystate.find("nodevars",Jali::NODE);
  ASSERT_NE(mystate.end(),itc);

  // Make sure the object we retrieved is identical to the one we put in

  Portage::StateVector myvec2_copy = *itc;

  ASSERT_EQ(myvec2.size(),myvec2_copy.size());
  for (int i = 0; i < myvec2.size(); ++i)
    ASSERT_EQ(myvec2[i],myvec2_copy[i]);
  

  // Retrieve state data through iterators and [] operators
 
  Portage::State::iterator it = mystate.begin();
  while (it != mystate.end()) {
    Portage::StateVector myvec4 = *it;

    ASSERT_TRUE((myvec4.name() == "cellvars" && myvec4.on_what() == Jali::CELL)
                ||
                (myvec4.name() == "cellvars2" && myvec4.on_what() == Jali::CELL)
                ||
                (myvec4.name() == "nodevars" && myvec4.on_what() == Jali::NODE));

    ++it;
  }
  

  myvec1_copy = mystate[0];
  ASSERT_TRUE(myvec1_copy.name() == "cellvars" && myvec1_copy.on_what() == Jali::CELL);

  myvec2_copy = mystate[1];
  ASSERT_TRUE(myvec2_copy.name() == "nodevars" && myvec2_copy.on_what() == Jali::NODE);

  // Print out state

  std::cout << mystate;

  
  // Finally make sure that we cannot find the vector that was on other mesh

  itc = mystate.find("cellvars2",Jali::CELL);
  ASSERT_EQ(mystate.end(),it);

}

  

  

