/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include "portage/wrappers/state/jali/jali_state.h"
#include "portage/wrappers/state/jali/jali_state_vector.h"

#include <iostream>

#include "gtest/gtest.h"
#include "mpi.h"

#include "Mesh.hh"
#include "MeshFactory.hh"

// Vector type for 2d doubles
struct Vec2d
{
  double x;
  double y;

  void set(double xvalue, double yvalue)
  {
    x = xvalue;  y = yvalue;
  }

  friend std::ostream &operator<<(std::ostream &output, const Vec2d &v)
  {
    output << "(" << v.x << ", " << v.y << ")";
    return output;
  }
};


TEST(Jali_State, DefineState) {

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  Jali::Mesh *mesh1 = mf(0.0,0.0,1.0,1.0,2,2);

  ASSERT_TRUE(mesh1 != NULL);

  // Define two state vectors

  std::vector<double> data1 = {1.0,3.0,2.5,4.5}; 
  Jali::StateVector<double> myvec1("cellvars",Jali::CELL,mesh1,&(data1[0]));

  std::vector<double> data2 = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0};
  Jali::StateVector<double> myvec2("nodevars",Jali::NODE,mesh1,&(data2[0]));

  // Define another mesh and another statevector on that mesh

  Jali::Mesh *mesh2 = mf(0.0,0.0,1.0,1.0,3,3);
  
  std::vector<double> data3 = {1.0,3.0,2.5,4.5,1.0,2.0}; 
  Jali::StateVector<double> myvec3("cellvars2",Jali::CELL,mesh2,&(data3[0]));
    

  // Create a state object and add the first two vectors to it

  Jali::State mystate(mesh1);

  int add_status;
  Jali::StateVector<double> &addvec1 = mystate.add(myvec1);
  ASSERT_EQ(addvec1.size(),myvec1.size());
  for (int i = 0; i < addvec1.size(); ++i)
    ASSERT_EQ(addvec1[i],myvec1[i]);

  Jali::StateVector<double> &addvec2 = mystate.add("nodevars",Jali::NODE,&(data2[0]));
  ASSERT_EQ(addvec2.size(),myvec2.size());
  for (int i = 0; i < addvec2.size(); ++i)
    ASSERT_EQ(addvec2[i],myvec2[i]);


  // Try to add the third vector (defined on a different mesh) to it - it 
  // should copy the data but be assigned to mesh1 instead of mesh2

  Jali::StateVector<double> &addvec3 = mystate.add(myvec3);
  ASSERT_NE(addvec3.mesh(),myvec3.mesh());


  // Now retrieve the state vectors from the state object in different ways

  Jali::State::const_iterator itc;
  
  // Make sure we can retrieve the object by name

  itc = mystate.find("cellvars",Jali::CELL);
  ASSERT_NE(mystate.end(),itc);

  // Make sure the object we retrieved is identical to the one we put in

  Jali::StateVector<double> myvec1_copy = *(std::static_pointer_cast<Jali::StateVector<double>>(*itc));
  
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

  Jali::StateVector<double> myvec2_copy = *(std::static_pointer_cast<Jali::StateVector<double>>(*itc));

  ASSERT_EQ(myvec2.size(),myvec2_copy.size());
  for (int i = 0; i < myvec2.size(); ++i)
    ASSERT_EQ(myvec2[i],myvec2_copy[i]);
  

  // Retrieve state data through iterators and [] operators
 
  Jali::State::iterator it = mystate.begin();
  while (it != mystate.end()) {
    Jali::StateVector<double> myvec4 = *(std::static_pointer_cast<Jali::StateVector<double>>(*it));

    ASSERT_TRUE((myvec4.name() == "cellvars" && myvec4.on_what() == Jali::CELL)
                ||
                (myvec4.name() == "cellvars2" && myvec4.on_what() == Jali::CELL)
                ||
                (myvec4.name() == "nodevars" && myvec4.on_what() == Jali::NODE));

    ++it;
  }
  

  myvec1_copy = *(std::static_pointer_cast<Jali::StateVector<double>>(mystate[0]));
  ASSERT_TRUE(myvec1_copy.name() == "cellvars" && myvec1_copy.on_what() == Jali::CELL);

  myvec2_copy = *(std::static_pointer_cast<Jali::StateVector<double>>(mystate[1]));
  ASSERT_TRUE(myvec2_copy.name() == "nodevars" && myvec2_copy.on_what() == Jali::NODE);

  // Print out state

  std::cout << mystate;

  
  // Finally make sure that we cannot find the vector that was on other mesh

  itc = mystate.find("cellvars2",Jali::CELL);
  ASSERT_EQ(mystate.end(),it);


  // Add state vectors of different data types

  int n_cells = 4;
  int n_nodes = 9;
  float ftest[] = {1.1, 2.2, 3.3, 4.4};
  int itest[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  Vec2d vtest[n_cells];
  for (unsigned int i=0; i<n_cells; i++) vtest[i].set(1.0*i, 2.0*i);

  Jali::MeshFactory factory(MPI_COMM_WORLD);

  Jali::Mesh* dataMesh = factory(0.0, 0.0, 1.0, 1.0, 2, 2);
  Jali::State dstate(dataMesh);

  dstate.add("f1", Jali::CELL, ftest);
  dstate.add("i1", Jali::NODE, itest);
  dstate.add("v1", Jali::CELL, vtest);

  // Iterate through all state vectors and count them

  int cnt = 0;
  for (Jali::State::iterator it = dstate.begin(); it != dstate.end(); it++) cnt++;
  ASSERT_EQ(cnt, 3);

  // Iterate through all cell state vectors and count them

  cnt = 0;
  for (Jali::State::permutation_type it = dstate.entity_begin(Jali::CELL); it != dstate.entity_end(Jali::CELL); it++) cnt++;
  ASSERT_EQ(cnt, 2);

  // Iterate through all node state vectors and count them

  cnt = 0;
  for (Jali::State::permutation_type it = dstate.entity_begin(Jali::NODE); it != dstate.entity_end(Jali::NODE); it++) cnt++;
  ASSERT_EQ(cnt, 1);

  // Iterate through all state vectors and get their type

  int testCnt = 0;
  for (Jali::State::iterator it = dstate.begin(); it != dstate.end(); it++)
  {
    if (typeid(float) == (*it)->get_type())      ASSERT_EQ(testCnt, 0);
    else if (typeid(int) == (*it)->get_type())   ASSERT_EQ(testCnt, 1);
    else if (typeid(Vec2d) == (*it)->get_type()) ASSERT_EQ(testCnt, 2);
    else                                         ASSERT_EQ(0, 1);        // This else should never be reached in this test
    testCnt++;
  }

}

  

  

