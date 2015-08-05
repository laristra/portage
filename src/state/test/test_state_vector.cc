/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include "../state_vector.h"

#include "gtest/gtest.h"
#include "mpi.h"

#include "Mesh.hh"
#include "MeshFactory.hh"

TEST(StateVector, DefineVectorOnCells) {

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  Jali::Mesh *mesh = mf(0.0,0.0,1.0,1.0,2,2);

  ASSERT_TRUE(mesh != NULL);

  std::vector<double> data1 = {1.0,3.0,2.5,4.5}; 
  Portage::StateVector myvec1("var1",Jali::CELL,mesh,&(data1[0]));

  int ncells = mesh->num_entities(Jali::CELL,Jali::ALL);
  ASSERT_EQ(ncells,myvec1.size());
  ASSERT_EQ(data1[0],myvec1[0]);
  ASSERT_EQ(data1[1],myvec1[1]);
  ASSERT_EQ(data1[2],myvec1[2]);
  ASSERT_EQ(data1[3],myvec1[3]);

}


TEST(StateVector, DefineVectorOnNodes) {

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  Jali::Mesh *mesh = mf(0.0,0.0,1.0,1.0,2,2);

  ASSERT_TRUE(mesh != NULL);

  std::vector<double> data1 = {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0}; 
  Portage::StateVector myvec1("var1",Jali::NODE,mesh,&(data1[0]));

  int nnodes = mesh->num_entities(Jali::NODE,Jali::ALL);
  ASSERT_EQ(nnodes,myvec1.size());
  for (int i = 0; i < nnodes; ++i)
    ASSERT_EQ(data1[i],myvec1[i]);
  
}


