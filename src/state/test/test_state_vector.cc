/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include "gtest/gtest.h"
#include "mpi.h"

#include "Mesh.hh"
#include "MeshFactory.hh"

#include "../state_vector.h"

TEST(STATE_VECTOR, DEFINE_VECTOR_ON_CELLS) {

  int mpi_initialized;
  MPI_Initialized(&mpi_initialized);

      
  int mpi_initialized_here = false;
  if (!mpi_initialized) {
    mpi_initialized_here = true;
    MPI_Init(NULL, NULL);
  }

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  Jali::Mesh *mesh = mf(0.0,0.0,1.0,1.0,2,2);

  ASSERT_TRUE(mesh != NULL);

  // purposely define data vector that is smaller than the number of cells 
  // so that we can test initialization

  std::vector<double> data1 = {1.0,3.0,2.5}; 
  NGC::Remap::StateVector myvec1("var1",Jali::CELL,mesh,&(data1[0]));

  int ncells = mesh->num_entities(Jali::CELL,Jali::ALL);
  ASSERT_EQ(ncells,myvec1.size());
  ASSERT_EQ(data1[0],myvec1[0]);
  ASSERT_EQ(data1[1],myvec1[1]);
  ASSERT_EQ(data1[2],myvec1[2]);
  ASSERT_EQ(0.0,myvec1[3]);
  
  if (mpi_initialized_here)
    MPI_Finalize();
}
