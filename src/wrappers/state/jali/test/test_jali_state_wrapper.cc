/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

//#include "portage/wrappers/state/jali/jali_state.h"
//#include "portage/wrappers/state/jali/jali_state_vector.h"
#include "portage/wrappers/state/jali/jali_state_wrapper.h"

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


TEST(Jali_State_Wrapper, DataTypes) {

  // Tests with multiple state vector types

  int n_cells = 4;
  int n_nodes = 9;
  float ftest[] = {1.1, 2.2, 3.3, 4.4};
  int itest[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  Vec2d vtest[n_cells];
  for (unsigned int i=0; i<n_cells; i++) vtest[i].set(1.0*i, 2.0*i);

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  Jali::Mesh* inputMesh = mf(0.0, 0.0, 1.0, 1.0, 2, 2);
  Jali::State source_state(inputMesh);
  Jali_State_Wrapper source_wrapper(source_state);

  source_state.add("f1", Jali::CELL, ftest);
  source_state.add("i1", Jali::NODE, itest);
  source_state.add("v1", Jali::CELL, vtest);

  float* fdata;
  source_wrapper.get_data(Jali::CELL, "f1", &fdata);
  std::cout << "Output f1 raw data using wrapper"<< std::endl;
  for (unsigned int i=0; i<n_cells; i++) std::cout << fdata[i] << "\t";
  std::cout << std::endl << std::endl;

  int* idata;
  source_wrapper.get_data(Jali::NODE, "i1", &idata);
  std::cout << "Output i1 raw data using wrapper" << std::endl;
  for (unsigned int i=0; i<n_nodes; i++) std::cout << idata[i] << "\t";
  std::cout << std::endl << std::endl;

  Vec2d* vdata;
  source_wrapper.get_data(Jali::CELL, "v1", &vdata);
  std::cout << "Output v1 raw data using wrapper" << std::endl;
  for (unsigned int i=0; i<n_cells; i++) std::cout << "(" << vdata[i].x << ", " << vdata[i].y << ")\t";
  std::cout << std::endl << std::endl;

  std::cout << "Iterate through all state vectors and print them" << std::endl;
  int cnt = 0;
  for (Jali::State::iterator it = source_state.begin(); it != source_state.end(); it++)
  {
    (*it)->print(std::cout);
    cnt++;
  }
  std::cout << std::endl;
  ASSERT_EQ(cnt, 3);

  std::cout << "Iterate through all cell state vectors and print them" << std::endl;
  cnt = 0;
  for (Jali::State::permutation_type it = source_state.entity_begin(Jali::CELL); it != source_state.entity_end(Jali::CELL); it++)
  {
    (*it)->print(std::cout);
    cnt++;
  }
  std::cout << std::endl;
  ASSERT_EQ(cnt, 2);

  std::cout << "Iterate through all node state vectors and print them" << std::endl;
  cnt = 0;
  for (Jali::State::permutation_type it = source_state.entity_begin(Jali::NODE); it != source_state.entity_end(Jali::NODE); it++)
  {
    (*it)->print(std::cout);
    cnt++;
  }
  std::cout << std::endl;
  ASSERT_EQ(cnt, 1);

}

  

  

