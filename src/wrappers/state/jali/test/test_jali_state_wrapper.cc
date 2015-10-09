/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include "portage/wrappers/state/jali/jali_state_wrapper.h"
#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"

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
  Jali_Mesh_Wrapper inputMeshWrapper(*inputMesh);
  Jali::State state(inputMesh);
  Jali_State_Wrapper wrapper(state);

  state.add("f1", Jali::CELL, ftest);
  state.add("i1", Jali::NODE, itest);
  state.add("v1", Jali::CELL, vtest);

  float* fdata;
  wrapper.get_data(Jali::CELL, "f1", &fdata);
  std::cout << "Output f1 raw data using wrapper"<< std::endl;
  for (unsigned int i=0; i<n_cells; i++) std::cout << fdata[i] << "\t";
  ASSERT_EQ(fdata[2], 3.3f);
  std::cout << std::endl << std::endl;

  int* idata;
  wrapper.get_data(Jali::NODE, "i1", &idata);
  std::cout << "Output i1 raw data using wrapper" << std::endl;
  for (unsigned int i=0; i<n_nodes; i++) std::cout << idata[i] << "\t";
  ASSERT_EQ(idata[6], 7);
  std::cout << std::endl << std::endl;

  Vec2d* vdata;
  wrapper.get_data(Jali::CELL, "v1", &vdata);
  std::cout << "Output v1 raw data using wrapper" << std::endl;
  for (unsigned int i=0; i<n_cells; i++) std::cout << "(" << vdata[i].x << ", " << vdata[i].y << ")\t";
  ASSERT_EQ(vdata[0].x, 0.0f);
  std::cout << std::endl << std::endl;

  std::cout << "Iterate through all state vectors and print them" << std::endl;
  int cnt = 0;
  for (Jali::State::iterator it = state.begin(); it != state.end(); it++)
  {
    (*it)->print(std::cout);
    cnt++;
  }
  std::cout << std::endl;
  ASSERT_EQ(cnt, 3);

  std::cout << "Iterate through all cell state vectors and print them" << std::endl;
  cnt = 0;
  for (Jali::State::permutation_type it = state.entity_begin(Jali::CELL); it != state.entity_end(Jali::CELL); it++)
  {
    (*it)->print(std::cout);
    cnt++;
  }
  std::cout << std::endl;
  ASSERT_EQ(cnt, 2);

  std::cout << "Iterate through all node state vectors and print them" << std::endl;
  cnt = 0;
  for (Jali::State::permutation_type it = state.entity_begin(Jali::NODE); it != state.entity_end(Jali::NODE); it++)
  {
    (*it)->print(std::cout);
    cnt++;
  }
  std::cout << std::endl;
  ASSERT_EQ(cnt, 1);

  std::cout << "Iterate through all state vectors and get their type" << std::endl;
  int testCnt = 0;
  for (Jali::State::iterator it = state.begin(); it != state.end(); it++)
  {
    if (typeid(float) == (*it)->get_type())
    {
      std::cout << (*it)->name() << " is a float" << std::endl;
      ASSERT_EQ(testCnt, 0);
    }
    else if (typeid(int) == (*it)->get_type())
    {
      std::cout << (*it)->name() << " is an int" << std::endl;
      ASSERT_EQ(testCnt, 1);
    }
    else if (typeid(Vec2d) == (*it)->get_type())
    {
      std::cout << (*it)->name() << " is a Vec2d" << std::endl;
      ASSERT_EQ(testCnt, 2);
    }
    else
    {
      std::cout << "Unknown data type" << std::endl;
      ASSERT_EQ(0, 1);
    }
    testCnt++;
  }
  std::cout << std::endl;

  std::cout << "Iterate through a vector of names" << std::endl;
  std::vector<std::string> fields;
  fields.push_back("v1");  fields.push_back("f1");  fields.push_back("i1");
  for (auto it = fields.begin(); it != fields.end(); it++)
  {
    std::cout << (*it) << std::endl;
    int on_what = wrapper.get_entity(*it);
    if (typeid(float) == wrapper.get_type(*it))
    {
      float* fdata;
      wrapper.get_data(on_what, *it, &fdata);
      std::cout << "Output " << (*it) << " raw data of float type using wrapper"<< std::endl;
      for (unsigned int i=0; i<inputMeshWrapper.num_entities(on_what); i++) std::cout << fdata[i] << "\t";
      ASSERT_EQ(fdata[0], 1.1f);
    }
    else if (typeid(int) == wrapper.get_type(*it))
    {
      int* idata;
      wrapper.get_data(on_what, *it, &idata);
      std::cout << "Output " << (*it) << " raw data of int type using wrapper"<< std::endl;
      for (unsigned int i=0; i<inputMeshWrapper.num_entities(on_what); i++) std::cout << idata[i] << "\t";
      ASSERT_EQ(idata[1], 2);
    }
    else if (typeid(Vec2d) == wrapper.get_type(*it))
    {
      Vec2d* vdata;
      wrapper.get_data(on_what, *it, &vdata);
      std::cout << "Output " << (*it) << " raw data of Vec2d type using wrapper"<< std::endl;
      for (unsigned int i=0; i<inputMeshWrapper.num_entities(on_what); i++) std::cout << vdata[i] << "\t";
      ASSERT_EQ(vdata[1].y, 2.0f);
    }
    else
    {
      std::cout << "Unknown data type" << std::endl;
      ASSERT_EQ(0, 1);
    }
    std::cout << std::endl << std::endl;
  }

  std::cout << "Iterate through all fields using the wrapper" << std::endl;
  for (auto it = wrapper.names_begin(); it != wrapper.names_end(); it++)
  {
    std::cout << (*it) << std::endl;
    int on_what = wrapper.get_entity(*it);
    if (typeid(float) == wrapper.get_type(*it))
    {
      float* fdata;
      wrapper.get_data(on_what, *it, &fdata);
      std::cout << "Output " << (*it) << " raw data of float type using wrapper"<< std::endl;
      for (unsigned int i=0; i<inputMeshWrapper.num_entities(on_what); i++) std::cout << fdata[i] << "\t";
      ASSERT_EQ(fdata[1], 2.2f);
    }
    else if (typeid(int) == wrapper.get_type(*it))
    {
      int* idata;
      wrapper.get_data(on_what, *it, &idata);
      std::cout << "Output " << (*it) << " raw data of int type using wrapper"<< std::endl;
      for (unsigned int i=0; i<inputMeshWrapper.num_entities(on_what); i++) std::cout << idata[i] << "\t";
      ASSERT_EQ(idata[0], 1);
    }
    else if (typeid(Vec2d) == wrapper.get_type(*it))
    {
      Vec2d* vdata;
      wrapper.get_data(on_what, *it, &vdata);
      std::cout << "Output " << (*it) << " raw data of Vec2d type using wrapper"<< std::endl;
      for (unsigned int i=0; i<inputMeshWrapper.num_entities(on_what); i++) std::cout << vdata[i] << "\t";
      ASSERT_EQ(vdata[2].x, 2.0f);
    }
    else
    {
      std::cout << "Unknown data type" << std::endl;
      ASSERT_EQ(0, 1);
    }
    std::cout << std::endl << std::endl;
  }

  std::cout << "Iterate through fields on cells only using the wrapper" << std::endl;
  for (auto it = wrapper.names_entity_begin(Jali::CELL); it != wrapper.names_entity_end(Jali::CELL); it++)
  {
    std::cout << (*it) << std::endl;
    int on_what = wrapper.get_entity(*it);
    if (typeid(float) == wrapper.get_type(*it))
    {
      float* fdata;
      wrapper.get_data(on_what, *it, &fdata);
      std::cout << "Output " << (*it) << " raw data of float type using wrapper"<< std::endl;
      for (unsigned int i=0; i<inputMeshWrapper.num_entities(on_what); i++) std::cout << fdata[i] << "\t";
      ASSERT_EQ(fdata[3], 4.4f);
    }
    else if (typeid(int) == wrapper.get_type(*it))
    {
      int* idata;
      wrapper.get_data(on_what, *it, &idata);
      std::cout << "Output " << (*it) << " raw data of int type using wrapper"<< std::endl;
      for (unsigned int i=0; i<inputMeshWrapper.num_entities(on_what); i++) std::cout << idata[i] << "\t";
      ASSERT_EQ(0, 1);
    }
    else if (typeid(Vec2d) == wrapper.get_type(*it))
    {
      Vec2d* vdata;
      wrapper.get_data(on_what, *it, &vdata);
      std::cout << "Output " << (*it) << " raw data of Vec2d type using wrapper"<< std::endl;
      for (unsigned int i=0; i<inputMeshWrapper.num_entities(on_what); i++) std::cout << vdata[i] << "\t";
      ASSERT_EQ(vdata[3].y, 6.0f);
    }
    else
    {
      std::cout << "Unknown data type" << std::endl;
      ASSERT_EQ(0, 1);
    }
    std::cout << std::endl << std::endl;
  }

}

  

  

