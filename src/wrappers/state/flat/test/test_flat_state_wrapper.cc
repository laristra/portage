/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include "portage/wrappers/state/jali/jali_state_wrapper.h"
#include "portage/wrappers/state/flat/flat_state_wrapper.h"

#include <iostream>

#include "gtest/gtest.h"
#include "mpi.h"

#include "Mesh.hh"
#include "MeshFactory.hh"

#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"

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


TEST(Flat_State_Wrapper, DataTypes) {

  // Add multiple state vector types
  int const n_cells = 4;
  int const n_nodes = 9;
  float ftest[] = {1.1, 2.2, 3.3, 4.4};
  int itest[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  Vec2d vtest[n_cells];
  for (unsigned int i=0; i<n_cells; i++) vtest[i].set(1.0*i, 2.0*i);

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  std::shared_ptr<Jali::Mesh> inputMesh = mf(0.0, 0.0, 1.0, 1.0, 2, 2);
  Portage::Jali_Mesh_Wrapper inputMeshWrapper(*inputMesh);
  Jali::State state(inputMesh);
  Portage::Jali_State_Wrapper wrapper(state);

  state.add("f1", inputMesh, Jali::Entity_kind::CELL,
            Jali::Entity_type::ALL, ftest);
  state.add("i1", inputMesh, Jali::Entity_kind::NODE,
            Jali::Entity_type::ALL, itest);
  state.add("v1", inputMesh, Jali::Entity_kind::CELL,
            Jali::Entity_type::ALL, vtest);

  Portage::Flat_State_Wrapper<> flat_state(wrapper, {"f1", "i1", "v1"});

  // Check the data using Jali as well as by the Flat_State_Wrapper wrapper

  // Get raw float data using wrapper
  float* fdata;
  wrapper.get_data(Portage::CELL, "f1", &fdata);
  for (unsigned int i=0; i<n_cells; i++) ASSERT_EQ(fdata[i], ftest[i]);
  flat_state.get_data(Portage::CELL, "f1", &fdata);
  for (unsigned int i=0; i<n_cells; i++) ASSERT_EQ(fdata[i], ftest[i]);

  // Get raw int data using wrapper
  int* idata;
  wrapper.get_data(Portage::NODE, "i1", &idata);
  for (unsigned int i=0; i<n_nodes; i++) ASSERT_EQ(idata[i], itest[i]);

  // Get raw Vec2d data using wrapper
  Vec2d* vdata;
  wrapper.get_data(Portage::CELL, "v1", &vdata);
  for (unsigned int i=0; i<n_cells; i++)
  {
    ASSERT_EQ(vdata[i].x, vtest[i].x);
    ASSERT_EQ(vdata[i].y, vtest[i].y);
  }
}
