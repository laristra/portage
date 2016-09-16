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


TEST(Flat_State_Wrapper, DataTypes2D) {

  // Add multiple state vector types
  int const n_cells = 4;
  double dtest1[] = {1.1, 2.2, 3.3, 4.4};
  double dtest2[] = {1.2, 2.2, 2.3, 2.4};
  double dtest3[] = {1.3, 3.2, 3.3, 3.4};

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  std::shared_ptr<Jali::Mesh> inputMesh = mf(0.0, 0.0, 1.0, 1.0, 2, 2);
  Portage::Jali_Mesh_Wrapper inputMeshWrapper(*inputMesh);
  Jali::State jali_state(inputMesh);
  Portage::Jali_State_Wrapper jali_state_wrapper(jali_state);

  jali_state.add("d1", inputMesh, Jali::Entity_kind::CELL,
            Jali::Entity_type::ALL, dtest1);
  jali_state.add("d2", inputMesh, Jali::Entity_kind::CELL,
            Jali::Entity_type::ALL, dtest2);
  jali_state.add("d3", inputMesh, Jali::Entity_kind::CELL,
            Jali::Entity_type::ALL, dtest3);

  Portage::Flat_State_Wrapper<> flat_state(jali_state_wrapper,
          {"d1", "d2", "d3"});

  // Check the data using Jali as well as by the Flat_State_Wrapper wrapper

  // Get raw float data using wrapper
  double* ddata = nullptr;
  jali_state_wrapper.get_data(Portage::CELL, "d1", &ddata);
  for (unsigned int i=0; i<n_cells; i++) ASSERT_EQ(ddata[i], dtest1[i]);

  // Get raw float data using the flat mesh wrapper (d1)
  ddata = nullptr;
  flat_state.get_data(Portage::CELL, "d1", &ddata);
  for (unsigned int i=0; i<n_cells; i++) ASSERT_EQ(ddata[i], dtest1[i]);

  // Get raw float data using the flat mesh wrapper (d2)
  ddata = nullptr;
  flat_state.get_data(Portage::CELL, "d2", &ddata);
  for (unsigned int i=0; i<n_cells; i++) ASSERT_EQ(ddata[i], dtest2[i]);

  // Get raw float data using the flat mesh wrapper (d3)
  ddata = nullptr;
  flat_state.get_data(Portage::CELL, "d3", &ddata);
  for (unsigned int i=0; i<n_cells; i++) ASSERT_EQ(ddata[i], dtest3[i]);
}

TEST(Flat_State_Wrapper, DataTypes3D) {

  // Add multiple state vector types
  int const n_cells = 8;
  double dtest1[] = {1.1, 2.2, 3.3, 4.4, 5, 6, 7, 8};
  double dtest2[] = {1.2, 2.2, 2.3, 2.4, 5.1, 6.1, 7.1, 8.1};
  double dtest3[] = {1.3, 3.2, 3.3, 3.4, 5.2, 6.2, 7.2, 8.2};

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  std::shared_ptr<Jali::Mesh> inputMesh = mf(0.0, 0.0, 0.0,
                                             1.0, 1.0, 1.0,
                                             2, 2, 2);
  Portage::Jali_Mesh_Wrapper inputMeshWrapper(*inputMesh);
  Jali::State jali_state(inputMesh);
  Portage::Jali_State_Wrapper jali_state_wrapper(jali_state);

  jali_state.add("d1", inputMesh, Jali::Entity_kind::CELL,
            Jali::Entity_type::ALL, dtest1);
  jali_state.add("d2", inputMesh, Jali::Entity_kind::CELL,
            Jali::Entity_type::ALL, dtest2);
  jali_state.add("d3", inputMesh, Jali::Entity_kind::CELL,
            Jali::Entity_type::ALL, dtest3);

  Portage::Flat_State_Wrapper<> flat_state(jali_state_wrapper,
          {"d1", "d2", "d3"});

  // Check the data using Jali as well as by the Flat_State_Wrapper wrapper

  // Get raw float data using wrapper
  double* ddata = nullptr;
  jali_state_wrapper.get_data(Portage::CELL, "d1", &ddata);
  for (unsigned int i=0; i<n_cells; i++) ASSERT_EQ(ddata[i], dtest1[i]);

  // Get raw float data using the flat mesh wrapper (d1)
  ddata = nullptr;
  flat_state.get_data(Portage::CELL, "d1", &ddata);
  for (unsigned int i=0; i<n_cells; i++) ASSERT_EQ(ddata[i], dtest1[i]);

  // Get raw float data using the flat mesh wrapper (d2)
  ddata = nullptr;
  flat_state.get_data(Portage::CELL, "d2", &ddata);
  for (unsigned int i=0; i<n_cells; i++) ASSERT_EQ(ddata[i], dtest2[i]);

  // Get raw float data using the flat mesh wrapper (d3)
  ddata = nullptr;
  flat_state.get_data(Portage::CELL, "d3", &ddata);
  for (unsigned int i=0; i<n_cells; i++) ASSERT_EQ(ddata[i], dtest3[i]);
}
