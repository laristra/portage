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


TEST(Flat_State_Wrapper, DataTypes) {

  // Add multiple state vector types
  int const n_cells = 4;
  double dtest[] = {1.1, 2.2, 3.3, 4.4};

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  std::shared_ptr<Jali::Mesh> inputMesh = mf(0.0, 0.0, 1.0, 1.0, 2, 2);
  Portage::Jali_Mesh_Wrapper inputMeshWrapper(*inputMesh);
  Jali::State state(inputMesh);
  Portage::Jali_State_Wrapper wrapper(state);

  state.add("d1", inputMesh, Jali::Entity_kind::CELL,
            Jali::Entity_type::ALL, dtest);

  Portage::Flat_State_Wrapper<> flat_state(wrapper, {"d1"});

  // Check the data using Jali as well as by the Flat_State_Wrapper wrapper

  // Get raw float data using wrapper
  double* ddata = nullptr;
  wrapper.get_data(Portage::CELL, "d1", &ddata);
  for (unsigned int i=0; i<n_cells; i++) ASSERT_EQ(ddata[i], dtest[i]);

  // Get raw float data using the flat mesh wrapper
  ddata = nullptr;
  flat_state.get_data(Portage::CELL, "d1", &ddata);
  for (unsigned int i=0; i<n_cells; i++) ASSERT_EQ(ddata[i], dtest[i]);
}
