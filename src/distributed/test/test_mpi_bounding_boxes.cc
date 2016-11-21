/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include "portage/support/Point.h"
#include "portage/distributed/mpi_bounding_boxes.h"
#include "portage/wrappers/state/jali/jali_state_wrapper.h"
#include "portage/wrappers/state/flat/flat_state_wrapper.h"
#include "portage/wrappers/mesh/flat/flat_mesh_wrapper.h"

#include <iostream>

#include "gtest/gtest.h"
#include "mpi.h"

#include "Mesh.hh"
#include "MeshFactory.hh"

#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"


TEST(MPI_Bounding_Boxes, SimpleTest) {

  // Add multiple state vector types
  double dtest1[] = {1.1, 2.2, 3.3, 4.4, 5., 6., 7., 8.};
  double dtest2[] = {1.2, 2.3, 3.4, 4.4, 5., 6., 9., 7.};

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  // Source mesh
  std::shared_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
        2, 2, 2);
  Portage::Jali_Mesh_Wrapper inputMeshWrapper(*source_mesh);
  Jali::State state(source_mesh);
  Portage::Jali_State_Wrapper wrapper(state);

  state.add("d1", source_mesh, Jali::Entity_kind::CELL,
            Jali::Entity_type::ALL, dtest1);
  state.add("d2", source_mesh, Jali::Entity_kind::CELL,
            Jali::Entity_type::ALL, dtest2);

  Portage::Flat_Mesh_Wrapper<> source_mesh_flat;
  source_mesh_flat.initialize(8, inputMeshWrapper);
  Portage::Flat_State_Wrapper<> source_state_flat;
  source_state_flat.initialize(wrapper, {"d1", "d2"});

  // Target mesh
  std::shared_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
        3, 3, 3);
  Portage::Jali_Mesh_Wrapper target_mesh_(*target_mesh);
  Jali::State target_state(target_mesh);
  Portage::Jali_State_Wrapper target_state_(target_state);

  // Use a bounding box distributor to send the source cells to the target
  // paritions where they are needed
  Portage::MPI_Bounding_Boxes distributor;
  distributor.distribute(source_mesh_flat, source_state_flat, target_mesh_,
                         target_state_);
}
