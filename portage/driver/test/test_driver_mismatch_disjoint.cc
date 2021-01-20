/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/


#include <iostream>
#include <memory>

#include "gtest/gtest.h"
#include "mpi.h"

#include "wonton/support/wonton.h"
#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"

#include "portage/support/portage.h"
#include "portage/driver/mmdriver.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliState.h"
#include "JaliStateVector.h"


// Remap between two completely disjoint meshes
// Make sure we don't barf or get stuck and only get 0 for the target

double TOL = 1e-12;

TEST(Test_Mismatch_Fixup, Test_Methods) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  if (Jali::framework_available(Jali::MSTK))
    mf.framework(Jali::MSTK);
  std::shared_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 1.0, 1.0, 3, 3);
  std::shared_ptr<Jali::Mesh> target_mesh = mf(2.0, 2.0, 4.0, 4.0, 2, 2);

  const int ncells_source =
      source_mesh->num_entities(Jali::Entity_kind::CELL,
                                Jali::Entity_type::PARALLEL_OWNED);
  const int ncells_target =
      target_mesh->num_entities(Jali::Entity_kind::CELL,
                                Jali::Entity_type::PARALLEL_OWNED);

  // Create state objects for source and target mesh

  std::shared_ptr<Jali::State> source_state(Jali::State::create(source_mesh));
  std::shared_ptr<Jali::State> target_state(Jali::State::create(target_mesh));

  // Add a state vector on the source cell with a constant value

  Jali::UniStateVector<double> sourcevec("cellvars", source_mesh, nullptr,
                                      Jali::Entity_kind::CELL,
                                      Jali::Entity_type::PARALLEL_OWNED);
  for (int c = 0; c < ncells_source; ++c) {
    sourcevec[c] = 1.0;
  }
  source_state->add(sourcevec);

  // Add zero value state vectors on target cells and nodes - once with
  // the old name and once with the new name

  Jali::UniStateVector<double, Jali::Mesh>& targetvec =
      target_state->add<double, Jali::Mesh, Jali::UniStateVector>("cellvars",
                                  target_mesh,
                                  Jali::Entity_kind::CELL,
                                  Jali::Entity_type::PARALLEL_OWNED);

  // Wrappers for interfacing with the underlying mesh data structures.

  Wonton::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Wonton::Jali_Mesh_Wrapper targetMeshWrapper(*target_mesh);

  // Wrappers for the source and target state managers

  Wonton::Jali_State_Wrapper sourceStateWrapper(*source_state);
  Wonton::Jali_State_Wrapper targetStateWrapper(*target_state);

  // Specify the fields to be remapped

  std::vector<std::string> source_var_names;
  source_var_names.push_back("cellvars");

  std::vector<std::string> target_var_names;
  target_var_names.push_back("cellvars");

  // Build the main driver object to test for default fixup options
  // Default fixup options are GLOBALLY_CONSERVATIVE, so even if the
  // meshes are disjoint the algorithm will attempt to distribute the
  // integral values of the source field on to the target

  Portage::MMDriver<Portage::SearchKDTree,
                  Portage::IntersectRnD,
                  Portage::Interpolate_1stOrder,
                  2,
                  Wonton::Jali_Mesh_Wrapper,
                  Wonton::Jali_State_Wrapper> remapper(sourceMeshWrapper,
                                                         sourceStateWrapper,
                                                         targetMeshWrapper,
                                                         targetStateWrapper);

  remapper.set_remap_var_names(source_var_names, target_var_names);

  // Execute remapper with default fixup

  remapper.run();

  // Verify that we got the fields we wanted (GLOBALLY CONSERVATIVE
  // means the total source value of 3.0 will be distributed equally
  // over the two target cells)
  for (int c = 0; c < ncells_target; c++) {
    ASSERT_NEAR(1.5, targetvec[c], TOL);
  }


  // Setup various partial fixup options and test

  remapper.set_partial_fixup_type(Portage::Partial_fixup_type::CONSTANT);
  remapper.set_empty_fixup_type(Portage::Empty_fixup_type::EXTRAPOLATE);

  remapper.run();

  for (int c = 0; c < ncells_target; c++) {
    ASSERT_NEAR(0.0, targetvec[c], TOL);
  }



  remapper.set_partial_fixup_type(Portage::Partial_fixup_type::LOCALLY_CONSERVATIVE);
  remapper.set_empty_fixup_type(Portage::Empty_fixup_type::EXTRAPOLATE);

  remapper.run();

  for (int c = 0; c < ncells_target; c++) {
    ASSERT_NEAR(0.0, targetvec[c], TOL);
  }



  remapper.set_partial_fixup_type(Portage::Partial_fixup_type::GLOBALLY_CONSERVATIVE);
  remapper.set_empty_fixup_type(Portage::Empty_fixup_type::EXTRAPOLATE);

  remapper.run();

  // With GLOBALLY_CONSERVATIVE, the integral of the source field which
  // is 3.0 will be spread out about the target cells
  for (int c = 0; c < ncells_target; c++) {
    ASSERT_NEAR(1.5, targetvec[c], TOL);
  }



  remapper.set_partial_fixup_type(Portage::Partial_fixup_type::CONSTANT);
  remapper.set_empty_fixup_type(Portage::Empty_fixup_type::LEAVE_EMPTY);

  remapper.run();

  for (int c = 0; c < ncells_target; c++) {
    ASSERT_NEAR(0.0, targetvec[c], TOL);
  }



  remapper.set_partial_fixup_type(Portage::Partial_fixup_type::LOCALLY_CONSERVATIVE);
  remapper.set_empty_fixup_type(Portage::Empty_fixup_type::LEAVE_EMPTY);

  remapper.run();

  for (int c = 0; c < ncells_target; c++) {
    ASSERT_NEAR(0.0, targetvec[c], TOL);
  }



  remapper.set_partial_fixup_type(Portage::Partial_fixup_type::GLOBALLY_CONSERVATIVE);
  remapper.set_empty_fixup_type(Portage::Empty_fixup_type::LEAVE_EMPTY);

  remapper.run();

  // Verify that we got the fields we wanted. Even though the integral
  // value of the source field is 3.0, we cannot put it anywhere
  // because all the target cells are empty and we requested that they
  // be left empty
  for (int c = 0; c < ncells_target; c++) {
    ASSERT_NEAR(0.0, targetvec[c], TOL);
  }

}
