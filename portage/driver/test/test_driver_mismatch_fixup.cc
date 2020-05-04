/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/


#include <iostream>
#include <memory>

#include "gtest/gtest.h"
#include "mpi.h"

#include "portage/support/portage.h"

#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"
#include "portage/driver/mmdriver.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliState.h"
#include "JaliStateVector.h"


double TOL = 1e-12;

TEST(Test_Mismatch_Fixup, Test_Methods) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  if (Jali::framework_available(Jali::MSTK))
    mf.framework(Jali::MSTK);
  std::shared_ptr<Jali::Mesh> source_mesh = mf(-0.8, 0.0, 0.4, 1.0, 1, 1);
  std::shared_ptr<Jali::Mesh> target_mesh = mf( 0.0, 0.0, 2.0, 1.0, 2, 1);

  // exact results for Partial_fixup_type:
  //	CONSTANT (C)
  //    LOCALLY_CONSERVATIVE (L)
  //    SHIFTED_CONSERVATIVE (S)
  // and Empty_fixup_type:
  // 	EXTRAPOLATE (E)
  //    LEAVE_EMPTY (L)

  double exact_C_E[2] = {1.0, 1.0};
  double exact_L_E[2] = {0.4, 0.4};
  double exact_S_E[2] = {0.6, 0.6};
  double exact_C_L[2] = {1.0, 0.0};
  double exact_L_L[2] = {0.4, 0.0};
  double exact_S_L[2] = {1.2, 0.0};

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

  Portage::MMDriver<Portage::SearchKDTree,
                  Portage::IntersectR2D,
                  Portage::Interpolate_1stOrder,
                  2,
                  Wonton::Jali_Mesh_Wrapper,
                  Wonton::Jali_State_Wrapper> remapper(sourceMeshWrapper,
                                                         sourceStateWrapper,
                                                         targetMeshWrapper,
                                                         targetStateWrapper);

  remapper.set_remap_var_names(source_var_names, target_var_names);

  // Execute remapper (No arguments implies serial execution)

  remapper.run();

  // Verify that we got the fields we wanted
  for (int c = 0; c < ncells_target; c++) {
    ASSERT_NEAR(exact_S_E[c], targetvec[c], TOL);
  }



  // Set fixup types

  remapper.set_partial_fixup_type(Portage::Partial_fixup_type::CONSTANT);
  remapper.set_empty_fixup_type(Portage::Empty_fixup_type::EXTRAPOLATE);

  // Execute remapper (No arguments implies serial execution)

  remapper.run();

  // Verify that we got the fields we wanted
  for (int c = 0; c < ncells_target; c++) {
    ASSERT_NEAR(exact_C_E[c], targetvec[c], TOL);
  }



  // Set fixup types

  remapper.set_partial_fixup_type(Portage::Partial_fixup_type::LOCALLY_CONSERVATIVE);
  remapper.set_empty_fixup_type(Portage::Empty_fixup_type::EXTRAPOLATE);

  // Execute remapper (No arguments implies serial execution)

  remapper.run();

  // Verify that we got the fields we wanted
  for (int c = 0; c < ncells_target; c++) {
    ASSERT_NEAR(exact_L_E[c], targetvec[c], TOL);
  }



  // Set fixup types

  remapper.set_partial_fixup_type(Portage::Partial_fixup_type::SHIFTED_CONSERVATIVE);
  remapper.set_empty_fixup_type(Portage::Empty_fixup_type::EXTRAPOLATE);

  // Execute remapper (No arguments implies serial execution)

  remapper.run();

  // Verify that we got the fields we wanted
  for (int c = 0; c < ncells_target; c++) {
    ASSERT_NEAR(exact_S_E[c], targetvec[c], TOL);
  }



  // Set fixup types

  remapper.set_partial_fixup_type(Portage::Partial_fixup_type::CONSTANT);
  remapper.set_empty_fixup_type(Portage::Empty_fixup_type::LEAVE_EMPTY);

  // Execute remapper (No arguments implies serial execution)

  remapper.run();

  // Verify that we got the fields we wanted
  for (int c = 0; c < ncells_target; c++) {
    ASSERT_NEAR(exact_C_L[c], targetvec[c], TOL);
  }



  // Set fixup types

  remapper.set_partial_fixup_type(Portage::Partial_fixup_type::LOCALLY_CONSERVATIVE);
  remapper.set_empty_fixup_type(Portage::Empty_fixup_type::LEAVE_EMPTY);

  // Execute remapper (No arguments implies serial execution)

  remapper.run();

  // Verify that we got the fields we wanted
  for (int c = 0; c < ncells_target; c++) {
    ASSERT_NEAR(exact_L_L[c], targetvec[c], TOL);
  }



  // Set fixup types

  remapper.set_partial_fixup_type(Portage::Partial_fixup_type::SHIFTED_CONSERVATIVE);
  remapper.set_empty_fixup_type(Portage::Empty_fixup_type::LEAVE_EMPTY);

  // Execute remapper (No arguments implies serial execution)

  remapper.run();

  // Verify that we got the fields we wanted
  for (int c = 0; c < ncells_target; c++) {
    ASSERT_NEAR(exact_S_L[c], targetvec[c], TOL);
  }

}
