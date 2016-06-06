/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include <iostream>
#include <memory>

#include "gtest/gtest.h"
#include "mpi.h"

#include "portage/support/portage.h"
#include "portage/support/Point.h"
#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"
#include "portage/wrappers/state/jali/jali_state_wrapper.h"
#include "portage/driver/driver.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "FrameworkTraits.hh"
#include "MeshFramework.hh"
#include "JaliState.h"
#include "JaliStateVector.h"

TEST(Test_MultiVar_Remap, Test1) {

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);
  mf.included_entities({Jali::Entity_kind::CORNER, Jali::Entity_kind::WEDGE});

  std::shared_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4);
  std::shared_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 1.0, 1.0, 5, 5);

  const int ncells_target =
      target_mesh->num_entities(Jali::Entity_kind::CELL,
                                Jali::Parallel_type::OWNED);
  const int nnodes_target =
      target_mesh->num_entities(Jali::Entity_kind::NODE,
                                Jali::Parallel_type::OWNED);

  // Create state objects for source and target mesh

  Jali::State source_state(source_mesh);
  Jali::State target_state(target_mesh);

  // Add a constant value state vector on source cells

  double Constant1 = 1.25;
  Jali::StateVector<double> myvec1("srccellvars1", source_mesh,
                                   Jali::Entity_kind::CELL,
                                   Jali::Parallel_type::OWNED,
                                   Constant1);
  source_state.add(myvec1);

  // Add another constant value state vector on source cells

  double Constant2 = -91.5;
  Jali::StateVector<double> myvec2("srccellvars2", source_mesh,
                                   Jali::Entity_kind::CELL,
                                   Jali::Parallel_type::OWNED,
                                   Constant2);
  source_state.add(myvec2);

  // Add a constant value state vector on source nodes

  double Constant3 = 3.14;
  Jali::StateVector<double> myvec3("srcnodevars", source_mesh,
                                   Jali::Entity_kind::NODE,
                                   Jali::Parallel_type::OWNED,
                                   Constant3);
  source_state.add(myvec3);


  // Add zero value state vectors on target cells and nodes - once with
  // the old name and once with the new name

  Jali::StateVector<double> myvec4("trgcellvars1", target_mesh,
                                   Jali::Entity_kind::CELL,
                                   Jali::Parallel_type::OWNED);

  target_state.add(myvec4);
  Jali::StateVector<double> myvec5("srccellvars1", target_mesh,
                                   Jali::Entity_kind::CELL,
                                   Jali::Parallel_type::OWNED);
  target_state.add(myvec5);

  Jali::StateVector<double> myvec6("trgcellvars2", target_mesh,
                                   Jali::Entity_kind::CELL,
                                   Jali::Parallel_type::OWNED);

  target_state.add(myvec6);
  Jali::StateVector<double> myvec7("srccellvars2", target_mesh,
                                   Jali::Entity_kind::CELL,
                                   Jali::Parallel_type::OWNED);
  target_state.add(myvec7);

  std::vector<double> zerodata2(nnodes_target, 0.0);
  Jali::StateVector<double> myvec8("trgnodevars", target_mesh,
                                   Jali::Entity_kind::NODE,
                                   Jali::Parallel_type::OWNED);

  target_state.add(myvec8);
  Jali::StateVector<double> myvec9("srcnodevars", target_mesh,
                                   Jali::Entity_kind::NODE,
                                   Jali::Parallel_type::OWNED);
  target_state.add(myvec9);

  // Wrappers for interfacing with the underlying mesh data structures.

  Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Portage::Jali_Mesh_Wrapper targetMeshWrapper(*target_mesh);

  // Wrappers for the source and target state managers

  Portage::Jali_State_Wrapper sourceStateWrapper(source_state);
  Portage::Jali_State_Wrapper targetStateWrapper(target_state);

  // Build the main driver object

  Portage::Driver<Portage::Jali_Mesh_Wrapper,
                  Portage::Jali_State_Wrapper> remapper(sourceMeshWrapper,
                                                        sourceStateWrapper,
                                                        targetMeshWrapper,
                                                        targetStateWrapper);
                                                      
  // Specify the fields to be remapped

  std::vector<std::string> source_var_names;
  source_var_names.push_back("srccellvars1");
  source_var_names.push_back("srccellvars2");
  source_var_names.push_back("srcnodevars");

  std::vector<std::string> target_var_names;
  target_var_names.push_back("trgcellvars1");
  target_var_names.push_back("trgcellvars2");
  target_var_names.push_back("trgnodevars");

  remapper.set_remap_var_names(source_var_names, target_var_names);

  // 1st order interpolation

  remapper.set_interpolation_order(1);

  // Execute remapper

  remapper.run();

  // Verify that we got the fields we wanted

  double *outcellvec1;
  targetStateWrapper.get_data(Portage::CELL, "trgcellvars1", &outcellvec1);

  for (int i = 0; i < ncells_target; i++)
    ASSERT_DOUBLE_EQ(Constant1, outcellvec1[i]);

  double *outcellvec2;
  targetStateWrapper.get_data(Portage::CELL, "trgcellvars2", &outcellvec2);

  for (int i = 0; i < ncells_target; i++)
    ASSERT_DOUBLE_EQ(Constant2, outcellvec2[i]);

  double *outnodevec;
  targetStateWrapper.get_data(Portage::NODE, "trgnodevars", &outnodevec);

  for (int i = 0; i < nnodes_target; i++)
    ASSERT_DOUBLE_EQ(Constant3, outnodevec[i]);



  // Remap between same name variables

  remapper.set_remap_var_names(source_var_names);

  // Execute remapper

  remapper.run();

  // Verify that we got the fields we wanted

  targetStateWrapper.get_data(Portage::CELL, "srccellvars1", &outcellvec1);

  for (int i = 0; i < ncells_target; i++)
    ASSERT_DOUBLE_EQ(Constant1, outcellvec1[i]);

  targetStateWrapper.get_data(Portage::CELL, "srccellvars2", &outcellvec2);

  for (int i = 0; i < ncells_target; i++)
    ASSERT_DOUBLE_EQ(Constant2, outcellvec2[i]);

  targetStateWrapper.get_data(Portage::NODE, "srcnodevars", &outnodevec);
  for (int i = 0; i < ncells_target; i++)
    ASSERT_DOUBLE_EQ(Constant3, outnodevec[i]);

}
