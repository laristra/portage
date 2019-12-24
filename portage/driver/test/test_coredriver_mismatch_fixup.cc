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
#include "portage/driver/coredriver.h"
#include "portage/driver/fix_mismatch.h"
#include "portage/search/search_kdtree.h"
#include "portage/intersect/intersect_rNd.h"
#include "portage/interpolate/interpolate_nth_order.h"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliState.h"
#include "JaliStateVector.h"


double TOL = 1e-12;

TEST(Test_Mismatch_Fixup, Test_Methods) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  if (Jali::framework_available(Jali::MSTK))
    mf.framework(Jali::MSTK);
    
  // Create meshes
  std::shared_ptr<Jali::Mesh> source_mesh = mf(-0.8, 0.0, 0.4, 1.0, 1, 1);
  std::shared_ptr<Jali::Mesh> target_mesh = mf( 0.0, 0.0, 2.0, 1.0, 2, 1);

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

  // Create mesh wrappers
  Wonton::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Wonton::Jali_Mesh_Wrapper targetMeshWrapper(*target_mesh);

  // Create state wrappers
  Wonton::Jali_State_Wrapper sourceStateWrapper(*source_state);
  Wonton::Jali_State_Wrapper targetStateWrapper(*target_state);


  
  // Create the driver
  Portage::CoreDriver<2,Wonton::Entity_kind::CELL,
                       Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
                       Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper>
      driver(sourceMeshWrapper, sourceStateWrapper,
        targetMeshWrapper, targetStateWrapper);

  auto candidates = driver.search<Portage::SearchKDTree>();
  auto source_weights = driver.intersect_meshes<Portage::IntersectR2D>(candidates);
  auto gradients = driver.compute_source_gradient("cellvars");

  // Create the mismatch fixer
  Portage::MismatchFixer<2, Portage::Entity_kind::CELL, Wonton::Jali_Mesh_Wrapper, 
    Wonton::Jali_State_Wrapper, Wonton::Jali_Mesh_Wrapper,  Wonton::Jali_State_Wrapper> 
    fixer(sourceMeshWrapper, sourceStateWrapper, targetMeshWrapper, targetStateWrapper, source_weights, nullptr);
      

  //-------------------------------------------------------------------
  // Expected Results
  //-------------------------------------------------------------------  
  
  // exact results for Partial_fixup_type:
  //	CONSTANT (C)
  //  LOCALLY_CONSERVATIVE (L)
  //  SHIFTED_CONSERVATIVE (S)
  // and Empty_fixup_type:
  // 	EXTRAPOLATE (E)
  //  LEAVE_EMPTY (L)

  double exact_C_E[2] = {1.0, 1.0};
  double exact_L_E[2] = {0.4, 0.4};
  double exact_S_E[2] = {0.6, 0.6};
  double exact_C_L[2] = {1.0, 0.0};
  double exact_L_L[2] = {0.4, 0.0};
  double exact_S_L[2] = {1.2, 0.0};


  //-------------------------------------------------------------------
  // Run Tests
  // Note: We need to run the interpolate_mesh_var every time, because 
  // mismatch fixup modifies the target state in place.
  // Interpolation no longer does any repair
  //-------------------------------------------------------------------      
  
  double dblmin = -std::numeric_limits<double>::max();
  double dblmax =  std::numeric_limits<double>::max();

  driver.interpolate_mesh_var<double,Portage::Interpolate_1stOrder>(
    "cellvars","cellvars", source_weights, &gradients);

  // test (C,E)
  if (fixer.has_mismatch())
    fixer.fix_mismatch("cellvars","cellvars", 0.0, dblmax, Portage::DEFAULT_CONSERVATION_TOL, 
      Portage::DEFAULT_MAX_FIXUP_ITER, Portage::CONSTANT, Portage::EXTRAPOLATE);

  for (int c = 0; c < ncells_target; c++) {
    ASSERT_NEAR(exact_C_E[c], targetvec[c], TOL);
  }


  driver.interpolate_mesh_var<double,Portage::Interpolate_1stOrder>(
    "cellvars","cellvars", source_weights, &gradients);

  // test (L,E)
  if (fixer.has_mismatch())
    fixer.fix_mismatch("cellvars","cellvars", 0.0, dblmax, Portage::DEFAULT_CONSERVATION_TOL, 
      Portage::DEFAULT_MAX_FIXUP_ITER, Portage::LOCALLY_CONSERVATIVE, Portage::EXTRAPOLATE);

  for (int c = 0; c < ncells_target; c++) {
    ASSERT_NEAR(exact_L_E[c], targetvec[c], TOL);
  }

 
  driver.interpolate_mesh_var<double,Portage::Interpolate_1stOrder>(
    "cellvars","cellvars", source_weights, &gradients);

  // test (S,E)
  if (fixer.has_mismatch())
    fixer.fix_mismatch("cellvars","cellvars", 0.0, dblmax, Portage::DEFAULT_CONSERVATION_TOL, 
      Portage::DEFAULT_MAX_FIXUP_ITER, Portage::SHIFTED_CONSERVATIVE, Portage::EXTRAPOLATE);

  for (int c = 0; c < ncells_target; c++) {
    ASSERT_NEAR(exact_S_E[c], targetvec[c], TOL);
  }


  driver.interpolate_mesh_var<double,Portage::Interpolate_1stOrder>(
    "cellvars","cellvars", source_weights, &gradients);

  // test (C,L)
  if (fixer.has_mismatch())
    fixer.fix_mismatch("cellvars","cellvars", 0.0, dblmax, Portage::DEFAULT_CONSERVATION_TOL, 
      Portage::DEFAULT_MAX_FIXUP_ITER, Portage::CONSTANT, Portage::LEAVE_EMPTY);

  for (int c = 0; c < ncells_target; c++) {
    ASSERT_NEAR(exact_C_L[c], targetvec[c], TOL);
  }

 
  driver.interpolate_mesh_var<double,Portage::Interpolate_1stOrder>(
    "cellvars","cellvars", source_weights, &gradients);

  // test (L,L)
  if (fixer.has_mismatch())
    fixer.fix_mismatch("cellvars","cellvars", 0.0, dblmax, Portage::DEFAULT_CONSERVATION_TOL, 
      Portage::DEFAULT_MAX_FIXUP_ITER, Portage::LOCALLY_CONSERVATIVE, Portage::LEAVE_EMPTY);

  for (int c = 0; c < ncells_target; c++) {
    ASSERT_NEAR(exact_L_L[c], targetvec[c], TOL);
  }

 
  driver.interpolate_mesh_var<double,Portage::Interpolate_1stOrder>(
    "cellvars","cellvars", source_weights, &gradients);

  // test (S,L)
  if (fixer.has_mismatch())
    fixer.fix_mismatch("cellvars","cellvars", 0.0, dblmax, Portage::DEFAULT_CONSERVATION_TOL, 
      Portage::DEFAULT_MAX_FIXUP_ITER, Portage::SHIFTED_CONSERVATIVE, Portage::LEAVE_EMPTY);

  for (int c = 0; c < ncells_target; c++) {
    ASSERT_NEAR(exact_S_L[c], targetvec[c], TOL);
  }

}
