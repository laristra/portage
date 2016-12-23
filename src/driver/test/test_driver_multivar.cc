/*
Copyright (c) 2016, Los Alamos National Security, LLC
All rights reserved.

Copyright 2016. Los Alamos National Security, LLC. This software was produced
under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National
Laboratory (LANL), which is operated by Los Alamos National Security, LLC for
the U.S. Department of Energy. The U.S. Government has rights to use,
reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS
NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
derivative works, such modified software should be clearly marked, so as not to
confuse it with the version available from LANL.

Additionally, redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. Neither the name of Los Alamos National Security, LLC, Los Alamos
   National Laboratory, LANL, the U.S. Government, nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL
SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/



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


double TOL = 1e-12;

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
                                Jali::Entity_type::PARALLEL_OWNED);
  const int nnodes_target =
      target_mesh->num_entities(Jali::Entity_kind::NODE,
                                Jali::Entity_type::PARALLEL_OWNED);

  // Create state objects for source and target mesh

  Jali::State source_state(source_mesh);
  Jali::State target_state(target_mesh);

  // Add a constant value state vector on source cells

  double Constant1 = 1.25;
  Jali::StateVector<double> myvec1("srccellvars1", source_mesh,
                                   Jali::Entity_kind::CELL,
                                   Jali::Entity_type::PARALLEL_OWNED,
                                   Constant1);
  source_state.add(myvec1);

  // Add another constant value state vector on source cells

  double Constant2 = -91.5;
  Jali::StateVector<double> myvec2("srccellvars2", source_mesh,
                                   Jali::Entity_kind::CELL,
                                   Jali::Entity_type::PARALLEL_OWNED,
                                   Constant2);
  source_state.add(myvec2);

  // Add a constant value state vector on source nodes

  double Constant3 = 3.14;
  Jali::StateVector<double> myvec3("srcnodevars", source_mesh,
                                   Jali::Entity_kind::NODE,
                                   Jali::Entity_type::PARALLEL_OWNED,
                                   Constant3);
  source_state.add(myvec3);


  // Add zero value state vectors on target cells and nodes - once with
  // the old name and once with the new name

  Jali::StateVector<double> myvec4("trgcellvars1", target_mesh,
                                   Jali::Entity_kind::CELL,
                                   Jali::Entity_type::PARALLEL_OWNED);

  target_state.add(myvec4);
  Jali::StateVector<double> myvec5("srccellvars1", target_mesh,
                                   Jali::Entity_kind::CELL,
                                   Jali::Entity_type::PARALLEL_OWNED);
  target_state.add(myvec5);

  Jali::StateVector<double> myvec6("trgcellvars2", target_mesh,
                                   Jali::Entity_kind::CELL,
                                   Jali::Entity_type::PARALLEL_OWNED);

  target_state.add(myvec6);
  Jali::StateVector<double> myvec7("srccellvars2", target_mesh,
                                   Jali::Entity_kind::CELL,
                                   Jali::Entity_type::PARALLEL_OWNED);
  target_state.add(myvec7);

  std::vector<double> zerodata2(nnodes_target, 0.0);
  Jali::StateVector<double> myvec8("trgnodevars", target_mesh,
                                   Jali::Entity_kind::NODE,
                                   Jali::Entity_type::PARALLEL_OWNED);

  target_state.add(myvec8);
  Jali::StateVector<double> myvec9("srcnodevars", target_mesh,
                                   Jali::Entity_kind::NODE,
                                   Jali::Entity_type::PARALLEL_OWNED);
  target_state.add(myvec9);

  // Wrappers for interfacing with the underlying mesh data structures.

  Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Portage::Jali_Mesh_Wrapper targetMeshWrapper(*target_mesh);

  // Wrappers for the source and target state managers

  Portage::Jali_State_Wrapper sourceStateWrapper(source_state);
  Portage::Jali_State_Wrapper targetStateWrapper(target_state);

  // Build the main driver object

  /////////

  Portage::Driver<Portage::SearchKDTree, 
      Portage::IntersectR2D, 
          Portage::Interpolate_1stOrder,
          2,
          Portage::Jali_Mesh_Wrapper, 
          Portage::Jali_State_Wrapper>  
          remapper(sourceMeshWrapper, sourceStateWrapper, targetMeshWrapper,
                   targetStateWrapper);
  /////////

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

  // Execute remapper

  remapper.run(false);

  // Verify that we got the fields we wanted

  double *outcellvec1;
  targetStateWrapper.get_data(Portage::CELL, "trgcellvars1", &outcellvec1);

  for (int i = 0; i < ncells_target; i++)
    ASSERT_NEAR(Constant1, outcellvec1[i], TOL);

  double *outcellvec2;
  targetStateWrapper.get_data(Portage::CELL, "trgcellvars2", &outcellvec2);

  for (int i = 0; i < ncells_target; i++)
    ASSERT_NEAR(Constant2, outcellvec2[i], TOL);

  // double *outnodevec;
  // targetStateWrapper.get_data(Portage::NODE, "trgnodevars", &outnodevec);

  // for (int i = 0; i < nnodes_target; i++)
  //   ASSERT_NEAR(Constant3, outnodevec[i], TOL);


  // Remap between same name variables

  remapper.set_remap_var_names(source_var_names);

  // Execute remapper

  remapper.run(false);

  // Verify that we got the fields we wanted

  targetStateWrapper.get_data(Portage::CELL, "srccellvars1", &outcellvec1);

  for (int i = 0; i < ncells_target; i++)
    ASSERT_NEAR(Constant1, outcellvec1[i], TOL);

  targetStateWrapper.get_data(Portage::CELL, "srccellvars2", &outcellvec2);

  for (int i = 0; i < ncells_target; i++)
    ASSERT_NEAR(Constant2, outcellvec2[i], TOL);

  // targetStateWrapper.get_data(Portage::NODE, "srcnodevars", &outnodevec);
  // for (int i = 0; i < ncells_target; i++)
  //   ASSERT_NEAR(Constant3, outnodevec[i], TOL);
}


TEST(Test_MultiVar_Remap, Nested_Meshes) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);

  std::shared_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 1.0, 1.0, 2, 2);
  std::shared_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4);

  const int ncells_source =
      source_mesh->num_entities(Jali::Entity_kind::CELL,
                                Jali::Entity_type::PARALLEL_OWNED);
  const int ncells_target =
      target_mesh->num_entities(Jali::Entity_kind::CELL,
                                Jali::Entity_type::PARALLEL_OWNED);

  // Create state objects for source and target mesh

  Jali::State source_state(source_mesh);
  Jali::State target_state(target_mesh);

  // Add a state vector on source cells with values dependent on the
  // centroid of each cell

  Jali::StateVector<double> sourcevec("cellvars", source_mesh,
                                      Jali::Entity_kind::CELL,
                                      Jali::Entity_type::PARALLEL_OWNED);
  for (int c = 0; c < ncells_source; ++c) {
    JaliGeometry::Point ccen = source_mesh->cell_centroid(c);
    sourcevec[c] = ccen[0] + ccen[1];
  }
  source_state.add(sourcevec);

  // Add zero value state vectors on target cells and nodes - once with
  // the old name and once with the new name

  Jali::StateVector<double>& targetvec =
      target_state.add<double>("cellvars", target_mesh,
                               Jali::Entity_kind::CELL,
                               Jali::Entity_type::PARALLEL_OWNED);
  
  // Wrappers for interfacing with the underlying mesh data structures.

  Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*source_mesh);
  Portage::Jali_Mesh_Wrapper targetMeshWrapper(*target_mesh);

  // Wrappers for the source and target state managers

  Portage::Jali_State_Wrapper sourceStateWrapper(source_state);
  Portage::Jali_State_Wrapper targetStateWrapper(target_state);

  // Build the main driver object

  Portage::Driver<Portage::SearchKDTree,
                  Portage::IntersectR2D,
                  Portage::Interpolate_1stOrder,
                  2,
                  Portage::Jali_Mesh_Wrapper,
                  Portage::Jali_State_Wrapper> remapper1(sourceMeshWrapper,
                                                         sourceStateWrapper,
                                                         targetMeshWrapper,
                                                         targetStateWrapper);

  // Specify the fields to be remapped

  std::vector<std::string> source_var_names;
  source_var_names.push_back("cellvars");

  std::vector<std::string> target_var_names;
  target_var_names.push_back("cellvars");

  remapper1.set_remap_var_names(source_var_names, target_var_names);

  // Execute remapper (distributed=false)

  remapper1.run(false);

  // Verify that we got the fields we wanted
  for (int c = 0; c < ncells_target; c++) {
    JaliGeometry::Point ccen = target_mesh->cell_centroid(c);
    double x, y;
    if (ccen[0] < 0.5) x = 0.25; else x = 0.75;
    if (ccen[1] < 0.5) y = 0.25; else y = 0.75;
    double expval = x + y;

    ASSERT_NEAR(expval, targetvec[c], TOL);
  }


  // Build the main driver object

  Portage::Driver<Portage::SearchKDTree,
                  Portage::IntersectR2D,
                  Portage::Interpolate_2ndOrder,
                  2,
                  Portage::Jali_Mesh_Wrapper,
                  Portage::Jali_State_Wrapper> remapper2(sourceMeshWrapper,
                                                         sourceStateWrapper,
                                                         targetMeshWrapper,
                                                         targetStateWrapper);


  remapper2.set_remap_var_names(source_var_names, target_var_names);

  // Execute remapper (distributed=false)

  remapper2.run(false);

  // Verify that we got the fields we wanted
  for (int c = 0; c < ncells_target; c++) {
    JaliGeometry::Point ccen = target_mesh->cell_centroid(c);
    double expval = ccen[0] + ccen[1];

    ASSERT_NEAR(expval, targetvec[c], TOL);
  }


}

