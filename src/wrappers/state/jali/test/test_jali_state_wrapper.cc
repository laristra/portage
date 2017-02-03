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



#include "portage/wrappers/state/jali/jali_state_wrapper.h"

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


TEST(Jali_State_Wrapper, DataTypes) {

  // Add multiple state vector types
  int const n_cells = 4;
  int const n_nodes = 9;
  float ftest[] = {1.1, 2.2, 3.3, 4.4};
  double dtest[] = {62., 78., 43., 22.};
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

  // Get raw float data using wrapper
  float* fdata;
  wrapper.get_data(Portage::CELL, "f1", &fdata);
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

  // add data through wrapper
  wrapper.add_data(inputMesh, Portage::CELL, "d1", dtest);
  wrapper.add_data(inputMesh, Portage::CELL, "d2", 123.456);

  double *ddata;
  wrapper.get_data(Portage::CELL, "d1", &ddata);
  for (unsigned int i=0; i<n_cells; i++) ASSERT_EQ(ddata[i], dtest[i]);
  wrapper.get_data(Portage::CELL, "d2", &ddata);
  for (unsigned int i=0; i<n_cells; i++) ASSERT_EQ(ddata[i], 123.456);

  // Iterate through a vector of names
  std::vector<std::string> fields;
  fields.push_back("v1");  fields.push_back("f1");  fields.push_back("i1");
#if 0
  for (auto it = fields.begin(); it != fields.end(); it++)
  {
    Portage::Entity_kind on_what = wrapper.get_entity(*it);
    if (typeid(float) == wrapper.get_type(*it))
    {
      float* fdata;
      wrapper.get_data(on_what, *it, &fdata);
      for (unsigned int i=0; i<inputMeshWrapper.num_entities(on_what); i++) ASSERT_EQ(fdata[i], ftest[i]);
    }
    else if (typeid(int) == wrapper.get_type(*it))
    {
      int* idata;
      wrapper.get_data(on_what, *it, &idata);
      for (unsigned int i=0; i<inputMeshWrapper.num_entities(on_what); i++) ASSERT_EQ(idata[i], itest[i]);
    }
    else if (typeid(Vec2d) == wrapper.get_type(*it))
    {
      Vec2d* vdata;
      wrapper.get_data(on_what, *it, &vdata);
      for (unsigned int i=0; i<inputMeshWrapper.num_entities(on_what); i++)
      {
        ASSERT_EQ(vdata[i].x, vtest[i].x);
        ASSERT_EQ(vdata[i].y, vtest[i].y);
      }
    }
    else
    {
      ASSERT_EQ(0, 1);    // This else should never be reached in this test
    }
  }
#endif

#if 0
  // Iterate through all fields using the wrapper
  for (auto it = wrapper.names_begin(); it != wrapper.names_end(); it++)
  {
    Portage::Entity_kind on_what = wrapper.get_entity(*it);
    if (typeid(float) == wrapper.get_type(*it))
    {
      float* fdata;
      wrapper.get_data(on_what, *it, &fdata);
      for (unsigned int i=0; i<inputMeshWrapper.num_entities(on_what); i++) ASSERT_EQ(fdata[i], ftest[i]);
    }
    else if (typeid(int) == wrapper.get_type(*it))
    {
      int* idata;
      wrapper.get_data(on_what, *it, &idata);
      for (unsigned int i=0; i<inputMeshWrapper.num_entities(on_what); i++) ASSERT_EQ(idata[i], itest[i]);
    }
    else if (typeid(Vec2d) == wrapper.get_type(*it))
    {
      Vec2d* vdata;
      wrapper.get_data(on_what, *it, &vdata);
      for (unsigned int i=0; i<inputMeshWrapper.num_entities(on_what); i++)
      {
        ASSERT_EQ(vdata[i].x, vtest[i].x);
        ASSERT_EQ(vdata[i].y, vtest[i].y);
      }
    }
    else
    {
      ASSERT_EQ(0, 1);    // This else should never be reached in this test
    }
  }
#endif

#if 0
  // Iterate through fields on cells only using the wrapper
  for (auto it = wrapper.names_entity_begin(Portage::CELL); it != wrapper.names_entity_end(Portage::CELL); it++)
  {
    Portage::Entity_kind on_what = wrapper.get_entity(*it);
    if (typeid(float) == wrapper.get_type(*it))
    {
      float* fdata;
      wrapper.get_data(on_what, *it, &fdata);
      for (unsigned int i=0; i<inputMeshWrapper.num_entities(on_what); i++) ASSERT_EQ(fdata[i], ftest[i]);
    }
    else if (typeid(int) == wrapper.get_type(*it))
    {
      ASSERT_EQ(0, 1);   // This else should never be reached in this test
    }
    else if (typeid(Vec2d) == wrapper.get_type(*it))
    {
      Vec2d* vdata;
      wrapper.get_data(on_what, *it, &vdata);
      for (unsigned int i=0; i<inputMeshWrapper.num_entities(on_what); i++)
      {
        ASSERT_EQ(vdata[i].x, vtest[i].x);
        ASSERT_EQ(vdata[i].y, vtest[i].y);
      }
    }
    else
    {
      ASSERT_EQ(0, 1);   // This else should never be reached in this test
    }
  }
#endif

}

  

  

