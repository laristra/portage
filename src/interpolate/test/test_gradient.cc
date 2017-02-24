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



#include "portage/interpolate/gradient.h"

#include <iostream>

#include "gtest/gtest.h"
#include "mpi.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "FrameworkTraits.hh"
#include "MeshFramework.hh"
#include "JaliState.h"
#include "JaliStateVector.h"

#include "portage/support/portage.h"
#include "portage/support/Vector.h"
#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"
#include "portage/wrappers/state/jali/jali_state_wrapper.h"

/// Test gradient computation for cell centered fields

TEST(Gradient, Fields_Cell_Ctr) {

  // Make a 4x4 mesh

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);
  mf.included_entities({Jali::Entity_kind::EDGE,
                        Jali::Entity_kind::FACE,
                        Jali::Entity_kind::WEDGE,
                        Jali::Entity_kind::CORNER});

  std::shared_ptr<Jali::Mesh> mesh1 = mf(0.0, 0.0, 1.0, 1.0, 4, 4);
  ASSERT_TRUE(mesh1 != nullptr);

  // Create a state object and add the first two vectors to it

  Jali::State mystate(mesh1);

  // Define three state vectors, one with constant value and the other
  // with a linear function that is x+2y

  int nc1 = mesh1->num_entities(Jali::Entity_kind::CELL,
                                Jali::Entity_type::PARALLEL_OWNED);
  std::vector<double> data1(nc1, 1.25);
  Jali::StateVector<double> myvec1("cellvars1", mesh1,
                                   Jali::Entity_kind::CELL,
                                   Jali::Entity_type::PARALLEL_OWNED,
                                   &(data1[0]));
  Jali::StateVector<double> &addvec1 = mystate.add(myvec1);


  std::vector<double> data2(nc1);
  for (int c = 0; c < nc1; c++) {
    JaliGeometry::Point ccen = mesh1->cell_centroid(c);
    data2[c] = ccen[0]+2*ccen[1];
  }

  Jali::StateVector<double> myvec2("cellvars2", mesh1,
                                   Jali::Entity_kind::CELL,
                                   Jali::Entity_type::PARALLEL_OWNED,
                                   &(data2[0]));
  Jali::StateVector<double> &addvec2 = mystate.add(myvec2);

  Portage::Jali_Mesh_Wrapper meshwrapper(*mesh1);
  Portage::Jali_State_Wrapper statewrapper(mystate);

  // Create Gradient objects

  Portage::Limited_Gradient<Portage::Jali_Mesh_Wrapper,
                            Portage::Jali_State_Wrapper,
                            Portage::CELL, 2>
      gradcalc1(meshwrapper, statewrapper, "cellvars1", Portage::NOLIMITER);
  Portage::Limited_Gradient<Portage::Jali_Mesh_Wrapper,
                            Portage::Jali_State_Wrapper,
                            Portage::CELL, 2>
      gradcalc2(meshwrapper, statewrapper, "cellvars2", Portage::NOLIMITER);
  Portage::Limited_Gradient<Portage::Jali_Mesh_Wrapper,
                            Portage::Jali_State_Wrapper,
                            Portage::CELL, 2>
      gradcalc3(meshwrapper, statewrapper, "cellvars1",
                Portage::BARTH_JESPERSEN);
  Portage::Limited_Gradient<Portage::Jali_Mesh_Wrapper,
                            Portage::Jali_State_Wrapper,
                            Portage::CELL, 2>
      gradcalc4(meshwrapper, statewrapper, "cellvars2",
                Portage::BARTH_JESPERSEN);


  // Compute the gradient for each of these fields

  Portage::Vector<2> grad;

  // Verify the gradient values
  // For field 1 (constant), it is is 0,0
  // For field 2 (x+2y), it is (1,2)

  for (int c = 0; c < nc1; ++c) {
    // unlimited gradient of constant function

    grad = gradcalc1(c);
    ASSERT_NEAR(0.0, grad[0], 1.0e-10);
    ASSERT_NEAR(0.0, grad[1], 1.0e-10);

    // unlimited gradient of linear function

    grad = gradcalc2(c);
    ASSERT_NEAR(1.0, grad[0], 1.0e-10);
    ASSERT_NEAR(2.0, grad[1], 1.0e-10);

    // limited gradient of constant function

    grad = gradcalc3(c);
    ASSERT_NEAR(0.0, grad[0], 1.0e-10);
    ASSERT_NEAR(0.0, grad[1], 1.0e-10);

    // limited gradient of linear function
    //
    // For now, the limiter does not know anything about boundary
    // conditions and therefore, it is a little unpredictable what the
    // gradient will be limited to on boundary cells (which are not
    // completely surrounded by other cells). So check only interior
    // cells.

    grad = gradcalc4(c);

    bool boundary_cell = false;
    std::vector<int> cfaces;
    mesh1->cell_get_faces(c, &cfaces);
    for (auto f : cfaces) {
      std::vector<int> fcells;
      mesh1->face_get_cells(f, Jali::Entity_type::ALL, &fcells);
      if (fcells.size() == 1) {
        boundary_cell = true;
        break;
      }
    }

    if (!boundary_cell) {
      ASSERT_NEAR(1.0, grad[0], 1.0e-10);
      ASSERT_NEAR(2.0, grad[1], 1.0e-10);
    }
  }
}

/// Test gradient computation with node centered fields

TEST(Gradient, Fields_Node_Ctr) {

  // Make a 3x3 mesh

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  Jali::FrameworkPreference pref;
  pref.push_back(Jali::MSTK);
  if (Jali::framework_available(Jali::MSTK))
    mf.preference(pref);
  mf.included_entities({Jali::Entity_kind::EDGE,
                        Jali::Entity_kind::FACE,
                        Jali::Entity_kind::WEDGE,
                        Jali::Entity_kind::CORNER});

  std::shared_ptr<Jali::Mesh> mesh1 = mf(0.0, 0.0, 1.0, 1.0, 3, 3);
  ASSERT_TRUE(mesh1 != nullptr);

  // Create a state object and add the first two vectors to it

  Jali::State mystate(mesh1);

  // Define three state vectors, one with constant value, the other
  // with a linear function

  int nn1 = mesh1->num_entities(Jali::Entity_kind::NODE,
                                Jali::Entity_type::PARALLEL_OWNED);

  std::vector<double> data1(nn1, 1.5);

  Jali::StateVector<double> myvec1("nodevars1", mesh1,
                                   Jali::Entity_kind::NODE,
                                   Jali::Entity_type::PARALLEL_OWNED,
                                   &(data1[0]));
  Jali::StateVector<double> &addvec1 = mystate.add(myvec1);

  std::vector<double> data2(nn1);
  for (int n = 0; n < nn1; ++n) {
    JaliGeometry::Point nodexy;
    mesh1->node_get_coordinates(n, &nodexy);
    data2[n] = 3*nodexy[0]+nodexy[1];
  }
  Jali::StateVector<double> myvec2("nodevars2", mesh1,
                                   Jali::Entity_kind::NODE,
                                   Jali::Entity_type::PARALLEL_OWNED,
                                   &(data2[0]));
  Jali::StateVector<double> &addvec2 = mystate.add(myvec2);

  // Create Gradient calculater objects

  Portage::Jali_Mesh_Wrapper meshwrapper(*mesh1);
  Portage::Jali_State_Wrapper statewrapper(mystate);

  Portage::Limited_Gradient<Portage::Jali_Mesh_Wrapper,
                            Portage::Jali_State_Wrapper,
                            Portage::NODE, 2>
      gradcalc1(meshwrapper, statewrapper, "nodevars1", Portage::NOLIMITER);
  Portage::Limited_Gradient<Portage::Jali_Mesh_Wrapper,
                            Portage::Jali_State_Wrapper,
                            Portage::NODE, 2>
      gradcalc2(meshwrapper, statewrapper, "nodevars2", Portage::NOLIMITER);
  Portage::Limited_Gradient<Portage::Jali_Mesh_Wrapper,
                            Portage::Jali_State_Wrapper,
                            Portage::NODE, 2>
      gradcalc3(meshwrapper, statewrapper, "nodevars1",
                Portage::BARTH_JESPERSEN);
  Portage::Limited_Gradient<Portage::Jali_Mesh_Wrapper,
                            Portage::Jali_State_Wrapper,
                            Portage::NODE, 2>
      gradcalc4(meshwrapper, statewrapper, "nodevars2",
                Portage::BARTH_JESPERSEN);


  // Make sure we retrieved the correct gradient value for each node
  // For field 1, it is a constant
  // For field 2, it is a linear function

  Portage::Vector<2> grad;

  for (int n = 0; n < nn1; ++n) {
    // unlimited gradient of constant function

    grad = gradcalc1(n);
    ASSERT_NEAR(0.0, grad[0], 1.0e-10);
    ASSERT_NEAR(0.0, grad[1], 1.0e-10);

    // unlimited gradient of linear function

    grad = gradcalc2(n);
    ASSERT_NEAR(3.0, grad[0], 1.0e-10);
    ASSERT_NEAR(1.0, grad[1], 1.0e-10);

    // limited gradient of constant function

    grad = gradcalc3(n);
    ASSERT_NEAR(0.0, grad[0], 1.0e-10);
    ASSERT_NEAR(0.0, grad[1], 1.0e-10);

    // unlimited gradient of linear function
    //
    // For now, the limiter does not know anything about boundary
    // conditions and therefore, it is a little unpredictable what the
    // limiting will be in boundary nodes/dual-cells (which are not
    // completely surrounded by other dual cells). So check only
    // interior nodes.

    grad = gradcalc4(n);

    // use the primary mesh to check if the node (dual cell) is a
    // boundary node and don't check answer if its a boundary node

    bool boundary_node = false;
    std::vector<int> nodecells;
    mesh1->node_get_cells(n, Jali::Entity_type::ALL, &nodecells);
    for (auto nc : nodecells) {
      std::vector<int> cfaces;
      mesh1->cell_get_faces(nc, &cfaces);
      for (auto f : cfaces) {
        std::vector<int> fcells;
        mesh1->face_get_cells(f, Jali::Entity_type::ALL, &fcells);
        if (fcells.size() == 1) {
          boundary_node = true;
          break;
        }
      }
      if (boundary_node)
        break;
    }

    if (!boundary_node) {
      ASSERT_NEAR(3.0, grad[0], 1.0e-10);
      ASSERT_NEAR(1.0, grad[1], 1.0e-10);
    }
  }
}


