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
#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "gtest/gtest.h"
#ifdef ENABLE_MPI
#include "mpi.h"
#endif

#include "portage/driver/driver.h"
#include "portage/driver/driver_mesh_swarm_mesh.h"
#include "portage/simple_mesh/simple_mesh.h"
#include "portage/simple_mesh/simple_state.h"
#include "portage/intersect/intersect_r2d.h"
#include "portage/intersect/intersect_r3d.h"
#include "portage/interpolate/interpolate_1st_order.h"
#include "portage/interpolate/interpolate_2nd_order.h"
#include "portage/search/search_points_by_cells.h"
#include "portage/accumulate/accumulate.h"
#include "portage/estimate/estimate.h"
#include "portage/support/Point.h"
#include "portage/wonton/mesh/simple_mesh/simple_mesh_wrapper.h"
#include "portage/wonton/state/simple_state/simple_state_wrapper.h"
#include "portage/wonton/mesh/flat/flat_mesh_wrapper.h"
// amh: TODO--change to simple mesh?
namespace {

double TOL = 1e-6;

// This is a set of integration tests based off of main.cc.  There
// will be at least one test corresponding to each case found in
// main.cc.  This is a test fixture and must be derived from the
// ::testing::Test class.  Specializations of this class, such as
// 2D/3D coincident and non-coincident remaps should be derived from
// this.
class MSMDriverTest : public ::testing::Test {
 protected:
  // Source and target meshes
  std::shared_ptr<Portage::Simple_Mesh> sourceMesh;
  std::shared_ptr<Portage::Simple_Mesh> targetMesh;
  // Source and target mesh state
  Portage::Simple_State sourceState;
  Portage::Simple_State targetState;
  Portage::Simple_State targetState2;
  //  Wrappers for interfacing with the underlying mesh data structures
  Wonton::Simple_Mesh_Wrapper sourceMeshWrapper;
  Wonton::Simple_Mesh_Wrapper targetMeshWrapper;
  Wonton::Simple_State_Wrapper sourceStateWrapper;
  Wonton::Simple_State_Wrapper targetStateWrapper;
  Wonton::Simple_State_Wrapper targetStateWrapper2;

  //  This is the basic test method to be called for each unit test. It will work
  //  for 2-D and 3-D, coincident and non-coincident cell-centered remaps.
  template <
    template<Portage::Entity_kind, class, class, class,
    template<class, int, class, class> class> class Intersect,
    template<int, Portage::Entity_kind, class, class, class,
    template<class, int, class, class> class> class Interpolate,
    template <int, class, class> class SwarmSearch,
    int Dimension = 3
  >
  void unitTest(double compute_initial_field(Portage::Point<Dimension> centroid),
                double smoothing_factor, Portage::Meshfree::Basis::Type basis, 
		Portage::Meshfree::WeightCenter center=Portage::Meshfree::Gather)
  {
    if (Dimension != 3) {
      throw std::runtime_error("2D not available yet");
    }

    //  Fill the source state data with the specified profile
    const int nsrccells = sourceMeshWrapper.num_owned_cells();
    std::vector<double> sourceData(nsrccells);
    const int nsrcnodes = sourceMeshWrapper.num_owned_nodes();
    std::vector<double> sourceDataNode(nsrcnodes);

    // Create the source data for given function
    Wonton::Flat_Mesh_Wrapper<double> sourceFlatMesh;
    sourceFlatMesh.initialize(sourceMeshWrapper);
    for (unsigned int c = 0; c < nsrccells; ++c) {
      Portage::Point<Dimension> cen;
      sourceFlatMesh.cell_centroid(c, &cen);
      sourceData[c] = compute_initial_field(cen);
    }
    Portage::Simple_State::vec &sourceVec(sourceState.add("celldata",
      Portage::Entity_kind::CELL, &(sourceData[0])));

    for (unsigned int c = 0; c < nsrcnodes; ++c) {
      Portage::Point<Dimension> cen;
      sourceFlatMesh.node_get_coordinates(c, &cen);
      sourceDataNode[c] = compute_initial_field(cen);
    }
    Portage::Simple_State::vec &sourceVecNode(sourceState.add("nodedata",
      Portage::Entity_kind::NODE, &(sourceDataNode[0])));

    // Build the target state storage
    const int ntarcells = targetMeshWrapper.num_owned_cells();
    const int ntarnodes = targetMeshWrapper.num_owned_nodes();
    std::vector<double> targetData(ntarcells), targetData2(ntarcells);
    std::vector<double> targetDataNode(ntarnodes), targetData2Node(ntarnodes);
    Portage::Simple_State::vec &targetVec(targetState.add("celldata",
      Portage::Entity_kind::CELL, &(targetData[0])));
    Portage::Simple_State::vec &targetVec2(targetState2.add("celldata",
      Portage::Entity_kind::CELL, &(targetData2[0])));
    Portage::Simple_State::vec &targetVecNode(targetState.add("nodedata",
      Portage::Entity_kind::NODE, &(targetDataNode[0])));
    Portage::Simple_State::vec &targetVec2Node(targetState2.add("nodedata",
      Portage::Entity_kind::NODE, &(targetData2Node[0])));

    //  Register the variable name and interpolation order with the driver
    std::vector<std::string> remap_fields;
    remap_fields.push_back("celldata");
    remap_fields.push_back("nodedata");

    //  Build the mesh-mesh driver data for this mesh type
    Portage::Driver<Portage::SearchKDTree,
    Intersect,
    Interpolate,
    Dimension,
    Wonton::Simple_Mesh_Wrapper, Wonton::Simple_State_Wrapper,
    Wonton::Simple_Mesh_Wrapper, Wonton::Simple_State_Wrapper>
    mmdriver(sourceMeshWrapper, sourceStateWrapper,
             targetMeshWrapper, targetStateWrapper);
    mmdriver.set_remap_var_names(remap_fields);
    // run on one processor
    mmdriver.run(false);

    //  Build the mesh-swarm-mesh driver data for this mesh type
    Portage::MSM_Driver<
      SwarmSearch,
      Portage::Meshfree::Accumulate,
      Portage::Meshfree::Estimate,
      Dimension,
      Wonton::Simple_Mesh_Wrapper, Wonton::Simple_State_Wrapper,
      Wonton::Simple_Mesh_Wrapper, Wonton::Simple_State_Wrapper
    >
    msmdriver(sourceMeshWrapper, sourceStateWrapper,
              targetMeshWrapper, targetStateWrapper2,
              smoothing_factor, basis, 
	      Portage::Meshfree::LocalRegression, 
	      Portage::Meshfree::Weight::TENSOR, 
	      Portage::Meshfree::Weight::B4, 
	      center);
    msmdriver.set_remap_var_names(remap_fields);
    // run on one processor
    msmdriver.run(false);

    // Check the answer
    double stdval, err;
    double toterr = 0.;

    Portage::Simple_State::vec &cellvecout(targetState.get("celldata",
                                                           Portage::CELL));
    Portage::Simple_State::vec &cellvecout2(targetState2.get("celldata",
                                                             Portage::CELL));
    Portage::Simple_State::vec &nodevecout(targetState.get("nodedata",
                                                           Portage::NODE));
    Portage::Simple_State::vec &nodevecout2(targetState2.get("nodedata",
                                                             Portage::NODE));

    Wonton::Flat_Mesh_Wrapper<double> targetFlatMesh;
    targetFlatMesh.initialize(targetMeshWrapper);
    for (int c = 0; c < ntarcells; ++c) {
      Portage::Point<Dimension> ccen;
      targetFlatMesh.cell_centroid(c, &ccen);
      double value = compute_initial_field(ccen);
      double merror, serror;
      merror = cellvecout[c] - value;
      serror = cellvecout2[c] - value;
      //  dump diagnostics for each cell
      std::printf("Cell=% 4d Centroid=(% 5.3lf,% 5.3lf,% 5.3lf)", c,
                  ccen[0], ccen[1], ccen[2]);
      std::printf(" Val=% 10.6lf MM=% 10.6lf Err=% 10.6lf MSM=% 10.6lf Err=% lf\n",
                  value, cellvecout[c], merror, cellvecout2[c], serror);
      toterr = std::max(toterr, std::fabs(serror));
    }

    std::printf("\n\nLinf NORM OF MSM CELL ERROR = %lf\n\n", toterr);
    ASSERT_LT(toterr, TOL);

    toterr = 0.;
    for (int n = 0; n < ntarnodes; ++n) {
      Portage::Point<Dimension> node;
      targetFlatMesh.node_get_coordinates(n, &node);
      double value = compute_initial_field(node);
      double merror, serror;
      merror = nodevecout[n] - value;
      serror = nodevecout2[n] - value;
      //  dump diagnostics for each node
      std::printf("Node=% 4d Coords=(% 5.3lf,% 5.3lf,% 5.3lf)", n,
                  node[0], node[1], node[2]);
      std::printf(" Val=% 10.6lf MM=% 10.6lf Err=% 10.6lf MSM=% 10.6lf Err=% lf\n",
                  value, nodevecout[n], merror, nodevecout2[n], serror);
      toterr = std::max(toterr, std::fabs(serror));
    }

    std::printf("\n\nLinf NORM OF MSM NODE ERROR = %lf\n\n", toterr);
    ASSERT_LT(toterr, TOL);
  }

  // Constructor for Driver test
  MSMDriverTest(std::shared_ptr<Portage::Simple_Mesh> s,
                std::shared_ptr<Portage::Simple_Mesh> t) :
    sourceMesh(s), targetMesh(t),
    sourceState(sourceMesh), targetState(targetMesh), targetState2(targetMesh),
    sourceMeshWrapper(*sourceMesh), targetMeshWrapper(*targetMesh),
    sourceStateWrapper(sourceState), targetStateWrapper(targetState),
    targetStateWrapper2(targetState2)
  {}

};

// Class which constructs a pair of simple 3-D meshes, target
// contained in source
struct MSMDriverTest3D : MSMDriverTest {
  MSMDriverTest3D(): MSMDriverTest(
      std::make_shared<Portage::Simple_Mesh>(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 10, 10, 10),
      std::make_shared<Portage::Simple_Mesh>(0.3, 0.3, 0.3, 0.7, 0.7, 0.7,  4,  4,  4)) {}
};

// Methods for computing initial field values

double compute_linear_field_2d(Portage::Point<2> centroid) {
  return centroid[0]+centroid[1];
}
double compute_quadratic_field_2d(Portage::Point<2> centroid) {
  return centroid[0]*centroid[0] + centroid[1]*centroid[1] +
      centroid[0]*centroid[1];
}

// Methods for computing initial field values
double compute_linear_field_3d(Portage::Point<3> centroid) {
  return centroid[0]+centroid[1]+centroid[2];
}
double compute_quadratic_field_3d(Portage::Point<3> centroid) {
  return centroid[0]*centroid[0] + centroid[1]*centroid[1] +
      centroid[2]*centroid[2] + centroid[0]*centroid[1] +
      centroid[1]*centroid[2] + centroid[2]*centroid[0];
}

// Test cases: these are constructed by calling TEST_F with the name
// of the test class you want to use.  The unit test method is then
// called inside each test with the appropriate template arguments.
// Google test will pick up each test and run it as part of the larger
// test fixture.  If any one of these fails the whole test_driver
// fails.

TEST_F(MSMDriverTest3D, 3D1stOrderLinear) {
  unitTest<Portage::IntersectR3D,
           Portage::Interpolate_1stOrder,
           Portage::SearchPointsByCells, 3>
    (&compute_linear_field_3d, .75, Portage::Meshfree::Basis::Linear);
}

TEST_F(MSMDriverTest3D, 3D2ndOrderQuadratic) {
  unitTest<Portage::IntersectR3D,
           Portage::Interpolate_2ndOrder,
           Portage::SearchPointsByCells, 3>
    (&compute_quadratic_field_3d, 1.5, Portage::Meshfree::Basis::Quadratic);
}

TEST_F(MSMDriverTest3D, 3D1stOrderLinearScatter) {
  unitTest<Portage::IntersectR3D,
           Portage::Interpolate_1stOrder,
           Portage::SearchPointsByCells, 3>
    (&compute_linear_field_3d, .75, Portage::Meshfree::Basis::Linear, 
     Portage::Meshfree::Scatter);
}

TEST_F(MSMDriverTest3D, 3D2ndOrderQuadraticScatter) {
  unitTest<Portage::IntersectR3D,
           Portage::Interpolate_2ndOrder,
           Portage::SearchPointsByCells, 3>
    (&compute_quadratic_field_3d, 1.5, Portage::Meshfree::Basis::Quadratic, 
     Portage::Meshfree::Scatter);
}

}  // namespace
