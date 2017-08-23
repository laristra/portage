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

#include "gtest/gtest.h"
#include "mpi.h"

#include "portage/driver/driver.h"
#include "portage/driver/driver_mesh_swarm_mesh.h"
#include "portage/wonton/mesh/jali/jali_mesh_wrapper.h"
#include "portage/wonton/state/jali/jali_state_wrapper.h"
#include "portage/intersect/intersect_r2d.h"
#include "portage/intersect/intersect_r3d.h"
#include "portage/interpolate/interpolate_1st_order.h"
#include "portage/interpolate/interpolate_2nd_order.h"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliStateVector.h"
#include "JaliState.h"
#include "portage/search/search_points_by_cells.h"
#include "portage/accumulate/accumulate.h"
#include "portage/estimate/estimate.h"
// amh: TODO--change to simple mesh?
namespace{

double TOL = 1e-6;

//This is a set of integration tests based off of main.cc.
//There will be at least one test corresponding to each case found in main.cc.
//This is a test fixture and  must be derived from the ::testing::Test class.
//Specializations of this class, such as 2D/3D coincident and non-coincident remaps
//should be derived from this.
class DriverTest : public ::testing::Test {
 protected:
  //Source and target meshes
  std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::Mesh> targetMesh;
  //Source and target mesh state
  Jali::State sourceState;
  Jali::State targetState;
  Jali::State targetState2;
  // Wrappers for interfacing with the underlying mesh data structures
  Portage::Jali_Mesh_Wrapper sourceMeshWrapper;
  Portage::Jali_Mesh_Wrapper targetMeshWrapper;
  Portage::Jali_State_Wrapper sourceStateWrapper;
  Portage::Jali_State_Wrapper targetStateWrapper;
  Portage::Jali_State_Wrapper targetStateWrapper2;

  // This is the basic test method to be called for each unit test. It will work
  // for 2-D and 3-D, coincident and non-coincident cell-centered remaps.
  template <
    template<class, class> class Intersect,
    template<class, class, class, Portage::Entity_kind, long> class Interpolate,
    template <int, class, class> class SwarmSearch,
    int Dimension
  >
  void unitTest(double compute_initial_field(JaliGeometry::Point centroid),
                double smoothing_factor, Portage::Meshfree::Basis::Type basis)
  {

    // Fill the source state data with the specified profile
    const int nsrccells = sourceMeshWrapper.num_owned_cells() +
        sourceMeshWrapper.num_ghost_cells();
    std::vector<double> sourceData(nsrccells);

    //Create the source data for given function
    for (unsigned int c = 0; c < nsrccells; ++c) {
      JaliGeometry::Point cen = sourceMesh->cell_centroid(c);
      sourceData[c] = compute_initial_field(cen);
    }
    sourceState.add("celldata", sourceMesh, Jali::Entity_kind::CELL,
                    Jali::Entity_type::ALL, &(sourceData[0]));

    //Build the target state storage
    const int ntarcells = targetMeshWrapper.num_owned_cells();
    std::vector<double> targetData(ntarcells, 0.0);
    targetState.add("celldata", targetMesh, Jali::Entity_kind::CELL,
                    Jali::Entity_type::ALL, &(targetData[0]));
    targetState2.add("celldata", targetMesh, Jali::Entity_kind::CELL,
                    Jali::Entity_type::ALL, &(targetData[0]));

    // Register the variable name and interpolation order with the driver
    std::vector<std::string> remap_fields;
    remap_fields.push_back("celldata");

    // Build the mesh-mesh driver data for this mesh type
    Portage::Driver<Portage::SearchKDTree,
    Intersect,
    Interpolate,
    Dimension,
    Portage::Jali_Mesh_Wrapper, Portage::Jali_State_Wrapper,
    Portage::Jali_Mesh_Wrapper, Portage::Jali_State_Wrapper>
    mmdriver(sourceMeshWrapper, sourceStateWrapper,
             targetMeshWrapper, targetStateWrapper);
    mmdriver.set_remap_var_names(remap_fields);
    //run on one processor
    mmdriver.run(false);

    // Build the mesh-swarm-mesh driver data for this mesh type
    Portage::MSM_Driver<
      SwarmSearch,
      Portage::Meshfree::Accumulate,
      Portage::Meshfree::Estimate,
      Dimension,
      Portage::Jali_Mesh_Wrapper, Portage::Jali_State_Wrapper,
      Portage::Jali_Mesh_Wrapper, Portage::Jali_State_Wrapper
    >
    msmdriver(sourceMeshWrapper, sourceStateWrapper,
              targetMeshWrapper, targetStateWrapper2,
	      smoothing_factor, basis);
    msmdriver.set_remap_var_names(remap_fields);
    //run on one processor
    msmdriver.run(false);

    //Check the answer
    Portage::Point<Dimension> nodexy;
    const int ntarnodes = targetMeshWrapper.num_owned_nodes();
    double stdval, err;
    double toterr=0.;

    Jali::StateVector<double, Jali::Mesh> cellvecout;
    bool found = targetState.get<double, Jali::Mesh>("celldata", targetMesh,
                                                     Jali::Entity_kind::CELL,
                                                     Jali::Entity_type::ALL,
                                                     &cellvecout);
    Jali::StateVector<double, Jali::Mesh> cellvecout2;
    bool found2 = targetState2.get<double, Jali::Mesh>("celldata", targetMesh,
                                                     Jali::Entity_kind::CELL,
                                                     Jali::Entity_type::ALL,
                                                     &cellvecout2);
    ASSERT_TRUE(found);
    ASSERT_TRUE(found2);

    for (int c = 0; c < ntarcells; ++c) {
      JaliGeometry::Point ccen = targetMesh->cell_centroid(c);
      double value = compute_initial_field(ccen);
      double error;
      error = cellvecout[c] - cellvecout2[c];
      // dump diagnostics for each cell
      std::printf("Cell=% 4d Centroid=(% 5.3lf,% 5.3lf)", c,
                  ccen[0], ccen[1]);
      std::printf(" Val=% 10.6lf MM=% 10.6lf MSM=% 10.6lf Err=% lf\n",
                  value, cellvecout[c], cellvecout2[c], error);
      toterr += std::max(toterr, std::fabs(error));
    }

    std::printf("\n\nLinf NORM OF ERROR = %lf\n\n", toterr);
    ASSERT_LT(toterr, TOL);
  }

  //Constructor for Driver test
  DriverTest(std::shared_ptr<Jali::Mesh> s,std::shared_ptr<Jali::Mesh> t) :
    sourceMesh(s), targetMesh(t), 
    sourceState(sourceMesh), targetState(targetMesh), targetState2(targetMesh),
    sourceMeshWrapper(*sourceMesh), targetMeshWrapper(*targetMesh),
    sourceStateWrapper(sourceState), targetStateWrapper(targetState), 
    targetStateWrapper2(targetState2)
  {}

};


//Class which constructs a pair of simple 2-D meshes, target contained in source
struct DriverTest2D : DriverTest {
  DriverTest2D() : DriverTest(
  Jali::MeshFactory(MPI_COMM_WORLD) (0.0, 0.0, 1.0, 1.0, 11, 11),
  Jali::MeshFactory(MPI_COMM_WORLD) (0.3, 0.3, 0.7, 0.7,  3,  3)) {}
};

//Class which constructs a pair of simple 3-D meshes, target contained in source
struct DriverTest3D : DriverTest {
  DriverTest3D(): DriverTest(
  Jali::MeshFactory(MPI_COMM_WORLD) (0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 11, 11, 11),
  Jali::MeshFactory(MPI_COMM_WORLD) (0.3, 0.3, 0.3, 0.7, 0.7, 0.7,  3,  3,  3)) {}
};

//Methods for computing initial field values
double compute_linear_field(JaliGeometry::Point centroid){
  return centroid[0]+centroid[1];
}
double compute_quadratic_field(JaliGeometry::Point centroid){
  return centroid[0]*centroid[0]+centroid[1]*centroid[1]+centroid[0]*centroid[1];
}

//Methods for computing initial field values
double compute_linear_field_3d(JaliGeometry::Point centroid){
  return centroid[0]+centroid[1]+centroid[2];
}
double compute_quadratic_field_3d(JaliGeometry::Point centroid){
  return centroid[0]*centroid[0]+centroid[1]*centroid[1]+centroid[2]*centroid[2]+
         centroid[0]*centroid[1]+centroid[1]*centroid[2]+centroid[2]*centroid[0];
}

//Test cases:  these are constructed by calling TEST_F with the name of the test
//class you want to use.  The unit test method is then called inside each test
//with the appropriate template arguments.  Google test will pick up each test
//and run it as part of the larger test fixture.  If any one of these fails the
//whole test_driver fails.

// Not working yet.

// TEST_F(DriverTest2D, 2D1stOrderLinear){
//   unitTest<Portage::IntersectR2D, 
// 	   Portage::Interpolate_1stOrder, 
//            Portage::SearchPointsByCells, 2>
//     (&compute_linear_field, 1.0, Portage::Meshfree::Basis::Linear);
// }
// TEST_F(DriverTest2D, 2D1stOrderQuadratic){
//   unitTest<Portage::IntersectR2D, 
// 	   Portage::Interpolate_1stOrder, 
//            Portage::SearchPointsByCells, 2>
//     (&compute_quadratic_field, 1.0, Portage::Meshfree::Basis::Linear);
// }

// TEST_F(DriverTest2D, 2D2ndOrderLinear){
//   unitTest<Portage::IntersectR2D, 
// 	   Portage::Interpolate_2ndOrder, 
//            Portage::SearchPointsByCells, 2>
//     (&compute_linear_field, 1.0, Portage::Meshfree::Basis::Quadratic);
// }
// TEST_F(DriverTest2D, 2D2ndOrderQuadratic){
//   unitTest<Portage::IntersectR2D, 
// 	   Portage::Interpolate_2ndOrder, 
//            Portage::SearchPointsByCells, 2>
//     (&compute_quadratic_field, 1.0, Portage::Meshfree::Basis::Quadratic);
// }

}
