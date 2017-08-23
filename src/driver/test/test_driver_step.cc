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

#include "portage/driver/driver.h"
#include "portage/wonton/mesh/jali/jali_mesh_wrapper.h"
#include "portage/wonton/state/jali/jali_state_wrapper.h"
#include "portage/intersect/intersect_r2d.h"
#include "portage/interpolate/interpolate_1st_order.h"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliStateVector.h"
#include "JaliState.h"
//amh: TODO--change to simple mesh?
namespace{

double TOL = 1e-6;

//This is a set of integration tests based off of main.cc.    There will be at
//least one test corresponding to each case found in main.cc.  This is a test
//fixture and  must be derived from the ::testing::Test class.  Specializations
//of this class, such as 2D/3D coincident and non-coincident remaps should be
//derived from this.
class DriverTest : public ::testing::Test {
 protected:
  //Source and target meshes
  std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::Mesh> targetMesh;
  //Source and target mesh state
  Jali::State sourceState;
  Jali::State targetState;
  // Wrappers for interfacing with the underlying mesh data structures
  Portage::Jali_Mesh_Wrapper sourceMeshWrapper;
  Portage::Jali_Mesh_Wrapper targetMeshWrapper;
  Portage::Jali_State_Wrapper sourceStateWrapper;
  Portage::Jali_State_Wrapper targetStateWrapper;

  //This is the basic test method to be called for each unit test.  It will
  //work for 2-D and 3-D, coincident and non-coincident cell-centered remaps.
  template <template<class, class>class Intersect,
  template<class, class, class, Portage::Entity_kind, long>
  class Interpolate, int Dimension>
  void unitTest(double compute_initial_field(JaliGeometry::Point centroid),
                double expected_answer,
                Portage::LimiterType limiter=Portage::NOLIMITER)
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

    // Build the main driver data for this mesh type
    // Register the variable name and interpolation order with the driver
    std::vector<std::string> remap_fields;
    remap_fields.push_back("celldata");

    Portage::Driver<Portage::SearchKDTree,
    Intersect,
    Interpolate, Dimension,
    Portage::Jali_Mesh_Wrapper, Portage::Jali_State_Wrapper,
    Portage::Jali_Mesh_Wrapper, Portage::Jali_State_Wrapper>
    d(sourceMeshWrapper, sourceStateWrapper,targetMeshWrapper,
      targetStateWrapper);
    d.set_remap_var_names(remap_fields, limiter);
    //run on one processor
    d.run(false);

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
    ASSERT_TRUE(found);

    std::printf("DATA_BEGIN\n");
    for (int c = 0; c < ntarcells; ++c) {
      JaliGeometry::Point ccen = targetMesh->cell_centroid(c);
      double error;
      error = compute_initial_field(ccen) - cellvecout[c];
      // dump diagnostics for each cell
      std::printf("Cell=% 4d Centroid = (% 5.3lf,% 5.3lf)", c,
                  ccen[0], ccen[1]);
      std::printf("  Value = % 10.6lf  Err = % lf\n",
                  cellvecout[c], error);
      toterr += error*error;
    }
    std::printf("DATA_END\n");

    //amh: FIXME!!  Compare individual, per-node  values/ error norms here
    std::printf("\n\nL2 NORM OF ERROR = %lf\n\n", sqrt(toterr));
    ASSERT_NEAR(expected_answer, sqrt(toterr), TOL);
  }

  //Constructor for Driver test
  DriverTest(std::shared_ptr<Jali::Mesh> s,std::shared_ptr<Jali::Mesh> t) :
    sourceMesh(s), targetMesh(t), sourceState(sourceMesh), targetState(targetMesh),
    sourceMeshWrapper(*sourceMesh), targetMeshWrapper(*targetMesh),
    sourceStateWrapper(sourceState), targetStateWrapper(targetState)
  {}

};

//Class which constructs a pair simple 2-D coincident meshes for remaps
struct DriverTest2D : DriverTest {
  DriverTest2D() : DriverTest(Jali::MeshFactory(MPI_COMM_WORLD)
  (0.0, 0.0, 1.0, 1.0, 11, 11),Jali::MeshFactory(MPI_COMM_WORLD)
  (0.0, 0.0, 1.0, 1.0, 11, 11)) {}
};

//Class which constructs a pair of simple 2-D non-coincident meshes for remaps
struct DriverTest2DNonCoincident : DriverTest {
  DriverTest2DNonCoincident() : DriverTest(Jali::MeshFactory(MPI_COMM_WORLD)
  (0.0, 0.0, 1.0, 1.0, 11, 11),Jali::MeshFactory(MPI_COMM_WORLD)
  (0.0, 0.0, 1.0, 1.0, 18, 18)) {}
};

//Class which constructs a pair simple 3-D coincident meshes for remaps
struct DriverTest3D : DriverTest {
  DriverTest3D(): DriverTest(Jali::MeshFactory(MPI_COMM_WORLD)
  (0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 11, 11, 11), Jali::MeshFactory(MPI_COMM_WORLD)
  (0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 11, 11, 11)) {}
};

//Class which constructs a pair simple 3-D non-coincident meshes for remaps
struct DriverTest3DNonCoincident : DriverTest {
  DriverTest3DNonCoincident(): DriverTest(Jali::MeshFactory(MPI_COMM_WORLD)
  (0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 11, 11, 11), Jali::MeshFactory(MPI_COMM_WORLD)
  (0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 18, 18, 18)) {}
};

double compute_step_field_2d(JaliGeometry::Point centroid) {
  const double norm = 0.8192319205190405;
  const double normal[2] = {1.0*norm, .7*norm};
  const double offset = .5*normal[0]+.5*normal[1];

  assert(centroid.dim()==2);
  double value=-1.0;
  if (normal[0]*centroid[0] + normal[1]*centroid[1] >= offset) value = 1.0;
  return value;
}

double compute_step_field_3d(JaliGeometry::Point centroid) {
  const double norm = 0.5607721540920443;
  const double normal[3] = {1.0*norm, .7*norm, 1.3*norm};
  const double offset = .5*normal[0]+.5*normal[1]+.5*normal[2];

  assert(centroid.dim()==3);
  double value=-1.0;
  if (normal[0]*centroid[0] + normal[1]*centroid[1] +
      normal[2]*centroid[2] >= offset) value = 1.0;
  return value;
}

//Test cases:  these are constructed by calling TEST_F with the name of the
//test class you want to use.  The unit test method is then called inside each
//test with the appropriate template arguments.  Google test will pick up each
//test and run it as part of the larger test fixture.  If any one of these
//fails the whole test_driver fails.

TEST_F(DriverTest2D, 2D_1stOrderStepCellCntrCoincident1proc){
  unitTest<Portage::IntersectR2D, Portage::Interpolate_1stOrder, 2>
  (compute_step_field_2d, 0.0);
}
TEST_F(DriverTest2D, 2D_2ndOrderStepCellCntrCoincident1proc){
  unitTest<Portage::IntersectR2D, Portage::Interpolate_2ndOrder, 2>
  (compute_step_field_2d, 0.0);
}
TEST_F(DriverTest2D, 2D_2ndOrderStepCellCntrCoincident1procBJ){
  unitTest<Portage::IntersectR2D, Portage::Interpolate_2ndOrder, 2>
  (compute_step_field_2d, 0.0, Portage::BARTH_JESPERSEN);
}
TEST_F(DriverTest2DNonCoincident, 2D_1stOrderStepCellCntrNonCoincident1proc){
  unitTest<Portage::IntersectR2D, Portage::Interpolate_1stOrder, 2>
  (compute_step_field_2d, 5.022128);
}
TEST_F(DriverTest2DNonCoincident, 2D_2ndOrderStepCellCntrNonCoincident1proc){
  unitTest<Portage::IntersectR2D, Portage::Interpolate_2ndOrder, 2>
  (compute_step_field_2d, 4.779890 );
}
TEST_F(DriverTest2DNonCoincident, 2D_2ndOrderStepCellCntrNonCoincident1procBJ){
  unitTest<Portage::IntersectR2D, Portage::Interpolate_2ndOrder, 2>
  (compute_step_field_2d, 5.022128, Portage::BARTH_JESPERSEN);
}
TEST_F(DriverTest3D, 3D_1stOrderStepCellCntrCoincident1proc){
  unitTest<Portage::IntersectR3D, Portage::Interpolate_1stOrder, 3>
  (compute_step_field_3d, 0.0);
}
TEST_F(DriverTest3D, 3D_2ndOrderStepCellCntrCoincident1proc){
  unitTest<Portage::IntersectR3D, Portage::Interpolate_2ndOrder, 3>
  (compute_step_field_3d, 0.0);
}
TEST_F(DriverTest3D, 3D_2ndOrderStepCellCntrCoincident1procBJ){
  unitTest<Portage::IntersectR3D, Portage::Interpolate_2ndOrder, 3>
  (compute_step_field_3d, 0.0, Portage::BARTH_JESPERSEN);
}
TEST_F(DriverTest3DNonCoincident, 3D_1stOrderStepCellCntrNonCoincident1proc){
  unitTest<Portage::IntersectR3D, Portage::Interpolate_1stOrder, 3>
  (compute_step_field_3d,  18.626843);
}
TEST_F(DriverTest3DNonCoincident, 3D_2ndOrderStepCellCntrNonCoincident1proc){
  unitTest<Portage::IntersectR3D, Portage::Interpolate_2ndOrder, 3>
  (compute_step_field_3d,  18.238559);
}
TEST_F(DriverTest3DNonCoincident, 3D_2ndOrderStepCellCntrNonCoincident1procBJ){
  unitTest<Portage::IntersectR3D, Portage::Interpolate_2ndOrder, 3>
  (compute_step_field_3d,  18.626843, Portage::BARTH_JESPERSEN);
}
}
