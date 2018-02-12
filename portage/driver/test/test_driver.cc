/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
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
// amh: TODO--change to simple mesh?
namespace {

double TOL = 1e-6;

// This is a set of integration tests based off of main.cc.  There
// will be at least one test corresponding to each case found in
// main.cc.  This is a test fixture and must be derived from the
// ::testing::Test class.  Specializations of this class, such as
// 2D/3D coincident and non-coincident remaps should be derived from
// this.
class DriverTest : public ::testing::Test {
 protected:
  // Source and target meshes
  std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::Mesh> targetMesh;
  // Source and target mesh state
  Jali::State sourceState;
  Jali::State targetState;
  //  Wrappers for interfacing with the underlying mesh data structures
  Wonton::Jali_Mesh_Wrapper sourceMeshWrapper;
  Wonton::Jali_Mesh_Wrapper targetMeshWrapper;
  Wonton::Jali_State_Wrapper sourceStateWrapper;
  Wonton::Jali_State_Wrapper targetStateWrapper;

  // This is the basic test method to be called for each unit test.
  //  It will work for 2-D and 3-D, coincident and non-coincident
  //  cell-centered remaps.
  template <
    template<Portage::Entity_kind, class, class> class Intersect,
    template<int, Portage::Entity_kind, class, class, class> class Interpolate,
    int Dimension
  >
  void unitTest(double compute_initial_field(JaliGeometry::Point centroid),
                double expected_answer) {

    //  Fill the source state data with the specified profile
    const int nsrccells = sourceMeshWrapper.num_owned_cells() +
        sourceMeshWrapper.num_ghost_cells();
    std::vector<double> sourceData(nsrccells);

    // Create the source data for given function
    for (unsigned int c = 0; c < nsrccells; ++c) {
      JaliGeometry::Point cen = sourceMesh->cell_centroid(c);
      sourceData[c] = compute_initial_field(cen);
    }
    sourceState.add("celldata", sourceMesh, Jali::Entity_kind::CELL,
                    Jali::Entity_type::ALL, &(sourceData[0]));

    // Build the target state storage
    const int ntarcells = targetMeshWrapper.num_owned_cells();
    std::vector<double> targetData(ntarcells, 0.0);
    targetState.add("celldata", targetMesh, Jali::Entity_kind::CELL,
                    Jali::Entity_type::ALL, &(targetData[0]));

    //  Build the main driver data for this mesh type
    //  Register the variable name and interpolation order with the driver
    std::vector<std::string> remap_fields;
    remap_fields.push_back("celldata");

    Portage::Driver<Portage::SearchKDTree,
    Intersect,
    Interpolate, Dimension,
    Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
    Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper>
    d(sourceMeshWrapper, sourceStateWrapper, targetMeshWrapper,
      targetStateWrapper);
    d.set_remap_var_names(remap_fields);
    // run on one processor
    d.run(false);

    // Check the answer
    Portage::Point<Dimension> nodexy;
    const int ntarnodes = targetMeshWrapper.num_owned_nodes();
    double stdval, err;
    double toterr = 0.;

    Jali::StateVector<double, Jali::Mesh> cellvecout;
    bool found = targetState.get<double, Jali::Mesh>("celldata", targetMesh,
                                                     Jali::Entity_kind::CELL,
                                                     Jali::Entity_type::ALL,
                                                     &cellvecout);
    ASSERT_TRUE(found);

    for (int c = 0; c < ntarcells; ++c) {
      JaliGeometry::Point ccen = targetMesh->cell_centroid(c);
      double error;
      error = compute_initial_field(ccen) - cellvecout[c];
      //  dump diagnostics for each cell
      std::printf("Cell=% 4d Centroid = (% 5.3lf,% 5.3lf)", c,
                  ccen[0], ccen[1]);
      std::printf("  Value = % 10.6lf  Err = % lf\n",
                  cellvecout[c], error);
      toterr += error*error;
    }

    // amh: FIXME!!  Compare individual, per-node  values/ error norms here
    std::printf("\n\nL2 NORM OF ERROR = %lf\n\n", sqrt(toterr));
    ASSERT_NEAR(expected_answer, sqrt(toterr), TOL);
  }
  // Constructor for Driver test
  DriverTest(std::shared_ptr<Jali::Mesh> s, std::shared_ptr<Jali::Mesh> t) :
    sourceMesh(s), targetMesh(t), sourceState(sourceMesh),
    targetState(targetMesh),
    sourceMeshWrapper(*sourceMesh), targetMeshWrapper(*targetMesh),
    sourceStateWrapper(sourceState), targetStateWrapper(targetState){
  }

};

// Class which constructs a pair simple 2-D coincident meshes for remaps
struct DriverTest2D : DriverTest {
  DriverTest2D() : DriverTest(Jali::MeshFactory(MPI_COMM_WORLD)
  (0.0, 0.0, 1.0, 1.0, 11, 11),
  Jali::MeshFactory(MPI_COMM_WORLD) (0.0, 0.0, 1.0, 1.0, 3, 3)) {}
};

// Class which constructs a pair of simple 2-D non-coincident meshes for remaps
struct DriverTest2DNonCoincident : DriverTest {
  DriverTest2DNonCoincident() : DriverTest(Jali::MeshFactory(MPI_COMM_WORLD)
  (0.0, 0.0, 1.0, 1.0, 11, 11),
  Jali::MeshFactory(MPI_COMM_WORLD) (0.0, 0.0, 1.0+1.5*(1/3.), 1.0, 3, 3)) {}
};

// Class which constructs a pair simple 3-D coincident meshes for remaps
struct DriverTest3D : DriverTest {
  DriverTest3D(): DriverTest(Jali::MeshFactory(MPI_COMM_WORLD)
  (0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 11, 11, 11),
  Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 3, 3, 3)) {}
};

// Class which constructs a pair simple 3-D non-coincident meshes for remaps
struct DriverTest3DNonCoincident : DriverTest {
  DriverTest3DNonCoincident(): DriverTest(Jali::MeshFactory(MPI_COMM_WORLD)
  (0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 11, 11, 11),
  Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 0.0, 1.0+1.5*1/3.,
      1.0+1.5*1/3., 1.0+1.5*1/3., 3, 3, 3)) {}
};

// Methods for computing initial field values
double compute_linear_field(JaliGeometry::Point centroid) {
  return centroid[0]+centroid[1];
}
double compute_quadratic_field(JaliGeometry::Point centroid) {
  return centroid[0]*centroid[0]+centroid[1]*centroid[1];
}
double compute_quadratic_field_3d(JaliGeometry::Point centroid) {
  return centroid[0]*centroid[0] + centroid[1]*centroid[1] +
      centroid[2]*centroid[2];
}

// Test cases: these are constructed by calling TEST_F with the name
// of the test class you want to use.  The unit test method is then
// called inside each test with the appropriate template arguments.
// Google test will pick up each test and run it as part of the larger
// test fixture.  If any one of these fails the whole test_driver
// fails.

// Example 0
TEST_F(DriverTest2D, 2D_1stOrderLinearCellCntrCoincident1proc) {
  unitTest<Portage::IntersectR2D, Portage::Interpolate_1stOrder, 2>
  (compute_linear_field, 0.0095429796560267122);
}
// Example 1
TEST_F(DriverTest2D, 2D_2ndOrderLinearCellCntrCoincident1proc) {
  unitTest<Portage::IntersectR2D, Portage::Interpolate_2ndOrder, 2>
  (compute_linear_field, 0.);
}
// Example 2
TEST_F(DriverTest2DNonCoincident, 2D_1stOrderLinearCellCntrNonCoincident1proc) {
  unitTest<Portage::IntersectR2D, Portage::Interpolate_1stOrder, 2>
  (compute_linear_field, 3.067536);
}
// Example 3
TEST_F(DriverTest2DNonCoincident, 2D_2ndOrderLinearCellCntrNonCoincident1proc) {
  unitTest<Portage::IntersectR2D, Portage::Interpolate_2ndOrder, 2>
  (compute_linear_field, 3.067527);
}
// Example 4
TEST_F(DriverTest2D, 2D_1stOrderQuadraticCellCntrCoincident1proc) {
  unitTest<Portage::IntersectR2D, Portage::Interpolate_1stOrder, 2>
  (compute_quadratic_field, 0.052627);
}
// Example 5
TEST_F(DriverTest2D, 2D_2ndOrderQuadraticCellCntrCoincident1proc) {
  unitTest<Portage::IntersectR2D, Portage::Interpolate_2ndOrder, 2>
  (compute_quadratic_field, 0.051424);
}
// Example 6
TEST_F(DriverTest2DNonCoincident, 2D_1stOrderQuadraticCellCntrNonCoincident1proc) {
  unitTest<Portage::IntersectR2D, Portage::Interpolate_1stOrder, 2>
  (compute_quadratic_field, 3.303476);
}
// Example 7
TEST_F(DriverTest2DNonCoincident, 2D_2ndOrderQuadraticCellCntrNonCoincident1proc) {
  unitTest<Portage::IntersectR2D, Portage::Interpolate_2ndOrder, 2>
  (compute_quadratic_field, 3.303466);
}
// Example 8
TEST_F(DriverTest3D, 3D_1stOrderQuadraticCellCntrCoincident1proc) {
  unitTest<Portage::IntersectR3D, Portage::Interpolate_1stOrder, 3>
  (compute_quadratic_field_3d, .135694);
}
// Example 9
TEST_F(DriverTest3D, 3D_2ndOrderQuadraticCellCntrCoincident1proc) {
  unitTest<Portage::IntersectR3D, Portage::Interpolate_2ndOrder, 3>
  (compute_quadratic_field_3d, .133602);
}
// Example 10
TEST_F(DriverTest3DNonCoincident, 3D_1stOrderQuadraticCellCntrNonCoincident1proc) {
  unitTest<Portage::IntersectR3D, Portage::Interpolate_1stOrder, 3>
  (compute_quadratic_field_3d, 12.336827);
}
// Example 11
TEST_F(DriverTest3DNonCoincident, 3D_2ndOrderQuadraticCellCntrNonCoincident1proc) {
  unitTest<Portage::IntersectR3D, Portage::Interpolate_2ndOrder, 3>
  (compute_quadratic_field_3d, 12.336822);
}
}
