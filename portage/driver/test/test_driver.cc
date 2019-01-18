/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/


#include <iostream>
#include <memory>

#include "gtest/gtest.h"
#include "mpi.h"

#include "portage/driver/mmdriver.h"
#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"
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
  std::shared_ptr<Jali::State> sourceState;
  std::shared_ptr<Jali::State> targetState;

  //  Wrappers for interfacing with the underlying mesh data structures
  Wonton::Jali_Mesh_Wrapper sourceMeshWrapper;
  Wonton::Jali_Mesh_Wrapper targetMeshWrapper;
  Wonton::Jali_State_Wrapper sourceStateWrapper;
  Wonton::Jali_State_Wrapper targetStateWrapper;

  // This is the basic test method to be called for each unit test.
  //  It will work for 2-D and 3-D, coincident and non-coincident
  //  cell-centered remaps.
  template <
    template<Portage::Entity_kind, class, class, class,
    template<class, int, class, class> class,
    class, class> class Intersect,
    template<int, Portage::Entity_kind, class, class, class,
    template<class, int, class, class> class,
    class, class> class Interpolate,
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
    sourceState->add("celldata", sourceMesh, Jali::Entity_kind::CELL,
                     Jali::Entity_type::ALL, &(sourceData[0]));

    // Build the target state storage
    const int ntarcells = targetMeshWrapper.num_owned_cells();
    std::vector<double> targetData(ntarcells, 0.0);
    targetState->add("celldata", targetMesh, Jali::Entity_kind::CELL,
                     Jali::Entity_type::ALL, &(targetData[0]));

    //  Build the main driver data for this mesh type
    //  Register the variable name and interpolation order with the driver
    std::vector<std::string> remap_fields;
    remap_fields.push_back("celldata");

    Portage::MMDriver<Portage::SearchKDTree,
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

    Jali::UniStateVector<double, Jali::Mesh> cellvecout;
    bool found = targetState->get<double, Jali::Mesh>("celldata", targetMesh,
                                                      Jali::Entity_kind::CELL,
                                                      Jali::Entity_type::ALL,
                                                      &cellvecout);
    ASSERT_TRUE(found);

    double source_integral = 0.0;
    for (int c = 0; c < nsrccells; ++c) {
      double cellvol = sourceMesh->cell_volume(c);
      source_integral += sourceData[c]*cellvol;
    }

    
    double field_err2 = 0., target_integral = 0.;
    for (int c = 0; c < ntarcells; ++c) {
      JaliGeometry::Point ccen = targetMesh->cell_centroid(c);

      double error = compute_initial_field(ccen) - cellvecout[c];
      field_err2 += error*error;

      double cellvol = targetMesh->cell_volume(c);
      target_integral += cellvecout[c]*cellvol;
    }

    double field_err = sqrt(field_err2);
    double conservation_err = source_integral - target_integral;
    
    std::printf("\n\nL2 NORM OF ERROR = %lf\n\n", field_err);
    ASSERT_NEAR(expected_answer, field_err, TOL);

    std::printf("\n\nConservation Error = %18.13lf\n", conservation_err);
    ASSERT_NEAR(0.0, conservation_err, TOL);
  }
  // Constructor for Driver test
  DriverTest(std::shared_ptr<Jali::Mesh> s, std::shared_ptr<Jali::Mesh> t) :
      sourceMesh(s), targetMesh(t),
      sourceState(Jali::State::create(sourceMesh)),
      targetState(Jali::State::create(targetMesh)),
    sourceMeshWrapper(*sourceMesh), targetMeshWrapper(*targetMesh),
    sourceStateWrapper(*sourceState), targetStateWrapper(*targetState) {
  }

};

// Class which constructs a pair simple 2-D coincident meshes for remaps
struct DriverTest2D : DriverTest {
  DriverTest2D() : DriverTest(Jali::MeshFactory(MPI_COMM_WORLD)
  (0.0, 0.0, 1.0, 1.0, 5, 5),
  Jali::MeshFactory(MPI_COMM_WORLD) (0.0, 0.0, 1.0, 1.0, 4, 4)) {}
};

// Class which constructs a pair of simple 2-D non-coincident meshes for remaps
struct DriverTest2DNonCoincident : DriverTest {
  DriverTest2DNonCoincident() : DriverTest(Jali::MeshFactory(MPI_COMM_WORLD)
  (0.0, 0.0, 1.0, 1.0, 5, 5),
  Jali::MeshFactory(MPI_COMM_WORLD) (0.2, 0.2, 1.2, 1.2, 4, 4)) {}
};

// Class which constructs a pair simple 3-D coincident meshes for remaps
struct DriverTest3D : DriverTest {
  DriverTest3D(): DriverTest(Jali::MeshFactory(MPI_COMM_WORLD)
  (0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 5, 5, 5),
  Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 4, 4)) {}
};

// Class which constructs a pair simple 3-D non-coincident meshes for remaps
struct DriverTest3DNonCoincident : DriverTest {
  DriverTest3DNonCoincident(): DriverTest(Jali::MeshFactory(MPI_COMM_WORLD)
  (0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 5, 5, 5),
  Jali::MeshFactory(MPI_COMM_WORLD)(0.2, 0.2, 0.2, 1.2, 1.2, 1.2, 4, 4, 4)) {}
};

// Methods for computing initial field values

double compute_constant_field(JaliGeometry::Point centroid) {
  return 10000.0;
}

double compute_linear_field(JaliGeometry::Point centroid) {
  return 100*(centroid[0]+centroid[1]);
}
double compute_linear_field_3d(JaliGeometry::Point centroid) {
  return 100*(centroid[0]+centroid[1]+centroid[2]);
}
double compute_quadratic_field(JaliGeometry::Point centroid) {
  return 3*3*centroid[0]*centroid[0]+40*40*centroid[1]*centroid[1];
}
double compute_quadratic_field_3d(JaliGeometry::Point centroid) {
  return 3*3*centroid[0]*centroid[0] + 40*40*centroid[1]*centroid[1] +
      500*500*centroid[2]*centroid[2];
}

// Test cases: these are constructed by calling TEST_F with the name
// of the test class you want to use.  The unit test method is then
// called inside each test with the appropriate template arguments.
// Google test will pick up each test and run it as part of the larger
// test fixture.  If any one of these fails the whole test_driver
// fails.

// Cell-centered Remap on 2D Meshes with Coincident Boundaries
// Example 0
TEST_F(DriverTest2D, 2D_Constant_1stOrderCellCntr_Coincident) {
  unitTest<Portage::IntersectR2D, Portage::Interpolate_1stOrder, 2>
  (compute_constant_field, 0.0);
}
// Example 1
TEST_F(DriverTest2D, 2D_Linear_1stOrderCellCntr_Coincident) {
  unitTest<Portage::IntersectR2D, Portage::Interpolate_1stOrder, 2>
  (compute_linear_field, 6.3245553203367573);
}
// Example 2
TEST_F(DriverTest2D, 2D_Linear_2ndOrderCellCntr_Coincident) {
  unitTest<Portage::IntersectR2D, Portage::Interpolate_2ndOrder, 2>
  (compute_linear_field, 0.0);
}
// Example 3
TEST_F(DriverTest2D, 2D_Quadratic_1stOrderCellCntr_Coincident) {
  unitTest<Portage::IntersectR2D, Portage::Interpolate_1stOrder, 2>
  (compute_quadratic_field, 83.644493615838698);
}
// Example 4
TEST_F(DriverTest2D, 2D_Quadratic_2ndOrderCellCntr_Coincident) {
  unitTest<Portage::IntersectR2D, Portage::Interpolate_2ndOrder, 2>
  (compute_quadratic_field, 12.484585705981344);
}

// Cell-centered Remap on 2D Meshes with Non-coincident (mismatched)
// Boundaries

// Example 5 - Since the domains have the same volume but are merely
// offset the same constant should be reproduced on the target mesh
// and the error should be zero
TEST_F(DriverTest2DNonCoincident, 2D_Constant_1stOrderCellCntr_NonCoincident) {
  unitTest<Portage::IntersectR2D, Portage::Interpolate_1stOrder, 2>
  (compute_constant_field, 0.0);
}
// Example 6
TEST_F(DriverTest2DNonCoincident, 2D_Linear_1stOrderCellCntr_NonCoincident) {
  unitTest<Portage::IntersectR2D, Portage::Interpolate_1stOrder, 2>
  (compute_linear_field, 166.01204775557699);
}
// Example 7
TEST_F(DriverTest2DNonCoincident, 2D_Linear_2ndOrderCellCntr_NonCoincident) {
  unitTest<Portage::IntersectR2D, Portage::Interpolate_2ndOrder, 2>
  (compute_linear_field, 161.86414056238655);
}
// Example 8
TEST_F(DriverTest2DNonCoincident, 2D_Quadratic_1stOrderCellCntr_NonCoincident) {
  unitTest<Portage::IntersectR2D, Portage::Interpolate_1stOrder, 2>
  (compute_quadratic_field, 1864.3054259128905);
}
// Example 9
TEST_F(DriverTest2DNonCoincident, 2D_Quadratic_2ndOrderCellCntrNonCoincident) {
  unitTest<Portage::IntersectR2D, Portage::Interpolate_2ndOrder, 2>
  (compute_quadratic_field, 1724.7774423007857);
}

// Cell-centered Remap on 3D Meshes with Coincident Boundaries
// Example 10
TEST_F(DriverTest3D, 3D_Constant_1stOrderCellCntr_Coincident) {
  unitTest<Portage::IntersectR3D, Portage::Interpolate_1stOrder, 3>
  (compute_constant_field, 0.0);
}
// Example 11
TEST_F(DriverTest3D, 3D_Linear_1stOrderCellCntr_Coincident) {
  unitTest<Portage::IntersectR3D, Portage::Interpolate_1stOrder, 3>
  (compute_linear_field_3d, 15.491933384829508);
}
// Example 12
TEST_F(DriverTest3D, 3D_Linear_2ndOrderCellCntr_Coincident) {
  unitTest<Portage::IntersectR3D, Portage::Interpolate_2ndOrder, 3>
  (compute_linear_field_3d, 0.0);
}
// Example 13
TEST_F(DriverTest3D, 3D_Quadratic_1stOrderCellCntr_Coincident) {
  unitTest<Portage::IntersectR3D, Portage::Interpolate_1stOrder, 3>
  (compute_quadratic_field_3d, 26139.462467794965);
}
// Example 14
TEST_F(DriverTest3D, 3D_Quadratic_2ndOrderCellCntr_Coincident) {
  unitTest<Portage::IntersectR3D, Portage::Interpolate_2ndOrder, 3>
  (compute_quadratic_field_3d, 3904.3739523156646);
}

// Cell-centered Remap on 3D Meshes with Non-coincident (mismatched) Boundaries

// Example 15 - Since the domains have the same volume but are merely
// offset the same constant should be reproduced on the target mesh
// and the error should be zero
TEST_F(DriverTest3DNonCoincident, 3D_Constant_1stOrderCellCntr_NonCoincident) {
  unitTest<Portage::IntersectR3D, Portage::Interpolate_1stOrder, 3>
  (compute_constant_field, 0.0);
}
// Example 16
TEST_F(DriverTest3DNonCoincident, 3D_Linear_1stOrderCellCntr_NonCoincident) {
  unitTest<Portage::IntersectR3D, Portage::Interpolate_1stOrder, 3>
  (compute_linear_field_3d, 492.09755130461593);
}
// Example 17
TEST_F(DriverTest3DNonCoincident, 3D_Linear_2ndOrderCellCntr_NonCoincident) {
  unitTest<Portage::IntersectR3D, Portage::Interpolate_2ndOrder, 3>
  (compute_linear_field_3d, 483.73546489791454);
}
// Example 18
TEST_F(DriverTest3DNonCoincident, 3D_Quadratic_1stOrderCellCntr_NonCoincident) {
  unitTest<Portage::IntersectR3D, Portage::Interpolate_1stOrder, 3>
  (compute_quadratic_field_3d, 582831.37671617966);
}
// Example 19
TEST_F(DriverTest3DNonCoincident, 3D_Quadratic_2ndOrderCellCntr_NonCoincident) {
  unitTest<Portage::IntersectR3D, Portage::Interpolate_2ndOrder, 3>
  (compute_quadratic_field_3d, 538971.52673063509);
}
}
