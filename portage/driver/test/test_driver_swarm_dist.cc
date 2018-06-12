/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/


#include <iostream>
#include <memory>
#include <vector>

#include "portage/accumulate/accumulate.h"
#include "portage/distributed/mpi_particle_distribute.h"
#include "portage/driver/driver_swarm.h"
#include "portage/estimate/estimate.h"
#include "portage/support/Point.h"
#include "portage/support/portage.h"
#include "portage/swarm/swarm.h"
#include "portage/swarm/swarm_state.h"
#include "portage/wonton/mesh/flat/flat_mesh_wrapper.h"
#include "portage/wonton/mesh/jali/jali_mesh_wrapper.h"
#include "portage/search/search_points_by_cells.h"

#include "gtest/gtest.h"
#include "mpi.h"
#include "Mesh.hh"
#include "MeshFactory.hh"

namespace {

using std::shared_ptr;
using std::make_shared;
using Portage::Meshfree::SwarmFactory;


double TOL = 1e-6;

// This is a set of integration tests for the swarm-swarm remap driver.
// There will be at least one test corresponding to each case found in
// main.cc. This is a test fixture and must be derived from the
// ::testing::Test class. Specializations of this class, such as 2D/3D
// coincident and non-coincident remaps should be derived from this.

template<size_t dim>
class DriverTest : public ::testing::Test {
 protected:
  
  using SmoothingLengthPtr = shared_ptr<Portage::vector<std::vector<std::vector<double>>>>;

  // Source and target swarms
  shared_ptr<Portage::Meshfree::Swarm<dim>> sourceSwarm;
  shared_ptr<Portage::Meshfree::Swarm<dim>> targetSwarm;

  // Source and target mesh state
  shared_ptr<Portage::Meshfree::SwarmState<dim>> sourceState;
  shared_ptr<Portage::Meshfree::SwarmState<dim>> targetState;

  SmoothingLengthPtr smoothing_lengths_;

  Portage::Meshfree::WeightCenter center_;

  // Constructor for Driver test
  DriverTest(shared_ptr<Jali::Mesh> source_mesh,
             shared_ptr<Jali::Mesh> target_mesh) :
    smoothing_lengths_(nullptr),
    center_(Portage::Meshfree::Gather) {

    Wonton::Jali_Mesh_Wrapper jali_smesh_wrapper(*source_mesh);
    Wonton::Jali_Mesh_Wrapper jali_tmesh_wrapper(*target_mesh);

    //Create flat mesh wrappers for source/target jali meshes
    Wonton::Flat_Mesh_Wrapper<> source_mesh_flat;
    source_mesh_flat.initialize(jali_smesh_wrapper);

    Wonton::Flat_Mesh_Wrapper<> target_mesh_flat;
    target_mesh_flat.initialize(jali_tmesh_wrapper);

    // Source and target swarms
    sourceSwarm = make_shared<Portage::Meshfree::Swarm<dim>>(source_mesh_flat, Portage::Entity_kind::CELL);
    targetSwarm = make_shared<Portage::Meshfree::Swarm<dim>>(target_mesh_flat, Portage::Entity_kind::CELL);

    sourceState = make_shared<Portage::Meshfree::SwarmState<dim>>(*sourceSwarm);
    targetState = make_shared<Portage::Meshfree::SwarmState<dim>>(*targetSwarm);
  }

  void set_smoothing_lengths(shared_ptr<Portage::vector<std::vector<std::vector<double>>>>
                             smoothing_lengths,
			     Portage::Meshfree::WeightCenter center=Portage::Meshfree::Gather) {
    smoothing_lengths_ = smoothing_lengths;
    center_ = center;
  }

  // This is the basic test method to be called for each unit test.
  // It will work for 1, 2-D and 3-D swarms
  //
  template <template<int, class, class> class Search,
            Portage::Meshfree::Basis::Type basis>
  void unitTest(double compute_initial_field(Portage::Point<dim> coord),
                double expected_answer) {

    // Fill the source state data with the specified profile
    const int nsrcpts = sourceSwarm->num_owned_particles();
    typename Portage::Meshfree::SwarmState<dim>::DblVecPtr sourceData = 
        make_shared<typename Portage::Meshfree::SwarmState<dim>::DblVec>(nsrcpts);

    // Create the source data for given function
    for (unsigned int p = 0; p < nsrcpts; ++p) {
      Portage::Point<dim> coord =
          sourceSwarm->get_particle_coordinates(p);
      (*sourceData)[p] = compute_initial_field(coord);
    }
    sourceState->add_field("particledata", sourceData);

    // Build the target state storage
    const int ntarpts = targetSwarm->num_owned_particles();
    typename Portage::Meshfree::SwarmState<dim>::DblVecPtr targetData = 
        make_shared<typename Portage::Meshfree::SwarmState<dim>::DblVec>(ntarpts, 0.0);
    targetState->add_field("particledata", targetData);

    // Build the main driver data for this mesh type
    // Register the variable name and interpolation order with the driver
    std::vector<std::string> remap_fields;
    remap_fields.push_back("particledata");

    Portage::Meshfree::SwarmDriver<Search,
                                   Portage::Meshfree::Accumulate,
                                   Portage::Meshfree::Estimate,
                                   dim,
                                   Portage::Meshfree::Swarm<dim>,
                                   Portage::Meshfree::SwarmState<dim>,
                                   Portage::Meshfree::Swarm<dim>,
                                   Portage::Meshfree::SwarmState<dim>>
        d(*sourceSwarm, *sourceState, *targetSwarm, *targetState, *smoothing_lengths_,
	  Portage::Meshfree::Weight::B4, Portage::Meshfree::Weight::ELLIPTIC, center_);
    d.set_remap_var_names(remap_fields, remap_fields,
                          Portage::Meshfree::LocalRegression,
                          basis);
    // run on multiple processors
    d.run(true);

    // Check the answer
    Portage::Point<dim> nodexy;
    double stdval, err;
    double toterr=0.;

    typename Portage::Meshfree::SwarmState<dim>::DblVecPtr vecout;
    targetState->get_field("particledata", vecout);
    ASSERT_NE(nullptr, vecout);

    for (int p = 0; p < ntarpts; ++p) {
      Portage::Point<dim> coord = targetSwarm->get_particle_coordinates(p);
      double error;
      error = compute_initial_field(coord) - (*vecout)[p];
      // dump diagnostics for each particle
      if (dim == 1)
        std::printf("Particle=% 4d Coord = (% 5.3lf)", p, coord[0]);
      else if (dim == 2)
        std::printf("Particle=% 4d Coord = (% 5.3lf,% 5.3lf)", p, coord[0],
                    coord[1]);
      else if (dim == 3)
        std::printf("Particle=% 4d Coord = (% 5.3lf,% 5.3lf,% 5.3lf)", p,
                    coord[0], coord[1], coord[2]);
      double dummy= (*vecout)[p];
      std::printf("  Value = % 10.6lf  Err = % lf\n", dummy, error);
      toterr += error*error;
    }
    
    std::printf("\n\nL2 NORM OF ERROR = %lf\n\n", sqrt(toterr));
    ASSERT_NEAR(expected_answer, sqrt(toterr), TOL);
  }

};

// Class which constructs a pair of 2-D swarms (from jali) for remaps
struct DriverTest2D : DriverTest<2> {
  DriverTest2D() : DriverTest(Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 1.0, 1.0, 16, 16),
                              Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 1.0, 1.0, 8, 8)) 
  {
    size_t ntarpts = targetSwarm->num_particles(Portage::Entity_type::PARALLEL_OWNED); 
    auto smoothing_lengths = make_shared<Portage::vector<std::vector<std::vector<double>>>>(ntarpts,
                   std::vector<std::vector<double>>(1, std::vector<double>(2, 2.0/4)));
    DriverTest::set_smoothing_lengths(smoothing_lengths);
  }
};

struct DriverTest2DScatter : DriverTest<2> {
  DriverTest2DScatter() : DriverTest(Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 1.0, 1.0, 8, 8),
                              Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 1.0, 1.0, 4, 4)) 
  {
    size_t nsrcpts = sourceSwarm->num_particles(Portage::Entity_type::PARALLEL_OWNED); 
    auto smoothing_lengths = make_shared<Portage::vector<std::vector<std::vector<double>>>>(nsrcpts,
                   std::vector<std::vector<double>>(1, std::vector<double>(2, 2.0/8)));
    DriverTest::set_smoothing_lengths(smoothing_lengths, Portage::Meshfree::Scatter);
  }
};

// Class which constructs a pair of 3-D swarms (from jali) for remaps
struct DriverTest3D : DriverTest<3> {
  DriverTest3D(): DriverTest(Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                          16, 16, 16),
                             Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                          8, 8, 8)) 
  {
    size_t ntarpts = targetSwarm->num_particles(Portage::Entity_type::PARALLEL_OWNED); 
    auto smoothing_lengths = make_shared<Portage::vector<std::vector<std::vector<double>>>>(ntarpts,
                   std::vector<std::vector<double>>(1, std::vector<double>(3, 2.0/4)));
    DriverTest::set_smoothing_lengths(smoothing_lengths);
  }
};

struct DriverTest3DScatter : DriverTest<3> {
  DriverTest3DScatter() : DriverTest(Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 8, 8, 8),
                              Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 4, 4)) 
  {
    size_t nsrcpts = sourceSwarm->num_particles(Portage::Entity_type::PARALLEL_OWNED); 
    auto smoothing_lengths = make_shared<Portage::vector<std::vector<std::vector<double>>>>(nsrcpts,
                   std::vector<std::vector<double>>(1, std::vector<double>(3, 2.0/8)));
    DriverTest::set_smoothing_lengths(smoothing_lengths, Portage::Meshfree::Scatter);
  }
};
// Methods for computing initial field values
template<size_t Dimension>
double compute_constant_field(Portage::Point<Dimension> coord) {
  return 25.0;
}

template<size_t Dimension>
double compute_linear_field(Portage::Point<Dimension> coord) {
  double val = 0.0;
  for (size_t i = 0; i < Dimension; i++) val += coord[i];
  return val;
}

template<size_t Dimension>
double compute_quadratic_field(Portage::Point<Dimension> coord) {
  double val = 0.0;
  for (size_t i = 0; i < Dimension; i++) val += coord[i]*coord[i];
  return val;
}

template<size_t Dimension>
double compute_cubic_field(Portage::Point<Dimension> coord) {
  double val = 0.0;
  for (size_t i = 0; i < Dimension; i++) val += coord[i]*coord[i]*coord[i];
  return val;
}


// Test cases: these are constructed by calling TEST_F with the name
// of the test class you want to use.  The unit test method is then
// called inside each test with the appropriate template arguments.
// Google test will pick up each test and run it as part of the larger
// test fixture.  If any one of these fails the whole test_driver
// fails.

TEST_F(DriverTest2D, 2D_ConstantFieldUnitaryBasis) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Unitary>
      (compute_constant_field<2>, 0.0);
}

TEST_F(DriverTest2D, 2D_LinearFieldLinearBasis) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Linear>
      (compute_linear_field<2>, 0.0);
}

TEST_F(DriverTest2D, 2D_QuadraticFieldQuadraticBasis) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Quadratic>
      (compute_quadratic_field<2>, 0.0);
}

TEST_F(DriverTest2DScatter, 2D_ConstantFieldUnitaryBasisScatter) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Unitary>
      (compute_constant_field<2>, 0.0);
}

TEST_F(DriverTest2DScatter, 2D_LinearFieldLinearBasisScatter) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Linear>
      (compute_linear_field<2>, 0.0);
}

TEST_F(DriverTest2DScatter, 2D_QuadraticFieldQuadraticBasisScatter) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Quadratic>
      (compute_quadratic_field<2>, 0.0);
}

TEST_F(DriverTest3D, 3D_ConstantFieldUnitaryBasis) {
   unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Unitary>
       (compute_constant_field<3>, 0.0);
}

TEST_F(DriverTest3D, 3D_LinearFieldLinearBasis) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Linear>
      (compute_linear_field<3>, 0.0);
}

TEST_F(DriverTest3D, 3D_QuadraticFieldQuadraticBasis) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Quadratic>
      (compute_quadratic_field<3>, 0.0);
}

TEST_F(DriverTest3DScatter, 3D_ConstantFieldUnitaryBasisScatter) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Unitary>
      (compute_constant_field<3>, 0.0);
}

TEST_F(DriverTest3DScatter, 3D_LinearFieldLinearBasisScatter) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Linear>
      (compute_linear_field<3>, 0.0);
}

TEST_F(DriverTest3DScatter, 3D_QuadraticFieldQuadraticBasisScatter) {
  unitTest<Portage::SearchPointsByCells, Portage::Meshfree::Basis::Quadratic>
      (compute_quadratic_field<3>, 0.0);
}
}  // end namespace
