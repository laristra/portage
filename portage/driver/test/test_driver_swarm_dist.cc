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

#include "portage/support/portage.h"
#include "portage/swarm/swarm.h"
#include "portage/swarm/swarm_state.h"
#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/support/Point.h"
#include "portage/search/search_points_by_cells.h"

#include "gtest/gtest.h"
#include "mpi.h"
#include "Mesh.hh"
#include "MeshFactory.hh"

namespace {
// avoid long namespaces
using namespace Portage::Meshfree;
// default numerical tolerance
double const epsilon = 1e-6;
// default MPI communicator
MPI_Comm comm = MPI_COMM_WORLD;

/**
 * @class BaseTest
 *
 * Set of unitary tests for the swarm-swarm remap driver.
 * There will be at least one test corresponding to each case found in
 * main.cc. This is a test fixture and specialized testcases, such as 2D/3D
 * coincident and non-coincident remaps should be derived from this.
 *
 * @tparam dim: dimension of the problem.
 */
template<int dim>
class BaseTest : public ::testing::Test {
  // skip long type names
  using Mesh = Wonton::Jali_Mesh_Wrapper;

public:

  /**
   * @brief Set swarms and states from distributed meshes
   *
   * @param source_mesh: a shared pointer to the source mesh
   * @param target_mesh: a shared pointer to the target mesh
   */
  BaseTest(std::shared_ptr<Jali::Mesh> source_mesh,
           std::shared_ptr<Jali::Mesh> target_mesh)
    : source_swarm(Mesh(*source_mesh), Wonton::CELL),
      target_swarm(Mesh(*target_mesh), Wonton::CELL),
      source_state(source_swarm),
      target_state(target_swarm) {}

  /**
   * @brief Set the smoothing lengths matrix.
   *
   * @param n: matrix dimensions.
   * @param h: the h-factor.
   * @param center: the weight center (gather, scatter).
   */
  void set_smoothing_lengths(const int* n, double h, WeightCenter center = Gather) {
    assert(n != nullptr);
    std::vector<std::vector<double>> const default_lengths(n[1], std::vector<double>(n[2], h));
    smoothing_lengths_.resize(n[0], default_lengths);
    center_ = center;
  }

  /**
   * @brief Test main method.
   *
   * @tparam Search: the particle search kernel to use.
   * @tparam basis: the basis type.
   * @param compute_initial_field: a function to compute a field to source swarm.
   * @param expected_answer: the expected value to compare against.
   */
  template <template<int, class, class> class Search,
            Basis::Type basis>
  void unitTest(double compute_initial_field(Wonton::Point<dim> const& coord),
                double expected_answer) {

    using Remapper = SwarmDriver<Search, Accumulate, Estimate, dim,
                                 Swarm<dim>, SwarmState<dim>,
                                 Swarm<dim>, SwarmState<dim>>;

    // Fill the source state data with the specified profile
    int const nb_source = source_swarm.num_owned_particles();
    int const nb_target = target_swarm.num_owned_particles();

    // Create the source data for given function
    Portage::vector<double> source_data(nb_source);

    for (int p = 0; p < nb_source; ++p) {
      auto coord = source_swarm.get_particle_coordinates(p);
      source_data[p] = compute_initial_field(coord);
    }

    source_state.add_field("particledata", source_data);
    target_state.add_field("particledata", 0.0);

    // Build the main driver data for this mesh type
    // Register the variable name and interpolation order with the driver
    std::vector<std::string> remap_fields = { "particledata" };
    Wonton::MPIExecutor_type executor(comm);

    Remapper remapper(source_swarm, source_state,
                      target_swarm, target_state,
                      smoothing_lengths_, Weight::B4, Weight::ELLIPTIC, center_);

    remapper.set_remap_var_names(remap_fields, remap_fields, LocalRegression, basis);
    remapper.run(&executor);

    // Check the answer
    double total_error = 0.;
    auto const& target_data = target_state.get_field("particledata");

    for (int i = 0; i < nb_target; ++i) {
      auto p = target_swarm.get_particle_coordinates(i);
      double error = compute_initial_field(p) - target_data[i];
      // dump diagnostics for each particle
      switch (dim) {
        case 1: std::printf("particle: %4d coord: (%5.3lf)", i, p[0]); break;
        case 2: std::printf("particle: %4d coord: (%5.3lf, %5.3lf)", i, p[0], p[1]); break;
        case 3: std::printf("particle: %4d coord: (%5.3lf, %5.3lf, %5.3lf)", i, p[0], p[1], p[2]); break;
        default: break;
      }
      double val = target_data[i];
      std::printf(" value: %10.6lf, error: %lf\n", val, error);
      total_error += error * error;
    }

    total_error = std::sqrt(total_error);
    std::printf("\n\nL2 NORM OF ERROR = %lf\n\n", total_error);
    ASSERT_NEAR(expected_answer, total_error, epsilon);
  }

protected:
  // Source and target swarms
  Swarm<dim> source_swarm;
  Swarm<dim> target_swarm;
  SwarmState<dim> source_state;
  SwarmState<dim> target_state;

  // smoothing lengths matrix and weight center type
  Portage::vector<std::vector<std::vector<double>>> smoothing_lengths_ {};
  WeightCenter center_ = Gather;
};

/**
 * @brief Build 2D swarms from distributed meshes and remap using gather scheme.
 *
 */
class DriverTest2D : public BaseTest<2> {
public:
  DriverTest2D()
    : BaseTest(Jali::MeshFactory(comm)(0.0, 0.0, 1.0, 1.0, 16, 16),
               Jali::MeshFactory(comm)(0.0, 0.0, 1.0, 1.0, 8, 8))
  {
    int const nb_target = target_swarm.num_particles(Wonton::PARALLEL_OWNED);
    int const dim[] = { nb_target, 1, 2 };
    BaseTest::set_smoothing_lengths(dim, 0.5);
  }
};

/**
 * @brief Build 2D swarms from distributed meshes and remap using scatter scheme.
 *
 */
class DriverTest2DScatter : public BaseTest<2> {
public:
  DriverTest2DScatter()
    : BaseTest(Jali::MeshFactory(comm)(0.0, 0.0, 1.0, 1.0, 8, 8),
               Jali::MeshFactory(comm)(0.0, 0.0, 1.0, 1.0, 4, 4))
  {
    int const nb_source = source_swarm.num_particles(Wonton::PARALLEL_OWNED);
    int const dim[] = { nb_source, 1, 2 };
    BaseTest::set_smoothing_lengths(dim, 0.25, Scatter);
  }
};

/**
 * @brief Build 3D swarms from distributed meshes and remap using gather scheme.
 *
 */
class DriverTest3D : public BaseTest<3> {
public:
  DriverTest3D()
    : BaseTest(Jali::MeshFactory(comm)(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,16, 16, 16),
               Jali::MeshFactory(comm)(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 8, 8, 8))
  {
    int const nb_target = target_swarm.num_particles(Wonton::PARALLEL_OWNED);
    int const dim[] = { nb_target, 1, 3 };
    BaseTest::set_smoothing_lengths(dim, 0.5);
  }
};

/**
 * @brief Build 3D swarms from distributed meshes and remap using scatter scheme.
 *
 */
class DriverTest3DScatter : public BaseTest<3> {
public:
  DriverTest3DScatter() : BaseTest(Jali::MeshFactory(comm)(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 8, 8, 8),
                                   Jali::MeshFactory(comm)(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 4, 4))
  {
    int const nb_source = source_swarm.num_particles(Wonton::PARALLEL_OWNED);
    int const dim[] = { nb_source, 1, 3 };
    BaseTest::set_smoothing_lengths(dim, 0.25, Scatter);
  }
};

template<int dim>
double compute_constant_field(Wonton::Point<dim> const& coord) {
  return 25.0;
}

template<int dim>
double compute_linear_field(Wonton::Point<dim> const& coord) {
  double val = 0.0;
  for (int i = 0; i < dim; i++)
    val += coord[i];
  return val;
}

template<int dim>
double compute_quadratic_field(Wonton::Point<dim> const& coord) {
  double val = 0.0;
  for (int i = 0; i < dim; i++)
    for (int j = i; j < dim; j++)
      val += coord[i] * coord[j];
  return val;
}

// Test cases: these are constructed by calling TEST_F with the name
// of the test class you want to use.  The unit test method is then
// called inside each test with the appropriate template arguments.
// Google test will pick up each test and run it as part of the larger
// test fixture.  If any one of these fails the whole test_driver
// fails.

TEST_F(DriverTest2D, 2D_ConstantFieldUnitaryBasis) {
  unitTest<Portage::SearchPointsByCells, Basis::Unitary>
      (compute_constant_field<2>, 0.0);
}

TEST_F(DriverTest2D, 2D_LinearFieldLinearBasis) {
  unitTest<Portage::SearchPointsByCells, Basis::Linear>
      (compute_linear_field<2>, 0.0);
}

TEST_F(DriverTest2D, 2D_QuadraticFieldQuadraticBasis) {
  unitTest<Portage::SearchPointsByCells, Basis::Quadratic>
      (compute_quadratic_field<2>, 0.0);
}

TEST_F(DriverTest2DScatter, 2D_ConstantFieldUnitaryBasisScatter) {
  unitTest<Portage::SearchPointsByCells, Basis::Unitary>
      (compute_constant_field<2>, 0.0);
}

TEST_F(DriverTest2DScatter, 2D_LinearFieldLinearBasisScatter) {
  unitTest<Portage::SearchPointsByCells, Basis::Linear>
      (compute_linear_field<2>, 0.0);
}

TEST_F(DriverTest2DScatter, 2D_QuadraticFieldQuadraticBasisScatter) {
  unitTest<Portage::SearchPointsByCells, Basis::Quadratic>
      (compute_quadratic_field<2>, 0.0);
}

TEST_F(DriverTest3D, 3D_ConstantFieldUnitaryBasis) {
   unitTest<Portage::SearchPointsByCells, Basis::Unitary>
       (compute_constant_field<3>, 0.0);
}

TEST_F(DriverTest3D, 3D_LinearFieldLinearBasis) {
  unitTest<Portage::SearchPointsByCells, Basis::Linear>
      (compute_linear_field<3>, 0.0);
}

TEST_F(DriverTest3D, 3D_QuadraticFieldQuadraticBasis) {
  unitTest<Portage::SearchPointsByCells, Basis::Quadratic>
      (compute_quadratic_field<3>, 0.0);
}

TEST_F(DriverTest3DScatter, 3D_ConstantFieldUnitaryBasisScatter) {
  unitTest<Portage::SearchPointsByCells, Basis::Unitary>
      (compute_constant_field<3>, 0.0);
}

TEST_F(DriverTest3DScatter, 3D_LinearFieldLinearBasisScatter) {
  unitTest<Portage::SearchPointsByCells, Basis::Linear>
      (compute_linear_field<3>, 0.0);
}

TEST_F(DriverTest3DScatter, 3D_QuadraticFieldQuadraticBasisScatter) {
  unitTest<Portage::SearchPointsByCells, Basis::Quadratic>
      (compute_quadratic_field<3>, 0.0);
}
}  // end namespace
