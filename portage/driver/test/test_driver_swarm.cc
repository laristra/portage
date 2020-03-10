/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/


#include <iostream>
#include <memory>
#include <cassert>
#include <cmath>
#include <string>

#include "gtest/gtest.h"
#include "mpi.h"

#include "portage/driver/driver_swarm.h"
#include "portage/swarm/swarm.h"
#include "portage/swarm/swarm_state.h"
#include "portage/accumulate/accumulate.h"
#include "portage/estimate/estimate.h"
#include "portage/support/operator.h"
#include "portage/support/faceted_setup.h"

#include "portage/support/portage.h"
#include "portage/search/search_points_by_cells.h"
#include "wonton/support/Point.h"
#include "wonton/mesh/simple/simple_mesh.h"
#include "wonton/mesh/simple/simple_mesh_wrapper.h"
#include "wonton/state/simple/simple_state.h"
#include "wonton/state/simple/simple_state_wrapper.h"

namespace {
// avoid long namespaces
using namespace Portage::Meshfree;
// default numerical tolerance
double const epsilon = 1e-6;

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
public:
  /**
   * @brief Set swarms and states for tests, initialize the operator to use.
   *
   * @param nb_source: number of source particles
   * @param nb_target: number of target particles
   * @param distrib: the distribution to use when generating random particles
   * @param x_min: lower bound on particle coordinates for x-axis
   * @param x_max: upper bound on particle coordinates for x-axis
   * @param y_min: lower bound on particle coordinates for y-axis (for 2D/3D)
   * @param y_max: upper bound on particle coordinates for x-axis (for 2D/3D)
   * @param z_min: lower bound on particle coordinates for z-axis (for 3D)
   * @param z_max: upper bound on particle coordinates for z-axis (for 3D)
   * @param op: the operator to use (volume/surface integral, last operator)
   */
  BaseTest(int nb_source, int nb_target, int distrib,
           double x_min = 0.0, double x_max = 0.0,
           double y_min = 0.0, double y_max = 0.0,
           double z_min = 0.0, double z_max = 0.0,
           oper::Type op = oper::LastOperator)
    : source_swarm(nb_source, distrib, 0, x_min, x_max, y_min, y_max, z_min, z_max),
      target_swarm(nb_target, distrib, 0, x_min, x_max, y_min, y_max, z_min, z_max),
      source_state(source_swarm),
      target_state(target_swarm),
      center_(Gather),
      operator_(op)
  {

    if (op != oper::LastOperator) {

      int const num_owned_target = target_swarm.num_owned_particles();
      int const npoints[] = { 2, 4, 8 };

      domains_ = Portage::vector<oper::Domain>(num_owned_target);
      operator_data_.resize(num_owned_target, std::vector<Wonton::Point<dim>>(npoints[dim-1]));

      int npdim = static_cast<int>(std::pow(1.001 * num_owned_target, 1./dim));
      int nptot = 1;
      for (int m = 0; m < dim; m++)
        nptot *= npdim;

      assert(nptot == num_owned_target);
      assert(npdim > 1);

      int i = 0;
      int j = 0;
      int k = 0;
      int ij[npoints[dim-1]];
      std::vector<Wonton::Point<dim>> points(npoints[dim-1]);

      // Create points determining the integration volume.
      // Assumes target swarm is created from SwarmFactory and represents a perfect cubic array of points.
      for (int n = 0; n < nptot; n++) {
        switch(dim) {
          case 1:
            i = n;
            ij[0] = i;
            ij[1] = (i+1);
            break;
          case 2:
            i = n / npdim;
            j = n - i * npdim;
            ij[0] = i * npdim + j;
            ij[1] = (i+1) * npdim + j;
            ij[2] = (i+1) * npdim + j + 1;
            ij[3] = i * npdim + j + 1;
            break;
          case 3:
            i = n / (npdim * npdim);
            j = (n - i * npdim * npdim) / npdim;
            k = n - i * npdim * npdim - j * npdim;
            ij[0] = npdim * (i * npdim + j) + k;
            ij[1] = npdim * ((i+1) * npdim + j) + k;
            ij[2] = npdim * ((i+1) * npdim + j + 1) + k;
            ij[3] = npdim * (i * npdim + j + 1) + k;
            ij[4] = npdim * (i * npdim + j) + k + 1;
            ij[5] = npdim * ((i+1) * npdim + j) + k + 1;
            ij[6] = npdim * ((i+1) * npdim + j + 1) + k + 1;
            ij[7] = npdim * (i * npdim + j + 1) + k + 1;
            break;
          default: break;
        }

        for (int m = 0; m < npoints[dim-1]; m++) {
          bool on_boundary = (i == npdim-1 or j == npdim-1 or k == npdim-1);
          // particles on the upper boundaries get an integration volume of zero
          points[m] = (on_boundary ? target_swarm.get_particle_coordinates(ij[0])
                                   : target_swarm.get_particle_coordinates(ij[m]));
        }
        operator_data_[n] = points;
        domains_[n] = oper::domain_from_points<dim>(points);
      }
    }
  }

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
            basis::Type basis>
  void unitTest(double compute_initial_field(Wonton::Point<dim> const& coord),
                double expected_answer) {

    using Remapper = SwarmDriver<Search, Accumulate, Estimate, dim,
                                 Swarm<dim>, SwarmState<dim>,
                                 Swarm<dim>, SwarmState<dim>>;

    // Fill the source state data with the specified profile
    const int nb_source = source_swarm.num_owned_particles();
    const int nb_target = target_swarm.num_owned_particles();

    Portage::vector<double> source_data(nb_source);
    Portage::vector<double> target_data(nb_target, 0.0);

    // Create the source data for given function
    for (int i = 0; i < nb_source; ++i) {
      auto coord = source_swarm.get_particle_coordinates(i);
      source_data[i] = compute_initial_field(coord);
    }

    source_state.add_field("particledata", source_data);
    target_state.add_field("particledata", target_data);

    // Build the main driver data for this mesh type
    Remapper remapper(source_swarm, source_state,
                      target_swarm, target_state,
                      smoothing_lengths_,Weight::B4, Weight::ELLIPTIC, center_);

    EstimateType estimator = LocalRegression;
    if (operator_ != oper::LastOperator)
      estimator = OperatorRegression;

    // Register the variable name with the driver
    std::vector<std::string> remap_fields = { "particledata" };
    remapper.set_remap_var_names(remap_fields, remap_fields,
                                 estimator, basis,
                                 operator_, domains_, operator_data_);

    // run on one processor (no argument implies serial run)
    remapper.run();

    // Check the answer
    double total_error = 0.;
    auto& remapped_field = target_state.get_field("particledata");

    if (operator_ == oper::LastOperator) {
      for (int i = 0; i < nb_target; ++i) {
        auto p = target_swarm.get_particle_coordinates(i);
        double error = compute_initial_field(p) - remapped_field[i];
        // dump diagnostics for each particle
        switch (dim) {
          case 1: std::printf("particle: %4d coord: (%5.3lf)", i, p[0]); break;
          case 2: std::printf("particle: %4d coord: (%5.3lf, %5.3lf)", i, p[0], p[1]); break;
          case 3: std::printf("particle: %4d coord: (%5.3lf, %5.3lf, %5.3lf)", i, p[0], p[1], p[2]); break;
          default: break;
        }
	      double val = remapped_field[i];
	      std::printf(" value: %10.6lf, error: %lf\n", val, error);
        total_error += error * error;
      }

      total_error = std::sqrt(total_error);
      std::printf("\n\nL2 NORM OF ERROR = %lf\n\n", total_error);
      ASSERT_NEAR(expected_answer, total_error, epsilon);

    } else if (operator_ == oper::VolumeIntegral) {
      double total = 0.;
      for (int i = 0; i < nb_target; ++i) {
        total += remapped_field[i];
      }
      ASSERT_NEAR(expected_answer, total, epsilon);
    }
  }

  /**
   * @brief Test alternate method using a more detailed constructor.
   *
   * @tparam Search: the particle search kernel to use.
   * @tparam basis: the basis type.
   * @param compute_initial_field: a function to compute a field to source swarm.
   * @param expected_answer: the expected value to compare against.
   */
  template <template<int, class, class> class Search,
            basis::Type basis>
  void unitTestAlt(double compute_initial_field(Wonton::Point<dim> const& coord),
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

    // Set kernels and geometries
    if (center_ == Gather) {
      kernels_.resize(nb_target, Weight::B4);
      geometries_.resize(nb_target, Weight::ELLIPTIC);
    } else if (center_ == Scatter) {
      kernels_.resize(nb_source, Weight::B4);
      geometries_.resize(nb_source, Weight::ELLIPTIC);
    }

    // Build the main driver data for this mesh type
    Remapper remapper(source_swarm, source_state,
                      target_swarm, target_state,
                      smoothing_lengths_, kernels_, geometries_, center_);

    EstimateType estimator = LocalRegression;
    if (operator_ != oper::LastOperator)
      estimator = OperatorRegression;

    // Register the variable name and interpolation order with the driver
    std::vector<std::string> remap_fields = { "particledata" };
    remapper.set_remap_var_names(remap_fields, remap_fields,
                                 estimator, basis,
                                 operator_, domains_, operator_data_);

    // run on one processor (no argument implies serial run)
    remapper.run();

    // Check the answer
    double total_error = 0.;
    auto const& target_data = target_state.get_field("particledata");

    if (operator_ == oper::LastOperator) {
      for (int i = 0; i < nb_target; ++i) {
        auto p = target_swarm.get_particle_coordinates(i);
        auto error = compute_initial_field(p) - target_data[i];
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

      //total_error = std::sqrt(total_error);
      std::printf("\n\nL2 NORM OF ERROR = %lf\n\n", std::sqrt(total_error));
      ASSERT_NEAR(expected_answer, std::sqrt(total_error/nb_target), epsilon);

    } else if (operator_ == oper::VolumeIntegral) {
      double total = 0.;
      for (int p = 0; p < nb_target; ++p) {
        total += target_data[p];
      }
      ASSERT_NEAR(expected_answer, total, epsilon);
    }
  }

protected:
  // swarms and states
  Swarm<dim> source_swarm;
  Swarm<dim> target_swarm;
  SwarmState<dim> source_state;
  SwarmState<dim> target_state;

  // smoothing lengths
  Portage::vector<std::vector<std::vector<double>>> smoothing_lengths_;

  // kernel and geometry specifications
  Portage::vector<Weight::Kernel> kernels_ {};
  Portage::vector<Weight::Geometry> geometries_ {};

  // operator info
  WeightCenter center_;
  oper::Type operator_;
  Portage::vector<oper::Domain> domains_;
  Portage::vector<std::vector<Wonton::Point<dim>>> operator_data_;
};

/**
 * @brief Remap 1D random swarm using a gather scheme.
 *
 */
class DriverTest1DGather : public BaseTest<1> {
public:
  DriverTest1DGather() : BaseTest<1>(7, 5, 2, 0.0, 1.0) {
    int const dim[] = { 5, 1, 1 };
    BaseTest<1>::set_smoothing_lengths(dim, 0.5, Gather);
  }
};

/**
 * @brief Remap 2D random swarm using a gather scheme.
 *
 */
class DriverTest2DGather : public BaseTest<2> {
public:
  DriverTest2DGather() : BaseTest<2>(7*7, 5*5, 2, 0.0, 1.0, 0.0, 1.0) {
    int const dim[] = { 5*5, 1, 2 };
    BaseTest<2>::set_smoothing_lengths(dim, 0.5, Gather);
  }
};

/**
 * @brief Remap 3D random swarm using a gather scheme.
 *
 */
class DriverTest3DGather : public BaseTest<3> {
public:
  DriverTest3DGather() : BaseTest<3>(7*7*7, 5*5*5, 2, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0) {
    int const dim[] = { 5*5*5, 1, 3 };
    BaseTest<3>::set_smoothing_lengths(dim, 0.5, Gather);
  }
};

/**
 * @brief Remap 1D random swarm using a scatter scheme.
 *
 */
class DriverTest1DScatter : public BaseTest<1> {
public:
  DriverTest1DScatter() : BaseTest<1>(7, 5, 2, 0.0, 1.0) {
    int const dim[] = { 7, 1, 1 };
    BaseTest<1>::set_smoothing_lengths(dim, 0.5, Scatter);
  }
};

/**
 * @brief Remap 2D random swarm using a scatter scheme.
 *
 */
class DriverTest2DScatter : public BaseTest<2> {
public:
  DriverTest2DScatter() : BaseTest<2>(7*7, 5*5, 2, 0.0, 1.0, 0.0, 1.0) {
    int const dim[] = { 7*7*7, 1, 2 };
    BaseTest<2>::set_smoothing_lengths(dim, 0.5, Scatter);
  }
};

/**
 * @brief Remap 3D random swarm using a scatter scheme.
 *
 */
class DriverTest3DScatter : public BaseTest<3> {
public:
  DriverTest3DScatter() : BaseTest<3>(7*7*7, 5*5*5, 2, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0) {
    int const dim[] = { 7*7*7, 1, 3 };
    BaseTest<3>::set_smoothing_lengths(dim, 0.5, Scatter);
  }
};

/**
 * @brief Remap 1D ordered swarm and integrate.
 *
 */
class IntegrateDriverTest1D : public BaseTest<1> {
public:
  IntegrateDriverTest1D() : BaseTest<1>(7, 5, 1, 0.0, 1.0,
                                        0.0, 1.0, 0.0, 1.0, oper::VolumeIntegral) {
    int const dim[] = { 5, 1, 1 };
    BaseTest<1>::set_smoothing_lengths(dim, 0.5);
  }
};

/**
 * @brief Remap 2D ordered swarm and integrate.
 *
 */
class IntegrateDriverTest2D : public BaseTest<2> {
public:
  IntegrateDriverTest2D() : BaseTest<2>(7*7, 5*5, 1, 0.0, 1.0,
                                        0.0, 1.0, 0.0, 1.0, oper::VolumeIntegral) {
    int const dim[] = { 5*5, 1, 2 };
    BaseTest<2>::set_smoothing_lengths(dim, 0.5);
  }
};

/**
 * @brief Remap 3D ordered swarm and integrate.
 *
 */
class IntegrateDriverTest3D : public BaseTest<3> {
public:
  IntegrateDriverTest3D() : BaseTest<3>(7*7*7, 5*5*5, 1, 0.0, 1.0,
                                        0.0, 1.0, 0.0, 1.0, oper::VolumeIntegral) {
    int const dim[] = { 5*5*5, 1, 3 };
    BaseTest<3>::set_smoothing_lengths(dim, 0.5);
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

/**
 * @brief Check if a particle is strictly negative or positive.
 *
 * @tparam dim: particle dimension
 * @param p: particle coordinates
 * @return true if strictly negative or positive, false otherwise.
 */
template<int dim>
bool not_zero(Wonton::Point<dim> const& p) {
  switch (dim) {
    case 1: return p[0] < 0. or p[0] > 0.;
    case 2: return (p[0] < 0. and p[1] < 0.) or (p[0] > 0. and p[1] > 0.);
    case 3: return (p[0] < 0. and p[1] < 0. and p[2] < 0.)
                or (p[0] > 0. and p[1] > 0. and p[2] > 0.);
    default: return false;
  }
}

// Test cases: these are constructed by calling TEST_F with the name
// of the test class you want to use.  The unit test method is then
// called inside each test with the appropriate template arguments.
// Google test will pick up each test and run it as part of the larger
// test fixture.  If any one of these fails the whole test_driver
// fails.
//
TEST_F(DriverTest1DGather, 1D_ConstantFieldUnitaryBasis) {
   unitTest<Portage::SearchPointsByCells, basis::Unitary>
       (compute_constant_field<1>, 0.0);
}

TEST_F(DriverTest1DGather, 1D_LinearFieldLinearBasis) {
  unitTest<Portage::SearchPointsByCells, basis::Linear>
      (compute_linear_field<1>, 0.0);
}

TEST_F(DriverTest1DGather, 1D_QuadraticFieldQuadraticBasis) {
  unitTest<Portage::SearchPointsByCells, basis::Quadratic>
      (compute_quadratic_field<1>, 0.0);
}

TEST_F(DriverTest1DScatter, 1D_QuadraticFieldQuadraticBasisScatter) {
  unitTest<Portage::SearchPointsByCells, basis::Quadratic>
      (compute_quadratic_field<1>, 0.0);
}

TEST_F(DriverTest2DGather, 2D_ConstantFieldUnitaryBasis) {
  unitTest<Portage::SearchPointsByCells, basis::Unitary>
      (compute_constant_field<2>, 0.0);
}

TEST_F(DriverTest2DGather, 2D_LinearFieldLinearBasis) {
  unitTest<Portage::SearchPointsByCells, basis::Linear>
      (compute_linear_field<2>, 0.0);
}

TEST_F(DriverTest2DGather, 2D_LinearFieldLinearBasisAlt) {
  unitTestAlt<Portage::SearchPointsByCells, basis::Linear>
      (compute_linear_field<2>, 0.0);
}

TEST_F(DriverTest2DGather, 2D_QuadraticFieldQuadraticBasis) {
  unitTest<Portage::SearchPointsByCells, basis::Quadratic>
      (compute_quadratic_field<2>, 0.0);
}

TEST_F(DriverTest2DScatter, 2D_QuadraticFieldQuadraticBasisScatter) {
  unitTest<Portage::SearchPointsByCells, basis::Quadratic>
      (compute_quadratic_field<2>, 0.0);
}

TEST_F(DriverTest2DScatter, 2D_QuadraticFieldQuadraticBasisScatterAlt) {
  unitTestAlt<Portage::SearchPointsByCells, basis::Quadratic>
      (compute_quadratic_field<2>, 0.0);
}

TEST_F(DriverTest3DGather, 3D_ConstantFieldUnitaryBasis) {
   unitTest<Portage::SearchPointsByCells, basis::Unitary>
       (compute_constant_field<3>, 0.0);
}

TEST_F(DriverTest3DGather, 3D_LinearFieldLinearBasis) {
  unitTest<Portage::SearchPointsByCells, basis::Linear>
      (compute_linear_field<3>, 0.0);
}

TEST_F(DriverTest3DGather, 3D_LinearFieldLinearBasisAlt) {
  unitTestAlt<Portage::SearchPointsByCells, basis::Linear>
      (compute_linear_field<3>, 0.0);
}

TEST_F(DriverTest3DGather, 3D_QuadraticFieldQuadraticBasis) {
  unitTest<Portage::SearchPointsByCells, basis::Quadratic>
      (compute_quadratic_field<3>, 0.0);
}

TEST_F(DriverTest3DScatter, 3D_QuadraticFieldQuadraticBasisScatter) {
  unitTest<Portage::SearchPointsByCells, basis::Quadratic>
      (compute_quadratic_field<3>, 0.0);
}

TEST_F(DriverTest3DScatter, 3D_QuadraticFieldQuadraticBasisScatterAlt) {
  unitTestAlt<Portage::SearchPointsByCells, basis::Quadratic>
      (compute_quadratic_field<3>, 0.0);
}

TEST_F(IntegrateDriverTest1D, 1D_LinearFieldLinearBasis) {
  unitTest<Portage::SearchPointsByCells, basis::Linear>
      (compute_linear_field<1>, 1./2.);
}

TEST_F(IntegrateDriverTest1D, 1D_QuadraticFieldQuadraticBasis) {
  unitTest<Portage::SearchPointsByCells, basis::Quadratic>
      (compute_quadratic_field<1>, 1./3.);
}

TEST_F(IntegrateDriverTest2D, 2D_LinearFieldLinearBasis) {
  unitTest<Portage::SearchPointsByCells, basis::Linear>
      (compute_linear_field<2>, 1.0);
}

TEST_F(IntegrateDriverTest2D, 2D_QuadraticFieldQuadraticBasis) {
  unitTest<Portage::SearchPointsByCells, basis::Quadratic>
    (compute_quadratic_field<2>, 2./3. + 1./4.); // 0.9166666666666666
}

TEST_F(IntegrateDriverTest3D, 3D_LinearFieldLinearBasis) {
  unitTest<Portage::SearchPointsByCells, basis::Linear>
    (compute_linear_field<3>, 3./2.);
}

TEST(Part, 2D) {

  using Remapper = SwarmDriver<Portage::SearchPointsByCells,
                               Accumulate, Estimate, 2,
                               Swarm<2>, SwarmState<2>,
                               Swarm<2>, SwarmState<2>>;

  int const nb_per_axis = 4;
  Wonton::Simple_Mesh mesh(-1.0, -1.0, 1.0, 1.0, nb_per_axis, nb_per_axis);
  Wonton::Simple_Mesh_Wrapper  mesh_wrapper(mesh);

  int const nb_target = std::pow(2 * nb_per_axis + 2, 2);
  int const nb_source = mesh.num_entities(Wonton::CELL, Wonton::PARALLEL_OWNED);
  assert(nb_source == nb_per_axis * nb_per_axis);

  std::vector<double> values(nb_source, 1.0);

  for (int i = 0; i < nb_source; i++) {
    Wonton::Point<2> p;
    mesh_wrapper.cell_centroid(i, &p);
    if (not_zero(p))
      values[i] = 2.0;
  }

  double* source_data = values.data();

  Wonton::Simple_State state(std::make_shared<Wonton::Simple_Mesh>(mesh));
  auto& added = state.add("indicate", Wonton::CELL, source_data);
  Wonton::Simple_State_Wrapper state_wrapper(state);

  // create source and target swarms and states
  Swarm<2> source_swarm(mesh_wrapper, Wonton::CELL);
  Swarm<2> target_swarm(nb_target, 1, 0, -1.0, 1.0, -1.0, 1.0);
  SwarmState<2> source_state(state_wrapper, Wonton::CELL);
  SwarmState<2> target_state(target_swarm);

  // keep all target swarm points from overlying any source points
  Portage::vector<double> target_data(nb_target, 0.0);
  target_state.add_field("indicate", target_data);

  Portage::vector<std::vector<std::vector<double>>> smoothing;
  Portage::vector<Wonton::Point<2>> extents;
  Portage::vector<Wonton::Point<2>> dummy;
  Weight::faceted_setup_cell<2, Wonton::Simple_Mesh_Wrapper>(mesh_wrapper, smoothing,
                                                          extents, 1.0, 2.0);

  Remapper remapper(source_swarm, source_state,
                    target_swarm, target_state,
                    smoothing, extents, dummy, Scatter);

  std::vector<std::string> const fields_names = {"indicate" };
  Portage::vector<std::vector<std::vector<double>>> psmoothing;
  Weight::faceted_setup_cell<2, Wonton::Simple_Mesh_Wrapper>(mesh_wrapper, psmoothing,
                                                            dummy, 0.25, 1.0);

  remapper.set_remap_var_names(fields_names, fields_names,
                               LocalRegression,
                               basis::Unitary,
                               oper::LastOperator,
                               Portage::vector<oper::Domain>(0),
                               Portage::vector<std::vector<Point<2>>>(0,std::vector<Point<2>>(0)),
                               "indicate", 0.25, psmoothing);

  remapper.run();

  auto& indicator = target_state.get_field("indicate");

  for (int i = 0; i < nb_target; i++) {
    auto const p = target_swarm.get_particle_coordinates(i);
    auto const expected = not_zero(p) ? 2.0 : 1.0;
    ASSERT_NEAR(expected, indicator[i], 1.E-12);
  }
}

TEST(Part, 3D) {

  using Remapper = SwarmDriver<Portage::SearchPointsByCells,
                               Accumulate, Estimate, 3,
                               Swarm<3>, SwarmState<3>,
                               Swarm<3>, SwarmState<3>>;

  int const nb_per_axis = 4;
  Wonton::Simple_Mesh mesh(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0,
                           nb_per_axis, nb_per_axis, nb_per_axis);
  Wonton::Simple_Mesh_Wrapper mesh_wrapper(mesh);

  int const nb_target = std::pow(2 * nb_per_axis + 2, 3);
  int const nb_source = mesh.num_entities(Wonton::CELL, Wonton::PARALLEL_OWNED);
  assert(nb_source == nb_per_axis * nb_per_axis * nb_per_axis);

  std::vector<double> values(nb_source, 1.0);

  for (int i = 0; i < nb_source; i++) {
    Wonton::Point<3> p;
    mesh_wrapper.cell_centroid(i, &p);
    if (not_zero(p))
      values[i] = 2.0;
  }

  double* source_data = values.data();

  Wonton::Simple_State state(std::make_shared<Wonton::Simple_Mesh>(mesh));
  auto& added = state.add("indicate", Wonton::CELL, source_data);
  Wonton::Simple_State_Wrapper state_wrapper(state);

  // create source and target swarms and states
  Swarm<3> source_swarm(mesh_wrapper, Wonton::CELL);
  Swarm<3> target_swarm(nb_target, 1, 0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
  SwarmState<3> source_state(state_wrapper, Wonton::CELL);
  SwarmState<3> target_state(target_swarm);

  // keep all target swarm points from overlying any source points
  Portage::vector<double> target_data(nb_target, 0.0);
  target_state.add_field("indicate", target_data);

  Portage::vector<std::vector<std::vector<double>>> smoothing;
  Portage::vector<Wonton::Point<3>> extents;
  Portage::vector<Wonton::Point<3>> dummy;
  Weight::faceted_setup_cell<3, Wonton::Simple_Mesh_Wrapper>(mesh_wrapper, smoothing,
                                                            extents, 1.0, 2.0);

  Remapper remapper(source_swarm, source_state,
                    target_swarm, target_state,
                    smoothing, extents, dummy, Scatter);

  std::vector<std::string> const fields_names = {"indicate" };
  Portage::vector<std::vector<std::vector<double>>> psmoothing;
  Weight::faceted_setup_cell<3, Wonton::Simple_Mesh_Wrapper>(mesh_wrapper, psmoothing,
                                                             dummy, 0.25, 1.0);

  remapper.set_remap_var_names(fields_names, fields_names,
                               LocalRegression,
                               basis::Unitary,
                               oper::LastOperator,
                               Portage::vector<oper::Domain>(0),
                               Portage::vector<std::vector<Point<3>>>(0,std::vector<Point<3>>(0)),
                               "indicate", 0.25, psmoothing);

  remapper.run();

  auto& indicator = target_state.get_field("indicate");

  for (int i = 0; i < nb_target; i++) {
    auto const p = target_swarm.get_particle_coordinates(i);
    auto const expected = not_zero(p) ? 2.0 : 1.0;
    ASSERT_NEAR(expected, indicator[i], 1.E-12);
  }
}

}  // end namespace
