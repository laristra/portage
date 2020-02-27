/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/
#include <cmath>
#include <vector>
#include <ostream>
#include "gtest/gtest.h"

#include "portage/estimate/estimate.h"

#include "portage/swarm/swarm.h"
#include "portage/accumulate/accumulate.h"
#include "portage/support/portage.h"
#include "wonton/support/Point.h"

using Wonton::Point;
using namespace Portage::Meshfree;

template<int dim>
void test_estimate(EstimateType etype, Basis::Type btype, WeightCenter center) {

  // useful aliases
  using Accumulator = Accumulate<dim, Swarm<dim>, Swarm<dim>>;
  using Estimator = Estimate<dim, SwarmState<dim>>;

  // create the source and target swarms input data
  int const nb_sides = 3;
  double const jitter = 0.;  // 0.25*smoothing;
  int const nb_source = static_cast<int>(std::pow(nb_sides + 1, dim));
  int const nb_target = static_cast<int>(std::pow(nb_sides + 3, dim));
  double const source_smoothing = 1. / nb_sides;
  double const target_smoothing = 1. / (nb_sides + 2);

  Portage::vector<Point<dim>> source_points(nb_source);
  Portage::vector<Point<dim>> target_points(nb_target);
  Portage::vector<Point<dim>> source_extents(nb_source);
  Portage::vector<Point<dim>> target_extents(nb_target);
  Portage::vector<Point<dim>> extents(nb_target);

  // set the random engine and generator
  std::random_device device;
  std::mt19937 engine { device() };
  std::uniform_real_distribution<double> generator(0.0, 0.5);

  for (int i = 0; i < nb_source; i++) {
    int offset = 0;
    for (int k = 0; k < dim; k++) {
      int index = (i - offset)/powl(nb_sides + 1, dim - k - 1);
      offset += index * powl(nb_sides + 1, dim - k - 1);

      double const delta = 2. * generator(engine) * jitter * 1.5;
      Point<dim> p = source_points[i];
      p[k] = -0.25 + 1.5 * index * source_smoothing + delta;
      source_points[i] = p;

      p = source_extents[i];
      p[k] = 2.0 * source_smoothing;
       if (index == 0 or index == nb_sides)
         p[k] *= 2;
       source_extents[i] = p;
    }
    extents = source_extents;
  }

  for (int i = 0; i < nb_target; i++) {
    int offset = 0;
    for (int k = 0; k < dim; k++) {
      int index = (i - offset)/powl(nb_sides + 3, dim - k - 1);
      offset += index*powl(nb_sides + 3, dim - k - 1);

      double const delta = 2. * generator(engine) * jitter * 1.5;
      Point<dim> p = target_points[i];
      p[k] = index * target_smoothing + delta;
      target_points[i] = p;

      p = target_extents[i];
      p[k] = 2.0 * target_smoothing;
      if (index == 0 or index == nb_sides + 2)
        p[k] *= 2;
      target_extents[i] = p;
    }
    extents = target_extents;
  }

  // create source+target swarms, kernels, geometries, and smoothing lengths
  Swarm<dim> source_swarm(source_points);
  Swarm<dim> target_swarm(target_points);
  SwarmState<dim> source_state(source_swarm);
  SwarmState<dim> target_state(target_swarm);

  int const nb_kernels = (center == Gather ? nb_target : nb_source);

  Portage::vector<Weight::Kernel> kernels(nb_kernels, Weight::B4);
  Portage::vector<Weight::Geometry> geometries(nb_kernels, Weight::TENSOR);
  Portage::vector<std::vector<std::vector<double>>> smoothing_lengths
    (nb_kernels, std::vector<std::vector<double>>(1, std::vector<double>(dim)));

  for (int i = 0; i < nb_kernels; i++) {
    for (int j = 0; j < dim; j++) {
      std::vector<std::vector<double>> temp = smoothing_lengths[i];
      Point<dim> p = extents[i];
      temp[0][j] = p[j];
      smoothing_lengths[i] = temp;
    }
  }

  // create the accumulator
  Accumulator accumulate(source_swarm, target_swarm,
                         etype, center, kernels,
                         geometries, smoothing_lengths, btype);

  // check sizes and allocate test array
  auto bsize = Basis::function_size<dim>(btype);
  auto jsize = Basis::jet_size<dim>(btype);
  ASSERT_EQ(jsize[0], bsize);
  ASSERT_EQ(jsize[1], bsize);

  // list of src swarm particles (indices)
  std::vector<unsigned> source_particles(nb_source);
  for (int i = 0; i < nb_source; i++)
    source_particles[i] = i;

  // add fields
  std::string cnums[15] = { "00","01","02","03","04",
                            "05","06","07","08","09",
                            "10","11","12","13","14" };
  std::string field_names[15];
  int const nb_basis = jsize[0];

  for (int i = 0; i < nb_basis; i++)
    field_names[i] = "field" + cnums[i];

  for (int i = 0; i < nb_basis; i++) {
    Portage::vector<double> source_field(nb_source);
    Portage::vector<double> target_field(nb_target, 0.);

    // fill source field
    for (int k = 0; k < nb_source; k++) {
      auto y = source_swarm.get_particle_coordinates(k);
      auto basis_y = Basis::function<dim>(btype, y);
      source_field[k] = basis_y[i];
    }

    // add fields to particle states
    source_state.add_field(field_names[i], source_field);
    target_state.add_field(field_names[i], target_field);
  }

  // make the estimator
  Estimator estimate(source_state);
    
  // Loop through target particles
  for (int i = 0; i < nb_target; i++) {
    auto x = target_swarm.get_particle_coordinates(i);
    auto jetx = Basis::jet<dim>(btype,x);
    auto sources_and_mults = accumulate(i, source_particles);

    // count actual neighbors
    int nb_neigh = 0;
    for (int k = 0; k < nb_source; k++) {
      if (accumulate.weight(i, k) != 0.)
        nb_neigh++;
    }

    if (nb_neigh < nb_basis) {
      std::cout << "number of neighbors " << nb_neigh;
      std::cout << " is too small at " << i << std::endl;
    }

    // Loop through fields
    for (int j = 0; j < nb_basis; j++) {
      // get the target field
      auto& field = target_state.get_field(field_names[j]);
      // Loop through derivatives
      for (int k = 0; k < nb_basis; k++) {
        // tell estimate what to do
        estimate.set_variable(field_names[j], k);
        // do the estimate
        double result = estimate(i, sources_and_mults);
        // save value
        if (k == 0)
          field[i] = result;
        // check the estimate, verifying reproducing property
        ASSERT_NEAR(result, jetx[j][k], 1.e-11);
      }
    }
  }
}

TEST(estimate, 1d_RUG) {
  test_estimate<1>(EstimateType::LocalRegression,
                   Basis::Type::Unitary,
                   WeightCenter::Gather);
}

TEST(estimate, 2d_RUG) {
  test_estimate<2>(EstimateType::LocalRegression,
                   Basis::Type::Unitary,
                   WeightCenter::Gather);
}

TEST(estimate, 3d_RUG) {
  test_estimate<3>(EstimateType::LocalRegression,
                   Basis::Type::Unitary,
                   WeightCenter::Gather);
}

TEST(estimate, 1d_RUS) {
  test_estimate<1>(EstimateType::LocalRegression,
                   Basis::Type::Unitary,
                   WeightCenter::Scatter);
}

TEST(estimate, 2d_RUS) {
  test_estimate<2>(EstimateType::LocalRegression,
                   Basis::Type::Unitary,
                   WeightCenter::Scatter);
}

TEST(estimate, 3d_RUS) {
  test_estimate<3>(EstimateType::LocalRegression,
                   Basis::Type::Unitary,
                   WeightCenter::Scatter);
}

TEST(estimate, 1d_RLG) {
  test_estimate<1>(EstimateType::LocalRegression,
                   Basis::Type::Linear,
                   WeightCenter::Gather);
}

TEST(estimate, 2d_RLG) {
  test_estimate<2>(EstimateType::LocalRegression,
                   Basis::Type::Linear,
                   WeightCenter::Gather);
}

TEST(estimate, 3d_RLG) {
  test_estimate<3>(EstimateType::LocalRegression,
                   Basis::Type::Linear,
                   WeightCenter::Gather);
}

TEST(estimate, 1d_RLS) {
  test_estimate<1>(EstimateType::LocalRegression,
                   Basis::Type::Linear,
                   WeightCenter::Scatter);
}

TEST(estimate, 2d_RLS) {
  test_estimate<2>(EstimateType::LocalRegression,
                   Basis::Type::Linear,
                   WeightCenter::Scatter);
}

TEST(estimate, 3d_RLS) {
  test_estimate<3>(EstimateType::LocalRegression,
                   Basis::Type::Linear,
                   WeightCenter::Scatter);
}

TEST(estimate, 1d_RQG) {
  test_estimate<1>(EstimateType::LocalRegression,
                   Basis::Type::Quadratic,
                   WeightCenter::Gather);
}

TEST(estimate, 2d_RQG) {
  test_estimate<2>(EstimateType::LocalRegression,
                   Basis::Type::Quadratic,
                    WeightCenter::Gather);
}

TEST(estimate, 3d_RQG) {
  test_estimate<3>(EstimateType::LocalRegression,
                   Basis::Type::Quadratic,
                   WeightCenter::Gather);
}

TEST(estimate, 1d_RQS) {
  test_estimate<1>(EstimateType::LocalRegression,
                   Basis::Type::Quadratic,
                   WeightCenter::Scatter);
}

TEST(estimate, 2d_RQS) {
  test_estimate<2>(EstimateType::LocalRegression,
                   Basis::Type::Quadratic,
                   WeightCenter::Scatter);
}

TEST(estimate, 3d_RQS) {
  test_estimate<3>(EstimateType::LocalRegression,
                   Basis::Type::Quadratic,
                   WeightCenter::Scatter);
}




