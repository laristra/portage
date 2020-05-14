/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <iostream>

#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"
#include "wonton/mesh/jali/jali_mesh_wrapper.h"

#include "portage/accumulate/accumulate.h"
#include "portage/distributed/mpi_particle_distribute.h"

#include "portage/support/portage.h"
#include "portage/swarm/swarm.h"
#include "portage/swarm/swarm_state.h"
#include "portage/support/faceted_setup.h"

#include "gtest/gtest.h"
#include "mpi.h"

#include "Mesh.hh"
#include "MeshFactory.hh"

TEST(MPI_Particle_Distribute, SimpleTest2DGather) {

  using Wonton::Point;
  using namespace Portage::Meshfree;

  // set MPI info
  int rank = 0;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  Wonton::MPIExecutor_type executor(comm);

  // Create a distributed jali source/target mesh
  Jali::MeshFactory mf(comm);
  mf.partitioner(Jali::Partitioner_type::BLOCK);
  auto source_mesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4);
  auto target_mesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4);

  Wonton::Jali_Mesh_Wrapper source_mesh_wrapper(*source_mesh);
  Wonton::Jali_Mesh_Wrapper target_mesh_wrapper(*target_mesh);

  // Source and target swarms 
  Swarm<2> source_swarm(source_mesh_wrapper, Wonton::CELL);
  Swarm<2> target_swarm(target_mesh_wrapper, Wonton::CELL);
  SwarmState<2> source_state(source_swarm);
  SwarmState<2> target_state(target_swarm);

  int const nb_source = source_swarm.num_particles(Wonton::ALL);
  int const nb_target = target_swarm.num_particles(Wonton::ALL);

  // Create an integer source data for given function
  Wonton::vector<int> source_data_int(nb_source);

  // Fill the source state data with the specified profile
  for (int i = 0; i < nb_source; ++i) {
    auto p = source_swarm.get_particle_coordinates(i);
    source_data_int[i] = static_cast<int>(p[0] * p[1]) * 100;
  }

  // Create a double source data for given function
  Wonton::vector<double> source_data_dbl(nb_source);

  // Fill the source state data with the specified profile
  for (int i = 0; i < nb_source; ++i) {
    auto p = source_swarm.get_particle_coordinates(i);
    source_data_dbl[i] = p[0] * p[1];
  }

  source_state.add_field("intdata", source_data_int);
  source_state.add_field("dbldata", source_data_dbl);

  // Build the target state storage
  Wonton::vector<int> target_data_int(nb_target, 0);
  Wonton::vector<double> target_data_dbl(nb_target, 0.0);

  target_state.add_field("intdata", target_data_int);
  target_state.add_field("dbldata", target_data_dbl);

  //Set smoothing lengths
  double const one_third = 1./3.;
  Wonton::Point<2> const default_point(one_third, one_third);
  std::vector<std::vector<double>> const default_lengths(1, std::vector<double>(2, one_third));

  Wonton::vector<std::vector<std::vector<double>>> smoothing_lengths(nb_target, default_lengths);
  Wonton::vector<Wonton::Point<2>> source_extents(1, default_point);
  Wonton::vector<Wonton::Point<2>> target_extents(nb_target, default_point);
  Wonton::vector<Weight::Kernel> kernel_types(nb_target, Weight::B4);
  Wonton::vector<Weight::Geometry> geom_types(nb_target, Weight::ELLIPTIC);

  // Distribute
  Portage::MPI_Particle_Distribute<2> distributor(&executor);
  distributor.distribute(source_swarm, source_state,
                         target_swarm, target_state,
                         smoothing_lengths, source_extents, target_extents,
                         kernel_types, geom_types, WeightCenter::Gather);

  MPI_Barrier(comm);
  // Check number of particles received
  int nb_source_after = source_swarm.num_particles(Wonton::ALL);
  ASSERT_EQ(nb_source_after, nb_source + 5);

  // check coordinates
  auto& source_data_after_int = source_state.get_field_int("intdata");
  auto& source_data_after_dbl = source_state.get_field_dbl("dbldata");

  for (int i = 0 ; i < nb_source_after; ++i) {
    auto coords = source_swarm.get_particle_coordinates(i);
    double const value = coords[0] * coords[1];
    ASSERT_EQ(source_data_after_int[i], static_cast<int>(value) * 100);
    ASSERT_DOUBLE_EQ(source_data_after_dbl[i], value);
  }
}

TEST(MPI_Particle_Distribute, SimpleTest2DScatter) {

  using namespace Portage::Meshfree;

  // set MPI info
  int rank = 0;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  Wonton::MPIExecutor_type executor(comm);

  // Create a distributed jali source/target mesh
  Jali::MeshFactory mf(comm);
  mf.partitioner(Jali::Partitioner_type::BLOCK);
  auto source_mesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4);
  auto target_mesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4);

  Wonton::Jali_Mesh_Wrapper source_mesh_wrapper(*source_mesh);
  Wonton::Jali_Mesh_Wrapper target_mesh_wrapper(*target_mesh);

  int const nb_source = source_mesh_wrapper.num_owned_cells();

  //Set smoothing lengths
  double const one_third = 1./3.;
  Wonton::Point<2> const default_point(one_third, one_third);
  std::vector<std::vector<double>> const default_lengths(1, std::vector<double>(2, one_third));

  Wonton::vector<std::vector<std::vector<double>>> smoothing_lengths(nb_source, default_lengths);
  Wonton::vector<Wonton::Point<2>> source_extents(nb_source, default_point);
  Wonton::vector<Wonton::Point<2>> target_extents(1, default_point);
  Wonton::vector<Weight::Kernel> kernel_types(nb_source, Weight::B4);
  Wonton::vector<Weight::Geometry> geom_types(nb_source, Weight::ELLIPTIC);

  // Source and target swarms
  Swarm<2> source_swarm(source_mesh_wrapper, Wonton::CELL);
  Swarm<2> target_swarm(target_mesh_wrapper, Wonton::CELL);
  SwarmState<2> source_state(source_swarm);
  SwarmState<2> target_state(target_swarm);

  // Create an integer source data for given function
  Wonton::vector<int> source_data_int(nb_source);
  Wonton::vector<double> source_data_dbl(nb_source);

  // Fill the source state data with the specified profile
  for (int i = 0; i < nb_source; ++i) {
    auto p = source_swarm.get_particle_coordinates(i);
    source_data_int[i] = static_cast<int>(p[0] * p[1]) * 100;
  }

  for (int i = 0; i < nb_source; ++i) {
    auto p = source_swarm.get_particle_coordinates(i);
    source_data_dbl[i] = p[0] * p[1];
  }

  source_state.add_field("intdata", source_data_int);
  source_state.add_field("dbldata", source_data_dbl);

  // Build the target state storage
  int const nb_target = target_swarm.num_particles(Wonton::ALL);
  Wonton::vector<int> target_data_int(nb_target);
  Wonton::vector<double> target_data_dbl(nb_target);

  target_state.add_field("intdata", target_data_int);
  target_state.add_field("dbldata", target_data_dbl);

  // Distribute
  Portage::MPI_Particle_Distribute<2> distributor(&executor);
  distributor.distribute(source_swarm, source_state, target_swarm,
                         target_state, smoothing_lengths,
                         source_extents, target_extents, kernel_types,
                         geom_types, WeightCenter::Scatter);

  // Check number of particles received
  int const nb_source_after = source_swarm.num_particles(Wonton::ALL);
  ASSERT_EQ(nb_source_after, nb_source + 5);

  // check coordinates
  auto& source_data_after_int = source_state.get_field_int("intdata");
  auto& source_data_after_dbl = source_state.get_field_dbl("dbldata");

  for (int i = 0 ; i < nb_source_after; ++i) {
    auto coords = source_swarm.get_particle_coordinates(i);
    double const value = coords[0] * coords[1];
    ASSERT_EQ(source_data_after_int[i], static_cast<int>(value) * 100);
    ASSERT_DOUBLE_EQ(source_data_after_dbl[i], value);
  }
}

TEST(MPI_Particle_Distribute, SimpleTest3DGather) {

  using Wonton::Point;
  using namespace Portage::Meshfree;

  // set MPI info
  int rank = 0;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  Wonton::MPIExecutor_type executor(comm);

  // Create a distributed jali source/target mesh
  Jali::MeshFactory mf(comm);
  mf.partitioner(Jali::Partitioner_type::BLOCK);
  auto source_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 4, 4);
  auto target_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 4, 4);

  Wonton::Jali_Mesh_Wrapper source_mesh_wrapper(*source_mesh);
  Wonton::Jali_Mesh_Wrapper target_mesh_wrapper(*target_mesh);

  // Source and target swarms
  Swarm<3> source_swarm(source_mesh_wrapper, Wonton::CELL);
  Swarm<3> target_swarm(target_mesh_wrapper, Wonton::CELL);
  SwarmState<3> source_state(source_swarm);
  SwarmState<3> target_state(target_swarm);

  int const nb_source = source_swarm.num_particles(Wonton::ALL);
  int const nb_target = target_swarm.num_particles(Wonton::ALL);

  // Create an integer source data for given function
  Wonton::vector<int> source_data_int(nb_source);

  // Fill the source state data with the specified profile
  for (int i = 0; i < nb_source; ++i) {
    auto p = source_swarm.get_particle_coordinates(i);
    source_data_int[i] = static_cast<int>(p[0] * p[1] * p[2]) * 1000;
  }

  // Create a double source data for given function
  Wonton::vector<double> source_data_dbl(nb_source);

  // Fill the source state data with the specified profile
  for (int i = 0; i < nb_source; ++i) {
    auto p = source_swarm.get_particle_coordinates(i);
    source_data_dbl[i] = p[0] * p[1] * p[2];
  }

  source_state.add_field("intdata", source_data_int);
  source_state.add_field("dbldata", source_data_dbl);

  // Build the target state storage
  Wonton::vector<int> target_data_int(nb_target, 0);
  Wonton::vector<double> target_data_dbl(nb_target, 0.0);

  target_state.add_field("intdata", target_data_int);
  target_state.add_field("dbldata", target_data_dbl);

  //Set smoothing lengths
  double const one_third = 1./3.;
  Wonton::Point<3> const default_point(one_third, one_third, one_third);
  std::vector<std::vector<double>> const default_lengths(1, std::vector<double>(3, one_third));

  Wonton::vector<std::vector<std::vector<double>>> smoothing_lengths(nb_target, default_lengths);
  Wonton::vector<Wonton::Point<3>> source_extents(1, default_point);
  Wonton::vector<Wonton::Point<3>> target_extents(nb_target, default_point);
  Wonton::vector<Weight::Kernel> kernel_types(nb_target, Weight::B4);
  Wonton::vector<Weight::Geometry> geom_types(nb_target, Weight::ELLIPTIC);

  // Distribute
  Portage::MPI_Particle_Distribute<3> distributor(&executor);
  distributor.distribute(source_swarm, source_state,
                         target_swarm, target_state,
                         smoothing_lengths, source_extents, target_extents,
                         kernel_types, geom_types, WeightCenter::Gather);

  // Check number of particles received
  int const nb_source_after = source_swarm.num_particles(Wonton::ALL);
  ASSERT_EQ(nb_source_after, nb_source + 20);

  // check coordinates
  auto& source_data_after_int = source_state.get_field_int("intdata");
  auto& source_data_after_dbl = source_state.get_field_dbl("dbldata");

  for (int i = 0 ; i < nb_source_after; ++i) {
    auto coords = source_swarm.get_particle_coordinates(i);
    double const value = coords[0] * coords[1] * coords[2];
    ASSERT_EQ(source_data_after_int[i], static_cast<int>(value) * 1000);
    ASSERT_DOUBLE_EQ(source_data_after_dbl[i], value);
  }
}

TEST(MPI_Particle_Distribute, SimpleTest3DScatter) {

  using namespace Portage::Meshfree;

  // set MPI info
  int rank = 0;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  Wonton::MPIExecutor_type executor(comm);

  // Create a distributed jali source/target mesh
  Jali::MeshFactory mf(comm);
  mf.partitioner(Jali::Partitioner_type::BLOCK);
  auto source_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 4, 4);
  auto target_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 4, 4);

  Wonton::Jali_Mesh_Wrapper source_mesh_wrapper(*source_mesh);
  Wonton::Jali_Mesh_Wrapper target_mesh_wrapper(*target_mesh);

  int const nb_source = source_mesh_wrapper.num_owned_cells();

  //Set smoothing lengths
  double const one_third = 1./3.;
  Wonton::Point<3> const default_point(one_third, one_third, one_third);
  std::vector<std::vector<double>> const default_lengths(1, std::vector<double>(3, one_third));

  Wonton::vector<std::vector<std::vector<double>>> smoothing_lengths(nb_source, default_lengths);
  Wonton::vector<Wonton::Point<3>> source_extents(nb_source, default_point);
  Wonton::vector<Wonton::Point<3>> target_extents(1, default_point);
  Wonton::vector<Weight::Kernel> kernel_types(nb_source, Weight::B4);
  Wonton::vector<Weight::Geometry> geom_types(nb_source, Weight::ELLIPTIC);

  // Source and target swarms
  Swarm<3> source_swarm(source_mesh_wrapper, Wonton::CELL);
  Swarm<3> target_swarm(target_mesh_wrapper, Wonton::CELL);
  SwarmState<3> source_state(source_swarm);
  SwarmState<3> target_state(target_swarm);

  // Create an integer source data for given function
  Wonton::vector<int> source_data_int(nb_source);
  Wonton::vector<double> source_data_dbl(nb_source);

  // Fill the source state data with the specified profile
  for (int i = 0; i < nb_source; ++i) {
    auto p = source_swarm.get_particle_coordinates(i);
    source_data_int[i] = static_cast<int>(p[0] * p[1] * p[2]) * 1000;
  }

  for (int i = 0; i < nb_source; ++i) {
    auto p = source_swarm.get_particle_coordinates(i);
    source_data_dbl[i] = p[0] * p[1] * p[2];
  }

  source_state.add_field("intdata", source_data_int);
  source_state.add_field("dbldata", source_data_dbl);

  // Build the target state storage
  int const nb_target = target_swarm.num_particles(Wonton::ALL);
  Wonton::vector<int> target_data_int(nb_target);
  Wonton::vector<double> target_data_dbl(nb_target);

  target_state.add_field("intdata", target_data_int);
  target_state.add_field("dbldata", target_data_dbl);

  // Distribute
  Portage::MPI_Particle_Distribute<3> distributor(&executor);
  distributor.distribute(source_swarm, source_state, target_swarm,
                         target_state, smoothing_lengths,
                         source_extents, target_extents, kernel_types,
                         geom_types, WeightCenter::Scatter);

  // Check number of particles received
  int const nb_source_after = source_swarm.num_particles(Wonton::ALL);
  ASSERT_EQ(nb_source_after, nb_source + 20);

  // check coordinates
  auto& source_data_after_int = source_state.get_field_int("intdata");
  auto& source_data_after_dbl = source_state.get_field_dbl("dbldata");

  for (int i = 0 ; i < nb_source_after; ++i) {
    auto coords = source_swarm.get_particle_coordinates(i);
    double const value = coords[0] * coords[1] * coords[2];
    ASSERT_EQ(source_data_after_int[i], static_cast<int>(value) * 1000);
    ASSERT_DOUBLE_EQ(source_data_after_dbl[i], value);
  }
}

TEST(MPI_Particle_Distribute, SimpleTest2DFaceted) {

  using namespace Portage::Meshfree;

  // set MPI info
  int rank = 0;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  Wonton::MPIExecutor_type executor(comm);

  // Create a distributed jali source/target mesh
  Jali::MeshFactory mf(comm);
  mf.partitioner(Jali::Partitioner_type::BLOCK);
  auto source_mesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4);
  auto target_mesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4);

  Wonton::Jali_Mesh_Wrapper source_mesh_wrapper(*source_mesh);
  Wonton::Jali_Mesh_Wrapper target_mesh_wrapper(*target_mesh);

  int const nb_source = source_mesh_wrapper.num_owned_cells();

  //Set smoothing lengths, etc
  double const one_third = 1./3.;
  Wonton::Point<2> const default_point(one_third, one_third);
  std::vector<std::vector<double>> facets {{ 0,-1, one_third}, { 1, 0, one_third},
                                           { 0, 1, one_third}, {-1, 0, one_third}};

  Wonton::vector<std::vector<std::vector<double>>> smoothing_lengths(nb_source, facets);
  Wonton::vector<Wonton::Point<2>> source_extents(nb_source, default_point);
  Wonton::vector<Wonton::Point<2>> target_extents(1, default_point);
  Wonton::vector<Weight::Kernel> kernel_types(nb_source, Weight::POLYRAMP);
  Wonton::vector<Weight::Geometry> geom_types(nb_source, Weight::FACETED);

  // Source and target swarms
  Swarm<2> source_swarm(source_mesh_wrapper, Wonton::CELL);
  Swarm<2> target_swarm(target_mesh_wrapper, Wonton::CELL);
  SwarmState<2> source_state(source_swarm);
  SwarmState<2> target_state(target_swarm);

  // Create an integer source data for given function
  Wonton::vector<int> source_data_int(nb_source);
  Wonton::vector<double> source_data_dbl(nb_source);

  // Fill the source state data with the specified profile
  for (int i = 0; i < nb_source; ++i) {
    auto p = source_swarm.get_particle_coordinates(i);
    source_data_int[i] = static_cast<int>(p[0] * p[1]) * 100;
  }

  for (int i = 0; i < nb_source; ++i) {
    auto p = source_swarm.get_particle_coordinates(i);
    source_data_dbl[i] = p[0] * p[1];
  }

  source_state.add_field("intdata", source_data_int);
  source_state.add_field("dbldata", source_data_dbl);

  // Build the target state storage
  int const nb_target = target_swarm.num_particles(Wonton::ALL);
  Wonton::vector<int> target_data_int(nb_target);
  Wonton::vector<double> target_data_dbl(nb_target);

  target_state.add_field("intdata", target_data_int);
  target_state.add_field("dbldata", target_data_dbl);

  // Distribute
  Portage::MPI_Particle_Distribute<2> distributor(&executor);
  distributor.distribute(source_swarm, source_state, target_swarm,
                         target_state, smoothing_lengths,
                         source_extents, target_extents, kernel_types,
                         geom_types, WeightCenter::Scatter);

  // Check number of particles received
  int const nb_source_after = source_swarm.num_particles(Wonton::ALL);
  ASSERT_EQ(nb_source_after, nb_source + 5);

  // Check coordinates
  auto& source_data_after_int = source_state.get_field_int("intdata");
  auto& source_data_after_dbl = source_state.get_field_dbl("dbldata");

  for (int i = 0 ; i < nb_source_after; ++i) {
    auto coords = source_swarm.get_particle_coordinates(i);
    std::vector<std::vector<double>> lengths = smoothing_lengths[i];

    for (int k = 0; k < 4; k++) {
      for (int d = 0 ; d < 3; ++d) {
        ASSERT_EQ(lengths[k][d], facets[k][d]);
      }
    }

    Wonton::Point<2> ext = source_extents[i];
    ASSERT_EQ(ext[0], one_third);
    ASSERT_EQ(ext[1], one_third);
    ASSERT_EQ(kernel_types[i], Weight::POLYRAMP);
    ASSERT_EQ(geom_types[i], Weight::FACETED);

    double const value = coords[0] * coords[1];
    ASSERT_EQ(source_data_after_int[i], static_cast<int>(value) * 100);
    ASSERT_DOUBLE_EQ(source_data_after_dbl[i], value);
  }
}

TEST(MPI_Particle_Distribute, SimpleTest3DFaceted) {

  using namespace Portage::Meshfree;

  // set MPI info
  int rank = 0;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  Wonton::MPIExecutor_type executor(comm);

  // Create a distributed jali source/target mesh
  Jali::MeshFactory mf(comm);
  mf.partitioner(Jali::Partitioner_type::BLOCK);
  auto source_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 4, 4);
  auto target_mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 4, 4);

  Wonton::Jali_Mesh_Wrapper source_mesh_wrapper(*source_mesh);
  Wonton::Jali_Mesh_Wrapper target_mesh_wrapper(*target_mesh);

  int const nb_source = source_mesh_wrapper.num_owned_cells();

  //Set smoothing lengths, etc
  double const one_third = 1./3.;
  Wonton::Point<3> const default_point(one_third, one_third, one_third);
  std::vector<std::vector<double>> facets {{ 0,-1, 0, one_third}, { 1, 0, 0, one_third},
                                           { 0, 1, 0, one_third}, {-1, 0, 0, one_third},
                                           { 0, 0,-1, one_third}, { 0, 0, 1, one_third}};

  Wonton::vector<std::vector<std::vector<double>>> smoothing_lengths(nb_source, facets);
  Wonton::vector<Wonton::Point<3>> source_extents(nb_source, default_point);
  Wonton::vector<Wonton::Point<3>> target_extents(1, default_point);
  Wonton::vector<Weight::Kernel> kernel_types(nb_source, Weight::POLYRAMP);
  Wonton::vector<Weight::Geometry> geom_types(nb_source, Weight::FACETED);

  // Source and target swarms
  Swarm<3> source_swarm(source_mesh_wrapper, Wonton::CELL);
  Swarm<3> target_swarm(target_mesh_wrapper, Wonton::CELL);
  SwarmState<3> source_state(source_swarm);
  SwarmState<3> target_state(target_swarm);

  // Create an integer source data for given function
  Wonton::vector<int> source_data_int(nb_source);
  Wonton::vector<double> source_data_dbl(nb_source);

  // Fill the source state data with the specified profile
  for (int i = 0; i < nb_source; ++i) {
    auto p = source_swarm.get_particle_coordinates(i);
    source_data_int[i] = static_cast<int>(p[0] * p[1] * p[2]) * 100;
  }

  for (int i = 0; i < nb_source; ++i) {
    auto p = source_swarm.get_particle_coordinates(i);
    source_data_dbl[i] = p[0] * p[1] * p[2];
  }

  source_state.add_field("intdata", source_data_int);
  source_state.add_field("dbldata", source_data_dbl);

  // Build the target state storage
  int const nb_target = target_swarm.num_particles(Wonton::ALL);
  Wonton::vector<int> target_data_int(nb_target);
  Wonton::vector<double> target_data_dbl(nb_target);

  target_state.add_field("intdata", target_data_int);
  target_state.add_field("dbldata", target_data_dbl);

  // Distribute
  Portage::MPI_Particle_Distribute<3> distributor(&executor);
  distributor.distribute(source_swarm, source_state, target_swarm,
                         target_state, smoothing_lengths,
                         source_extents, target_extents, kernel_types,
                         geom_types, WeightCenter::Scatter);

  // Check number of particles received
  int const nb_source_after = source_swarm.num_particles(Wonton::ALL);
  ASSERT_EQ(nb_source_after, nb_source + 20);

  // Check coordinates
  auto& source_data_after_int = source_state.get_field_int("intdata");
  auto& source_data_after_dbl = source_state.get_field_dbl("dbldata");

  for (int i = 0 ; i < nb_source_after; ++i) {
    auto coords = source_swarm.get_particle_coordinates(i);
    std::vector<std::vector<double>> lengths = smoothing_lengths[i];

    for (int k = 0; k < 6; k++) {
      for (int d = 0 ; d < 4; ++d) {
        ASSERT_EQ(lengths[k][d], facets[k][d]);
      }
    }

    Wonton::Point<3> ext = source_extents[i];
    ASSERT_EQ(ext[0], one_third);
    ASSERT_EQ(ext[1], one_third);
    ASSERT_EQ(ext[2], one_third);
    ASSERT_EQ(kernel_types[i], Weight::POLYRAMP);
    ASSERT_EQ(geom_types[i], Weight::FACETED);

    double const value = coords[0] * coords[1] * coords[2];
    ASSERT_EQ(source_data_after_int[i], static_cast<int>(value) * 100);
    ASSERT_DOUBLE_EQ(source_data_after_dbl[i], value);
  }
}
