/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <iostream>
#include "portage/accumulate/accumulate.h"
#include "portage/distributed/mpi_particle_distribute.h"
#include "portage/driver/driver_swarm.h"
#include "portage/support/Point.h"
#include "portage/swarm/swarm.h"
#include "portage/swarm/swarm_state.h"
#include "portage/wonton/mesh/flat/flat_mesh_wrapper.h"
#include "portage/wonton/mesh/jali/jali_mesh_wrapper.h"

#include "gtest/gtest.h"
#include "mpi.h"

#include "Mesh.hh"
#include "MeshFactory.hh"

double TOL = 1e-6;

TEST(SwarmDriver, Test2D) {

  int commRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &commRank);

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  // Create a distributed jali source/target mesh 
  std::shared_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 1.0, 1.0, 8, 8);
  Wonton::Jali_Mesh_Wrapper jali_smesh_wrapper(*source_mesh);

  std::shared_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4);
  Wonton::Jali_Mesh_Wrapper jali_tmesh_wrapper(*target_mesh);

  //Create flat mesh wrappers for source/target jali meshes
  Wonton::Flat_Mesh_Wrapper<> source_mesh_flat;
  source_mesh_flat.initialize(jali_smesh_wrapper);

  Wonton::Flat_Mesh_Wrapper<> target_mesh_flat;
  target_mesh_flat.initialize(jali_tmesh_wrapper);

  // Source and target swarms 
  Portage::Meshfree::Swarm<2> source_swarm(source_mesh_flat, Portage::Entity_kind::CELL);
  Portage::Meshfree::Swarm<2> target_swarm(target_mesh_flat, Portage::Entity_kind::CELL);

  // Source and target mesh state
  std::shared_ptr<Portage::Meshfree::SwarmState<2>> source_state;
  std::shared_ptr<Portage::Meshfree::SwarmState<2>> target_state;
  source_state = std::make_shared<Portage::Meshfree::SwarmState<2>>(source_swarm);
  target_state = std::make_shared<Portage::Meshfree::SwarmState<2>>(target_swarm);

  // Create a double source data for given function
  const int nsrcpts = source_swarm.num_particles(Portage::Entity_type::PARALLEL_OWNED);
  typename Portage::Meshfree::SwarmState<2>::DblVecPtr source_data_dbl = 
       std::make_shared<typename Portage::Meshfree::SwarmState<2>::DblVec>(nsrcpts);

  // Fill the source state data with the specified profile
  for (size_t p = 0; p < nsrcpts; ++p) {
   Portage::Point<2> coord =
          source_swarm.get_particle_coordinates(p);
    (*source_data_dbl)[p] = coord[0]*coord[0]+coord[1]*coord[1];
  }
  source_state->add_field("particledata", source_data_dbl);

  // Build the target state storage
  const int ntarpts = target_swarm.num_particles(Portage::Entity_type::PARALLEL_OWNED);
  typename Portage::Meshfree::SwarmState<2>::DblVecPtr target_data_dbl = 
        std::make_shared<typename Portage::Meshfree::SwarmState<2>::DblVec>(ntarpts, 0.0);
    target_state->add_field("particledata", target_data_dbl);

  //Set smoothing lengths 
  auto smoothing_lengths = std::vector<std::vector<std::vector<double>>>(ntarpts,
                   std::vector<std::vector<double>>(1, std::vector<double>(2, 2.0/4)));
 
  // Build the swarm driver
  // Register the variable name and interpolation order with the driver
  std::vector<std::string> remap_fields;
  remap_fields.push_back("particledata");

  Portage::Meshfree::SwarmDriver<Portage::SearchSimplePoints,
                                 Portage::Meshfree::Accumulate,
                                 Portage::Meshfree::Estimate,
                                 2,
                                 Portage::Meshfree::Swarm<2>,
                                 Portage::Meshfree::SwarmState<2>,
                                 Portage::Meshfree::Swarm<2>,
                                 Portage::Meshfree::SwarmState<2>>
        d(source_swarm, *source_state, target_swarm, *target_state, smoothing_lengths,
          Portage::Meshfree::Weight::B4, Portage::Meshfree::Weight::ELLIPTIC, Portage::Meshfree::WeightCenter::Gather);
  
  d.set_remap_var_names(remap_fields, remap_fields,
                          Portage::Meshfree::LocalRegression,
                          Portage::Meshfree::Basis::Quadratic);

  // run on multiple processor
  d.run(true);

  // Check the answer
  double stdval, err;
  double toterr=0.;

  typename Portage::Meshfree::SwarmState<2>::DblVecPtr vecout;
  target_state->get_field("particledata", vecout);
  ASSERT_NE(nullptr, vecout);

  for (int p = 0; p < ntarpts; ++p) {
    Portage::Point<2> coord = target_swarm.get_particle_coordinates(p);
    double error;
    error = coord[0]*coord[0]+coord[1]*coord[1] - (*vecout)[p];
    std::printf("Particle=% 4d Coord = (% 5.3lf,% 5.3lf)", p, coord[0],
                    coord[1]);
    std::printf("  Value = % 10.6lf  Err = % lf\n", (*vecout)[p], error);
    toterr += error*error;
   }

  std::printf("\n\nL2 NORM OF ERROR = %lf\n\n", sqrt(toterr));
  ASSERT_NEAR(0.0, sqrt(toterr), TOL);
  
}
