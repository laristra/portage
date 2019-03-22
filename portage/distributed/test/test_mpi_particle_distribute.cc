/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <iostream>
#include "portage/accumulate/accumulate.h"
#include "portage/distributed/mpi_particle_distribute.h"

#include "portage/support/portage.h"
#include "portage/swarm/swarm.h"
#include "portage/swarm/swarm_state.h"
#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/support/Point.h"

#include "gtest/gtest.h"
#include "mpi.h"

#include "Mesh.hh"
#include "MeshFactory.hh"



TEST(MPI_Particle_Distribute, SimpleTest2DGather) {
  using Wonton::Point;
  int commRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &commRank);

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  // Create a distributed jali source/target mesh 
  std::shared_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4);
  Wonton::Jali_Mesh_Wrapper jali_smesh_wrapper(*source_mesh);

  std::shared_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4);
  Wonton::Jali_Mesh_Wrapper jali_tmesh_wrapper(*target_mesh);

  // Source and target swarms 
  std::shared_ptr<Portage::Meshfree::Swarm<2>> source_swarm_ptr =
      Portage::Meshfree::SwarmFactory<2>(jali_smesh_wrapper, Portage::Entity_kind::CELL);
  std::shared_ptr<Portage::Meshfree::Swarm<2>> target_swarm_ptr =
      Portage::Meshfree::SwarmFactory<2>(jali_tmesh_wrapper, Portage::Entity_kind::CELL);
  Portage::Meshfree::Swarm<2> &source_swarm(*source_swarm_ptr);
  Portage::Meshfree::Swarm<2> &target_swarm(*target_swarm_ptr);

  // Source and target mesh state
  std::shared_ptr<Portage::Meshfree::SwarmState<2>> source_state;
  std::shared_ptr<Portage::Meshfree::SwarmState<2>> target_state;
  source_state = std::make_shared<Portage::Meshfree::SwarmState<2>>(source_swarm);
  target_state = std::make_shared<Portage::Meshfree::SwarmState<2>>(target_swarm);

  // Create an integer source data for given function
  const int nsrcpts = source_swarm.num_particles(Portage::Entity_type::ALL);
  typename Portage::Meshfree::SwarmState<2>::IntVecPtr source_data_int = 
      std::make_shared<typename Portage::Meshfree::SwarmState<2>::IntVec>(nsrcpts);

  // Fill the source state data with the specified profile
  for (size_t p = 0; p < nsrcpts; ++p) {
    Point<2> coords = source_swarm.get_particle_coordinates(p);
    (*source_data_int)[p] = (int)(coords[0]*coords[1]*100);
  }
  source_state->add_field("intdata", source_data_int);
  
  // Create a double source data for given function
  typename Portage::Meshfree::SwarmState<2>::DblVecPtr source_data_dbl = 
      std::make_shared<typename Portage::Meshfree::SwarmState<2>::DblVec>(nsrcpts);

  // Fill the source state data with the specified profile
  for (size_t p = 0; p < nsrcpts; ++p) {
    Point<2> coords = source_swarm.get_particle_coordinates(p);
    (*source_data_dbl)[p] = coords[0]*coords[1];
  }
  source_state->add_field("dbldata", source_data_dbl);

  // Build the target state storage
  const int ntarpts = target_swarm.num_particles(Portage::Entity_type::ALL);
  typename Portage::Meshfree::SwarmState<2>::IntVecPtr target_data_int = 
      std::make_shared<typename Portage::Meshfree::SwarmState<2>::IntVec>(ntarpts, 0);
  target_state->add_field("intdata", target_data_int);

  // Build the target state storage
  typename Portage::Meshfree::SwarmState<2>::DblVecPtr target_data_dbl = 
      std::make_shared<typename Portage::Meshfree::SwarmState<2>::DblVec>(ntarpts, 0.0);
  target_state->add_field("dbldata", target_data_dbl);

  //Set smoothing lengths 
  auto smoothing_lengths = std::make_shared<Portage::vector<std::vector<std::vector<double>>>>(ntarpts,
                                                                                               std::vector<std::vector<double>>(1, std::vector<double>(2, 1.0/3)));

  auto kernel_types = std::make_shared<Portage::vector<Portage::Meshfree::Weight::Kernel>>
      (ntarpts, Portage::Meshfree::Weight::B4);
   
  auto geom_types = std::make_shared<Portage::vector<Portage::Meshfree::Weight::Geometry>>
      (ntarpts, Portage::Meshfree::Weight::ELLIPTIC);

  // Distribute
  Wonton::MPIExecutor_type executor(MPI_COMM_WORLD);
  Portage::MPI_Particle_Distribute<2> distributor(&executor);
  distributor.distribute(source_swarm, *source_state, 
			 target_swarm, *target_state, 
                         *smoothing_lengths, *kernel_types,
			 *geom_types, Portage::Meshfree::WeightCenter::Gather);

  // Check number of particles received
  int nsrcpts_after = source_swarm.num_particles(Portage::Entity_type::ALL);
  typename Portage::Meshfree::SwarmState<2>::IntVecPtr sd_int_after = 
      std::make_shared<typename Portage::Meshfree::SwarmState<2>::IntVec>(nsrcpts_after);
  source_state->get_field("intdata",sd_int_after);

  typename Portage::Meshfree::SwarmState<2>::DblVecPtr sd_dbl_after = 
      std::make_shared<typename Portage::Meshfree::SwarmState<2>::DblVec>(nsrcpts_after);
  source_state->get_field("dbldata",sd_dbl_after);
   
  ASSERT_EQ(nsrcpts_after, nsrcpts+5);  

  for (size_t p = 0 ; p < nsrcpts_after; ++p)
  { 
    Point<2> coords = source_swarm.get_particle_coordinates(p);
    ASSERT_EQ((*sd_int_after)[p],(int)(coords[0]*coords[1]*100));  
    ASSERT_EQ((*sd_dbl_after)[p],(coords[0]*coords[1]));  
  }
}

TEST(MPI_Particle_Distribute, SimpleTest2DScatter) {
  using Wonton::Point;
  int commRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
  
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  
  // Create a distributed jali source/target mesh 
  std::shared_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4);
  Wonton::Jali_Mesh_Wrapper jali_smesh_wrapper(*source_mesh);
  
  std::shared_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4);
  Wonton::Jali_Mesh_Wrapper jali_tmesh_wrapper(*target_mesh);

  //Set smoothing lengths 
  int nsrcpts = jali_smesh_wrapper.num_owned_cells(); 
  auto smoothing_lengths = std::make_shared<Portage::vector<std::vector<std::vector<double>>>>(nsrcpts,
                                                                                               std::vector<std::vector<double>>(1, std::vector<double>(2, 1.0/3)));
  
  auto kernel_types = std::make_shared<Portage::vector<Portage::Meshfree::Weight::Kernel>>
      (nsrcpts, Portage::Meshfree::Weight::B4);
   
  auto geom_types = std::make_shared<Portage::vector<Portage::Meshfree::Weight::Geometry>>
      (nsrcpts, Portage::Meshfree::Weight::ELLIPTIC);

  // Source and target swarms 
  std::shared_ptr<Portage::Meshfree::Swarm<2>> source_swarm_ptr =
      Portage::Meshfree::SwarmFactory<2>(jali_smesh_wrapper, Portage::Entity_kind::CELL);
  std::shared_ptr<Portage::Meshfree::Swarm<2>> target_swarm_ptr =
      Portage::Meshfree::SwarmFactory<2>(jali_tmesh_wrapper, Portage::Entity_kind::CELL);
  Portage::Meshfree::Swarm<2> &source_swarm(*source_swarm_ptr);
  Portage::Meshfree::Swarm<2> &target_swarm(*target_swarm_ptr);

  // Source and target mesh state
  std::shared_ptr<Portage::Meshfree::SwarmState<2>> source_state;
  std::shared_ptr<Portage::Meshfree::SwarmState<2>> target_state;
  source_state = std::make_shared<Portage::Meshfree::SwarmState<2>>(source_swarm);
  target_state = std::make_shared<Portage::Meshfree::SwarmState<2>>(target_swarm);

  // Create an integer source data for given function
  typename Portage::Meshfree::SwarmState<2>::IntVecPtr source_data_int = 
      std::make_shared<typename Portage::Meshfree::SwarmState<2>::IntVec>(nsrcpts);

  // Fill the source state data with the specified profile
  for (size_t p = 0; p < nsrcpts; ++p) {
    Point<2> coords = source_swarm.get_particle_coordinates(p);
    (*source_data_int)[p] = (int)(coords[0]*coords[1]*100);
  }
  source_state->add_field("intdata", source_data_int);
  
  // Create a double source data for given function
  typename Portage::Meshfree::SwarmState<2>::DblVecPtr source_data_dbl = 
      std::make_shared<typename Portage::Meshfree::SwarmState<2>::DblVec>(nsrcpts);

  // Fill the source state data with the specified profile
  for (size_t p = 0; p < nsrcpts; ++p) {
    Point<2> coords = source_swarm.get_particle_coordinates(p);
    (*source_data_dbl)[p] = coords[0]*coords[1];
  }
  source_state->add_field("dbldata", source_data_dbl);

  // Build the target state storage
  const int ntarpts = target_swarm.num_particles(Portage::Entity_type::ALL);
  typename Portage::Meshfree::SwarmState<2>::IntVecPtr target_data_int = 
      std::make_shared<typename Portage::Meshfree::SwarmState<2>::IntVec>(ntarpts, 0);
  target_state->add_field("intdata", target_data_int);

  // Build the target state storage
  typename Portage::Meshfree::SwarmState<2>::DblVecPtr target_data_dbl = 
      std::make_shared<typename Portage::Meshfree::SwarmState<2>::DblVec>(ntarpts, 0.0);
  target_state->add_field("dbldata", target_data_dbl);
    
  // Distribute
  Wonton::MPIExecutor_type executor(MPI_COMM_WORLD);
  Portage::MPI_Particle_Distribute<2> distributor(&executor);
  distributor.distribute(source_swarm, *source_state, target_swarm,
                         *target_state, *smoothing_lengths, *kernel_types,
                         *geom_types, Portage::Meshfree::WeightCenter::Scatter);
   
  // Check number of particles received
  int nsrcpts_after = source_swarm.num_particles(Portage::Entity_type::ALL);
  typename Portage::Meshfree::SwarmState<2>::IntVecPtr sd_int_after = 
      std::make_shared<typename Portage::Meshfree::SwarmState<2>::IntVec>(nsrcpts_after);
  source_state->get_field("intdata",sd_int_after);

  typename Portage::Meshfree::SwarmState<2>::DblVecPtr sd_dbl_after = 
      std::make_shared<typename Portage::Meshfree::SwarmState<2>::DblVec>(nsrcpts_after);
  source_state->get_field("dbldata",sd_dbl_after);
   
  ASSERT_EQ(nsrcpts_after, nsrcpts+5);  

  for (size_t p = 0 ; p < nsrcpts_after; ++p)
  { 
    Point<2> coords = source_swarm.get_particle_coordinates(p);
    std::vector<std::vector<double>> smlen = (*smoothing_lengths)[p];
    for (size_t d = 0 ; d < 2; ++d)
      ASSERT_EQ(smlen[0][d],1.0/3);
    ASSERT_EQ((*sd_int_after)[p],(int)(coords[0]*coords[1]*100));  
    ASSERT_EQ((*sd_dbl_after)[p],(coords[0]*coords[1]));  
  }
}


TEST(MPI_Particle_Distribute, SimpleTest3DGather) {
  using Wonton::Point;
  int commRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &commRank);

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  // Create a distributed jali source/target mesh 
  std::shared_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 0.0, 
                                               1.0, 1.0, 1.0, 
                                               4, 4, 4);
  Wonton::Jali_Mesh_Wrapper jali_smesh_wrapper(*source_mesh);

  std::shared_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 0.0,
                                               1.0, 1.0, 1.0,
                                               4, 4, 4);
  Wonton::Jali_Mesh_Wrapper jali_tmesh_wrapper(*target_mesh);

  // Source and target swarms 
  std::shared_ptr<Portage::Meshfree::Swarm<3>> source_swarm_ptr =
      Portage::Meshfree::SwarmFactory<3>(jali_smesh_wrapper, Portage::Entity_kind::CELL);
  std::shared_ptr<Portage::Meshfree::Swarm<3>> target_swarm_ptr =
      Portage::Meshfree::SwarmFactory<3>(jali_tmesh_wrapper, Portage::Entity_kind::CELL);
  Portage::Meshfree::Swarm<3> &source_swarm(*source_swarm_ptr);
  Portage::Meshfree::Swarm<3> &target_swarm(*target_swarm_ptr);

  // Source and target mesh state
  std::shared_ptr<Portage::Meshfree::SwarmState<3>> source_state;
  std::shared_ptr<Portage::Meshfree::SwarmState<3>> target_state;
  source_state = std::make_shared<Portage::Meshfree::SwarmState<3>>(source_swarm);
  target_state = std::make_shared<Portage::Meshfree::SwarmState<3>>(target_swarm);

  // Create an integer source data for given function
  const int nsrcpts = source_swarm.num_particles(Portage::Entity_type::ALL);
  typename Portage::Meshfree::SwarmState<3>::IntVecPtr source_data_int = 
      std::make_shared<typename Portage::Meshfree::SwarmState<3>::IntVec>(nsrcpts);

  // Fill the source state data with the specified profile
  for (size_t p = 0; p < nsrcpts; ++p) {
    Point<3> coords = source_swarm.get_particle_coordinates(p);
    (*source_data_int)[p] = (int)(coords[0]*coords[1]*coords[2]*1000);
  }
  source_state->add_field("intdata", source_data_int);
  
  // Create a double source data for given function
  typename Portage::Meshfree::SwarmState<3>::DblVecPtr source_data_dbl = 
      std::make_shared<typename Portage::Meshfree::SwarmState<3>::DblVec>(nsrcpts);

  // Fill the source state data with the specified profile
  for (size_t p = 0; p < nsrcpts; ++p) {
    Point<3> coords = source_swarm.get_particle_coordinates(p);
    (*source_data_dbl)[p] = coords[0]*coords[1]*coords[2];
  }
  source_state->add_field("dbldata", source_data_dbl);

  // Build the target state storage
  const int ntarpts = target_swarm.num_particles(Portage::Entity_type::ALL);
  typename Portage::Meshfree::SwarmState<3>::IntVecPtr target_data_int = 
      std::make_shared<typename Portage::Meshfree::SwarmState<3>::IntVec>(ntarpts, 0);
  target_state->add_field("intdata", target_data_int);

  // Build the target state storage
  typename Portage::Meshfree::SwarmState<3>::DblVecPtr target_data_dbl = 
      std::make_shared<typename Portage::Meshfree::SwarmState<3>::DblVec>(ntarpts, 0.0);
  target_state->add_field("dbldata", target_data_dbl);

  //Set smoothing lengths 
  auto smoothing_lengths = std::make_shared<Portage::vector<std::vector<std::vector<double>>>>
      (ntarpts, std::vector<std::vector<double>>(1, std::vector<double>(3, 1.0/3)));

  auto kernel_types = std::make_shared<Portage::vector<Portage::Meshfree::Weight::Kernel>>
      (ntarpts, Portage::Meshfree::Weight::B4);
   
  auto geom_types = std::make_shared<Portage::vector<Portage::Meshfree::Weight::Geometry>>
      (ntarpts, Portage::Meshfree::Weight::ELLIPTIC);
 
  // Distribute
  Wonton::MPIExecutor_type executor(MPI_COMM_WORLD);
  Portage::MPI_Particle_Distribute<3> distributor(&executor);
  distributor.distribute(source_swarm, *source_state, 
			 target_swarm, *target_state, *smoothing_lengths,
                         *kernel_types, *geom_types, 
                         Portage::Meshfree::WeightCenter::Gather);

  // Check number of particles received
  int nsrcpts_after = source_swarm.num_particles(Portage::Entity_type::ALL);
  typename Portage::Meshfree::SwarmState<3>::IntVecPtr sd_int_after = 
      std::make_shared<typename Portage::Meshfree::SwarmState<3>::IntVec>(nsrcpts_after);
  source_state->get_field("intdata",sd_int_after);

  typename Portage::Meshfree::SwarmState<3>::DblVecPtr sd_dbl_after = 
      std::make_shared<typename Portage::Meshfree::SwarmState<3>::DblVec>(nsrcpts_after);
  source_state->get_field("dbldata",sd_dbl_after);
   
  ASSERT_EQ(nsrcpts_after, nsrcpts+20);  

  for (size_t p = 0 ; p < nsrcpts_after; ++p)
  { 
    Point<3> coords = source_swarm.get_particle_coordinates(p);
    ASSERT_EQ((*sd_int_after)[p],(int)(coords[0]*coords[1]*coords[2]*1000));  
    ASSERT_EQ((*sd_dbl_after)[p],(coords[0]*coords[1]*coords[2]));  
  }
}

TEST(MPI_Particle_Distribute, SimpleTest3DScatter) {
  using Wonton::Point;
  int commRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &commRank);

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  // Create a distributed jali source/target mesh 
  std::shared_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 0.0, 
                                               1.0, 1.0, 1.0, 
                                               4, 4, 4);
  Wonton::Jali_Mesh_Wrapper jali_smesh_wrapper(*source_mesh);

  std::shared_ptr<Jali::Mesh> target_mesh = mf(0.0, 0.0, 0.0,
                                               1.0, 1.0, 1.0,
                                               4, 4, 4);
  Wonton::Jali_Mesh_Wrapper jali_tmesh_wrapper(*target_mesh);
  
  //Set smoothing lengths 
  int nsrcpts = jali_smesh_wrapper.num_owned_cells();
  auto smoothing_lengths = std::make_shared<Portage::vector<std::vector<std::vector<double>>>>(nsrcpts,
                                                                                               std::vector<std::vector<double>>(1, std::vector<double>(3, 1.0/3)));

  auto kernel_types = std::make_shared<Portage::vector<Portage::Meshfree::Weight::Kernel>>
      (nsrcpts, Portage::Meshfree::Weight::B4);
   
  auto geom_types = std::make_shared<Portage::vector<Portage::Meshfree::Weight::Geometry>>
      (nsrcpts, Portage::Meshfree::Weight::ELLIPTIC);

  // Source and target swarms 
  std::shared_ptr<Portage::Meshfree::Swarm<3>> source_swarm_ptr =
      Portage::Meshfree::SwarmFactory<3>(jali_smesh_wrapper, Portage::Entity_kind::CELL);
  std::shared_ptr<Portage::Meshfree::Swarm<3>> target_swarm_ptr =
      Portage::Meshfree::SwarmFactory<3>(jali_tmesh_wrapper, Portage::Entity_kind::CELL);
  Portage::Meshfree::Swarm<3> &source_swarm(*source_swarm_ptr);
  Portage::Meshfree::Swarm<3> &target_swarm(*target_swarm_ptr);

  // Source and target mesh state
  std::shared_ptr<Portage::Meshfree::SwarmState<3>> source_state;
  std::shared_ptr<Portage::Meshfree::SwarmState<3>> target_state;
  source_state = std::make_shared<Portage::Meshfree::SwarmState<3>>(source_swarm);
  target_state = std::make_shared<Portage::Meshfree::SwarmState<3>>(target_swarm);

  // Create an integer source data for given function
  typename Portage::Meshfree::SwarmState<3>::IntVecPtr source_data_int = 
      std::make_shared<typename Portage::Meshfree::SwarmState<3>::IntVec>(nsrcpts);

  // Fill the source state data with the specified profile
  for (size_t p = 0; p < nsrcpts; ++p) {
    Point<3> coords = source_swarm.get_particle_coordinates(p);
    (*source_data_int)[p] = (int)(coords[0]*coords[1]*coords[2]*1000);
  }
  source_state->add_field("intdata", source_data_int);
  
  // Create a double source data for given function
  typename Portage::Meshfree::SwarmState<3>::DblVecPtr source_data_dbl = 
      std::make_shared<typename Portage::Meshfree::SwarmState<3>::DblVec>(nsrcpts);

  // Fill the source state data with the specified profile
  for (size_t p = 0; p < nsrcpts; ++p) {
    Point<3> coords = source_swarm.get_particle_coordinates(p);
    (*source_data_dbl)[p] = coords[0]*coords[1]*coords[2];
  }
  source_state->add_field("dbldata", source_data_dbl);

  // Build the target state storage
  const int ntarpts = target_swarm.num_particles(Portage::Entity_type::ALL);
  typename Portage::Meshfree::SwarmState<3>::IntVecPtr target_data_int = 
      std::make_shared<typename Portage::Meshfree::SwarmState<3>::IntVec>(ntarpts, 0);
  target_state->add_field("intdata", target_data_int);

  // Build the target state storage
  typename Portage::Meshfree::SwarmState<3>::DblVecPtr target_data_dbl = 
      std::make_shared<typename Portage::Meshfree::SwarmState<3>::DblVec>(ntarpts, 0.0);
  target_state->add_field("dbldata", target_data_dbl);

  // Distribute
  Wonton::MPIExecutor_type executor(MPI_COMM_WORLD);
  Portage::MPI_Particle_Distribute<3> distributor(&executor);
  distributor.distribute(source_swarm, *source_state, target_swarm,
                         *target_state, *smoothing_lengths,
                         *kernel_types, *geom_types, 
                         Portage::Meshfree::WeightCenter::Scatter);

  // Check number of particles received
  int nsrcpts_after = source_swarm.num_particles(Portage::Entity_type::ALL);
  typename Portage::Meshfree::SwarmState<3>::IntVecPtr sd_int_after = 
      std::make_shared<typename Portage::Meshfree::SwarmState<3>::IntVec>(nsrcpts_after);
  source_state->get_field("intdata",sd_int_after);

  typename Portage::Meshfree::SwarmState<3>::DblVecPtr sd_dbl_after = 
      std::make_shared<typename Portage::Meshfree::SwarmState<3>::DblVec>(nsrcpts_after);
  source_state->get_field("dbldata",sd_dbl_after);
   
  ASSERT_EQ(nsrcpts_after, nsrcpts+20);  

  for (size_t p = 0 ; p < nsrcpts_after; ++p)
  { 
    Point<3> coords = source_swarm.get_particle_coordinates(p);
    std::vector<std::vector<double>> smlen = (*smoothing_lengths)[p];
    for (size_t d = 0 ; d < 3; ++d)
      ASSERT_EQ(smlen[0][d],1.0/3);
    ASSERT_EQ((*sd_int_after)[p],(int)(coords[0]*coords[1]*coords[2]*1000));  
    ASSERT_EQ((*sd_dbl_after)[p],(coords[0]*coords[1]*coords[2]));  
  }
}

