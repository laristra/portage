/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <iostream>
#include "portage/accumulate/accumulate.h"
#include "portage/distributed/mpi_particle_distribute.h"
#include "portage/support/Point.h"
#include "portage/swarm/swarm.h"
#include "portage/swarm/swarm_state.h"
#include "portage/wonton/mesh/flat/flat_mesh_wrapper.h"
#include "portage/wonton/mesh/jali/jali_mesh_wrapper.h"
#include "portage/wonton/state/flat/flat_state_wrapper.h"

#include "gtest/gtest.h"
#include "mpi.h"

#include "Mesh.hh"
#include "MeshFactory.hh"

TEST(MPI_Particle_Distribute, SimpleTest2D) {

  int commRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &commRank);

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  // Create a distributed jali source/target mesh 
  std::shared_ptr<Jali::Mesh> source_mesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4);
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

  // Fill the source state data with the specified profile
  const int nsrcpts = source_swarm.num_particles(Portage::Entity_type::ALL);
  typename Portage::Meshfree::SwarmState<2>::IntVecPtr source_data = 
       std::make_shared<typename Portage::Meshfree::SwarmState<2>::IntVec>(nsrcpts);

  // Create the source data for given function
  for (size_t p = 0; p < nsrcpts; ++p) {
    (*source_data)[p] = 1;
  }
  source_state->add_field("particledata", source_data);

  // Build the target state storage
  const int ntarpts = target_swarm.num_particles(Portage::Entity_type::ALL);
  typename Portage::Meshfree::SwarmState<2>::IntVecPtr target_data = 
        std::make_shared<typename Portage::Meshfree::SwarmState<2>::IntVec>(ntarpts, 0);
    target_state->add_field("particledata", target_data);

   //Print out source swarm before distribution
   typename Portage::Meshfree::SwarmState<2>::IntVecPtr sd_before = 
        std::make_shared<typename Portage::Meshfree::SwarmState<2>::IntVec>(nsrcpts);
   source_state->get_field("particledata",sd_before);

   for (size_t p = 0 ; p < nsrcpts; ++p)
   {
	Portage::Point<2> coords = source_swarm.get_particle_coordinates(p);
  //      std::cout<<"Rank-"<<commRank<<"::SOURCE-before::coords = [ "<<coords[0]<<", "<<coords[1]<<" ] "<<std::endl;
  //      std::cout<<"Rank-"<<commRank<<"::SOURCE-before::particledata = "<<(*sd_before)[p]<<std::endl;
   }

   typename Portage::Meshfree::SwarmState<2>::IntVecPtr td_before = 
        std::make_shared<typename Portage::Meshfree::SwarmState<2>::IntVec>(ntarpts);
   target_state->get_field("particledata",td_before);

   for (size_t p = 0 ; p < ntarpts; ++p)
   {
	Portage::Point<2> coords = target_swarm.get_particle_coordinates(p);
    //    std::cout<<"Rank-"<<commRank<<"::TARGET-before::coords = [ "<<coords[0]<<", "<<coords[1]<<" ] "<<std::endl;
    //    std::cout<<"Rank-"<<commRank<<"::TARGET-before::particledata = "<<(*td_before)[p]<<std::endl;
   }

  //Set smoothing lengths 
   auto smoothing_lengths = std::vector<std::vector<std::vector<double>>>(4*4,
                   std::vector<std::vector<double>>(1, std::vector<double>(2, 1.0/3)));

  //Distribute 
  Portage::MPI_Particle_Distribute<2> distributor;
  distributor.distribute(source_swarm, source_state, target_swarm,
                         target_state, smoothing_lengths, 
                         Portage::Meshfree::WeightCenter::Gather);

  // Check number of particles received
  
   //Print out source swarm before distribution
   int nsrcpts_after = source_swarm.num_particles(Portage::Entity_type::ALL);
   typename Portage::Meshfree::SwarmState<2>::IntVecPtr sd_after = 
        std::make_shared<typename Portage::Meshfree::SwarmState<2>::IntVec>(nsrcpts_after);
   source_state->get_field("particledata",sd_after);

   for (size_t p = 0 ; p < nsrcpts_after; ++p)
   {
	Portage::Point<2> coords = source_swarm.get_particle_coordinates(p);
      //  std::cout<<"Rank-"<<commRank<<"::SOURCE-after::coords = [ "<<coords[0]<<", "<<coords[1]<<" ] "<<std::endl;
      //  std::cout<<"Rank-"<<commRank<<"::SOURCE-after::particledata = "<<(*sd_after)[p]<<std::endl;
   }

   int ntarpts_after = target_swarm.num_particles(Portage::Entity_type::ALL);
   typename Portage::Meshfree::SwarmState<2>::IntVecPtr td_after = 
        std::make_shared<typename Portage::Meshfree::SwarmState<2>::IntVec>(ntarpts_after);
   target_state->get_field("particledata",td_after);

   for (size_t p = 0 ; p < ntarpts_after; ++p)
   {
	Portage::Point<2> coords = target_swarm.get_particle_coordinates(p);
      //  std::cout<<"Rank-"<<commRank<<"::TARGET-after::coords = [ "<<coords[0]<<", "<<coords[1]<<" ] "<<std::endl;
      //  std::cout<<"Rank-"<<commRank<<"::TARGET-after::particledata = "<<(*td_after)[p]<<std::endl;
   }
}

