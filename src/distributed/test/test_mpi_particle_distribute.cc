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
//#include "portage/wonton/state/flat/flat_state_wrapper.h"

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

  // Create an integer source data for given function
  const int nsrcpts = source_swarm.num_particles(Portage::Entity_type::ALL);
  typename Portage::Meshfree::SwarmState<2>::IntVecPtr source_data_int = 
       std::make_shared<typename Portage::Meshfree::SwarmState<2>::IntVec>(nsrcpts);

  // Fill the source state data with the specified profile
  for (size_t p = 0; p < nsrcpts; ++p) {
    (*source_data_int)[p] = p*p;
  }
  source_state->add_field("intdata", source_data_int);
  
  // Create a double source data for given function
  typename Portage::Meshfree::SwarmState<2>::DblVecPtr source_data_dbl = 
       std::make_shared<typename Portage::Meshfree::SwarmState<2>::DblVec>(nsrcpts);

  // Fill the source state data with the specified profile
  for (size_t p = 0; p < nsrcpts; ++p) {
    (*source_data_dbl)[p] = p*p;
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

   //Print out source swarm before distribution
   typename Portage::Meshfree::SwarmState<2>::IntVecPtr sd_int_before = 
        std::make_shared<typename Portage::Meshfree::SwarmState<2>::IntVec>(nsrcpts);
   source_state->get_field("intdata",sd_int_before);
   
   typename Portage::Meshfree::SwarmState<2>::DblVecPtr sd_dbl_before = 
        std::make_shared<typename Portage::Meshfree::SwarmState<2>::DblVec>(nsrcpts);
   source_state->get_field("dbldata",sd_dbl_before);

   typename Portage::Meshfree::SwarmState<2>::IntVecPtr td_int_before = 
        std::make_shared<typename Portage::Meshfree::SwarmState<2>::IntVec>(ntarpts);
   target_state->get_field("intdata",td_int_before);

   typename Portage::Meshfree::SwarmState<2>::DblVecPtr td_dbl_before = 
        std::make_shared<typename Portage::Meshfree::SwarmState<2>::DblVec>(ntarpts);
   target_state->get_field("dbldata",td_dbl_before);

   if (commRank == 1) {
   for (size_t p = 0 ; p < nsrcpts; ++p)
   {
	Portage::Point<2> coords = source_swarm.get_particle_coordinates(p);
        std::cout<<"Rank-"<<commRank<<"::SOURCE-before::coords = [ "<<coords[0]<<", "<<coords[1]<<" ] "<<std::endl;
        std::cout<<"Rank-"<<commRank<<"::SOURCE-before::intdata = "<<(*sd_int_before)[p]<<std::endl;
        std::cout<<"Rank-"<<commRank<<"::SOURCE-before::dbldata = "<<(*sd_dbl_before)[p]<<std::endl;
   }

   for (size_t p = 0 ; p < ntarpts; ++p)
   {
	Portage::Point<2> coords = target_swarm.get_particle_coordinates(p);
        std::cout<<"Rank-"<<commRank<<"::TARGET-before::coords = [ "<<coords[0]<<", "<<coords[1]<<" ] "<<std::endl;
        std::cout<<"Rank-"<<commRank<<"::TARGET-before::intdata = "<<(*td_int_before)[p]<<std::endl;
        std::cout<<"Rank-"<<commRank<<"::TARGET-before::dbldata = "<<(*td_dbl_before)[p]<<std::endl;
   }
  }
  //Set smoothing lengths 
   auto smoothing_lengths = std::vector<std::vector<std::vector<double>>>(4*4,
                   std::vector<std::vector<double>>(1, std::vector<double>(2, 1.0/3)));

  //Distribute 
  Portage::MPI_Particle_Distribute<2> distributor;
  distributor.distribute(source_swarm, *source_state, target_swarm,
                         *target_state, smoothing_lengths, 
                         Portage::Meshfree::WeightCenter::Gather);

  // Check number of particles received
  
   //Print out source swarm before distribution
   int nsrcpts_after = source_swarm.num_particles(Portage::Entity_type::ALL);
   typename Portage::Meshfree::SwarmState<2>::IntVecPtr sd_int_after = 
        std::make_shared<typename Portage::Meshfree::SwarmState<2>::IntVec>(nsrcpts_after);
   source_state->get_field("intdata",sd_int_after);

   typename Portage::Meshfree::SwarmState<2>::DblVecPtr sd_dbl_after = 
        std::make_shared<typename Portage::Meshfree::SwarmState<2>::DblVec>(nsrcpts_after);
   source_state->get_field("dbldata",sd_dbl_after);
   
   int ntarpts_after = target_swarm.num_particles(Portage::Entity_type::ALL);
   typename Portage::Meshfree::SwarmState<2>::IntVecPtr td_int_after = 
        std::make_shared<typename Portage::Meshfree::SwarmState<2>::IntVec>(ntarpts_after);
   target_state->get_field("intdata",td_int_after);

   typename Portage::Meshfree::SwarmState<2>::DblVecPtr td_dbl_after = 
        std::make_shared<typename Portage::Meshfree::SwarmState<2>::DblVec>(ntarpts_after);
   target_state->get_field("dbldata",td_dbl_after);
   
   if (commRank == 1 ){
   for (size_t p = 0 ; p < nsrcpts_after; ++p)
   {
	Portage::Point<2> coords = source_swarm.get_particle_coordinates(p);
        std::cout<<"Rank-"<<commRank<<"::SOURCE-after::coords = [ "<<coords[0]<<", "<<coords[1]<<" ] "<<std::endl;
        std::cout<<"Rank-"<<commRank<<"::SOURCE-after::intdata = "<<(*sd_int_after)[p]<<std::endl;
        std::cout<<"Rank-"<<commRank<<"::SOURCE-after::dbldata = "<<(*sd_dbl_after)[p]<<std::endl;
   }

   for (size_t p = 0 ; p < ntarpts_after; ++p)
   {
	Portage::Point<2> coords = target_swarm.get_particle_coordinates(p);
        std::cout<<"Rank-"<<commRank<<"::TARGET-after::coords = [ "<<coords[0]<<", "<<coords[1]<<" ] "<<std::endl;
        std::cout<<"Rank-"<<commRank<<"::TARGET-after::intdata = "<<(*td_int_after)[p]<<std::endl;
        std::cout<<"Rank-"<<commRank<<"::TARGET-after::dbldata = "<<(*td_dbl_after)[p]<<std::endl;
   }
  }
}

