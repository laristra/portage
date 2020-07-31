/*
 * This file is part of the Ristra portage project.
 * Please see the license file at the root of this repository, or at:
 * https://github.com/laristra/portage/blob/master/LICENSE
 */


#include "gtest/gtest.h"

// wonton
#include "wonton/support/wonton.h"
#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"

// portage
#include "portage/support/portage.h"
#include "portage/driver/coredriver.h"
#include "portage/distributed/mpi_update_ghosts.h"
#include "portage/search/search_kdtree.h"
#include "portage/intersect/intersect_r2d.h"
#include "portage/interpolate/interpolate_2nd_order.h"

// jali
#include "Mesh.hh"
#include "MeshFactory.hh"

TEST(GhostManager, CommunicationMatrices) {

  int rank = 0;
  int num_ranks = 1;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &num_ranks);

  ASSERT_GE(num_ranks, 1);

  auto jali_mesh  = Jali::MeshFactory(comm)(0.0, 0.0, 1.0, 1.0, 7, 6);
  auto jali_state = Jali::State::create(jali_mesh);
  Wonton::Jali_Mesh_Wrapper mesh(*jali_mesh);
  Wonton::Jali_State_Wrapper state(*jali_state);

  using GhostManager = Portage::MPI_GhostManager<Wonton::Jali_Mesh_Wrapper,
                                                 Wonton::Jali_State_Wrapper,
                                                 Wonton::CELL>;

  std::vector<std::vector<int>> owned[num_ranks]; /* owned cells that are ghost for some rank */
  std::vector<std::vector<int>> ghost[num_ranks]; /* ghost cells that are owned by some rank */
  std::vector<int> num_owned[num_ranks];          /* number of owned cells sent by each rank */
  std::vector<int> num_ghost[num_ranks];          /* number of ghost cells receive from each rank */
  std::vector<MPI_Request> requests;              /* list of asynchronous MPI requests */

  GhostManager ghost_manager(mesh, state, comm);
  owned[rank] = ghost_manager.owned_comm_matrix();
  ghost[rank] = ghost_manager.ghost_comm_matrix();

  for (int i = 0; i < num_ranks; ++i) {
    num_owned[rank].emplace_back(owned[rank][i].size());
    num_ghost[rank].emplace_back(ghost[rank][i].size());
    // store GID for owned cells to ease verification
    for (int& id : owned[rank][i]) {
      id = mesh.get_global_id(id, Wonton::CELL);
    }
  }

  if (rank > 0) {
    MPI_Request request;
    MPI_Isend(num_owned[rank].data(), num_ranks, MPI_INT, 0, rank, comm, &request);
    requests.emplace_back(request);
    MPI_Isend(num_ghost[rank].data(), num_ranks, MPI_INT, 0, rank, comm, &request);
    requests.emplace_back(request);
  } else {
    for (int i = 1; i < num_ranks; ++i) {
      MPI_Request request;
      // send
      num_owned[i].resize(num_ranks, 0);
      MPI_Irecv(num_owned[i].data(), num_ranks, MPI_INT, i, i, comm, &request);
      requests.emplace_back(request);
      // receive
      num_ghost[i].resize(num_ranks, 0);
      MPI_Irecv(num_ghost[i].data(), num_ranks, MPI_INT, i, i, comm, &request);
      requests.emplace_back(request);
    }
  }

  MPI_Waitall(requests.size(), requests.data(), MPI_STATUS_IGNORE);

  // check that the count of sent owned cells matches
  // the count of received ghost cells for each rank pair.
  if (rank == 0) {
    for (int i = 0; i < num_ranks; ++i) {
      for (int j = 0; j < num_ranks; ++j) {
        ASSERT_EQ(num_owned[i][j], num_ghost[j][i]);
      }
    }
  }

  requests.clear();

  if (rank > 0) {
    int const offset = num_ranks * num_ranks;
    for (int i = 0; i < num_ranks; ++i) {
      if (i != rank) {
        int tag = num_ranks * (rank - 1) + i;
        MPI_Request request;
        MPI_Isend(owned[rank][i].data(), owned[rank][i].size(), MPI_INT, 0, tag, comm, &request);
        requests.emplace_back(request);
        MPI_Isend(ghost[rank][i].data(), ghost[rank][i].size(), MPI_INT, 0, tag + offset, comm, &request);
        requests.emplace_back(request);
      }
    }
  } else {
    int const offset = num_ranks * num_ranks;
    for (int i = 1; i < num_ranks; ++i) {
      owned[i].resize(num_ranks);
      ghost[i].resize(num_ranks);
      for (int j = 0; j < num_ranks; ++j) {
        if (i != j) {
          int tag = num_ranks * (i - 1) + j;
          MPI_Request request;
          // send
          owned[i][j].resize(num_owned[i][j]);
          MPI_Irecv(owned[i][j].data(), num_owned[i][j], MPI_INT, i, tag, comm, &request);
          requests.emplace_back(request);
          // receive
          ghost[i][j].resize(num_ghost[i][j]);
          MPI_Irecv(ghost[i][j].data(), num_ghost[i][j], MPI_INT, i, tag + offset, comm, &request);
          requests.emplace_back(request);
        }
      }
    }
  }

  MPI_Waitall(requests.size(), requests.data(), MPI_STATUS_IGNORE);

  // check that we have a perfect matching on
  // sent and received cells for each rank pair.
  if (rank == 0) {
    for (int i = 0; i < num_ranks; ++i) {
      for (int j = 0; j < num_ranks; ++j) {
        if (i == j) { continue; }
        int const num_sent_owned = owned[i][j].size();
        int const num_recv_ghost = ghost[j][i].size();
        ASSERT_EQ(num_sent_owned, num_recv_ghost);
        for (int k = 0; k < num_sent_owned; ++k) {
          ASSERT_EQ(owned[i][j][k], ghost[j][i][k]);
        }
      }
    }
  }

  MPI_Barrier(comm);
}

TEST(GhostManager, UpdateValuesSingleMat) {

  using GhostManager = Portage::MPI_GhostManager<Wonton::Jali_Mesh_Wrapper,
                                                 Wonton::Jali_State_Wrapper,
                                                 Wonton::CELL>;

  using Remapper = Portage::CoreDriver<2, Wonton::CELL,
                                       Wonton::Jali_Mesh_Wrapper,
                                       Wonton::Jali_State_Wrapper>;

  int rank = 0;
  int num_ranks = 1;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &num_ranks);
  Wonton::MPIExecutor_type executor(comm);

  auto jali_source_mesh  = Jali::MeshFactory(comm)(0.0, 0.0, 1.0, 1.0, 5, 5);
  auto jali_target_mesh  = Jali::MeshFactory(comm)(0.0, 0.0, 1.0, 1.0, 7, 6);
  auto jali_source_state = Jali::State::create(jali_source_mesh);
  auto jali_target_state = Jali::State::create(jali_target_mesh);

  Wonton::Jali_Mesh_Wrapper  source_mesh(*jali_source_mesh);
  Wonton::Jali_Mesh_Wrapper  target_mesh(*jali_target_mesh);
  Wonton::Jali_State_Wrapper source_state(*jali_source_state);
  Wonton::Jali_State_Wrapper target_state(*jali_target_state);

  // set temperature fields on states
  int const num_source_cells = source_mesh.num_entities(Wonton::CELL, Wonton::ALL);

  double source_field[num_source_cells];
  for (int i = 0; i < num_source_cells; ++i) {
    Wonton::Point<2> centroid;
    source_mesh.cell_centroid(i, &centroid);
    source_field[i] = centroid[0] + 2 * centroid[1];
  }

  source_state.mesh_add_data<double>(Wonton::CELL, "temperature", source_field);
  target_state.mesh_add_data<double>(Wonton::CELL, "temperature", 0.0);

  // remap
  Remapper remapper(source_mesh, source_state, target_mesh, target_state, &executor);

  auto candidates = remapper.search<Portage::SearchKDTree>();
  auto weights    = remapper.intersect_meshes<Portage::IntersectR2D>(candidates);
  auto gradients  = remapper.compute_source_gradient("temperature");

  remapper.interpolate_mesh_var<double, Portage::Interpolate_2ndOrder>(
    "temperature", "temperature", weights, &gradients
  );

  // fill ghost values on target mesh
  GhostManager ghost_manager(target_mesh, target_state, comm);
  ghost_manager.update_ghost_values({"temperature"});

}

TEST(GhostManager, UpdateValuesMultiMat) {

}