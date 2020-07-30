/*
 * This file is part of the Ristra portage project.
 * Please see the license file at the root of this repository, or at:
 * https://github.com/laristra/portage/blob/master/LICENSE
 */

#include "gtest/gtest.h"

// wonton includes
#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"
#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"

// portage includes
#include "portage/support/portage.h"
#include "portage/distributed/mpi_update_ghosts.h"

// Jali includes
#include "Mesh.hh"
#include "MeshFactory.hh"

TEST(GhostManager, CommunicationMatrices) {

  int rank = 0;
  int num_ranks = 1;

  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &num_ranks);

  ASSERT_GE(num_ranks, 1);

  using GhostManager = Portage::MPI_GhostManager<Wonton::Jali_Mesh_Wrapper,
                                                 Wonton::Jali_State_Wrapper,
                                                 Wonton::CELL>;

  auto jali_mesh  = Jali::MeshFactory(comm)(0.0, 0.0, 1.0, 1.0, 7, 6);
  auto jali_state = Jali::State::create(jali_mesh);
  Wonton::Jali_Mesh_Wrapper mesh(*jali_mesh);
  Wonton::Jali_State_Wrapper state(*jali_state);

  GhostManager ghost_manager(mesh, state, comm);

  // step 1: retrieve all communication matrices of each rank to the master.

  std::vector<std::vector<int>> owned[num_ranks]; /* owned cells that are ghost for some rank */
  std::vector<std::vector<int>> ghost[num_ranks]; /* ghost cells that are owned by some rank */
  std::vector<int> num_owned[num_ranks];          /* number of owned cells sent by each rank */
  std::vector<int> num_ghost[num_ranks];          /* number of ghost cells receive from each rank */
  std::vector<MPI_Request> requests;              /* list of asynchronous MPI requests */

  owned[rank] = ghost_manager.owned_comm_matrix();
  ghost[rank] = ghost_manager.ghost_comm_matrix();

  for (int i = 0; i < num_ranks; ++i) {
    num_owned[rank].emplace_back(owned[rank][i].size());
    num_ghost[rank].emplace_back(ghost[rank][i].size());
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

#ifdef DEBUG
  if (rank == 0) {
    for (int i = 0; i < num_ranks; ++i) {
      for (int j = 0; j < num_ranks; ++j) {
        if (i != j) {
          std::cout << "owned: from ["<< i <<"] to ["<< j <<"]: "<< num_owned[i][j] << std::endl;
        }
      }
    }

    for (int i = 0; i < num_ranks; ++i) {
      for (int j = 0; j < num_ranks; ++j) {
        if (i != j) {
          std::cout << "ghost: to ["<< i <<"] from ["<< j <<"]: "<< num_ghost[i][j] << std::endl;
        }
      }
    }
  }
#endif

  requests.clear();

  if (rank > 0) {
    for (int i = 0; i < num_ranks; ++i) {
      if (i != rank) {
        int tag = rank * i;
        MPI_Request request;
        MPI_Isend(owned[rank][i].data(), owned[rank][i].size(), MPI_INT, 0, tag, comm, &request);
        requests.emplace_back(request);
        MPI_Isend(ghost[rank][i].data(), ghost[rank][i].size(), MPI_INT, 0, tag + 1, comm, &request);
        requests.emplace_back(request);
      }
    }
  } else {
    for (int i = 1; i < num_ranks; ++i) {
      owned[i].resize(num_ranks);
      ghost[i].resize(num_ranks);
      for (int j = 1; j < num_ranks; ++j) {
        if (i != j) {
          int tag = i * j;
          MPI_Request request;
          // send
          owned[i][j].resize(num_owned[i][j]);
          MPI_Irecv(owned[i][j].data(), num_owned[i][j], MPI_INT, i, tag, comm, &request);
          requests.emplace_back(request);
          // receive
          ghost[i][j].resize(num_ghost[i][j]);
          MPI_Irecv(ghost[i][j].data(), num_ghost[i][j], MPI_INT, i, tag + 1, comm, &request);
          requests.emplace_back(request);
        }
      }
    }
  }

  MPI_Waitall(requests.size(), requests.data(), MPI_STATUS_IGNORE);
  MPI_Barrier(comm);

  // check that we have a perfect matching between
  // sent and received ghost cells for each pair of ranks.


//  for (int i = 0; i < num_ranks; ++i) {
//    int const count = receives[rank][i].size();
//    int* buffer = num_receives[rank].data() + i;
//    MPI_Gather(&count, 1, MPI_INT, buffer, 1, MPI_INT, 0, comm);
//  }


}