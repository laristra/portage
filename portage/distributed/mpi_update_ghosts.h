/* This file is part of the Ristra portage project.
 * Please see the license file at the root of this repository, or at:
 * https://github.com/laristra/portage/blob/master/LICENSE
 */

#ifndef PORTAGE_MPI_UPDATE_GHOSTS_H
#define PORTAGE_MPI_UPDATE_GHOSTS_H

#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"

#ifdef WONTON_ENABLE_MPI
namespace Portage {

/**
 * @brief Check if the given entity is a ghost one.
 *
 * @tparam Mesh: the given mesh type.
 * @param mesh: the mesh instance.
 * @param entity: the entity kind.
 * @param i: the index of the entity.
 * @return true if a ghost entity, false otherwise.
 */
template<typename Mesh>
bool is_ghost(Mesh const& mesh, Wonton::Entity_kind entity, int i) {
  switch (entity) {
    case Wonton::CELL: return mesh.cell_get_type(i) == Wonton::PARALLEL_GHOST;
    case Wonton::NODE: return mesh.node_get_type(i) == Wonton::PARALLEL_GHOST;
    default: return false;
  }
}

/**
 * @brief Build the communication matrix required to update ghost cells values.
 *
 * @tparam Mesh: the given mesh type.
 * @param mesh: the mesh instance.
 * @param entity: the entity kind.
 * @param comm: the communication channel.
 * @return a sparse matrix that maps the owned entities to send to each rank.
 */
template<typename Mesh>
std::vector<std::vector<int>> build_ghost_comm_matrix(Mesh const& mesh,
                                                      Wonton::Entity_kind entity,
                                                      MPI_Comm comm) {
  int rank = 0;
  int num_ranks = 1;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &num_ranks);

  std::map<int, int> owned;
  std::vector<int> ghosts;
  int num_ghosts[num_ranks];
  int offsets[num_ranks];
  std::vector<int> received;
  std::vector<std::vector<int>> send(num_ranks);

  // step 1: retrieve ghost entities on current rank
  int const num_entities = mesh.num_entities(entity, Wonton::ALL);

  for (int i = 0; i < num_entities; ++i) {
    int const& gid = mesh.get_global_id(i, entity);
    if (is_ghost(mesh, entity, i)) {
      ghosts.emplace_back(gid);
    } else {
      owned[gid] = i;
    }
  }

  // step 2: gather number of ghosts for all ranks and deduce offsets
  num_ghosts[rank] = ghosts.size();
  MPI_Allgather(num_ghosts + rank, 1, MPI_INT, num_ghosts, num_ranks, MPI_INT, comm);

  int const total_ghosts = std::accumulate(num_ghosts, num_ghosts + num_ranks, 0);
  received.resize(total_ghosts);

  int index = 0;
  for (int i = 0; i < num_ranks; ++i) {
    offsets[i] = index;
    index += num_ghosts[i];
  }

  // step 3: gather all ghost entities on all ranks
  MPI_Allgatherv(ghosts.data(), num_ghosts[rank], MPI_INT, received.data(), num_ghosts, offsets, MPI_INT, comm);

  // step 4: check received ghost cells and build map
  for (int i = 0; i < num_ranks; ++i) {
    if (i != rank) {
      int const& start = offsets[i];
      int const& extent = (i < num_ranks - 1 ? offsets[i+1] : total_ghosts);
      for (int j = start; j < extent; ++j) {
        int const& gid = received[j];
        if (owned.count(gid)) {
          send[i].emplace_back(owned[gid]);
        }
      }
    }
  }

#ifdef DEBUG
  for (int i = 0; i < num_ranks; ++i) {
    std::cout << "["<< rank <<"]: received[" << i << "]: (";
    int const& start = offsets[i];
    int const& extent = (i < num_ranks - 1 ? offsets[i+1] : total_ghosts);
    for (int j = start; j < extent; ++j) {
      int const& current = received[j];
      std::cout << current;
      if (j < extent - 1) {
        std::cout << ", ";
      }
    }
    std::cout << ")" << std::endl;
  }

  std::cout << " ------------------ " << std::endl;

  for (int i = 0; i < num_ranks; ++i) {
    std::cout << "["<< rank <<"]: send[" << i << "]: (";
    int const num_to_send = send[i].size();
    for (int j = 0; j < num_to_send; ++j) {
      std::cout << mesh.get_global_id(send[i][j], entity);
      if (j < num_to_send - 1) {
        std::cout << ", ";
      }
    }
    std::cout << ")" << std::endl;
  }

  std::cout << " ------------------ " << std::endl;
#endif
  return send;
}

/**
 * @brief
 *
 * @tparam Mesh
 * @param mesh
 * @param comm_matrix
 */
template<typename Mesh,
         typename State,
         Wonton::Entity_kind onwhat>
bool fill_ghost_values(Mesh& mesh,
                       State& state,
                       std::vector<std::string> const& fields,
                       std::vector<std::vector<int>> const& comm_matrix,
                       MPI_Comm comm) {

  int rank = 0;
  int num_ranks = 1;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &num_ranks);

  if (num_ranks == 1) {
    std::cerr << "Warning: serial run, nothing to do." << std::endl;
    return false;
  }

  for (auto&& field : fields) {
#ifdef PORTAGE_HAS_TANGRAM
    auto field_type = state.field_type(onwhat, field);
    bool multimat = (onwhat == Wonton::CELL and field_type == Field_type::MULTIMATERIAL_FIELD);

    if (multimat) {
      // todo
    } else {
      std::vector<double> send[num_ranks];
      std::vector<double> recv[num_ranks];
      std::vector<int> gids[num_ranks];
      int count[num_ranks];

      // step 1: retrieve field values and populate data
      double* values = nullptr;
      state.mesh_get_data(onwhat, &values);
      assert(values != nullptr);

      for (int i = 0; i < num_ranks; ++i) {
        if (i != rank) {
          for (auto&& j : comm_matrix[i]) {  // 'j' local entity index
            assert(not is_ghost(mesh, onwhat, j));
            send[i].emplace_back(values[j]);
            gids[i].emplace_back(mesh.get_global_id(j, onwhat));
          }
          count[i] = send[i].size();
        }
      }

      // step 2: send them
      for (int i = 0; i < num_ranks; ++i) {
        if (i != rank) {
          MPI_Send(count + i, 1, MPI_INT, i, 0, comm);
          if (count[i] > 0) {
            MPI_Send(send[i].data(), count[i], MPI_DOUBLE, i, 1, comm);
            MPI_Send(gids[i].data(), count[i], MPI_INT, i, 2, comm);
          }
        }
      }

      // step 3: receive
      for (int i = 0; i < num_ranks; ++i) {
        if (i != rank) {
          MPI_Status status;
          MPI_Recv(count + i, 1, MPI_INT, i, 0, comm, &status);
          if (count[i] > 0) {
            recv[i].resize(count[i]);
            gids[i].resize(count[i], -1);  // reset and reuse
            MPI_Recv(recv[i].data(), count[i], MPI_DOUBLE, i, 1, comm, &status);
            MPI_Recv(gids[i].data(), count[i], MPI_INT, i, 2, comm, &status);
          }
        }
      }

      // step 4: store them
      std::map<int, int> ghosts;
      int const num_entities = mesh.num_entities(onwhat, Wonton::ALL);
      for (int i = 0; i < num_entities; ++i) {
        if (is_ghost(mesh, onwhat, i)) {
          int const& gid = mesh.get_global_id(i, onwhat);
          ghosts[gid] = i;
        }
      }

      for (int i = 0; i < num_ranks; ++i) {
        if (i != rank and count[i] > 0) {
          for (int j = 0; j < count[i]; ++j) {
            int const& gid = gids[i][j];
            int const& lid = ghosts[gid];
            values[lid] = recv[i][j];
          }
        }
      }
    }
#endif
  }

  return true;
}

} // namespace Portage
#endif
#endif //PORTAGE_MPI_UPDATE_GHOSTS_H
