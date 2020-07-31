/*
 * This file is part of the Ristra portage project.
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
 * @brief A communication manager to fill values on ghost cells.
 *
 * @tparam Mesh: the mesh type.
 * @tparam State: its state type.
 * @tparam entity: the entity kind.
 */
template<typename Mesh, typename State, Wonton::Entity_kind entity>
class MPI_GhostManager {

public:
  /**
   * @brief Create and initialize the manager.
   *
   * @param in_mesh: the given mesh.
   * @param in_state: its state.
   * @param in_comm: the MPI communicator.
   */
  MPI_GhostManager(Mesh const& in_mesh, State& in_state, MPI_Comm in_comm)
    : mesh(in_mesh), state(in_state), comm(in_comm)
  {
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &num_ranks);

    if (num_ranks > 1) {
      send.resize(num_ranks);
      receive.resize(num_ranks);
      count.resize(num_ranks);
      cache_comm_matrices();
    }
  }

  /**
   * @brief Delete the manager.
   *
   */
  ~MPI_GhostManager() = default;

  /**
   * @brief Send and store field values on ghost cells.
   *
   * @param fields: list of remapped fields.
   */
  void update_ghost_values(std::vector<std::string> const& fields) {

    for (auto&& field : fields) {
#ifdef PORTAGE_HAS_TANGRAM
      auto field_type = state.field_type(entity, field);
      bool multimat = (entity == Wonton::CELL and field_type == Field_type::MULTIMATERIAL_FIELD);

      if (multimat) {
        // todo
      } else /* single material */{
#endif
        std::vector<double> buffer[num_ranks];
        std::vector<MPI_Request> requests;

        // step 1: retrieve field values and send them
        double* values = nullptr;
        state.mesh_get_data(entity, field, &values);
        assert(values != nullptr);

        for (int i = 0; i < num_ranks; ++i) {
          if (i != rank and not send[i].empty()) {
            buffer[i].clear();
            for (auto&& j : send[i]) {  // 'j' local entity index
              assert(not is_ghost(j));
              buffer[i].emplace_back(values[j]);
            }
            MPI_Request request;
            MPI_Isend(buffer[i].data(), buffer[i].size(), MPI_DOUBLE, i, 0, comm, &request);
            requests.emplace_back(request);
          }
        }

        // step 2: receive field data
        for (int i = 0; i < num_ranks; ++i) {
          if (i != rank and count[i]) {
            buffer[i].resize(count[i]);
            MPI_Request request;
            MPI_Irecv(buffer[i].data(), count[i], MPI_DOUBLE, i, 0, comm, &request);
            requests.emplace_back(request);
          }
        }

        // step 3: update state
        MPI_Waitall(requests.size(), requests.data(), status);

        for (int i = 0; i < num_ranks; ++i) {
          if (i != rank and count[i]) {
            for (int j = 0; j < count[i]; ++j) {
              int const& gid = receive[i][j];
              int const& lid = ghost[gid];
              values[lid] = buffer[i][j];
            }
            buffer[i].clear();
          }
        }

#ifdef PORTAGE_HAS_TANGRAM
      } // if not multimat
#endif
      MPI_Barrier(comm);
    }
  }

  /**
   * @brief Retrieve the matrix of owned entities to send to each rank.
   *
   * @return the matrix of owned entities to send to each rank.
   */
  std::vector<std::vector<int>> const& owned_comm_matrix() { return send; }

  /**
   * @brief Retrieve the matrix of ghost entities to receive from each rank.
   *
   * @return the matrix of ghost entities to receive from each rank.
   */
  std::vector<std::vector<int>> const& ghost_comm_matrix() { return receive; }

private:
  /**
   * @brief Build and store communication matrices to identify
   *        which owned entities data should be sent to each rank,
   *        and which ghost entities data should be received from
   *        each rank.
   */
  void cache_comm_matrices() {

    std::vector<int> entities;
    int num_ghosts[num_ranks];
    int offsets[num_ranks];
    std::vector<int> received;
    std::vector<int> gids[num_ranks];

    // step 1: retrieve ghost entities on current rank
    int const num_entities = mesh.num_entities(entity, Wonton::ALL);

    for (int i = 0; i < num_entities; ++i) {
      int const& gid = mesh.get_global_id(i, entity);
      if (is_ghost(i)) {
        entities.emplace_back(gid);
        ghost[gid] = i;
      } else {
        owned[gid] = i;
      }
    }

    // step 2: gather number of ghosts for all ranks and deduce offsets
    num_ghosts[rank] = entities.size();
    MPI_Allgather(num_ghosts + rank, 1, MPI_INT, num_ghosts, 1, MPI_INT, comm);

    int const total_ghosts = std::accumulate(num_ghosts, num_ghosts + num_ranks, 0);
    received.resize(total_ghosts);

    int index = 0;
    for (int i = 0; i < num_ranks; ++i) {
      offsets[i] = index;
      index += num_ghosts[i];
    }

    // step 3: gather all ghost entities on all ranks
    MPI_Allgatherv(entities.data(), num_ghosts[rank], MPI_INT, received.data(), num_ghosts, offsets, MPI_INT, comm);

    // step 4: check received ghost cells, build send matrix and send entities count.
    int num_sent[num_ranks];
    std::vector<MPI_Request> requests;

    for (int i = 0; i < num_ranks; ++i) {
      if (i != rank) {
        int const& start = offsets[i];
        int const& extent = (i < num_ranks - 1 ? offsets[i+1] : total_ghosts);
        for (int j = start; j < extent; ++j) {
          int const& gid = received[j];
          if (owned.count(gid)) {
            send[i].emplace_back(owned[gid]);
            gids[i].emplace_back(gid);
          }
        }

        MPI_Request request;
        num_sent[i] = send[i].size();
        MPI_Isend(num_sent + i, 1, MPI_INT, i, 0, comm, &request);
        requests.emplace_back(request);
      }
    }

#ifdef DEBUG
  for (int i = 0; i < num_ranks; ++i) {
    if (i != rank) {
      std::cout << "[" << rank << "] -> [" << i << "]: (";
      int const num_to_send = send[i].size();
      for (int j = 0; j < num_to_send; ++j) {
        std::cout << mesh.get_global_id(send[i][j], entity);
        if (j < num_to_send - 1) {
          std::cout << ", ";
        }
      }
      std::cout << ")" << std::endl;
    }
  }
#endif

    // step 5: receive ghost count per rank
    for (int i = 0; i < num_ranks; ++i) {
      if (i != rank) {
        MPI_Request request;
        MPI_Irecv(count.data() + i, 1, MPI_INT, i, 0, comm, &request);
        requests.emplace_back(request);
      }
    }

    // step 6: send owned entities, receive expected ghosts.
    for (int i = 0; i < num_ranks; ++i) {
      if (i != rank and num_sent[i]) {
        MPI_Request request;
        MPI_Isend(gids[i].data(), num_sent[i], MPI_INT, i, 1, comm, &request);
        requests.emplace_back(request);
      }
    }

    MPI_Waitall(requests.size(), requests.data(), status);
    requests.clear();

    for (int i = 0; i < num_ranks; ++i) {
      if (i != rank and count[i]) {
        MPI_Request request;
        receive[i].resize(count[i]);
        MPI_Irecv(receive[i].data(), count[i], MPI_INT, i, 1, comm, &request);
        requests.emplace_back(request);
      }
    }

    MPI_Waitall(requests.size(), requests.data(), status);

#ifdef DEBUG
    for (int i = 0; i < num_ranks; ++i) {
      if (i != rank) {
        std::cout << "[" << rank << "] <- [" << i << "]: (";
        for (int j = 0; j < count[i]; ++j) {
          std::cout << receive[i][j];
          if (j < count[i] - 1) {
            std::cout << ", ";
          }
        }
        std::cout << ")" << std::endl;
      }
    }
#endif

    // verification
    int total_received = 0;
    for (int i = 0; i < num_ranks; ++i) {
      total_received += receive[i].size();
    }
    assert(total_received == num_ghosts[rank]);
    MPI_Barrier(comm);
  }

  /**
   * @brief Check if the given entity is a ghost one.
   *
   * @param i: the index of the entity.
   * @return true if a ghost entity, false otherwise.
   */
  bool is_ghost(int i) const {
    switch (entity) {
      case Wonton::CELL: return mesh.cell_get_type(i) == Wonton::PARALLEL_GHOST;
      case Wonton::NODE: return mesh.node_get_type(i) == Wonton::PARALLEL_GHOST;
      default: return false;
    }
  }

  /** mesh instance */
  Mesh const& mesh;
  /** mesh state */
  State& state;
  /** MPI communicator */
  MPI_Comm comm = MPI_COMM_NULL;
  /** MPI status */
  MPI_Status* status = MPI_STATUS_IGNORE;
  /** current MPI rank */
  int rank = 0;
  /** number of ranks */
  int num_ranks = 1;
  /** list of owned entities to send to each rank */
  std::vector<std::vector<int>> send {};
  /** list of ghost entities received from each rank */
  std::vector<std::vector<int>> receive {};
  /** count of received entities from each rank */
  std::vector<int> count {};
  /** owned entities indexed by their GID */
  std::map<int, int> owned {};
  /** ghost entities indexed by their GID */
  std::map<int, int> ghost {};
};

} // namespace Portage
#endif // ifdef WONTON_ENABLE_MPI
#endif // ifndef PORTAGE_MPI_UPDATE_GHOSTS_H
