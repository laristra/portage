/*
 * This file is part of the Ristra portage project.
 * Please see the license file at the root of this repository, or at:
 * https://github.com/laristra/portage/blob/master/LICENSE
 */

#ifndef PORTAGE_MPI_GHOST_MANAGER_H
#define PORTAGE_MPI_GHOST_MANAGER_H

#include <map>
#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"
#include "portage/support/portage.h"

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

protected:
  /**
   * @struct Data for MPI communications.
   *
   */
  struct Data {
    /** lookup table for local indices */
    std::map<int, int> lookup {};
    /** communication matrices per material */
    std::vector<std::vector<std::vector<int>>> matrix {};
    /** exchanged entities count per material */
    std::vector<std::vector<int>> count {};
    /** cached values per rank for tests */
    std::vector<std::vector<double>> values {};
  };

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
    cache_comm_matrices();
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
   * @param cache_values: whether to cache values or not for this field.
   */
  void update_ghost_values(std::string const& field, bool cache_values = false) {

    if (cache_values) {
      if (owned.values.empty()) { owned.values.resize(num_ranks); }
      if (ghost.values.empty()) { ghost.values.resize(num_ranks); }
    }

#ifdef PORTAGE_HAS_TANGRAM
    auto field_type = state.field_type(entity, field);
    bool multimat = (entity == Wonton::CELL and field_type == Field_type::MULTIMATERIAL_FIELD);

    if (multimat) {
      for (int m = 1; m < num_mats; ++m) {
        std::vector<int> material_cells;
        state.mat_get_cells(m, &material_cells);


      }

    } else /* single material */{
#endif
      std::vector<double> buffer[num_ranks];
      std::vector<MPI_Request> requests;
      int const m = 0;

      // step 1: retrieve field values and send them
      double* data = nullptr;
      state.mesh_get_data(entity, field, &data);
      assert(data != nullptr);

      for (int i = 0; i < num_ranks; ++i) {
        if (i != rank and not owned.matrix[m][i].empty()) {
          buffer[i].clear();
          for (auto&& j : owned.matrix[m][i]) {  // 'j' local entity index
            assert(not is_ghost(j));
            buffer[i].emplace_back(data[j]);
          }
          MPI_Request request;
          MPI_Isend(buffer[i].data(), buffer[i].size(), MPI_DOUBLE, i, 0, comm, &request);
          requests.emplace_back(request);

          if (cache_values) {
            owned.values[i].resize(buffer[i].size());
            std::copy(buffer[i].begin(), buffer[i].end(), owned.values[i].begin());
          }
        }
      }

      // step 2: receive field data
      for (int i = 0; i < num_ranks; ++i) {
        if (i != rank and ghost.count[m][i]) {
          buffer[i].resize(ghost.count[m][i]);
          MPI_Request request;
          MPI_Irecv(buffer[i].data(), ghost.count[m][i], MPI_DOUBLE, i, 0, comm, &request);
          requests.emplace_back(request);
        }
      }

      // step 3: update state
      MPI_Waitall(requests.size(), requests.data(), status);

      for (int i = 0; i < num_ranks; ++i) {
        if (i != rank and ghost.count[m][i]) {
          for (int j = 0; j < ghost.count[m][i]; ++j) {
            int const& gid = ghost.matrix[m][i][j];
            int const& lid = ghost.lookup[gid];
            data[lid] = buffer[i][j];
          }
          if (cache_values) {
            ghost.values[i].resize(ghost.count[m][i]);
            std::copy(buffer[i].begin(), buffer[i].end(), ghost.values[i].begin());
          }
          buffer[i].clear();
        }
      }
#ifdef PORTAGE_HAS_TANGRAM
    } // if not multimat
#endif
    MPI_Barrier(comm);
  }

  /**
   * @brief Retrieve the list of owned entities to send to each rank.
   *
   * @param m: material index.
   * @return matrix of communication.
   */
  std::vector<std::vector<int>> const& owned_matrix(int m = 0) const { return owned.matrix[m]; }

  /**
   * @brief Retrieve the list of ghost entities to receive from each rank.
   *
   * @param m: material index.
   * @return matrix of communication.
   */
  std::vector<std::vector<int>> const& ghost_matrix(int m = 0) const { return ghost.matrix[m]; }

  /**
   * @brief Retrieve values of owned entities sent to each rank.
   *
   * @return matrix of values.
   */
  std::vector<std::vector<double>> const& owned_values() const { return owned.values; }

  /**
   * @brief Retrieve values of ghost entities received from each rank.
   *
   * @return matrix of values.
   */
  std::vector<std::vector<double>> const& ghost_values() const { return ghost.values; }

private:
  /**
   * @brief Build and store communication matrices to identify
   *        which owned entities data should be sent to each rank,
   *        and which ghost entities data should be received from
   *        each rank.
   */
  void cache_comm_matrices() {

    // step 0: initialization
    num_mats = 1 + state.num_materials();
    owned.matrix.resize(num_mats);
    ghost.matrix.resize(num_mats);
    ghost.count.resize(num_mats);

    for (int m = 0; m < num_mats; ++m) {
      owned.matrix[m].resize(num_ranks);
      ghost.matrix[m].resize(num_ranks);
      ghost.count[m].resize(num_ranks);
    }

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
        ghost.lookup[gid] = i;
      } else {
        owned.lookup[gid] = i;
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
          if (owned.lookup.count(gid)) {
            owned.matrix[0][i].emplace_back(owned.lookup[gid]);
            gids[i].emplace_back(gid);
          }
        }

        MPI_Request request;
        num_sent[i] = owned.matrix[0][i].size();
        MPI_Isend(num_sent + i, 1, MPI_INT, i, 0, comm, &request);
        requests.emplace_back(request);
      }
    }

#ifdef DEBUG
  for (int i = 0; i < num_ranks; ++i) {
    if (i != rank) {
      std::cout << "[" << rank << "] -> [" << i << "]: (";
      int const num_to_send = owned.matrix[0][i].size();
      for (int j = 0; j < num_to_send; ++j) {
        std::cout << mesh.get_global_id(owned.matrix[0][i][j], entity);
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
        MPI_Irecv(ghost.count[0].data() + i, 1, MPI_INT, i, 0, comm, &request);
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
      if (i != rank and ghost.count[0][i]) {
        MPI_Request request;
        ghost.matrix[0][i].resize(ghost.count[0][i]);
        MPI_Irecv(ghost.matrix[0][i].data(), ghost.count[0][i], MPI_INT, i, 1, comm, &request);
        requests.emplace_back(request);
      }
    }

    MPI_Waitall(requests.size(), requests.data(), status);

#ifdef DEBUG
    for (int i = 0; i < num_ranks; ++i) {
      if (i != rank) {
        std::cout << "[" << rank << "] <- [" << i << "]: (";
        for (int j = 0; j < ghost.count[0][i]; ++j) {
          std::cout << ghost.matrix[0][i][j];
          if (j < ghost.count[0][i] - 1) {
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
      total_received += ghost.matrix[0][i].size();
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
  /** number of materials */
  int num_mats = 0;
  /** sent data */
  Data owned;
  /** received data */
  Data ghost;
};

} // namespace Portage
#endif // ifdef WONTON_ENABLE_MPI
#endif // ifndef PORTAGE_MPI_UPDATE_GHOSTS_H
