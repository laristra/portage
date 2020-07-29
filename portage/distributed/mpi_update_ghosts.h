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
  MPI_GhostManager(Mesh const& in_mesh, State& in_state, MPI_Comm in_comm) {
    mesh  = in_mesh;
    state = in_state;
    comm  = in_comm;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &num_ranks);

    if (num_ranks > 1) {
      send_matrix.resize(num_ranks);
      count.resize(num_ranks);
    }
  }

  /**
   * @brief Delete the manager.
   *
   */
  ~MPI_GhostManager() = default;

  /**
   * @brief Build and store the communication matrix to identify
   *        which owned entities data should be sent to each rank.
   */
  void build_comm_matrix() {

    std::vector<int> entities;
    int num_ghosts[num_ranks];
    int offsets[num_ranks];
    std::vector<int> received;

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
    MPI_Allgather(num_ghosts + rank, 1, MPI_INT, num_ghosts, num_ranks, MPI_INT, comm);

    int const total_ghosts = std::accumulate(num_ghosts, num_ghosts + num_ranks, 0);
    received.resize(total_ghosts);

    int index = 0;
    for (int i = 0; i < num_ranks; ++i) {
      offsets[i] = index;
      index += num_ghosts[i];
    }

    // step 3: gather all ghost entities on all ranks
    MPI_Allgatherv(entities.data(), num_ghosts[rank], MPI_INT, received.data(), num_ghosts, offsets, MPI_INT, comm);

    // step 4: check received ghost cells and build map
    for (int i = 0; i < num_ranks; ++i) {
      if (i != rank) {
        int const& start = offsets[i];
        int const& extent = (i < num_ranks - 1 ? offsets[i+1] : total_ghosts);
        for (int j = start; j < extent; ++j) {
          int const& gid = received[j];
          if (owned.count(gid)) {
            send_matrix[i].emplace_back(owned[gid]);
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
    int const num_to_send = send_matrix[i].size();
    for (int j = 0; j < num_to_send; ++j) {
      std::cout << mesh.get_global_id(send_matrix[i][j], entity);
      if (j < num_to_send - 1) {
        std::cout << ", ";
      }
    }
    std::cout << ")" << std::endl;
  }

  std::cout << " ------------------ " << std::endl;
#endif
  }

  /**
   * @brief Send and store field values on ghost cells.
   *
   * @param fields: list of remapped fields.
   */
  void fill_ghost_values(std::vector<std::string> const& fields) {

    for (auto&& field : fields) {
#ifdef PORTAGE_HAS_TANGRAM
      auto field_type = state.field_type(entity, field);
      bool multimat = (entity == Wonton::CELL and field_type == Field_type::MULTIMATERIAL_FIELD);

      if (multimat) {
        // todo
      } else /* single material */{
#endif
        std::vector<double> send[num_ranks];
        std::vector<double> recv[num_ranks];
        std::vector<int> gids[num_ranks];

        // step 1: retrieve field values and populate data
        double* values = nullptr;
        state.mesh_get_data(entity, &values);
        assert(values != nullptr);

        for (int i = 0; i < num_ranks; ++i) {
          if (i != rank) {
            for (auto&& j : send_matrix[i]) {  // 'j' local entity index
              assert(not is_ghost(j));
              send[i].emplace_back(values[j]);
              gids[i].emplace_back(mesh.get_global_id(j, entity));
            }
            count[i] = send[i].size();
          }
        }

        // step 2: send them
        for (int i = 0; i < num_ranks; ++i) {
          if (i != rank) {
            MPI_Send(count.data() + i, 1, MPI_INT, i, 0, comm);
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
            MPI_Recv(count.data() + i, 1, MPI_INT, i, 0, comm, &status);
            if (count[i] > 0) {
              recv[i].resize(count[i]);
              gids[i].resize(count[i], -1);  // reset and reuse
              MPI_Recv(recv[i].data(), count[i], MPI_DOUBLE, i, 1, comm, &status);
              MPI_Recv(gids[i].data(), count[i], MPI_INT, i, 2, comm, &status);
            }
          }
        }

        // step 4: store them
        for (int i = 0; i < num_ranks; ++i) {
          if (i != rank and count[i] > 0) {
            for (int j = 0; j < count[i]; ++j) {
              int const& gid = gids[i][j];
              int const& lid = ghost[gid];
              values[lid] = recv[i][j];
            }
          }
        }
#ifdef PORTAGE_HAS_TANGRAM
      } // if not multimat
#endif
    }
  }

private:
  /** mesh instance */
  Mesh const& mesh;
  /** mesh state */
  State& state;
  /** current MPI rank */
  int rank = 0;
  /** number of ranks */
  int num_ranks = 1;
  /** MPI communicator */
  MPI_Comm comm = MPI_COMM_NULL;
  /** list of owned entities to send to each rank */
  std::vector<std::vector<int>> send_matrix {};
  /** count of received entities from each rank */
  std::vector<int> count {};
  /** owned entities indexed by their GID */
  std::map<int, int> owned {};
  /** ghost entities indexed by their GID */
  std::map<int, int> ghost {};

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
};


} // namespace Portage
#endif
#endif //PORTAGE_MPI_UPDATE_GHOSTS_H
