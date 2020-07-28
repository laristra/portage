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

#ifndef DEBUG
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
}
#endif
#endif //PORTAGE_MPI_UPDATE_GHOSTS_H
