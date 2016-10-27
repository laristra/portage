/*----------------------------------------------------------------------*
 * Copyright (c) 2016 Los Alamos National Security, LLC
 * All rights reserved.
 *----------------------------------------------------------------------*/

#ifndef SIMPLE_MESH_WRAPPER_H_
#define SIMPLE_MESH_WRAPPER_H_

#include "portage/simple_mesh/simple_mesh.h"

#include <vector>
#include <algorithm>

#include "portage/wrappers/mesh/AuxMeshTopology.h"
#include "portage/support/portage.h"
#include "portage/support/Point.h"

namespace Portage {
class Simple_Mesh_Wrapper : public AuxMeshTopology<Simple_Mesh_Wrapper> {
 public:
  explicit Simple_Mesh_Wrapper(Simple_Mesh const & mesh,
                               bool request_sides = true,
                               bool request_wedges = true,
                               bool request_corners = true) :
  mesh_(mesh),
      AuxMeshTopology<Simple_Mesh_Wrapper>(request_sides, request_wedges,
                                           request_corners) {
    AuxMeshTopology<Simple_Mesh_Wrapper>::build_aux_entities();
  }  // explicit constructor

  Simple_Mesh_Wrapper(Simple_Mesh_Wrapper const & inmesh) = delete;

  Simple_Mesh_Wrapper & operator=(Simple_Mesh_Wrapper const & inmesh) = delete;

  ~Simple_Mesh_Wrapper() {}

  //////////////////////////////////////////////////////////////////////
  // The following methods are needed somewhere within AuxMeshTopology
  int space_dimension() const {
    return mesh_.space_dimension();
  }

  // number of OWNED data entities
  int num_owned_cells() const {
    return mesh_.num_entities(Entity_kind::CELL,
                              Entity_type::PARALLEL_OWNED);
  }

  int num_owned_faces() const {
    return mesh_.num_entities(Entity_kind::FACE,
                              Entity_type::PARALLEL_OWNED);
  }

  int num_owned_nodes() const {
    return mesh_.num_entities(Entity_kind::NODE,
                              Entity_type::PARALLEL_OWNED);
  }

  // number of ghost data entities
  int num_ghost_cells() const {
    return mesh_.num_entities(Entity_kind::CELL,
                              Entity_type::PARALLEL_GHOST);
  }

  int num_ghost_faces() const {
    return mesh_.num_entities(Entity_kind::FACE,
                              Entity_type::PARALLEL_GHOST);
  }

  int num_ghost_nodes() const {
    return mesh_.num_entities(Entity_kind::NODE,
                              Entity_type::PARALLEL_GHOST);
  }

  // cell type: Simple_Mesh only has OWNED
  Entity_type cell_get_type(int const cellid) const {
    return Entity_type::PARALLEL_OWNED;
  }

  // cell element type: Simple_Mesh only deals with HEX
  Element_type cell_get_element_type(int const cellid) const {
    return Element_type::HEX;
  }

  // Connectivity information
  void cell_get_faces_and_dirs(int const cellid, std::vector<int> *cfaces,
                               std::vector<int> *cfdirs) const {
    mesh_.cell_get_faces_and_dirs(cellid, cfaces, cfdirs);
  }

  void cell_get_nodes(int const cellid, std::vector<int> *nodes) const {
    mesh_.cell_get_nodes(cellid, nodes);
  }

  void face_get_nodes(int const faceid, std::vector<int> *nodes) const {
    mesh_.face_get_nodes(faceid, nodes);
  }

  template<long D=3>
      void node_get_coordinates(int const nodeid, Point<D>* pp) const {
    mesh_.node_get_coordinates(nodeid, pp);
  }

  void node_get_cell_adj_nodes(int const nodeid,
                               Entity_type const ptype,
                               std::vector<int> *adjnodes) const {
    adjnodes->clear();

    // Find the cells attached to this node
    std::vector<int> nodecells;
    mesh_.node_get_cells(nodeid, &nodecells);

    // Loop over these cells, and find their nodes; these are the ones we seek
    // but make sure we aren't duplicating them
    for (auto const& c : nodecells) {
      std::vector<int> cellnodes;
      cell_get_nodes(c, &cellnodes);

      for (auto const& n : cellnodes) {
        if (n == nodeid) continue;
        if (std::find(adjnodes->begin(), adjnodes->end(), n) == adjnodes->end())
          adjnodes->emplace_back(n);
      }
    }
  }  // node_get_cell_adj_nodes

  //////////////////////////////////////////////////////////////////////
  // The following methods are needed elsewhere in the source.
  void cell_get_node_adj_cells(int const cellid,
                               Entity_type const ptype,
                               std::vector<int> *adjcells) const {
    adjcells->clear();

    // Find the nodes attached to this cell
    std::vector<int> cellnodes;
    cell_get_nodes(cellid, &cellnodes);

    // Loop over these nodes and find the associated cells; these are the ones
    // we seek, but make sure there are not duplicates.
    for (auto const& n : cellnodes) {
      std::vector<int> nodecells;
      mesh_.node_get_cells(n, &nodecells);

      for (auto const& c : nodecells) {
        if (c == cellid) continue;
        if (std::find(adjcells->begin(), adjcells->end(), c) == adjcells->end())
          adjcells->emplace_back(c);
      }
    }
  }  // cell_get_node_adj_cells

 private:
  Simple_Mesh const & mesh_;
};  // class Simple_Mesh_Wrapper
}  // namespace Portage

#endif  // SIMPLE_MESH_WRAPPER_H_
