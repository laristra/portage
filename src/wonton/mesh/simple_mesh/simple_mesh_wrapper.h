/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/





#ifndef SIMPLE_MESH_WRAPPER_H_
#define SIMPLE_MESH_WRAPPER_H_

#include "portage/simple_mesh/simple_mesh.h"

#include <vector>
#include <algorithm>

#include "portage/wonton/mesh/AuxMeshTopology.h"
#include "portage/support/portage.h"
#include "portage/support/Point.h"

/*!
  @file simple_mesh_wrapper.h
  @brief Definitions for a wrapper to Simple_Mesh
 */

namespace Wonton {
  /*!
    @class Simple_Mesh_Wrapper simple_mesh_wrapper.h
    @brief A thin wrapper that implements mesh methods for Simple_Mesh

    The methods implemented are those required elsewhere in Portage to
    answer specific queries about the mesh.  Simple_Mesh_Wrapper derives
    from the AuxMeshTopology class, which helps to build further mesh
    entities and connectivities (e.g. CORNERS and WEDGES) from the existant
    mesh entities.  This uses a template pattern called the Curiously
    Recurring Template Pattern (CRTP) to allow the Simple_Mesh_Wrapper
    itself to answer queries within the AuxMeshTopology class.  See
    https://en.m.wikipedia.org/wiki/Curiously_recurring_template_pattern
   */
class Simple_Mesh_Wrapper : public AuxMeshTopology<Simple_Mesh_Wrapper> {
 public:
  /*!
    @brief Constructor for the mesh wrapper.
    @param[in] mesh The Simple_Mesh we wish to wrap.
    @param[in] request_sides Should the AuxMeshTopology class build side
    datastructures?
    @param[in] request_wedges Should the AuxMeshTopology class build wedge
    datastructures?
    @param[in] request_corners Should the AuxMeshToplogy class build corner
    datastructures?
   */
  explicit Simple_Mesh_Wrapper(Simple_Mesh const & mesh,
                               bool request_sides = true,
                               bool request_wedges = true,
                               bool request_corners = true) :
  mesh_(mesh),
      AuxMeshTopology<Simple_Mesh_Wrapper>(request_sides, request_wedges,
                                           request_corners) {
    AuxMeshTopology<Simple_Mesh_Wrapper>::build_aux_entities();
  }  // explicit constructor

  /// Copy constructor (disabled).
  Simple_Mesh_Wrapper(Simple_Mesh_Wrapper const & inmesh) = delete;

  /// Assignment operator (disabled).
  Simple_Mesh_Wrapper & operator=(Simple_Mesh_Wrapper const & inmesh) = delete;

  /// Destructor
  ~Simple_Mesh_Wrapper() {}

  //////////////////////////////////////////////////////////////////////
  // The following methods are needed somewhere within AuxMeshTopology

  /// The spatial dimension of the mesh.
  int space_dimension() const {
    return mesh_.space_dimension();
  }

  // The number of OWNED entities

  /// The number of OWNED cells in the mesh.
  int num_owned_cells() const {
    return mesh_.num_entities(Entity_kind::CELL,
                              Entity_type::PARALLEL_OWNED);
  }

  /// The number of OWNED faces in the mesh.
  int num_owned_faces() const {
    return mesh_.num_entities(Entity_kind::FACE,
                              Entity_type::PARALLEL_OWNED);
  }

  /// The number of ONWED nodes in the mesh.
  int num_owned_nodes() const {
    return mesh_.num_entities(Entity_kind::NODE,
                              Entity_type::PARALLEL_OWNED);
  }

  // number of ghost data entities

  /// The number of GHOST cells in the mesh.
  int num_ghost_cells() const {
    return mesh_.num_entities(Entity_kind::CELL,
                              Entity_type::PARALLEL_GHOST);
  }

  /// The number of GHOST faces in the mesh.
  int num_ghost_faces() const {
    return mesh_.num_entities(Entity_kind::FACE,
                              Entity_type::PARALLEL_GHOST);
  }

  /// The number of GHOST nodes in the mesh.
  int num_ghost_nodes() const {
    return mesh_.num_entities(Entity_kind::NODE,
                              Entity_type::PARALLEL_GHOST);
  }

  // cell type: Simple_Mesh only has OWNED
  /*!
    @brief Get the Entity_type (e.g. PARALLEL_OWNED) of a specific cell.
    @param[in] cellid The ID of the cell.
    @returns The Entity_type of this cell. @b NOTE: Simple_Mesh _only_ has
    PARALLEL_OWNED data as it is a serial implementation
   */
  Entity_type cell_get_type(int const cellid) const {
    return Entity_type::PARALLEL_OWNED;
  }

  // node type: Simple_Mesh only has OWNED
  /*!
    @brief Get the Entity_type (e.g. PARALLEL_OWNED) of a specific node.
    @param[in] nodeid The ID of the node.
    @returns The Entity_type of this node. @b NOTE: Simple_Mesh _only_ has
    PARALLEL_OWNED data as it is a serial implementation
   */
  Entity_type node_get_type(int const nodeid) const {
    return Entity_type::PARALLEL_OWNED;
  }

  // cell element type: Simple_Mesh only deals with HEX
  /*!
    @brief Get the Element_type (e.g. HEX) of a specific cell.
    @param[in] cellid The ID of the cell.
    @returns The Element_type of this cell. @b NOTE: Simple_Mesh _only_ has
    HEX cells.
   */
  Element_type cell_get_element_type(int const cellid) const {
    return Element_type::HEX;
  }

  // Connectivity information

  /*!
    @brief Get the list of face IDs and face normal directions for a specific
    cell.
    @param[in] cellid The ID of the cell.
    @param[out] cfaces The vector of face IDs for the faces that make up cell
    @c cellid.
    @param[out] cfdirs The vector of face normal directions for the faces that
    make up cell @c cellid.
  */
  void cell_get_faces_and_dirs(int const cellid, std::vector<int> *cfaces,
                               std::vector<int> *cfdirs) const {
    mesh_.cell_get_faces_and_dirs(cellid, cfaces, cfdirs);
  }

  /*!
    @brief Get the list of node IDs for a specific cell.
    @param[in] cellid The ID of the cell.
    @param[out] nodes The vector of node IDs for the nodes that make up cell
    @c cellid.
   */
  void cell_get_nodes(int const cellid, std::vector<int> *nodes) const {
    mesh_.cell_get_nodes(cellid, nodes);
  }

  /*!
    @brief Get the list of node IDs for a specific face.
    @param[in] faceid The ID of the face.
    @param[out] The list of node IDs that make up face @c faceid.
   */
  void face_get_nodes(int const faceid, std::vector<int> *nodes) const {
    mesh_.face_get_nodes(faceid, nodes);
  }

  /*!
    @brief Get the list of IDs of all cells of a particular parallel type attached to a node.
    @param[in] nodeid The ID of the node.
    @param[in] ptype The Entity_type (e.g. PARALLEL_OWNED); @b NOTE: Simple_Mesh
    only contains PARALLEL_OWNED data, no ghosts.
    @param[out] nodecells The list of IDs for cells of @c ptype attached to @c nodeid
   */
  void node_get_cells(int const nodeid,
                      Entity_type const ptype,
                      std::vector<int> *nodecells) const {
    mesh_.node_get_cells(nodeid, nodecells);
  }

  /// Get the global ID. @b NOTE: Simple_Mesh only has local IDs.
  int get_global_id(int const id, Entity_kind const kind) const {
    return id;
  }


  /*!
    @brief Get the coordinates of a specific node as a Portage::Point.
    @tparam D Dimensionality -- this is a specialization, as Simple_Mesh
    only supports 3D.
    @param[in] nodeid The ID of the node.
    @param[out] The Portage::Point containing the coordiantes of node @c nodeid.
   */
  template<long D=3>
  void node_get_coordinates(int const nodeid, Point<D>* pp) const {
    mesh_.node_get_coordinates(nodeid, pp);
  }

#ifdef HAVE_TANGRAM
  // TEMPORARY - until we pull WONTON out as a separate repository
  int get_global_id(int const id, Tangram::Entity_kind const kind) const {
    return get_global_id(id, static_cast<Portage::Entity_kind>(kind));
  }

  template<long D=3>
  void node_get_coordinates(int const nodeid, Tangram::Point<D>* tcoord) const {
    Point<D> pcoord;
    node_get_coordinates(nodeid, &pcoord);
    *tcoord = pcoord;
  }
#endif
    
 private:
  /// The mesh to wrap.
  Simple_Mesh const & mesh_;
};  // class Simple_Mesh_Wrapper
}  // namespace Wonton

#endif  // SIMPLE_MESH_WRAPPER_H_
