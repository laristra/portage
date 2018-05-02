/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/



// Copyright 2016, Los Alamos National Laboratory, NM, USA

#ifndef AUX_MESH_TOPOLOGY_H_
#define AUX_MESH_TOPOLOGY_H_

#include <vector>
#include <array>
#include <algorithm>
#include <utility>
#include <limits>

#include "portage/support/portage.h"
#include "portage/support/Point.h"


// Temporary - until we pull WONTON out as a separate repository
#ifdef HAVE_TANGRAM
#include "tangram/support/tangram.h"
#endif

namespace Wonton {

using namespace Portage;
// Some helper functions
//! Compute volume of 1D side
inline
double calc_side_volume(std::array<Point<1>, 2> const& sxyz) {
  return (sxyz[1][0]-sxyz[0][0]);
}

//! Compute volume of 2D side
inline
double calc_side_volume(std::array<Point<2>, 3> const& sxyz) {
  Vector<2> vec1 = sxyz[1]-sxyz[0];
  Vector<2> vec2 = sxyz[2]-sxyz[0];
  return 0.5*cross(vec1, vec2);
}

//! Compute volume of 3D side
inline
double calc_side_volume(std::array<Point<3>, 4> const& sxyz) {
  Vector<3> vec1 = sxyz[1]-sxyz[0];
  Vector<3> vec2 = sxyz[2]-sxyz[0];
  Vector<3> vec3 = sxyz[3]-sxyz[0];
  Vector<3> cpvec = cross(vec1, vec2);
  return dot(cpvec, vec3)/6.0;
}

template<typename BasicMesh> class AuxMeshTopology;
template<typename BasicMesh>
void build_sides_1D(AuxMeshTopology<BasicMesh>& mesh);
template<typename BasicMesh>
void build_sides_2D(AuxMeshTopology<BasicMesh>& mesh);
template<typename BasicMesh>
void build_sides_3D(AuxMeshTopology<BasicMesh>& mesh);


//! \class AuxMeshTopology.h
//! \brief A class to build enhanced mesh topology (mainly subcell
//! entities - sides, wedges and corners)
//!
//! @tparam BasicMesh Mesh wrapper class that provides methods for answering
//! questions about cells, nodes (and depending on the mesh dimension,
//! mesh faces and edges)
//!
//! A class to build enhanced mesh topological entities. For now this
//! class builds the following subcell entities - sides, corners and
//! wedges. The non-standard entities are required for more accurate
//! remapping when one or both meshes have elements with non-planar
//! faces
//!
//! 1D:
//! A side is a line segment from a node to the cell. Wedges and
//! corners are the same as sides.
//
//! 2D:
//! A side is a triangle formed by the two nodes of an edge/face and
//! the cell center. A wedge is half of a side formed by one node of
//! the edge, the edge center and the cell center. A corner is a
//! quadrilateral formed by the two wedges in a cell at a node
//!
//! 3D:
//! A side is a tet formed by the two nodes of an edge, a face center
//! and a cell center. A wedge is half a side, formed by a node of
//! the edge, the edge center, the face center and the cell center. A
//! corner is formed by all the wedges of a cell at a node.
//!
//!
//! The basic mesh class must support cells, faces and nodes and
//! adjacency queries between these entities (In 2D, faces are the
//! same as edges and in 1D, faces are the same as nodes). In
//! particular, the basic mesh class is expected to support the
//! following methods to successfully instantiate this class:
//!
//!~~~
//! int space_dimension() const;  // dimensionality of mesh points (1, 2, 3)
//!
//! int num_owned_cells() const;
//! int num_ghost_cells() const;
//! int num_owned_faces() const;
//! int num_ghost_faces() const;
//! int num_owned_nodes() const;
//! int num_ghost_nodes() const;
//! Portage::Entity_type cell_get_type(int const cellid) const;
//! Portage::Entity_type node_get_type(int const nodeid) const;
//!~~~
//!
//! NOTE: Entity_type is Portage::OWNED or Portage::GHOST
//!
//!~~~
//! Portage::Element_type cell_get_element_type(int const cellid) const;
//!~~~
//!
//! Can be Portage::UNKNOWN_TOPOLOGY, Portage::TRI, Portage::QUAD,
//! Portage::POLYGON, Portage::TET, Portage::PRISM, Portage::PYRAMID,
//! Portage::HEX, Portage::POLYHEDRON
//!
//!
//!~~~
//! void cell_get_faces_and_dirs(int const cellid, std::vector<int> *cfaces,
//!                              std::vector<int> *cfdirs) const;
//!~~~
//!
//! NOTE: The 'cfdirs' conveys the directions in which the faces are used by
//! the cell. If the natural normal of the face points out of the cell, its
//! direction should be returned as 1, if not, it should be returned as -1
//!
//!
//!~~~
//! void cell_get_nodes(int const cellid, std::vector<int> *cnodes) const;
//!
//! void face_get_nodes(int const faceid, std::vector<int> *fnodes) const;
//!
//! void face_get_cells(int const faceid, Portage::Entity_type etype,
//!                     std::vector<int> *fcells) const;
//!
//! void node_get_cells(int const nodeid, Portage::Entity_type etype,
//!                     std::vector<int> *ncells) const;
//!
//! int get_global_id(int const id, Entity_kind const kind) const;
//!
//! template<long D>
//! void node_get_coordinates(int const nodeid, Portage::Point<D> *pp) const;
//!
//!~~~
//! ******************************** NOTE ***********************************
//!
//! THIS IS AN INCOMPLETE CLASS DESIGNED TO BE USED IN A 'CRTP' (CURIOUSLY
//! RECURRING TEMPLATE PATTERN) DESIGN ALONG WITH THE BASICMESH CLASS TO
//! PROVIDE A COMPLETE MESH CLASS
//! (See https://en.m.wikipedia.org/wiki/Curiously_recurring_template_pattern)
//!
//! So, If one is writing a mesh wrapper class called @c MY_MESH_WRAPPER, it is
//! declared like so
//!
//!~~~
//! class MY_MESH_WRAPPER : public AuxMeshTopology<MY_MESH_WRAPPER>
//! {......}
//!~~~
//! and it will automatically have the methods of the AuxMeshTopology class.
//!
//! If MY_MESH_WRAPPER has equivalent classes for the ones in
//! AuxMeshTopology (possibly because they are more efficient), then
//! the ones from AuxMeshTopology are overridden
//!
//! NOTE THAT THIS CLASS IS NOT DESIGNED TO EVER BE INSTANTIATED DIRECTLY
//!
//!***************************************************************************



template<typename BasicMesh>
class AuxMeshTopology {
 public:

  //! @brief Constructor indicating which entities are wanted.
  //!
  //! It is possible to request none of the auxiliary entities so that
  //! one can instantiate a version of a derived class with the
  //! auxiliary entities or without (to save memory)

  AuxMeshTopology(bool request_sides = true,
                  bool request_wedges = true,
                  bool request_corners = true) :
      sides_requested_(request_sides || request_wedges || request_corners),
      wedges_requested_(request_wedges || request_corners),
      corners_requested_(request_corners) {

    // we really cannot do anything without sides, so we will turn
    // them on all the time. Don't want to change the interface at
    // this time though.
    
    sides_requested_ = true;

    // Part of the CRTP magic. Knowing that we will be deriving the
    // BasicMesh class from this class, we can cast a pointer to this
    // class type to a pointer of BasicMesh class type and call its
    // methods

    basicmesh_ptr_ = static_cast<BasicMesh*>(this);

    // Sadly we cannot trigger building of entities here because, this
    // routine reaches into the derived class which has not yet had a chance
    // to initialize some class data which will support the basic mesh queries
    // needed. So we have to call this routine from the derived class

    // build_aux_entities()
  }


  // //! A method expected to be found in the BasicMesh class but defined
  // //! here as pure virtual to prevent this class from ever being
  // //! instantiated directly

  // virtual int space_dimension() const = 0;


  //! Number of owned sides in the mesh

  int num_owned_sides() const {
    return num_sides_owned_;
  }


  //! Number of owned wedges in the mesh

  int num_owned_wedges() const {
    return num_wedges_owned_;
  }


  //! Number of owned corners in the mesh

  int num_owned_corners() const {
    return num_corners_owned_;
  }


  //! Number of ghost sides in the mesh

  int num_ghost_sides() const {
    return num_sides_ghost_;
  }


  //! Number of ghost wedges in the mesh

  int num_ghost_wedges() const {
    return num_wedges_ghost_;
  }


  //! Number of ghost corners in the mesh

  int num_ghost_corners() const {
    return num_corners_ghost_;
  }


#ifdef HAVE_TANGRAM
  // TEMPORARY - until we pull WONTON out as a separate repository
  int num_entities(Tangram::Entity_kind const entity,
                   Tangram::Entity_type const etype = Tangram::Entity_type::ALL)
      const {
    return num_entities(static_cast<Portage::Entity_kind>(entity),
                        static_cast<Portage::Entity_type>(etype));
  }
#endif

  //! Number of items of given entity
  int num_entities(Entity_kind const entity,
                   Entity_type const etype = Entity_type::ALL) const {
    switch (entity) {
      case CELL:
        switch (etype) {
          case PARALLEL_OWNED: return basicmesh_ptr_->num_owned_cells();
          case PARALLEL_GHOST: return basicmesh_ptr_->num_ghost_cells();
          case ALL: return (basicmesh_ptr_->num_owned_cells() +
                            basicmesh_ptr_->num_ghost_cells());
          default: return 0;
        }
      case FACE:
        switch (etype) {
          case PARALLEL_OWNED: return basicmesh_ptr_->num_owned_faces();
          case PARALLEL_GHOST: return basicmesh_ptr_->num_ghost_faces();
          case ALL: return (basicmesh_ptr_->num_owned_faces() +
                            basicmesh_ptr_->num_ghost_faces());
          default: return 0;
        }
      case NODE:
        switch (etype) {
          case PARALLEL_OWNED: return basicmesh_ptr_->num_owned_nodes();
          case PARALLEL_GHOST: return basicmesh_ptr_->num_ghost_nodes();
          case ALL: return (basicmesh_ptr_->num_owned_nodes() +
                            basicmesh_ptr_->num_ghost_nodes());
          default: return 0;
        }
      case SIDE:
        switch (etype) {
          case PARALLEL_OWNED: return num_owned_sides();
          case PARALLEL_GHOST: return num_ghost_sides();
          case ALL: return (num_owned_sides() + num_ghost_sides());
          default: return 0;
        }
      case WEDGE:
        switch (etype) {
          case PARALLEL_OWNED: return num_owned_wedges();
          case PARALLEL_GHOST: return num_ghost_wedges();
          case ALL: return (num_owned_wedges() + num_ghost_wedges());
          default: return 0;
        }
      case CORNER:
        switch (etype) {
          case PARALLEL_OWNED: return num_owned_corners();
          case PARALLEL_GHOST: return num_ghost_corners();
          case ALL: return (num_owned_corners() + num_ghost_corners());
          default: return 0;
        }
      default:
        return 0;
    }
  }

  //! Iterators on mesh entity - begin
  counting_iterator begin(Entity_kind const entity,
                          Entity_type const etype = Entity_type::ALL) const {
    int start_index = 0;
    return make_counting_iterator(start_index);
  }

  //! Iterator on mesh entity - end
  counting_iterator end(Entity_kind const entity,
                        Entity_type const etype = Entity_type::ALL) const {
    int start_index = 0;
    return (make_counting_iterator(start_index) + num_entities(entity, etype));
  }

#ifdef HAVE_TANGRAM
  // TEMPORARY - until we pull WONTON out as a separate repository
  //! Iterators on mesh entity - begin
  counting_iterator begin(Tangram::Entity_kind const entity,
                          Tangram::Entity_type const etype = Tangram::Entity_type::ALL) const {
    return begin(static_cast<Portage::Entity_kind>(entity),
                 static_cast<Portage::Entity_type>(etype));
  }

  //! Iterator on mesh entity - end
  counting_iterator end(Tangram::Entity_kind const entity,
                        Tangram::Entity_type const etype = Tangram::Entity_type::ALL) const {
    return end(static_cast<Portage::Entity_kind>(entity),
               static_cast<Portage::Entity_type>(etype));
  }
#endif
  /*!
    @brief Get the list of cell IDs for all cells attached to a specific
    cell through its nodes.
    @param[in] cellid The ID of the cell.
    @param[in] ptype The Entity_type (e.g. PARALLEL_OWNED)
    @param[out] adjcells The list of cell IDs for all cells attached to
    cell @c cellid through its nodes, excluding @c cellid.
   */
  void cell_get_node_adj_cells(int const cellid,
                               Entity_type const ptype,
                               std::vector<int> *adjcells) const {

    // If it does not interfere with the GPU implementation, we could
    // insert into a std::set (nlog(n) instead of n^2) and then copy
    // into the output vector

    adjcells->clear();

    // Find the nodes attached to this cell
    std::vector<int> cellnodes;
    basicmesh_ptr_->cell_get_nodes(cellid, &cellnodes);

    // Loop over these nodes and find the associated cells; these are the ones
    // we seek, but make sure there are not duplicates.
    for (auto const& n : cellnodes) {
      std::vector<int> nodecells;
      basicmesh_ptr_->node_get_cells(n, ptype, &nodecells);

      for (auto const& c : nodecells) {
        if (c == cellid) continue;
        if (std::find(adjcells->begin(), adjcells->end(), c) == adjcells->end())
          adjcells->push_back(c);
      }
    }
  }  // cell_get_node_adj_cells

  //! Get cells of given Entity_type connected to face (in no particular order)
  void face_get_cells(int const faceid, Entity_type const etype,
                      std::vector<int> *cells) const {
    cells->clear();
    int c0 = face_cell_ids_[faceid][0];
    if (c0 != -1) {
      Entity_type ctype0 = basicmesh_ptr_->cell_get_type(c0);
      if (etype == Entity_type::ALL || ctype0 == etype)
        cells->push_back(c0);
    }
    int c1 = face_cell_ids_[faceid][1];
    if (c1 != -1) {
      Entity_type ctype1 = basicmesh_ptr_->cell_get_type(c1);
      if (etype == Entity_type::ALL || ctype1 == etype)
        cells->push_back(c1);
    }
  }

  /*!
    @brief Get the list of node IDs for all nodes attached to all cells
    attached to a specific node.
    @param[in] nodeid The ID of the node.
    @param[in] ptype The Entity_type (e.g. PARALLEL_OWNED)
    @param[out] adjnodes The list of node IDs for all cells attached to
    @c nodeid, excluding @c nodeid.
   */
  void node_get_cell_adj_nodes(int const nodeid,
                               Entity_type const ptype,
                               std::vector<int> *adjnodes) const {
    adjnodes->clear();

    // If it does not interfere with the GPU implementation, we could
    // insert into a std::set (nlog(n) instead of n^2) and then copy
    // into the output vector

    // Find the cells attached to this node
    std::vector<int> nodecells;
    basicmesh_ptr_->node_get_cells(nodeid, Entity_type::ALL, &nodecells);

    // Loop over these cells, and find their nodes; these are the ones we seek
    // but make sure we aren't duplicating them
    for (auto const& c : nodecells) {
      std::vector<int> cellnodes;
      basicmesh_ptr_->cell_get_nodes(c, &cellnodes);

      for (auto const& n : cellnodes) {
        if (n == nodeid) continue;
        Entity_type ntype = basicmesh_ptr_->node_get_type(n);
        if (ptype == Entity_type::ALL || ntype == ptype) {
          if (std::find(adjnodes->begin(), adjnodes->end(), n) == adjnodes->end())
            adjnodes->push_back(n);
        }
      }
    }
  }  // node_get_cell_adj_nodes


  //! if entity is on exterior boundary
  bool on_exterior_boundary(Entity_kind const entity, int const entity_id) const {
    switch (entity) {
      case NODE:
        return node_on_exterior_boundary_[entity_id];
      case FACE:
        return face_on_exterior_boundary_[entity_id];
      case CELL:
        return cell_on_exterior_boundary_[entity_id];
      case SIDE:
        return cell_on_exterior_boundary_[side_cell_id_[entity_id]];
      case CORNER:
        return cell_on_exterior_boundary_[corner_cell_id_[entity_id]];
      default:
        return false;
    }
  }
        

  //! Coordinates of nodes of cell

  template <long D>
  void cell_get_coordinates(int const cellid,
                            std::vector<Point<D>> *pplist) const;

  //! Centroid of a cell (actually geometric center of cell nodes)

  template <long D>
  void cell_centroid(int const cellid, Point<D> *ccen) const {
#ifdef DEBUG
    assert(cellid < num_entities(CELL, ALL));
#endif
    for (int d = 0; d < D; ++d)
      (*ccen)[d] = cell_centroids_[cellid][d];    
  }


  //! Volume of a cell
  
  double cell_volume(int const cellid) const {
#ifdef DEBUG
    assert(cellid < num_entities(CELL, ALL));
#endif
    assert(sides_requested_);
    return cell_volumes_[cellid];
  }

  
  //! Centroid of a face (actually geometric center of face nodes)

  template <long D>
  void face_centroid(int const faceid, Point<D> *fcen) const {
#ifdef DEBUG
    assert(faceid < num_entities(FACE, ALL));
#endif
    for (int d = 0; d < D; ++d)
      (*fcen)[d] = face_centroids_[faceid][d];
  }


  //! Node of a side
  //!
  //! Each side is tied to two mesh nodes in 2D and 3D and inode = 0
  //! or 1 indicates which one to return. In 1D, the same node is
  //! returned whether inode = 0 or 1. In 2D and 3D, the node ordering
  //! is such that node 0, node 1 and the cell centroid form a positive
  //! area triangle and in 3D, node 0, node 1, the face centroid and cell
  //! centroid form a positive volume tet.

  int side_get_node(int const sideid, int const inode) const {
#ifdef DEBUG
    assert(sideid < num_entities(SIDE, ALL));
#endif
    assert(sides_requested_);
    assert(inode == 0 || inode == 1);
    return side_node_ids_[sideid][inode];
  }


  //! Cell of side

  int side_get_cell(int const sideid) const {
#ifdef DEBUG
    assert(sideid < num_entities(SIDE, ALL));
#endif
    assert(sides_requested_);
    return side_cell_id_[sideid];
  }


  //! Face of side

  int side_get_face(int const sideid) const {
#ifdef DEBUG
    assert(sideid < num_entities(SIDE, ALL));
#endif
    assert(sides_requested_);
    return side_face_id_[sideid];
  }


  //! Wedge of side
  //!
  //! Each side points to two wedges - iwedge (0, 1)
  //! indicates which one to return; the wedge returned will be
  //! consistent with the node returned by side_get_node.
  //! So, side_get_node(s,i) = wedge_get_node(side_get_wedge(s,i))

  int side_get_wedge(int const sideid, int iwedge) const {
#ifdef DEBUG
    assert(sideid < num_entities(SIDE, ALL));
#endif
    assert(wedges_requested_);
    return 2*sideid + iwedge;
  }


  //! Opposite side in neighboring cell of a side.
  //!
  //! The two sides share facet 0 of wedge comprised of nodes 0,1 of
  //! the common edge and center point of the common face in 3D, and
  //! nodes 0,1 of the common edge in 2D. At boundaries, this routine
  //! returns -1

  int side_get_opposite_side(int const sideid) const {
#ifdef DEBUG
    assert(sideid < num_entities(SIDE, ALL));
#endif
    assert(sides_requested_);
    return side_opp_side_id_[sideid];
  }


  //! Get all the sides of a cell

  void cell_get_sides(int const cellid, std::vector<int> *csides) const {
#ifdef DEBUG
    assert(cellid < num_entities(CELL, ALL));
#endif
    assert(sides_requested_);
    csides->resize(cell_side_ids_[cellid].size());
    std::copy(cell_side_ids_[cellid].begin(), cell_side_ids_[cellid].end(),
              csides->begin());
  }

  //! Get coordinates of side in 3D
  //!
  //! If posvol_order = true, then the coordinates will be returned in
  //! an order that will result in a positive volume (in 3D this
  //! assumes that the computation for volume is done as (V01 x
  //! V02).V03 where V0k is a vector from coordinate 0 to coordinate k
  //! of the tet). If posvol_order is false, the coordinates will be
  //! returned in a fixed order - in 3D, this is node point 0, node
  //! point 1, face center, cell center. By default the coordinates are
  //! returned in the natural order (posvol_order = false)

  void side_get_coordinates(int const sideid,
                            std::array<Point<3>, 4> *scoords,
                            bool posvol_order = false) const;

  //! Get coordinates of side in 2D
  //!
  //! If posvol_order = true, then the coordinates will be returned in
  //! an order that will result in a positive area (in 2D this assumes
  //! that the computation for area is done as (V01 x V02).V03 where
  //! V0k is a vector from coordinate 0 to coordinate k of the
  //! tri). If posvol_order is false, the coordinates will be returned
  //! in a fixed order - in 2D, this is node point 0, node point 1,
  //! cell center. By default the coordinates are returned in the
  //! natural order (posvol_order = false)

  void side_get_coordinates(int const sideid,
                            std::array<Point<2>, 3> *scoords,
                            bool posvol_order = false) const;

  //! Get coordinates of side in 1D
  //!
  //! If posvol_order = true, then the coordinates will be returned in
  //! an order that will result in a positive length. If posvol_order
  //! is false, the coordinates will be returned in a fixed order - in
  //! 1D, this is node point and cell center. By default the
  //! coordinates are returned in the natural order (posvol_order =
  //! false)

  void side_get_coordinates(int const sideid,
                            std::array<Point<1>, 2> *scoords,
                            bool posvol_order = false) const;


  //! Volume of a side

  double side_volume(int const sideid) const {
#ifdef DEBUG
    assert(sideid < num_entities(SIDE, ALL));
#endif
    assert(sides_requested_);
    return side_volumes_[sideid];
  }

  //! Side of wedge

  int wedge_get_side(int const wedgeid) const {
#ifdef DEBUG
    assert(wedgeid < num_entities(WEDGE, ALL));
#endif
    assert(wedges_requested_);
    return static_cast<int>(wedgeid/2);
  }


  //! Cell of wedge

  int wedge_get_cell(int const wedgeid) const {
#ifdef DEBUG
    assert(wedgeid < num_entities(WEDGE, ALL));
#endif
    assert(wedges_requested_);
    int sideid = static_cast<int>(wedgeid/2);
    return side_cell_id_[sideid];
  }


  //! Face of wedge

  int wedge_get_face(int const wedgeid) const {
#ifdef DEBUG
    assert(wedgeid < num_entities(WEDGE, ALL));
#endif
    assert(wedges_requested_);
    int sideid = static_cast<int>(wedgeid/2);
    return side_face_id_[sideid];
  }


  //! Corner of a wedge

  int wedge_get_corner(int const wedgeid) const {
#ifdef DEBUG
    assert(wedgeid < num_entities(WEDGE, ALL));
#endif
    assert(corners_requested_);
    return wedge_corner_id_[wedgeid];
  }


  //! node of a wedge

  int wedge_get_node(int const wedgeid) const {
#ifdef DEBUG
    assert(wedgeid < num_entities(WEDGE, ALL));
#endif
    assert(wedges_requested_);
    int sideid = static_cast<int>(wedgeid/2);
    return wedgeid%2 ? side_node_ids_[sideid][1] : side_node_ids_[sideid][0];
  }


  //! Opposite wedge in neighboring cell of a wedge.
  //!
  //! The two wedges share facet 0 of wedge comprised of the node,
  //! center point of the common edge and center point of the common
  //! face in 3D, and node and edge center in 2D. At boundaries, this
  //! routine returns -1

  int wedge_get_opposite_wedge(const int wedgeid) const {
#ifdef DEBUG
    assert(wedgeid < num_entities(WEDGE, ALL));
#endif
    assert(wedges_requested_);
    int sideid = static_cast<int>(wedgeid/2);
    int oppsideid = side_opp_side_id_[sideid];
    if (oppsideid >= 0)
      return (wedgeid%2 ? 2*oppsideid : 2*oppsideid + 1);
    else
      return -1;
  }


  //! adjacent wedge along edge in the same cell.
  //!
  //! The two wedges share facet 1 of wedge comprised of edge center,
  //! face center and zone center in 3D, and node and zone center in
  //! 2D

  int wedge_get_adjacent_wedge(const int wedgeid) const {
#ifdef DEBUG
    assert(wedgeid < num_entities(WEDGE, ALL));
#endif
    // Wedges come in pairs; their IDs are (2*sideid) and (2*sideid+1)
    // If the wedge ID is an odd number, then the adjacent wedge ID is
    // wedge ID minus one; If it is an even number, the adjacent wedge
    // ID is wedge ID plus one

    return (wedgeid%2 ? wedgeid - 1 : wedgeid + 1);
  }


  //! Volume of a wedge - half its side volume

  double wedge_volume(int const wedgeid) const {
#ifdef DEBUG
    assert(wedgeid < num_entities(WEDGE, ALL));
#endif
    int sideid = static_cast<int>(wedgeid/2);
    return 0.5*side_volumes_[sideid];
  }


  //! Get coordinates of wedge in 3D
  //!
  //! If posvol_order = true, then the coordinates will be returned in
  //! an order that will result in a positive volume (in 3D this
  //! assumes that the computation for volume is done as (V01 x
  //! V02).V03 where V0k is a vector from coordinate 0 to coordinate k
  //! of the tet). If posvol_order is false, the coordinates will be
  //! returned in a fixed order - in 3D, this is node point, edge
  //! center, face center, cell center. By default the coordinates are
  //! returned in the natural order (posvol_order = false)

  void wedge_get_coordinates(int const wedgeid,
                             std::array<Point<3>, 4> *wcoords,
                             bool posvol_order = false) const;


  //! Get coordinates of wedge in 2D
  //!
  //! If posvol_order = true, then the coordinates will be returned in
  //! an order that will result in a positive area (in 2D this assumes
  //! that the computation for area is done as (V01 x V02).V03 where
  //! V0k is a vector from coordinate 0 to coordinate k of the
  //! tri). If posvol_order is false, the coordinates will be returned
  //! in a fixed order - in 2D, this is node point, face center, cell
  //! center. By default the coordinates are returned in the natural
  //! order (posvol_order = false)

  void wedge_get_coordinates(int const wedgeid,
                             std::array<Point<2>, 3> *wcoords,
                             bool posvol_order = false) const;


  //! Get coordinates of wedge in 1D
  //!
  //! If posvol_order = true, then the coordinates will be returned in
  //! an order that will result in a positive length. If posvol_order
  //! is false, the coordinates will be returned in a fixed order - in
  //! 1D, this is node point and cell center. By default the
  //! coordinates are returned in the natural order (posvol_order =
  //! false)

  void wedge_get_coordinates(int const wedgeid,
                             std::array<Point<1>, 2> *wcoords,
                             bool posvol_order = false) const;


  //! Get all the wedges in a cell

  void cell_get_wedges(int const cellid, std::vector<int> *wedgeids) const {
#ifdef DEBUG
    assert(cellid < num_entities(CELL, ALL));
#endif
    int nsides = cell_side_ids_[cellid].size();
    wedgeids->resize(2*nsides);
    int iw = 0;
    for (auto const& s : cell_side_ids_[cellid]) {
      (*wedgeids)[iw++] = 2*s;
      (*wedgeids)[iw++] = 2*s + 1;
    }
  }


  //! Get wedges at a node

  void node_get_wedges(int const nodeid, Entity_type const type,
                       std::vector<int> *wedgeids) const {
#ifdef DEBUG
    assert(nodeid < num_entities(NODE, ALL));
#endif
    assert(wedges_requested_ && corners_requested_);
    wedgeids->clear();
    for (auto const& cn : node_corner_ids_[nodeid]) {
      std::vector<int> cnwedges = corner_wedge_ids_[cn];
      for (auto const& w : cnwedges) {
        int s = static_cast<int>(w/2);
        int c = side_cell_id_[s];

        if (type == ALL ||
            static_cast<int>(basicmesh_ptr_->cell_get_type(c)) == type)
          wedgeids->push_back(w);
      }
    }
  }


  //! Get node of corner

  int corner_get_node(const int cornerid) const {
#ifdef DEBUG
    assert(cornerid < num_entities(CORNER, ALL));
#endif
    assert(corners_requested_);
    return corner_node_id_[cornerid];
  }


  //! Get cell of corner

  int corner_get_cell(int const cornerid) const {
#ifdef DEBUG
    assert(cornerid < num_entities(CORNER, ALL));
#endif
    assert(corners_requested_);
    return corner_cell_id_[cornerid];
  }


  //! Get wedges of a corner

  void corner_get_wedges(int const cornerid, std::vector<int> *wedgeids) const {
#ifdef DEBUG
    assert(cornerid < num_entities(CORNER, ALL));
#endif
    assert(corners_requested_);
    wedgeids->resize(corner_wedge_ids_[cornerid].size());
    std::copy(corner_wedge_ids_[cornerid].begin(),
              corner_wedge_ids_[cornerid].end(), wedgeids->begin());
  }


  //! Get corners connected to a node

  void node_get_corners(int const nodeid, Entity_type const type,
                        std::vector<int> *cornerids) const {
#ifdef DEBUG
    assert(nodeid < num_entities(NODE, ALL));
#endif
    assert(corners_requested_);
    if (type == ALL) {
      cornerids->resize(node_corner_ids_[nodeid].size());
      std::copy(node_corner_ids_[nodeid].begin(),
                node_corner_ids_[nodeid].end(),
                cornerids->begin());
    } else {
      cornerids->clear();
      for (auto const& cn : node_corner_ids_[nodeid]) {
        int c = corner_cell_id_[cn];
        if (static_cast<int>(basicmesh_ptr_->cell_get_type(c)) == type)
          cornerids->push_back(cn);
      }
    }
  }


  //! Get corners in a cell

  void cell_get_corners(int const cellid, std::vector<int> *cornerids) const {
#ifdef DEBUG
    assert(cellid < num_entities(CELL, ALL));
#endif
    assert(corners_requested_);
    cornerids->resize(cell_corner_ids_[cellid].size());
    std::copy(cell_corner_ids_[cellid].begin(), cell_corner_ids_[cellid].end(),
              cornerids->begin());
  }

  //! Get a cell's corner at a particular node of the cell

  int cell_get_corner_at_node(int const cellid, int const nodeid) const {
#ifdef DEBUG
    assert(cellid < num_entities(CELL, ALL));
#endif
    assert(corners_requested_);
    for (auto const& cn : cell_corner_ids_[cellid])
      if (corner_get_node(cn) == nodeid)
        return cn;
  }

  //! Volume of a corner

  double corner_volume(int const cornerid) const {
#ifdef DEBUG
    assert(cornerid < num_entities(CORNER, ALL));
#endif
    assert(corners_requested_);
    double cnvol = 0.0;
    for (auto const& w : corner_wedge_ids_[cornerid])
      cnvol += wedge_volume(w);
    return cnvol;
  }

  //! Get a triangular facetization of polyhedral cell boundary
  void cell_get_facetization(int const cellid,
                             std::vector<std::vector<int>> *facetpoints,
                             std::vector<Point<3>> *points) const;


  //! Get a triangular facetization of boundary of dual cell (or in
  //! other words, control volume) corresponding to a node
  void dual_cell_get_facetization(int const nodeid,
                                  std::vector<std::vector<int>> *facetpoints,
                                  std::vector<Point<3>> *points) const;


  //! Get the simplest possible decomposition of a 3D cell into tets.

  void decompose_cell_into_tets(int cellid,
              std::vector<std::array<Portage::Point<3>, 4>> *tcoords,
                                const bool planar_hex) const {
#ifdef DEBUG
    assert(cellid < num_entities(CELL, ALL));
#endif
    if (planar_hex
            && basicmesh_ptr_->cell_get_element_type(cellid) == HEX) {
      // Decompose a hex into a 5-tet decomposition.

      // IMPORTANT: This only works well if the hex is planar. Otherwise we
      // need to implement some more sophisticated solution.
      std::vector<Point<3>> coords;
      basicmesh_ptr_->cell_get_coordinates(cellid, &coords);
      /*

      The hex is returned using the following numbering::

        7---6
       /|  /|
      4---5 |
      | 3-|-2
      |/  |/
      0---1

      Then `indexes` contains the 5-tet decomposition of this hex.

      The tet's vertices are ordered in the following way:

         3
       / | \
      /  |  \
      2--|---1
       \ | /
         0
      */
      const int indexes[5][4] = {
          {0, 1, 3, 4},
          {1, 4, 5, 6},
          {1, 3, 4, 6},
          {1, 6, 2, 3},
          {4, 7, 6, 3}
        };
      for (unsigned int t = 0; t < 5; t++) {
        std::array<Portage::Point<3>, 4> tmp;
        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 3; j++) {
            tmp[i][j] = coords[indexes[t][i]][j];
          }
        }
        tcoords->push_back(tmp);
      }
    } else {
      sides_get_coordinates(cellid, tcoords);
    }
  }


  //! 2D version of coords of nodes of a dual cell
  // Input is the node ID 'nodeid', and it returns the vertex coordinates of
  // the dual cell around this node in `pplist`. The vertices are ordered CCW.
  // For boundary node 'nodeid', the first vertex is the node itself, this
  // uniquely determines the 'pplist' vector. For node 'nodeid' not on a
  // boundary, the vector 'pplist' starts with a random vertex, but it is still
  // ordered CCW. Use the 'dual_cell_coordinates_canonical_rotation' to rotate
  // the 'pplist' into a canonical (unique) form.

  void dual_cell_get_coordinates(int const nodeid,
                    std::vector<Point<2>> *pplist) const {
#ifdef DEBUG
    assert(nodeid < num_entities(NODE, ALL));
#endif
    assert(basicmesh_ptr_->space_dimension() == 2);
    assert(wedges_requested_);

    int cornerid, wedgeid, wedgeid0;
    std::vector<int> cornerids, wedgeids;
    std::array<Point<2>, 3> wcoords;  // (node, edge midpoint, centroid)

    // Start with an arbitrary corner
    node_get_corners(nodeid, ALL, &cornerids);
    cornerid = cornerids[0];

    // Process this corner
    corner_get_wedges(cornerid, &wedgeids);
    order_wedges_ccw(&wedgeids);
    wedge_get_coordinates(wedgeids[0], &wcoords);
    pplist->emplace_back(wcoords[2]);  // centroid
    wedge_get_coordinates(wedgeids[1], &wcoords);
    pplist->emplace_back(wcoords[1]);  // edge midpoint

    wedgeid0 = wedgeids[0];

    // Process the rest of the corners in the CCW manner
    wedgeid = wedge_get_opposite_wedge(wedgeids[1]);
    // wedgeid == -1 means we are on the boundary, and wedgeid == wedgeid0
    // means we are not on the boundary and we finished the loop
    while (wedgeid != -1 && wedgeid != wedgeid0) {
        cornerid = wedge_get_corner(wedgeid);
        corner_get_wedges(cornerid, &wedgeids);
        order_wedges_ccw(&wedgeids);
        assert(wedgeids[0] == wedgeid);
        wedge_get_coordinates(wedgeids[0], &wcoords);
        pplist->push_back(wcoords[2]);  // centroid
        wedge_get_coordinates(wedgeids[1], &wcoords);
        pplist->push_back(wcoords[1]);  // edge midpoint
        wedgeid = wedge_get_opposite_wedge(wedgeids[1]);
    }

    if (wedgeid == -1) {
        // This is a boundary node, go the other way in a CW manner to get all
        // the coordinates and include the node (nodeid) itself
        wedge_get_coordinates(wedgeid0, &wcoords);
        // edge midpoint
        pplist->insert(pplist->begin(), wcoords[1]);

        wedgeid = wedge_get_opposite_wedge(wedgeid0);
        // We must encounter the other boundary, so we only test for wedgeid ==
        // -1
        while (wedgeid != -1) {
            cornerid = wedge_get_corner(wedgeid);
            corner_get_wedges(cornerid, &wedgeids);
            order_wedges_ccw(&wedgeids);
            assert(wedgeids[1] == wedgeid);
            wedge_get_coordinates(wedgeids[1], &wcoords);
            // centroid
            pplist->insert(pplist->begin(), wcoords[2]);
            wedge_get_coordinates(wedgeids[0], &wcoords);
            // edge midpoint
            pplist->insert(pplist->begin(), wcoords[1]);
            wedgeid = wedge_get_opposite_wedge(wedgeids[0]);
        }

        // Include the node itself
        pplist->insert(pplist->begin(), wcoords[0]);
    }
  }


  //! Order wedges around a node in ccw manner

  void order_wedges_ccw(std::vector<int> *wedgeids) const {
    assert(wedges_requested_);
    assert(wedgeids->size() == 2);
    std::array<Point<2>, 3> wcoords;
    wedge_get_coordinates((*wedgeids)[0], &wcoords);

    // Ensure (*wedgeids)[0] is the first wedge in a CCW direction
    if (!ccw(wcoords[0], wcoords[1], wcoords[2])) {
      std::swap((*wedgeids)[0], (*wedgeids)[1]);
    }
  }

  //! Returns true if the three 2D points (p1, p2, p3) are a counter-clockwise
  //! turn, otherwise returns false (corresponding to clockwise or collinear)

  bool ccw(Point<2> const& p1, Point<2> const& p2, Point<2> const& p3) const {
      return (p2[0] - p1[0]) * (p3[1] - p1[1]) -
             (p2[1] - p1[1]) * (p3[0] - p1[0]) > 0;
  }

  // Can go away when we stop using Clipper
  std::vector<Point<2>> cellToXY(int cellID) const {
    std::vector<Point<2>> cellPoints;
    basicmesh_ptr_->cell_get_coordinates(cellID, &cellPoints);
    return cellPoints;
  }

  //! Get coordinates of wedge in 3D

  void wedges_get_coordinates(int cellID,
      std::vector<std::array<Point<3>, 4>> *wcoords) const {
    assert(basicmesh_ptr_->space_dimension() == 3);
    assert(wedges_requested_);

    std::vector<int> wedges;
    cell_get_wedges(cellID, &wedges);
    int nwedges = wedges.size();

    wcoords->resize(nwedges);
    for (int i = 0; i < nwedges; ++i)
      wedge_get_coordinates(wedges[i], &((*wcoords)[i]), true);
  }

  //! Get coordinates of side in 3D

  void sides_get_coordinates(int cellID,
      std::vector<std::array<Point<3>, 4>> *scoords) const {
    assert(basicmesh_ptr_->space_dimension() == 3);
    assert(sides_requested_);

    std::vector<int> sides;
    cell_get_sides(cellID, &sides);
    int nsides = sides.size();

    scoords->resize(nsides);
    for (int i = 0; i < nsides; ++i)
      side_get_coordinates(sides[i], &((*scoords)[i]), true);
  }

  //! \brief Get adjacent "dual cells" of a given "dual cell"
  void dual_cell_get_node_adj_cells(int const nodeid,
                                    Entity_type const ptype,
                                    std::vector<int> *adjnodes) const {
#ifdef DEBUG
    assert(nodeid < num_entities(NODE, ALL));
#endif
    basicmesh_ptr_->node_get_cell_adj_nodes(nodeid, ptype, adjnodes);
  }


  //! 3D version of coords of nodes of a dual cell
  // Input is the node ID 'nodeid', and it returns the vertex coordinates of
  // the dual cell around this node in `xyzlist`.  The vertices are NOT ordered
  // in any particular way

  void dual_cell_get_coordinates(int const nodeid,
         std::vector<Point<3>> *pplist) const {
#ifdef DEBUG
    assert(nodeid < num_entities(NODE, ALL));
#endif
    assert(basicmesh_ptr_->space_dimension() == 3);
    assert(wedges_requested_);

    // If it does not interfere with the GPU implementation, we could
    // insert into a std::set (nlog(n) instead of n^2) and then copy
    // into the output vector

    // wedge_get_coordinates - (node, edge center, face centroid, cell centroid)
    std::array<Point<3>, 4> wcoords;

    std::vector<int> wedgeids;
    node_get_wedges(nodeid, ALL, &wedgeids);

    std::vector<int> face_list, cell_list;
    std::vector<std::pair<int, int>> edge_list;

    for (auto wedgeid : wedgeids) {
      wedge_get_coordinates(wedgeid, &wcoords);

      int sideid = wedge_get_side(wedgeid);
      int n0 = side_get_node(sideid, 0);
      int n1 = side_get_node(sideid, 1);
      if (n0 > n1) std::swap(n0, n1);
      std::pair<int, int> node_pair(n0, n1);
      if (std::find(edge_list.begin(), edge_list.end(), node_pair) ==
          edge_list.end()) {
        // This "edge" not encountered yet - put it in the edge list and add
        // the corresponding wedge point to the coordinate list

        edge_list.emplace_back(node_pair);
        pplist->emplace_back(wcoords[1]);
      }

      int faceid = wedge_get_face(wedgeid);
      if (std::find(face_list.begin(), face_list.end(), faceid) ==
          face_list.end()) {
        // This face not encountered yet - put it in the face list and add
        // the corresponding wedge point to the coordinate list

        face_list.push_back(faceid);
        pplist->emplace_back(wcoords[2]);
      }

      int cellid = wedge_get_cell(wedgeid);
      if (std::find(cell_list.begin(), cell_list.end(), cellid) ==
          cell_list.end()) {
        // This cell not encountered yet - put it in the cell list and add
        // the cooresponding wedge point to the coordinate list

        cell_list.push_back(cellid);
        pplist->emplace_back(wcoords[3]);
      }
    }
  }

  // Get the coordinates of the wedges of the dual mesh in 3D

  void dual_wedges_get_coordinates(int nodeid,
      std::vector<std::array<Point<3>, 4>> *wcoords) const {
#ifdef DEBUG
    assert(nodeid < num_entities(NODE, ALL));
#endif
    assert(wedges_requested_);
    std::vector<int> wedges;
    node_get_wedges(nodeid, ALL, &wedges);
    int nwedges = wedges.size();
    wcoords->resize(nwedges);

    for (int i = 0; i < nwedges; ++i)
      wedge_get_coordinates(wedges[i], &((*wcoords)[i]), true);
  }

  /// \brief Centroid of a dual cell
  //
  // Centroid of a dual cell.

  template <long D>
  void dual_cell_centroid(int nodeid, Point<D> *centroid) const {
#ifdef DEBUG
    assert(nodeid < num_entities(NODE, ALL));
#endif
    assert(wedges_requested_);
    
    bool posvol_order = true;
    std::vector<int> wedgeids;
    node_get_wedges(nodeid, ALL, &wedgeids);
    double vol = 0.0;
    for (int i = 0; i < D; i++) (*centroid)[i] = 0.0;
    Point<D> geom_cen;
    for (auto const wid : wedgeids) {
      double wvol = wedge_volume(wid);

      Point<D> wcen;
      std::array<Point<D>, D+1> wcoords;
      wedge_get_coordinates(wid, &wcoords, posvol_order);
      for (int i = 0; i < D+1; i++)
        wcen += wcoords[i];
      wcen /= D+1;

      geom_cen += wcen;
      *centroid += wvol*wcen;
      vol += wvol;
    }
    //If we have a degenerate dual cell of zero volume, we use the geometric center
    if (vol > std::numeric_limits<double>::epsilon()) *centroid /= vol;
    else *centroid = geom_cen/wedgeids.size();
  }


  //! Get the volume of dual cell by finding the corners that attach to the node
  double dual_cell_volume(int const nodeid) const {
#ifdef DEBUG
    assert(nodeid < num_entities(NODE, ALL));
#endif
    assert(corners_requested_);
    std::vector<int> cornerids;
    node_get_corners(nodeid, ALL, &cornerids);
    double vol = 0.0;
    for (auto const cid : cornerids)
      vol += corner_volume(cid);
    return vol;
  }


#ifdef HAVE_TANGRAM
  // TEMPORARY - Until we pull out WONTON into separate repository
  template <long D>
  void cell_get_coordinates(int const cellid,
                            std::vector<Tangram::Point<D>> *tplist) const {
    std::vector<Point<D>> pplist;
    cell_get_coordinates(cellid, &pplist);

    tplist->clear();
    for (auto const& pp : pplist)
      tplist->emplace_back(pp);
  }

  //! Get cells of given Entity_type connected to face (in no particular order)
  void face_get_cells(int const faceid, Tangram::Entity_type const etype,
                      std::vector<int> *cells) const {
    face_get_cells(faceid, static_cast<Entity_type>(etype), cells);
  }
#endif




 protected:
  void build_aux_entities() {
    int dim = basicmesh_ptr_->space_dimension();

    // We will use the simplicial decomposition of cells to compute
    // the real centroid of the cells. The sides of a cell provide
    // such a decomposition but we have a chicken-and-egg problem
    // since one of the points of the sides is centroid of the
    // cell. To solve this, we build the simplices by using the
    // topology of the sides but the geometric centers of the cells
    // instead of the centroids. Then we compute the real centroids as
    // the volume weighted average of the centroids of these temporary
    // simplices. 
    //
    // So the steps should be:
    //
    // compute_approximate_face_centroids
    // compute_approximate_cell_centroids
    //
    // build_sides (topology only)
    //
    // compute_face_centroids
    // compute_cell_centroids
    //
    // compute_side_volumes

    if (dim == 1) {
      compute_approximate_face_centroids<1>();
      compute_approximate_cell_centroids<1>();

      if (sides_requested_)
        build_sides_1D(*this);

      compute_face_centroids<1>();
      compute_cell_centroids<1>();

      if (sides_requested_)
        compute_side_volumes<1>();

    } else if (dim == 2) {
      compute_approximate_face_centroids<2>();
      compute_approximate_cell_centroids<2>();

      if (sides_requested_)
        build_sides_2D(*this);

      compute_face_centroids<2>();
      compute_cell_centroids<2>();

      if (sides_requested_)
        compute_side_volumes<2>();

    } else if (dim == 3) {
      compute_approximate_face_centroids<3>();
      compute_approximate_cell_centroids<3>();

      if (sides_requested_)
        build_sides_3D(*this);

      compute_face_centroids<3>();
      compute_cell_centroids<3>();

      if (sides_requested_)
        compute_side_volumes<3>();
    }

    compute_cell_volumes();  // needs side volumes

    if (wedges_requested_) build_wedges();
    if (corners_requested_) build_corners();

    build_face_to_cell_adjacency();
    flag_entities_on_exterior_boundary();
  }


 private:
  template<int dim> void compute_approximate_cell_centroids();
  template<int dim> void compute_approximate_face_centroids();
  template<int dim> void compute_cell_centroids();
  template<int dim> void compute_face_centroids();
  template<int dim> void compute_side_volumes();
  void compute_cell_volumes();

  void build_wedges();
  void build_corners();

  void build_face_to_cell_adjacency();
  void flag_entities_on_exterior_boundary();

  // External helper functions to build dimension-dependent side info
  // - forward declared at the top of the file so that only an
  // instantiation with the 'BasicMesh' template parameter of this
  // class not just any other 'BasicMesh'

  // WHEN WE TEMPLATE THIS CLASS ON DIMENSION, ALL OF THESE FUNCTIONS
  // WILL BE NAMED 'build_sides' BUT WILL BE DISTINGUISHED BY THE
  // DIFFERENT INTEGERS (1, 2 or 3) USED AS THE SECOND ARGUMENT FOR
  // AuxMeshTopology ARGUMENT

  friend
  void build_sides_1D<BasicMesh>(AuxMeshTopology<BasicMesh>& mesh);
  friend
  void build_sides_2D<BasicMesh>(AuxMeshTopology<BasicMesh>& mesh);
  friend
  void build_sides_3D<BasicMesh>(AuxMeshTopology<BasicMesh>& mesh);

  BasicMesh *basicmesh_ptr_ = nullptr;
  bool sides_requested_ = true;
  bool wedges_requested_ = true;
  bool corners_requested_ = true;

  std::vector<int> sideids_owned_, sideids_ghost_, sideids_all_;
  std::vector<int> wedgeids_owned_, wedgeids_ghost_, wedgeids_all_;
  std::vector<int> cornerids_owned_, cornerids_ghost_, cornerids_all_;

  int num_sides_owned_ = 0, num_sides_ghost_ = 0;
  int num_wedges_owned_ = 0, num_wedges_ghost_ = 0;
  int num_corners_owned_ = 0, num_corners_ghost_ = 0;

  // Cells
  std::vector<double> cell_volumes_;

  // If this class were templated on dimension D, we could make these
  // declarations std::vector<Portage::Point<D>>

  std::vector<std::vector<double>> cell_centroids_;

  // Faces
  std::vector<std::vector<double>> face_centroids_;

  // compute this adjacency and store it - while many mesh frameworks
  // have it some may not and in particular, the flat mesh wrapper we
  // use within Portage does not.
  std::vector<std::array<int, 2>> face_cell_ids_;

  // Sides
  std::vector<int> side_cell_id_;
  std::vector<int> side_face_id_;
  std::vector<std::array<int, 2>> side_node_ids_;
  std::vector<int> side_opp_side_id_;
  std::vector<double> side_volumes_;

  // Wedges - most wedge info is derived from sides
  std::vector<int> wedge_corner_id_;  // need only in 2D if we want to build
  //                                  // an oriented polygon out of corners
  //                                  // or wedges connected to a node

  // Corners
  std::vector<int> corner_node_id_;
  std::vector<int> corner_cell_id_;

  // some other one-many adjacencies - MAY NEED TO REWORK TO BE ABLE
  // TO SEND A MESH FROM ONE PROCESSOR TO ANOTHER
  std::vector<std::vector<int>> cell_side_ids_;
  std::vector<std::vector<int>> cell_corner_ids_;  // do we need this?
  std::vector<std::vector<int>> node_corner_ids_;
  std::vector<std::vector<int>> corner_wedge_ids_;

  // Flag indicating if entities (cells, faces, nodes) are on exterior boundary 
  std::vector<bool> cell_on_exterior_boundary_;
  std::vector<bool> face_on_exterior_boundary_;
  std::vector<bool> node_on_exterior_boundary_;

};  // class AuxMeshTopology


//! build face to cell adjacency
template<typename BasicMesh>
void AuxMeshTopology<BasicMesh>::build_face_to_cell_adjacency() {
  int nfaces = basicmesh_ptr_->num_entities(Entity_kind::FACE,
                                            Entity_type::ALL);
  //  face_cell_ids_.resize(nfaces, {-1, -1});  // I think intel 15 barfs
  //                                            // if I do this
  std::array<int, 2> iniarr = {-1, -1};
  face_cell_ids_.assign(nfaces, iniarr);
  
  int ncells = basicmesh_ptr_->num_entities(Entity_kind::CELL,
                                            Entity_type::ALL);
  for (int c = 0; c < ncells; c++) {
    std::vector<int> cfaces;
    std::vector<int> cfdirs;
    basicmesh_ptr_->cell_get_faces_and_dirs(c, &cfaces, &cfdirs);
    for (auto const& f : cfaces) {
      if (face_cell_ids_[f][0] == -1)
        face_cell_ids_[f][0] = c;
      else
        face_cell_ids_[f][1] = c;
    }
  }
}


// Flag entities on exterior boundaries
template <typename BasicMesh>
void AuxMeshTopology<BasicMesh>::flag_entities_on_exterior_boundary() {  

  // Check and flag all entities, owned or ghost. However, ghost faces that
  // are the outer faces of the ghost layer of cells will be incorrectly
  // tagged as being exterior faces because they will have only one (ghost)
  // cell attached. We expect that this won't be a problem since we won't
  // have to query these faces anyway but if it is, then we will have to
  // do an exchange of information across nodes

  int ncells = basicmesh_ptr_->num_entities(Entity_kind::CELL,
                                            Entity_type::ALL);
  int nfaces = basicmesh_ptr_->num_entities(Entity_kind::FACE,
                                            Entity_type::ALL);
  int nnodes = basicmesh_ptr_->num_entities(Entity_kind::NODE,
                                            Entity_type::ALL);
  
  cell_on_exterior_boundary_.assign(ncells, false);
  face_on_exterior_boundary_.assign(nfaces, false);
  node_on_exterior_boundary_.assign(nnodes, false);

  for (int f = 0; f < nfaces; f++) {
    std::vector<int> fcells;
    basicmesh_ptr_->face_get_cells(f, Entity_type::ALL, &fcells);

    if (fcells.size() == 1) {
      face_on_exterior_boundary_[f] = true;
      cell_on_exterior_boundary_[fcells[0]] = true;
      
      // if the face is on an exterior boundary, all its nodes are too
      std::vector<int> fnodes;
      basicmesh_ptr_->face_get_nodes(f, &fnodes);
      for (auto const &n : fnodes)
        node_on_exterior_boundary_[n] = true;
    }
  }
}

//! side coordinates in 1D
template<typename BasicMesh>
void AuxMeshTopology<BasicMesh>::side_get_coordinates(int const s,
                                            std::array<Point<1>, 2> *scoords,
                                            bool posvol_order) const {
  int c = side_cell_id_[s];
  int n = side_node_ids_[s][0];
  Point<1> nxyz;
  basicmesh_ptr_->node_get_coordinates(n, &nxyz);

  // There are two sides per cell in 1D, and going from left to right
  // (in a +ve direction), the side numbers are even, odd, even, odd
  // starting from 0. So for even sides, we need to return node
  // coordinate and cell centroid if we want the side coordinates to
  // be in positive volume order. For odd sides, we need to return
  // cell centroid and then the node coordinate

  if (posvol_order && s%2) {
    cell_centroid(c, &((*scoords)[0]));
    (*scoords)[1] = nxyz;
  } else {
    (*scoords)[0] = nxyz;
    cell_centroid(c, &((*scoords)[1]));
  }
}

//! side coordinates in 2D

template<typename BasicMesh>
void AuxMeshTopology<BasicMesh>::side_get_coordinates(int const s,
                                            std::array<Point<2>, 3> *scoords,
                                            bool posvol_order) const {
  int c = side_cell_id_[s];
  int n[2];
  n[0] = side_node_ids_[s][0];
  n[1] = side_node_ids_[s][1];

  // Due to the way we set up the sides, the natural ordering of side
  // coordinates also results in a positive side volume in 3D

  basicmesh_ptr_->node_get_coordinates(n[0], &((*scoords)[0]));
  basicmesh_ptr_->node_get_coordinates(n[1], &((*scoords)[1]));
  cell_centroid(c, &((*scoords)[2]));
}

//! side coordinates in 3D

template<typename BasicMesh>
void AuxMeshTopology<BasicMesh>::side_get_coordinates(int const s,
                                            std::array<Point<3>, 4> *scoords,
                                            bool posvol_order) const {
  int c = side_cell_id_[s];
  int f = side_face_id_[s];
  int n[2];
  n[0] = side_node_ids_[s][0];
  n[1] = side_node_ids_[s][1];

  // Due to the way we set up the sides, the natural ordering of side
  // coordinates also results in a positive side volume in 3D

  basicmesh_ptr_->node_get_coordinates(n[0], &((*scoords)[0]));
  basicmesh_ptr_->node_get_coordinates(n[1], &((*scoords)[1]));
  face_centroid(f, &((*scoords)[2]));
  cell_centroid(c, &((*scoords)[3]));
}


//! Wedge coordinates in 1D

template<typename BasicMesh>
void AuxMeshTopology<BasicMesh>::wedge_get_coordinates(int const w,
                                             std::array<Point<1>, 2> *wcoords,
                                             bool posvol_order) const {
  int s = w/2;
  int c = side_cell_id_[s];
  int iw = w%2;   // wedge index in side
  int n = side_node_ids_[s][iw];
  Point<1> nxyz;
  basicmesh_ptr_->node_get_coordinates(n, &nxyz);
  if (posvol_order && iw) {
    cell_centroid(c, &((*wcoords)[0]));
    (*wcoords)[1] = nxyz;
  } else {
    (*wcoords)[0] = nxyz;
    cell_centroid(c, &((*wcoords)[1]));
  }
}

//! Wedge coordinates in 2D

template<typename BasicMesh>
void AuxMeshTopology<BasicMesh>::wedge_get_coordinates(int const w,
                                             std::array<Point<2>, 3> *wcoords,
                                             bool posvol_order) const {
  int s = w/2;
  int c = side_cell_id_[s];
  int iw = w%2;
  int n[2];
  n[0] = side_node_ids_[s][0];
  n[1] = side_node_ids_[s][1];
  Point<2> nxyz[2];
  basicmesh_ptr_->node_get_coordinates(n[0], &nxyz[0]);
  basicmesh_ptr_->node_get_coordinates(n[1], &nxyz[1]);
  Point<2> exyz = (nxyz[0] + nxyz[1])/2;  // mid-point of edge
  if (posvol_order && iw) {  // wedge 1 of side
    (*wcoords)[0] = nxyz[iw];
    cell_centroid(c, &((*wcoords)[1]));
    (*wcoords)[2] = exyz;
  } else {
    (*wcoords)[0] = nxyz[iw];
    (*wcoords)[1] = exyz;
    cell_centroid(c, &((*wcoords)[2]));
  }
}

//! Wedge coordinates in 3D

template<typename BasicMesh>
void AuxMeshTopology<BasicMesh>::wedge_get_coordinates(int const w,
                                             std::array<Point<3>, 4> *wcoords,
                                             bool posvol_order) const {
  int s = w/2;
  int c = side_cell_id_[s];
  int f = side_face_id_[s];
  int iw = w%2;
  int n[2];
  n[0] = side_node_ids_[s][0];
  n[1] = side_node_ids_[s][1];
  Point<3> nxyz[2];
  basicmesh_ptr_->node_get_coordinates(n[0], &nxyz[0]);
  basicmesh_ptr_->node_get_coordinates(n[1], &nxyz[1]);
  Point<3> exyz = (nxyz[0] + nxyz[1])/2.0;  // mid-point of edge
  if (posvol_order && iw) {
    (*wcoords)[0] = nxyz[iw];
    face_centroid(f, &((*wcoords)[1]));
    (*wcoords)[2] = exyz;
    cell_centroid(c, &((*wcoords)[3]));
  } else {
    (*wcoords)[0] = nxyz[iw];
    (*wcoords)[1] = exyz;
    face_centroid(f, &((*wcoords)[2]));
    cell_centroid(c, &((*wcoords)[3]));
  }
}


// Build sides for 1D meshes (even though we currently don't support
// 1D remapping)

template<typename BasicMesh>
void build_sides_1D(AuxMeshTopology<BasicMesh>& mesh) {
  int ncells_owned = mesh.basicmesh_ptr_->num_owned_cells();
  int ncells_ghost = mesh.basicmesh_ptr_->num_ghost_cells();
  int ncells = ncells_owned + ncells_ghost;
  
  mesh.cell_side_ids_.clear();
  mesh.cell_side_ids_.resize(ncells);
  
  int nnodes_owned = mesh.basicmesh_ptr_->num_owned_nodes();
  int nnodes_ghost = mesh.basicmesh_ptr_->num_ghost_nodes();
  int nnodes = nnodes_owned + nnodes_ghost;
  
  int num_sides_all = 2*ncells;
  mesh.num_sides_owned_ = 2*ncells_owned;
  mesh.num_sides_ghost_ = 2*ncells_ghost;
  
  for (int c = 0; c < ncells; ++c)
    mesh.cell_side_ids_[c].reserve(2);

  mesh.sideids_owned_.resize(mesh.num_sides_owned_);
  mesh.sideids_ghost_.resize(mesh.num_sides_ghost_);
  mesh.sideids_all_.resize(num_sides_all);
  mesh.side_cell_id_.assign(num_sides_all, -1);
  mesh.side_face_id_.assign(num_sides_all, -1);
  mesh.side_opp_side_id_.assign(num_sides_all, -1);
  mesh.side_node_ids_.assign(num_sides_all, {{-1, -1}});

  int iall = 0, iown = 0, ighost = 0;
  for (int c = 0; c < ncells; ++c) {
    int sideid = 2*c;    // always 2 sides per cell
    
    std::vector<int> nodeids;
    mesh.basicmesh_ptr_->cell_get_nodes(c, &nodeids);
    
    mesh.cell_side_ids_[c].push_back(sideid);
    mesh.cell_side_ids_[c].push_back(sideid+1);
    mesh.sideids_all_[iall++] = sideid;
    mesh.sideids_all_[iall++] = sideid+1;
    if (mesh.basicmesh_ptr_->cell_get_type(c) == PARALLEL_OWNED) {
      mesh.sideids_owned_[iown++] = sideid;
      mesh.sideids_owned_[iown++] = sideid+1;
    } else if (mesh.basicmesh_ptr_->cell_get_type(c) == PARALLEL_GHOST) {
      mesh.sideids_ghost_[ighost++] = sideid;
      mesh.sideids_ghost_[ighost++] = sideid+1;
    }
    
    // Sides in 1D are degenerate - instead of two nodes of an edge
    // the sides point to the same node
    mesh.side_node_ids_[sideid][0] =  nodeids[0];
    mesh.side_node_ids_[sideid][1] =  nodeids[0];
    mesh.side_node_ids_[sideid+1][0] = nodeids[1];
    mesh.side_node_ids_[sideid+1][1] = nodeids[1];
    
    mesh.side_cell_id_[sideid  ] = c;
    mesh.side_cell_id_[sideid+1] = c;
    
    // face ids are same as node ids in 1D
    mesh.side_face_id_[sideid  ] = nodeids[0];
    mesh.side_face_id_[sideid+1] = nodeids[1];
    
    // across face boundaries
    mesh.side_opp_side_id_[sideid  ] = (sideid == 0) ? -1 : sideid-1;
    mesh.side_opp_side_id_[sideid+1] = (sideid+1 == num_sides_all-1) ?
        -1 : sideid+2;
  }  // for c = 0, ncells-1
}  // build_sides in 1D


// Build sides for 2D meshes

template<typename BasicMesh>
void build_sides_2D(AuxMeshTopology<BasicMesh>& mesh) {
  int ncells_owned = mesh.basicmesh_ptr_->num_owned_cells();
  int ncells_ghost = mesh.basicmesh_ptr_->num_ghost_cells();
  int ncells = ncells_owned + ncells_ghost;

  mesh.cell_side_ids_.clear();
  mesh.cell_side_ids_.resize(ncells);

  int nnodes_owned = mesh.basicmesh_ptr_->num_owned_nodes();
  int nnodes_ghost = mesh.basicmesh_ptr_->num_ghost_nodes();
  int nnodes = nnodes_owned + nnodes_ghost;

  int num_sides_all = 0;
  mesh.num_sides_owned_ = 0;
  mesh.num_sides_ghost_ = 0;

  for (int c = 0; c < ncells; ++c) {
    std::vector<int> cfaces, cfdirs;
    mesh.basicmesh_ptr_->cell_get_faces_and_dirs(c, &cfaces, &cfdirs);

    int numsides_in_cell = cfaces.size();
    num_sides_all += numsides_in_cell;
    if (mesh.basicmesh_ptr_->cell_get_type(c) == PARALLEL_OWNED)
      mesh.num_sides_owned_ += numsides_in_cell;
    else if (mesh.basicmesh_ptr_->cell_get_type(c) == PARALLEL_GHOST)
      mesh.num_sides_ghost_ += numsides_in_cell;

    mesh.cell_side_ids_[c].reserve(numsides_in_cell);
  }

  mesh.sideids_owned_.resize(mesh.num_sides_owned_);
  mesh.sideids_ghost_.resize(mesh.num_sides_ghost_);
  mesh.sideids_all_.resize(num_sides_all);
  mesh.side_cell_id_.assign(num_sides_all, -1);
  mesh.side_face_id_.assign(num_sides_all, -1);
  mesh.side_opp_side_id_.clear();
  mesh.side_opp_side_id_.resize(num_sides_all, -1);
  mesh.side_node_ids_.resize(num_sides_all, {{-1, -1}});

  std::vector<std::vector<int>> sides_of_node(nnodes);  // Temp. var.

  int sideid = 0;
  int iall = 0, iown = 0, ighost = 0;
  for (int c = 0; c < ncells; ++c) {
    std::vector<int> cfaces, cfdirs;
    mesh.basicmesh_ptr_->cell_get_faces_and_dirs(c, &cfaces, &cfdirs);

    int ncfaces = cfaces.size();
    for (int i = 0; i < ncfaces; ++i) {
      int f = cfaces[i];
      int fdir = cfdirs[i];
      std::vector<int> fnodes;
      mesh.basicmesh_ptr_->face_get_nodes(f, &fnodes);  // face=edge in 2D
      int n0 = (fdir == 1) ? fnodes[0] : fnodes[1];
      int n1 = (fdir == 1) ? fnodes[1] : fnodes[0];
      mesh.side_node_ids_[sideid][0] = n0;
      mesh.side_node_ids_[sideid][1] = n1;

      mesh.side_cell_id_[sideid] = c;
      mesh.cell_side_ids_[c].push_back(sideid);

      mesh.side_face_id_[sideid] = f;

      mesh.sideids_all_[iall++] = sideid;
      if (mesh.basicmesh_ptr_->cell_get_type(c) == PARALLEL_OWNED)
        mesh.sideids_owned_[iown++] = sideid;
      else if (mesh.basicmesh_ptr_->cell_get_type(c) == PARALLEL_GHOST)
        mesh.sideids_ghost_[ighost++] = sideid;

      // See if any of the other sides attached to the node
      // shares the same edge (pair of nodes) but is in the
      // adjacent cell. This is called the opposite side

      for (auto const& s2 : sides_of_node[n0]) {
        if (mesh.side_node_ids_[s2][0] == n1 && mesh.side_node_ids_[s2][1] == n0) {
          mesh.side_opp_side_id_[sideid] = s2;
          mesh.side_opp_side_id_[s2] = sideid;
          break;
        }
      }
      sides_of_node[n0].push_back(sideid);

      for (auto const& s2 : sides_of_node[n1]) {
        if (mesh.side_node_ids_[s2][0] == n1 && mesh.side_node_ids_[s2][1] == n0) {
          mesh.side_opp_side_id_[sideid] = s2;
          mesh.side_opp_side_id_[s2] = sideid;
          break;
        }
      }
      sides_of_node[n1].push_back(sideid);
      sideid++;
    }
  }  // for c = 0, ncells-1
}  // build sides for 2D meshes


// build sides for 3D meshes

template<typename BasicMesh>
void build_sides_3D(AuxMeshTopology<BasicMesh>& mesh) {
  int ncells_owned = mesh.basicmesh_ptr_->num_owned_cells();
  int ncells_ghost = mesh.basicmesh_ptr_->num_ghost_cells();
  int ncells = ncells_owned + ncells_ghost;

  mesh.cell_side_ids_.clear();
  mesh.cell_side_ids_.resize(ncells);

  int nnodes_owned = mesh.basicmesh_ptr_->num_owned_nodes();
  int nnodes_ghost = mesh.basicmesh_ptr_->num_ghost_nodes();
  int nnodes = nnodes_owned + nnodes_ghost;

  int num_sides_all = 0;
  mesh.num_sides_owned_ = 0;
  mesh.num_sides_ghost_ = 0;

  for (int c = 0; c < ncells; ++c) {
    std::vector<int> cfaces, cfdirs;
    mesh.basicmesh_ptr_->cell_get_faces_and_dirs(c, &cfaces, &cfdirs);

    int numsides_in_cell = 0;
    for (auto const & f : cfaces) {
      std::vector<int> fnodes;
      mesh.basicmesh_ptr_->face_get_nodes(f, &fnodes);

      int nfnodes = fnodes.size();
      num_sides_all += nfnodes;
      numsides_in_cell += nfnodes;

      if (mesh.basicmesh_ptr_->cell_get_type(c) == PARALLEL_OWNED)
        mesh.num_sides_owned_ += nfnodes;
      else if (mesh.basicmesh_ptr_->cell_get_type(c) == PARALLEL_GHOST)
        mesh.num_sides_ghost_ += nfnodes;
    }

    mesh.cell_side_ids_[c].reserve(numsides_in_cell);
  }

  mesh.sideids_owned_.resize(mesh.num_sides_owned_);
  mesh.sideids_ghost_.resize(mesh.num_sides_ghost_);
  mesh.sideids_all_.resize(num_sides_all);
  mesh.side_cell_id_.assign(num_sides_all, -1);
  mesh.side_face_id_.assign(num_sides_all, -1);
  mesh.side_opp_side_id_.clear();
  mesh.side_opp_side_id_.assign(num_sides_all, -1);
  mesh.side_node_ids_.assign(num_sides_all, {{-1, -1}});

  std::vector<std::vector<int>> sides_of_node(nnodes);  // Temporary variable

  int sideid = 0;
  int iall = 0, iown = 0, ighost = 0;
  for (int c = 0; c < ncells; ++c) {
    std::vector<int> cfaces;
    std::vector<int> cfdirs;
    mesh.basicmesh_ptr_->cell_get_faces_and_dirs(c, &cfaces, &cfdirs);
    
    std::vector<int>::iterator itf = cfaces.begin();
    std::vector<int>::iterator itfd = cfdirs.begin();
    while (itf != cfaces.end()) {
      int f = *itf;
      int fdir = *itfd;

      std::vector<int> fnodes;
      mesh.basicmesh_ptr_->face_get_nodes(f, &fnodes);

      // We want the facet formed by side point 0, side point 1 and
      // the face center to point into the cell - this makes the side
      // volume calculation from its natural order of coordinates
      // direct. So, if the face is being used by the cell in the +ve
      // sense, then reverse the node order

      if (fdir == 1)
        std::reverse(fnodes.begin(), fnodes.end());

      int nfnodes = fnodes.size();

      for (int i = 0; i < nfnodes; ++i) {
        mesh.side_node_ids_[sideid][0] = fnodes[i];
        mesh.side_node_ids_[sideid][1] = fnodes[(i+1)%nfnodes];

        mesh.side_cell_id_[sideid] = c;
        mesh.cell_side_ids_[c].push_back(sideid);

        mesh.side_face_id_[sideid] = f;

        mesh.sideids_all_[iall++] = sideid;
        if (mesh.basicmesh_ptr_->cell_get_type(c) == PARALLEL_OWNED)
          mesh.sideids_owned_[iown++] = sideid;
        else if (mesh.basicmesh_ptr_->cell_get_type(c) == PARALLEL_GHOST)
          mesh.sideids_ghost_[ighost++] = sideid;

        // See if any of the other sides attached to the edge (pair
        // of nodes) shares the same edge and face but is in the
        // adjacent cell. This is called the opposite side

        for (auto const& s2 : sides_of_node[fnodes[i]]) {
          if ((mesh.side_node_ids_[s2][0] == fnodes[(i+1)%nfnodes] &&
               mesh.side_node_ids_[s2][1] == fnodes[i]) ||
              (mesh.side_node_ids_[s2][0] == fnodes[i] &&
               mesh.side_node_ids_[s2][1] == fnodes[(i+1)%nfnodes])) {
            if (mesh.side_face_id_[sideid] == mesh.side_face_id_[s2] &&
                mesh.side_cell_id_[sideid] != mesh.side_cell_id_[s2]) {
              mesh.side_opp_side_id_[sideid] = s2;
              mesh.side_opp_side_id_[s2] = sideid;
              break;
            }
          }
        }
        sides_of_node[fnodes[i]].push_back(sideid);
        sideid++;
      }  // for (int i = 0; i < nfnodes; ++i)

      ++itf;
      ++itfd;
    }  // while (itf != cfaces.end())
  }  // for c = 0, ncells-1
}  // build sides for 3D meshes


// Build wedge information - wedges are half a side

template<typename BasicMesh>
void AuxMeshTopology<BasicMesh>::build_wedges() {
  int ncells_owned = basicmesh_ptr_->num_owned_cells();
  int ncells_ghost = basicmesh_ptr_->num_ghost_cells();
  int ncells = ncells_owned + ncells_ghost;

  int nsides_owned = num_sides_owned_;
  int nsides_ghost = num_sides_ghost_;
  int nsides_all = nsides_owned + nsides_ghost;

  int num_wedges_all = 2*nsides_all;
  num_wedges_owned_ = 2*nsides_owned;
  num_wedges_ghost_ = 2*nsides_ghost;

  wedgeids_owned_.resize(num_wedges_owned_);
  wedgeids_ghost_.resize(num_wedges_ghost_);
  wedge_corner_id_.assign(num_wedges_all, -1);  // filled when building corners

  int iown = 0, ighost = 0;
  for (int s = 0; s < nsides_all; ++s) {
    int wedgeid0 = 2*s;
    int wedgeid1 = 2*s + 1;
    int c = side_cell_id_[s];
    if (basicmesh_ptr_->cell_get_type(c) == PARALLEL_OWNED) {
      wedgeids_owned_[iown++] = wedgeid0;
      wedgeids_owned_[iown++] = wedgeid1;
    } else if (basicmesh_ptr_->cell_get_type(c) == PARALLEL_GHOST) {
      wedgeids_ghost_[ighost++] = wedgeid0;
      wedgeids_ghost_[ighost++] = wedgeid1;
    }
  }

  wedgeids_all_.reserve(num_wedges_all);
  wedgeids_all_ = wedgeids_owned_;  // copy
  wedgeids_all_.insert(wedgeids_all_.end(), wedgeids_ghost_.begin(),
              wedgeids_ghost_.end());

}  // build_wedges


// Build corner info

template <typename BasicMesh>
void AuxMeshTopology<BasicMesh>::build_corners() {
  int ncells_owned = basicmesh_ptr_->num_owned_cells();
  int ncells_ghost = basicmesh_ptr_->num_ghost_cells();
  int ncells = ncells_owned + ncells_ghost;

  int nnodes_owned = basicmesh_ptr_->num_owned_nodes();
  int nnodes_ghost = basicmesh_ptr_->num_ghost_nodes();
  int nnodes = nnodes_owned + nnodes_ghost;

  cell_corner_ids_.clear();
  cell_corner_ids_.resize(ncells);
  node_corner_ids_.clear();
  node_corner_ids_.resize(nnodes);

  int num_corners_all = 0;
  num_corners_owned_ = 0;
  num_corners_ghost_ = 0;

  for (int c = 0; c < ncells; ++c) {
    std::vector<int> cnodes;
    basicmesh_ptr_->cell_get_nodes(c, &cnodes);
    cell_corner_ids_[c].reserve(cnodes.size());

    num_corners_all += cnodes.size();  // as many corners as nodes in cell
    if (basicmesh_ptr_->cell_get_type(c) == Entity_type::PARALLEL_OWNED)
      num_corners_owned_ += cnodes.size();
    else if (basicmesh_ptr_->cell_get_type(c) == Entity_type::PARALLEL_GHOST)
      num_corners_ghost_ += cnodes.size();
  }

  cornerids_owned_.resize(num_corners_owned_);
  cornerids_ghost_.resize(num_corners_ghost_);
  corner_wedge_ids_.clear();
  corner_wedge_ids_.resize(num_corners_all);
  corner_cell_id_.resize(num_corners_all);
  corner_node_id_.resize(num_corners_all);

  int cornerid = 0;
  int iown = 0, ighost = 0, ibndry = 0;
  for (int c = 0; c < ncells; ++c) {
    std::vector<int> cnodes;
    basicmesh_ptr_->cell_get_nodes(c, &cnodes);

    std::vector<int> cwedges;
    cell_get_wedges(c, &cwedges);

    for (auto const& n : cnodes) {
      corner_cell_id_[cornerid] = c;
      cell_corner_ids_[c].push_back(cornerid);

      corner_node_id_[cornerid] = n;
      node_corner_ids_[n].push_back(cornerid);

      if (basicmesh_ptr_->cell_get_type(c) == Entity_type::PARALLEL_OWNED)
        cornerids_owned_[iown++] = cornerid;
      else if (basicmesh_ptr_->cell_get_type(c) == Entity_type::PARALLEL_GHOST)
        cornerids_ghost_[ighost++] = cornerid;

      for (auto const& w : cwedges) {
        int n2 = wedge_get_node(w);
        if (n == n2) {
          corner_wedge_ids_[cornerid].push_back(w);
          wedge_corner_id_[w] = cornerid;
        }
      }  // for (w : cwedges)

      ++cornerid;
    }  // for (n : cnodes)
  }  // for (c : cells())

  cornerids_all_.reserve(num_corners_all);
  cornerids_all_ = cornerids_owned_;  // list copy
  cornerids_all_.insert(cornerids_all_.end(), cornerids_ghost_.begin(),
                        cornerids_ghost_.end());

}  // build_corners


// Compute an approximate centroid as the geometric center of the face nodes
template<typename BasicMesh>
template<int dim>
void AuxMeshTopology<BasicMesh>::compute_approximate_face_centroids() {
  int nfaces = basicmesh_ptr_->num_owned_faces() +
      basicmesh_ptr_->num_ghost_faces();

  std::vector<double> pnt(dim, 0.0);
  face_centroids_.assign(nfaces, pnt);

  for (int f = 0; f < nfaces; f++) {
    std::vector<int> fnodes;
    basicmesh_ptr_->face_get_nodes(f, &fnodes);
    int nfnodes = fnodes.size();
    Point<dim> geomcen;
    
    Portage::Point<dim> ncoord;
    for (int n = 0; n < nfnodes; ++n) {
      basicmesh_ptr_->node_get_coordinates(fnodes[n], &ncoord);
      geomcen += ncoord;
    }
    geomcen /= nfnodes;
    for (int d = 0; d < dim; ++d)
      face_centroids_[f][d] = geomcen[d];
  }
}


// using approximate face centroids and triangulation of the face,
// compute the real centroid as the area weighted average of the
// centroids of the triangular facets

template<typename BasicMesh>
template<int dim>
void AuxMeshTopology<BasicMesh>::compute_face_centroids() {
  int nfaces = basicmesh_ptr_->num_owned_faces() +
      basicmesh_ptr_->num_ghost_faces();

  // Resize (just to be sure) but don't reinitialize since it contains
  // the approximate centroids (geometric centers)
  face_centroids_.resize(nfaces);

  std::vector<bool> face_processed(nfaces, false);

  int ncells = basicmesh_ptr_->num_owned_cells() +
      basicmesh_ptr_->num_ghost_cells();

  std::array<Point<dim>, dim+1> sxyz;
  bool posvol_order = true;

  for (int c = 0; c < ncells; c++) {
    std::vector<int> csides;
    cell_get_sides(c, &csides);

    std::vector<int> cfaces, cfdirs;
    basicmesh_ptr_->cell_get_faces_and_dirs(c, &cfaces, &cfdirs);
    
    std::vector<int> fsides;
    for (auto const& f : cfaces) {
      if (face_processed[f]) continue;

      // compute the real centroid as the area weighted average of the
      // facet centroids

      Point<dim> fcen;
      double farea = 0.0;
      for (auto const& s : csides) {
        if (side_get_face(s) == f) {
          side_get_coordinates(s, &sxyz, posvol_order);
          
          // Centroid of facet of side that lies on face (point in 1D,
          // line in 2D, triangle in 3D)
          Point<dim> fctcen;
          for (int j = 0; j < dim; j++)
            fctcen += sxyz[j];
          fctcen /= dim;

          double fctarea = 0.0;
          if (dim == 1)
            fctarea = 1.0;
          else if (dim == 2) {
            Vector<dim> evec = sxyz[1]-sxyz[0];
            fctarea = evec.norm();
          }
          else if (dim == 3) {
            Vector<3> vec0(sxyz[2][0]-sxyz[0][0], sxyz[2][1]-sxyz[0][1],
                           sxyz[2][2]-sxyz[0][2]);
            Vector<3> vec1(sxyz[1][0]-sxyz[0][0], sxyz[1][1]-sxyz[0][1],
                           sxyz[1][2]-sxyz[0][2]);
            Vector<3> cpvec = cross(vec0, vec1);
            fctarea = cpvec.norm();
          }
          farea += fctarea;

          fcen += fctcen*fctarea;
        }
      }
      //If we have a degenerate face of zero area, we use the geometric center
      if (farea > std::numeric_limits<double>::epsilon()) {
        fcen /= farea;
        for (int d = 0; d < dim; d++)
          face_centroids_[f][d] = fcen[d];
      }
      face_processed[f] = true;
    }
  }
}  // compute_face_centroids


// Compute an approximate centroid as the geometric mean of the cell nodes

template<typename BasicMesh>
template<int dim>
void AuxMeshTopology<BasicMesh>::compute_approximate_cell_centroids() {
  int ncells = basicmesh_ptr_->num_owned_cells() +
      basicmesh_ptr_->num_ghost_cells();

  std::vector<double> pnt(dim, 0.0);
  cell_centroids_.assign(ncells, pnt);

  for (int c = 0; c < ncells; ++c) {

    // temporarily set cell_centroid to be geometric center
    std::vector<int> cnodes;
    basicmesh_ptr_->cell_get_nodes(c, &cnodes);
    int ncnodes = cnodes.size();
    Point<dim> geomcen;
    Point<dim> ncoord;
    for (int n = 0; n < ncnodes; ++n) {
      basicmesh_ptr_->node_get_coordinates(cnodes[n], &ncoord);
      geomcen += ncoord;
    }
    geomcen /= ncnodes;
    for (int d = 0; d < dim; ++d)
      cell_centroids_[c][d] = geomcen[d];
  }
}

// Use approximate cell centroids and the geometry of sides (tets)
// using these approximate centroids to compute the real
// centroids. The real centroid is computed as the volume weighted
// average of the centroids of the sides (tets)

template<typename BasicMesh>
template<int dim>
void AuxMeshTopology<BasicMesh>::compute_cell_centroids() {
  int ncells = basicmesh_ptr_->num_owned_cells() +
      basicmesh_ptr_->num_ghost_cells();

  // Resize (just to be sure) but don't reinitialize since it contains
  // the approximate centroids (geometric centers)

  cell_centroids_.resize(ncells);

  std::array<Point<dim>, dim+1> sxyz;
  bool posvol_order = true;
  for (int c = 0; c < ncells; ++c) {
    // Now compute the real centroid as the volume weighted average of the
    // centroids of the sides with the temporary geometry

    std::vector<int> csides;
    basicmesh_ptr_->cell_get_sides(c, &csides);

    Point<dim> ccen;
    double cellvol = 0.0;
    for (auto const& s : csides) {
      side_get_coordinates(s, &sxyz, posvol_order);
      Point<dim> scen;
      for (int i = 0; i < dim+1; i++)
        scen += sxyz[i];
      scen /= (dim+1);

      double svol = calc_side_volume(sxyz);
      cellvol += svol;

      ccen += scen*svol;
    }

    //If we have a degenerate cell of zero volume, we use the geometric center
    if (cellvol > std::numeric_limits<double>::epsilon()) {    
      ccen /= cellvol;
      for (int d = 0; d < dim; d++)
        cell_centroids_[c][d] = ccen[d];  // true centroid
    }
  }
}


template<typename BasicMesh>
template<int dim>
void AuxMeshTopology<BasicMesh>::compute_side_volumes() {
  int num_sides_all = num_sides_owned_ + num_sides_ghost_;
  side_volumes_.resize(num_sides_all);

  std::array<Point<dim>, dim+1> sxyz;
  bool posvol_order = true;
  for (int s = 0; s < num_sides_all; s++) {
    side_get_coordinates(s, &sxyz, posvol_order);
    side_volumes_[s] = calc_side_volume(sxyz);
  }
}


template<typename BasicMesh>
void AuxMeshTopology<BasicMesh>::compute_cell_volumes() {
  int ncells = basicmesh_ptr_->num_owned_cells() +
      basicmesh_ptr_->num_ghost_cells();

  cell_volumes_.clear();
  cell_volumes_.assign(ncells, 0.0);

  for (int c = 0; c < ncells; ++c)
    for (auto s : cell_side_ids_[c])
      cell_volumes_[c] += side_volumes_[s];
}

//! coords of nodes of a cell
template <typename BasicMesh>
template<long D>
void AuxMeshTopology<BasicMesh>::cell_get_coordinates(int const cellid,
                          std::vector<Point<D>> *pplist) const {
  std::vector<int> cnodes;
  basicmesh_ptr_->cell_get_nodes(cellid, &cnodes);

  int ncnodes = cnodes.size();
  pplist->resize(ncnodes);
  for (int n = 0; n < ncnodes; ++n)
    basicmesh_ptr_->node_get_coordinates(cnodes[n], &((*pplist)[n]));
}

//! Get a faceted (all planar faces) representation of a general 3D
//! polyhedralcell
//
// If a face of the cell is not guaranteed to be planar, it is faceted
// by connecting its edges to a "central point" of the face. For now
// the "central point" is merely the geometric center of its nodes. If
// a client code wants to use a different point, it should furnish
// this routine in its wrapper
//
// WE COULD DO THIS USING THE SIDES BUT IT WOULD BE A BIT MORE WORK TO
// ELIMINATE DUPLICATE POINTS (SEE dual_cell_get_facetization). ALSO,
// IF WE WANT TO HAVE THE OPTION OF A SIMPLER FACETIZATION IF WE KNOW
// SOMETHING ABOUT THE CELL

template <typename BasicMesh>
void AuxMeshTopology<BasicMesh>::cell_get_facetization(int const cellid,
           std::vector<std::vector<int>> *facetpoints,
           std::vector<Point<3>> *points) const {
  facetpoints->clear();
  points->clear();

  std::vector<int> cnodes;
  basicmesh_ptr_->cell_get_nodes(cellid, &cnodes);
  int ncnodes = cnodes.size();

  points->resize(ncnodes);
  for (int n = 0; n < ncnodes; ++n)
    basicmesh_ptr_->node_get_coordinates(cnodes[n], &((*points)[n]));

  std::vector<int> cfaces, cfdirs;
  basicmesh_ptr_->cell_get_faces_and_dirs(cellid, &cfaces, &cfdirs);
  int ncfaces = cfaces.size();

  for (int f = 0; f < ncfaces; f++) {
    std::vector<int> fnodes;
    basicmesh_ptr_->face_get_nodes(cfaces[f], &fnodes);
    int nfnodes = fnodes.size();
    
    // Get the local indices (in the cell node list) of the face nodes
    std::vector<int> fnodes_local(nfnodes);
    for (int n = 0; n < nfnodes; n++) {
      bool found = false;
      int i = 0;
      while (!found && i < ncnodes) {
        found = (fnodes[n] == cnodes[i]);
        if (!found) i++;
      }
      assert(found);
      fnodes_local[n] = i;
    }

    if (nfnodes == 3) {  // Triangle; guaranteed to be planar - no need to split
      if (cfdirs[f] != 1) {  // reverse direction of nodes
        int tmp = fnodes_local[0];
        fnodes_local[0] = fnodes_local[1];
        fnodes_local[1] = tmp;
      }
      facetpoints->emplace_back(fnodes_local);
    } else {
      // quad or more general polygonal face which could be
      // curved. Facetize/Triangulate it by connecting each edge to a
      // central point in the face. This central point is computed as
      // the geometric center of the nodes of the face
      
      // Add centroid of face a new point to the point list
      Point<3> fcen;
      face_centroid(cfaces[f], &fcen);
      points->push_back(fcen);
      int icen = points->size() - 1;

      // Add the triangular facets formed using edges of face and centroid
      std::vector<int> fctpnts(3);
      for (int n = 0; n < nfnodes; n++) {
        if (cfdirs[f] == 1) {
          fctpnts[0] = fnodes_local[n];
          fctpnts[1] = fnodes_local[(n+1)%nfnodes];
        } else {
          fctpnts[0] = fnodes_local[(n+1)%nfnodes];
          fctpnts[1] = fnodes_local[n];
        }
        fctpnts[2] = icen;
        facetpoints->push_back(fctpnts);
      }
    }
  }  // for (f...)
}  // cell_get_facetization


//! Get a triangular facetization of boundary of dual cell (or in
//! other words, control volume) corresponding to a node
template <typename BasicMesh>
void AuxMeshTopology<BasicMesh>::dual_cell_get_facetization(int const nodeid,
                                                            std::vector<std::vector<int>> *facetpoints,
                                                            std::vector<Point<3>> *points) const {
  int64_t factor = 1e10;  // used to manufacture unique edge ID from 2 node IDs

  facetpoints->clear();
  points->clear();

#ifdef DEBUG
  assert(nodeid < num_entities(NODE, ALL));
#endif
  assert(wedges_requested_);
  std::vector<int> wedges;
  node_get_wedges(nodeid, ALL, &wedges);
  int nwedges = wedges.size();
  
  points->reserve(3*nwedges);
  facetpoints->reserve(3*nwedges);
  
  // Tracking wedge facet points using these ID/type pairs will help
  // us identify duplicates without doing a distance check between
  // coordinates
  std::vector<int64_t> pntentids;  // ID of entity that wedge point is on
  std::vector<int> pntenttypes;  // Type of entity that wedge point is on

  pntentids.reserve(3*nwedges);
  pntenttypes.reserve(3*nwedges);
  
  int np = pntentids.size();
  for (auto const &w : wedges) {
    int s = w/2;
    int c = side_cell_id_[s];
    int f = side_face_id_[s];
    int iw = w%2;
    int n[2];
    n[0] = side_node_ids_[s][0];
    n[1] = side_node_ids_[s][1];
    // Since we don't explicitly have edges, manufacture a
    // (hopefully) unique edge ID
    int64_t e = (n[0] < n[1]) ? n[0]*factor + n[1] : n[1]*factor + n[0];

    int64_t wpentid[2][3];
    int wpenttyp[2][3];
    int nfct = 0;
    if (iw) {
      // facet in the interior of cell c
      wpentid[0][0] = e;      wpentid[0][1] = c;   wpentid[0][2] = f;
      wpenttyp[0][0] = 1;     wpenttyp[0][1] = 3;  wpenttyp[0][2] = 2;
      nfct++;

      // If the face associated with the wedge is on the exterior boundary
      // we have to add the facet on the boundary of cell c
      if (on_exterior_boundary(FACE, f)) {
        wpentid[1][0] = n[iw];  wpentid[1][1] = e;   wpentid[1][2] = f;
        wpenttyp[1][0] = 0;     wpenttyp[1][1] = 1;  wpenttyp[1][2] = 2;
        nfct++;
      }
    } else {
      // facet in the interior of cell c
      wpentid[0][0] = e;      wpentid[0][1] = f;   wpentid[0][2] = c;
      wpenttyp[0][0] = 1;     wpenttyp[0][1] = 2;  wpenttyp[0][2] = 3;
      nfct++;

      // If the face associated with the wedge is on the exterior boundary
      // we have to add the facet on the boundary of cell c
      if (on_exterior_boundary(FACE, f)) {
        wpentid[1][0] = n[iw];  wpentid[1][1] = f;   wpentid[1][2] = e;
        wpenttyp[1][0] = 0;     wpenttyp[1][1] = 2;  wpenttyp[1][2] = 1;
        nfct++;
      }
    }
    
    // Get local IDs for wedge facet points
    for (int i = 0; i < nfct; i++) {
      std::vector<int> fctpnts(3);
      for (int j = 0; j < 3; j++) {
        bool found = false;
        if (wpenttyp[i][j] != FACE) {  // use id and type
          int k = 0;
          while (!found && k < np) {
            if (wpentid[i][j] == pntentids[k] &&
                wpenttyp[i][j] == pntenttypes[k]) {
              found = true;
              fctpnts[j] = k;
            } else
              k++;
          }
        } else {  // use coordinate comparison because in distributed case
          //      // flat mesh wrapper may have duplicate faces
          Point<3> fcen1;
          face_centroid(wpentid[i][j], &fcen1);
          int k = 0;
          while (!found && k < np) {
            if (pntenttypes[k] == FACE) {
              Point<3> fcen2;
              face_centroid(pntentids[k], &fcen2);
              if (approxEq(fcen1, fcen2, 1.0e-20)) {
                found = true;
                fctpnts[j] = k;
              }
            }
            k++;
          }
        }
        if (!found) {  // point not in list
          pntentids.push_back(wpentid[i][j]);
          pntenttypes.push_back(wpenttyp[i][j]);
          fctpnts[j] = np++;
        }
      }
      facetpoints->push_back(fctpnts);
    }
  }  // for (w : wedges)
  
  for (int i = 0; i < np; i++) {
    int64_t id = pntentids[i];
    int type = pntenttypes[i];
    Point<3> pxyz;
    switch (type) {
      case 0:
        basicmesh_ptr_->node_get_coordinates(id, &pxyz);
        break;
      case 1: {
        int n1 = id%factor;
        int n0 = (id-n1)/factor;
        Point<3> nxyz0, nxyz1;
        basicmesh_ptr_->node_get_coordinates(n0, &nxyz0);
        basicmesh_ptr_->node_get_coordinates(n1, &nxyz1);
        pxyz = (nxyz0 + nxyz1)/2.0;  // mid-point of edge
        break;
      }
      case 2:
        basicmesh_ptr_->face_centroid(id, &pxyz);
        break;
      case 3:
        basicmesh_ptr_->cell_centroid(id, &pxyz);
        break;
      default:
        std::cerr << "Unknown type" << std::endl;
    }
    points->push_back(pxyz);
  }  // Get coordinates of each unique point
}  // dual_cell_get_facetization


}  // namespace Wonton

#endif  // AUX_MESH_TOPOLOGY_H_
