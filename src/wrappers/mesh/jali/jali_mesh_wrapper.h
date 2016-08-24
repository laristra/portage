/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef JALI_MESH_WRAPPER_H_
#define JALI_MESH_WRAPPER_H_

#include <cassert>
#include <algorithm>
#include <array>

#include "Mesh.hh"                      // Jali mesh header

#include "portage/support/portage.h"
#include "portage/support/Point.h"


namespace Portage {

/*!
  \class Jali_Mesh_Wrapper jali_mesh_wrapper.h
  \brief Jali_Mesh_Wrapper implements mesh methods for Jali

  Jali_Mesh_Wrapper implements methods required for Portage mesh
  queries for the Jali mesh infrastructure
*/

namespace {  // unnamned
//! helper function:  convert Jali point to Portage point
template <long D>
Portage::Point<D> toPortagePoint(const JaliGeometry::Point& jp) {
  Portage::Point<D> pp;
  assert(jp.dim() == D);
  for (int d = 0; d < D; ++d)
    pp[d] = jp[d];
  return pp;
}

}  // namespace unnamed

class Jali_Mesh_Wrapper {
 public:

  //! Constructor
  Jali_Mesh_Wrapper(Jali::Mesh const & mesh) :
      jali_mesh_(mesh)
  {}

  //! Copy constructor
  Jali_Mesh_Wrapper(Jali_Mesh_Wrapper const & inmesh) :
      jali_mesh_(inmesh.jali_mesh_)
  {}

  //! Assignment operator (disabled) - don't know how to implement (RVG)
  Jali_Mesh_Wrapper & operator=(Jali_Mesh_Wrapper const &) = delete;

  //! Empty destructor
  ~Jali_Mesh_Wrapper() {};


  //! Cell area/volume
  double cell_volume(int cellID) const {
    return jali_mesh_.cell_volume(cellID);
  }

  //! Dimension of space or mesh points
  int space_dimension() const {
    return jali_mesh_.space_dimension();
  }

  //! Number of owned cells in the mesh
  int num_owned_cells() const {
    return jali_mesh_.num_entities(Jali::Entity_kind::CELL,
                                   Jali::Entity_type::PARALLEL_OWNED);
  }

  //! Number of owned nodes in the mesh
  int num_owned_nodes() const {
    return jali_mesh_.num_entities(Jali::Entity_kind::NODE,
                                   Jali::Entity_type::PARALLEL_OWNED);
  }

  //! Number of ghost cells in the mesh
  int num_ghost_cells() const {
    return jali_mesh_.num_entities(Jali::Entity_kind::CELL,
                                   Jali::Entity_type::PARALLEL_GHOST);
  }

  //! Number of ghost nodes in the mesh
  int num_ghost_nodes() const {
    return jali_mesh_.num_entities(Jali::Entity_kind::NODE,
                                   Jali::Entity_type::PARALLEL_GHOST);
  }

  //! Number of items of given entity
  int num_entities(Entity_kind const entity, Entity_type const etype=Entity_type::ALL) const {
    return jali_mesh_.num_entities((Jali::Entity_kind)entity,
                                   (Jali::Entity_type)etype);
  }

  //! Iterators on mesh entity - begin
  counting_iterator begin(Entity_kind const entity) const {
    int start_index = 0;
    return make_counting_iterator(start_index);
  }

  //! Iterator on mesh entity - end
  counting_iterator end(Entity_kind const entity, Entity_type const etype=Entity_type::ALL) const {
    int start_index = 0;
    return (make_counting_iterator(start_index) + num_entities(entity, etype));
  }

  //! Get list of nodes for a cell
  void cell_get_nodes(int cellid, std::vector<int> *nodes) const {
    jali_mesh_.cell_get_nodes(cellid, nodes);
  }


  //! Get node connected neighbors of cell
  void cell_get_node_adj_cells(int const cellid,
                               Entity_type const ptype,
                               std::vector<int> *adjcells) const {
    jali_mesh_.cell_get_node_adj_cells(cellid, (Jali::Entity_type) ptype,
                                       adjcells);
  }

  //! \brief Get "adjacent" nodes of given node
  //!
  //! Get "adjacent" nodes of given node - nodes that share a common
  //! cell with given node
  void node_get_cell_adj_nodes(int const nodeid,
                               Entity_type const ptype,
                               std::vector<int> *adjnodes) const {
    adjnodes->clear();

    Jali::Entity_ID_List nodecells;
    jali_mesh_.node_get_cells(nodeid, (Jali::Entity_type) ptype, &nodecells);

    for (auto const& c : nodecells) {
      Jali::Entity_ID_List cellnodes;
      jali_mesh_.cell_get_nodes(c, &cellnodes);

      for (auto const& n : cellnodes) {
        if (n == nodeid) continue;
        if (std::find(adjnodes->begin(), adjnodes->end(), n) == adjnodes->end()) 
          adjnodes->emplace_back(n);
      }
    }
  }

  //! Get the volume of dual cell by finding the corners that attach to the node
  double dual_cell_volume(int const nodeid) const {
    Jali::Entity_ID_List cornerids;
    jali_mesh_.node_get_corners(nodeid, Jali::Entity_type::ALL, &cornerids);
    double vol = 0.0;
    for (auto const cid : cornerids)
      vol += jali_mesh_.corner_volume(cid);
    return vol;
  }

  //! \brief Get adjacent "dual cells" of a given "dual cell"
  void dual_cell_get_node_adj_cells(int const nodeid,
                                    Entity_type const ptype,
                                    std::vector<int> *adjnodes) const {
    node_get_cell_adj_nodes(nodeid,ptype,adjnodes);
  }

  //! coords of a node
  template <long D>
  void node_get_coordinates(int const nodeid, Portage::Point<D>* pp) const {
    JaliGeometry::Point jp;
    jali_mesh_.node_get_coordinates(nodeid, &jp);
    assert(jp.dim() == D);
    *pp = toPortagePoint<D>(jp);
  }

  //! coords of nodes of a cell
  template<long D>
  void cell_get_coordinates(int const cellid,
                            std::vector<Portage::Point<D>> *pplist) const {
    assert(jali_mesh_.space_dimension() == D);

    std::vector<JaliGeometry::Point> jplist;
    jali_mesh_.cell_get_coordinates(cellid, &jplist);

    pplist->resize(jplist.size());
    // This cast appears necessary for proper template deduction
    std::transform(jplist.begin(), jplist.end(), pplist->begin(),
                   (Point<D>(*)(const JaliGeometry::Point&))toPortagePoint<D>);
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
                    std::vector<Portage::Point<2>> *pplist) const {
    assert(jali_mesh_.space_dimension() == 2);

    Jali::Entity_ID cornerid, wedgeid, wedgeid0;
    Jali::Entity_ID_List cornerids, wedgeids;
    std::vector<JaliGeometry::Point> wcoords; // (node, edge midpoint, centroid)

    // Start with an arbitrary corner
    jali_mesh_.node_get_corners(nodeid, Jali::Entity_type::ALL,
                                &cornerids);
    cornerid = cornerids[0];

    // Process this corner
    jali_mesh_.corner_get_wedges(cornerid, &wedgeids);
    order_wedges_ccw(&wedgeids);
    jali_mesh_.wedge_get_coordinates(wedgeids[0], &wcoords);
    pplist->push_back({wcoords[2].x(), wcoords[2].y()}); // centroid
    jali_mesh_.wedge_get_coordinates(wedgeids[1], &wcoords);
    pplist->push_back({wcoords[1].x(), wcoords[1].y()}); // edge midpoint

    wedgeid0 = wedgeids[0];

    // Process the rest of the corners in the CCW manner
    wedgeid = jali_mesh_.wedge_get_opposite_wedge(wedgeids[1]);
    // wedgeid == -1 means we are on the boundary, and wedgeid == wedgeid0
    // means we are not on the boundary and we finished the loop
    while (wedgeid != -1 and wedgeid != wedgeid0) {
        cornerid = jali_mesh_.wedge_get_corner(wedgeid);
        jali_mesh_.corner_get_wedges(cornerid, &wedgeids);
        order_wedges_ccw(&wedgeids);
        assert(wedgeids[0] == wedgeid);
        jali_mesh_.wedge_get_coordinates(wedgeids[0], &wcoords);
        pplist->push_back({wcoords[2].x(), wcoords[2].y()}); // centroid
        jali_mesh_.wedge_get_coordinates(wedgeids[1], &wcoords);
        pplist->push_back({wcoords[1].x(), wcoords[1].y()}); // edge midpoint
        wedgeid = jali_mesh_.wedge_get_opposite_wedge(wedgeids[1]);
    }

    if (wedgeid == -1) {
        // This is a boundary node, go the other way in a CW manner to get all
        // the coordinates and include the node (nodid) itself
        jali_mesh_.wedge_get_coordinates(wedgeid0, &wcoords);
        // edge midpoint
        pplist->insert(pplist->begin(), {wcoords[1].x(), wcoords[1].y()});

        wedgeid = jali_mesh_.wedge_get_opposite_wedge(wedgeid0);
        // We must encounter the other boundary, so we only test for wedgeid ==
        // -1
        while (wedgeid != -1) {
            cornerid = jali_mesh_.wedge_get_corner(wedgeid);
            jali_mesh_.corner_get_wedges(cornerid, &wedgeids);
            order_wedges_ccw(&wedgeids);
            assert(wedgeids[1] == wedgeid);
            jali_mesh_.wedge_get_coordinates(wedgeids[1], &wcoords);
            // centroid
            pplist->insert(pplist->begin(), {wcoords[2].x(), wcoords[2].y()});
            jali_mesh_.wedge_get_coordinates(wedgeids[0], &wcoords);
            // edge midpoint
            pplist->insert(pplist->begin(), {wcoords[1].x(), wcoords[1].y()});
            wedgeid = jali_mesh_.wedge_get_opposite_wedge(wedgeids[0]);
        }

        // Include the node itself
        pplist->insert(pplist->begin(), {wcoords[0].x(), wcoords[0].y()});
    }
  }

  void order_wedges_ccw(Jali::Entity_ID_List *wedgeids) const {
    assert(wedgeids->size() == 2);
    std::vector<JaliGeometry::Point> wcoords;
    jali_mesh_.wedge_get_coordinates((*wedgeids)[0], &wcoords);

    // Ensure (*wedgeids)[0] is the first wedge in a CCW direction
    if (not ccw(wcoords[0], wcoords[1], wcoords[2])) {
      std::swap((*wedgeids)[0], (*wedgeids)[1]);
    }
  }

  // Returns true if the three 2D points (p1, p2, p3) are a counter-clockwise
  // turn, otherwise returns false (corresponding to clockwise or collinear)
  bool ccw(const JaliGeometry::Point& p1,
           const JaliGeometry::Point& p2,
           const JaliGeometry::Point& p3) const {
      return (p2[0] - p1[0]) * (p3[1] - p1[1]) -
             (p2[1] - p1[1]) * (p3[0] - p1[0]) > 0;
  }

  std::vector<Portage::Point<2>>
      cellToXY(Jali::Entity_ID cellID) const {
    std::vector<Portage::Point<2>> cellPoints;
    cell_get_coordinates(cellID, &cellPoints);
    return cellPoints;
  }


  void wedges_get_coordinates(Jali::Entity_ID cellID,
      std::vector<std::array<Portage::Point<3>, 4>> *wcoords) const {
    assert(jali_mesh_.space_dimension() == 3);
    std::vector<Jali::Entity_ID> wedges;
    jali_mesh_.cell_get_wedges(cellID, &wedges);
    for (const auto &wedge : wedges) {
      std::vector<JaliGeometry::Point> coords;
      jali_mesh_.wedge_get_coordinates(wedge, &coords, true);
      std::array<Portage::Point<3>, 4> tmp;
      for (int i=0; i<4; i++)
        tmp[i] = toPortagePoint<3>(coords[i]);
      wcoords->push_back(tmp);
    }
  }

  void sides_get_coordinates(Jali::Entity_ID cellID,
      std::vector<std::array<Portage::Point<3>, 4>> *scoords) const {
    assert(jali_mesh_.space_dimension() == 3);
    std::vector<Jali::Entity_ID> sides;
    jali_mesh_.cell_get_sides(cellID, &sides);
    for (const auto &side : sides) {
      std::vector<JaliGeometry::Point> coords;
      jali_mesh_.side_get_coordinates(side, &coords, true);
      std::array<Portage::Point<3>, 4> tmp;
      for (int i=0; i<4; i++)
        tmp[i] = toPortagePoint<3>(coords[i]);
      scoords->push_back(tmp);
    }
  }

  // Get the simplest possible decomposition of a 3D cell into tets.
  // For this mesh type, that means returning a list of sides.
  void decompose_cell_into_tets(Jali::Entity_ID cellID,
      std::vector<std::array<Portage::Point<3>, 4>> *tcoords) const {
    sides_get_coordinates(cellID, tcoords);
  }


  //! 3D version of coords of nodes of a dual cell
  // Input is the node ID 'nodeid', and it returns the vertex coordinates of
  // the dual cell around this node in `xyzlist`.  The vertices are NOT ordered
  // in any particular way
  void dual_cell_get_coordinates(int const nodeid,
         std::vector<Portage::Point<3>> *pplist) const {
    assert(jali_mesh_.space_dimension() == 3);

    Jali::Entity_ID wedgeid;
    Jali::Entity_ID_List wedgeids;

    // wedge_get_coordinates - (node, edge center, face centroid, cell centroid)
    std::vector<JaliGeometry::Point> wcoords;

    jali_mesh_.node_get_wedges(nodeid, Jali::Entity_type::ALL,
                               &wedgeids);

    std::vector<int> edge_list, face_list, cell_list;

    for (auto wedgeid : wedgeids) {
      jali_mesh_.wedge_get_coordinates(wedgeid,&wcoords);

      int edgeid = jali_mesh_.wedge_get_edge(wedgeid);
      if (std::find(edge_list.begin(),edge_list.end(),edgeid) ==
          edge_list.end()) {
        // This edge not encountered yet - put it in the edge list and add
        // the corresponding wedge point to the coordinate list

        edge_list.push_back(edgeid);
        pplist->emplace_back(wcoords[1][0],wcoords[1][1],wcoords[1][2]);
      }

      int faceid = jali_mesh_.wedge_get_face(wedgeid);
      if (std::find(face_list.begin(),face_list.end(),faceid) ==
          face_list.end()) {
        // This face not encountered yet - put it in the face list and add
        // the corresponding wedge point to the coordinate list

        face_list.push_back(faceid);
        pplist->emplace_back(wcoords[2][0],wcoords[2][1],wcoords[2][2]);
      }

      int cellid = jali_mesh_.wedge_get_cell(wedgeid);
      if (std::find(cell_list.begin(),cell_list.end(),cellid) ==
          cell_list.end()) {
        // This cell not encountered yet - put it in the cell list and add
        // the cooresponding wedge point to the coordinate list

        cell_list.push_back(cellid);
        pplist->emplace_back(wcoords[3][0],wcoords[3][1],wcoords[3][2]);
      }
    }
  }

  // Get the coordinates of the wedges of the dual mesh
  void dual_wedges_get_coordinates(Jali::Entity_ID nodeID,
      std::vector<std::array<Portage::Point<3>, 4>> *wcoords) const {
    std::vector<Jali::Entity_ID> wedges;
    jali_mesh_.node_get_wedges(nodeID, Jali::Entity_type::ALL,
                               &wedges);
    for (const auto &wedge : wedges) {
      std::vector<JaliGeometry::Point> coords;
      jali_mesh_.wedge_get_coordinates(wedge, &coords, true);
      std::array<Portage::Point<3>, 4> tmp;
      for (int i=0; i<4; i++)
        tmp[i] = toPortagePoint<3>(coords[i]);
      wcoords->push_back(tmp);
    }
  }

  /// \brief Centroid of a cell
  //
  // Return the centroid of a cell - THIS ROUTINE IS VIOLATING THE
  // CONVENTION THAT NODE_GET_COORDINATES AND CELL_GET_COORDINATES
  // USES FOR THE VARIABLE TYPE OF THE RETURN COORDINATES BECAUSE
  // BUILDING A GRADIENT OPERATOR WITH DIFFERENT TYPES FOR 2D
  // COORDINATES AND 3D COORDINATES IS VERY CONVOLUTED

  void cell_centroid(Jali::Entity_ID cellid,
                     std::vector<double> *centroid) const {
    JaliGeometry::Point ccen = jali_mesh_.cell_centroid(cellid);
    int dim = ccen.dim();
    centroid->resize(dim);
    for (int i = 0; i < dim; ++i)
      (*centroid)[i] = ccen[i];
  }

  /// \brief Centroid of a dual cell
  //
  // Centroid of a dual cell.

  //! \todo NOTE: THIS IS ASSUMED TO BE THE NODE COORDINATE BECAUSE
  //! THE NODAL VARIABLES LIVE THERE, BUT FOR DISTORTED GRIDS, THE
  //! NODE COORDINATED MAY NOT BE THE CENTROID OF THE DUAL CELL

  void dual_cell_centroid(Jali::Entity_ID nodeid,
                          std::vector<double> *centroid) const {
    JaliGeometry::Point nodepnt;
    jali_mesh_.node_get_coordinates(nodeid, &nodepnt);
    int dim = nodepnt.dim();
    centroid->resize(dim);
    for (int i = 0; i < dim; ++i)
      (*centroid)[i] = nodepnt[i];
  }

 private:
  Jali::Mesh const & jali_mesh_;

}; // class Jali_Mesh_Wrapper


} // end namespace Portage

#endif // JALI_MESH_WRAPPER_H_
