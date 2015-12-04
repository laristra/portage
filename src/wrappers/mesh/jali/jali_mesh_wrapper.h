/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef JALI_MESH_WRAPPER_H_
#define JALI_MESH_WRAPPER_H_

namespace std
{
    typedef decltype(nullptr) nullptr_t;
}


#include <cassert>
#include <algorithm>

#include "Mesh.hh"                      // Jali mesh header

#include "portage/support/portage.h"

/*!
  \class Jali_Mesh_Wrapper jali_mesh_wrapper.h
  \brief Jali_Mesh_Wrapper implements mesh methods for Jali

  Jali_Mesh_Wrapper implements methods required for Portage mesh
  queries for the Jali mesh infrastructure
*/

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
  
  //! Number of owned cells in the mesh
  int num_owned_cells() const {
    return jali_mesh_.num_entities(Jali::CELL, Jali::OWNED);
  }

  //! Number of owned nodes in the mesh
  int num_owned_nodes() const {
    return jali_mesh_.num_entities(Jali::NODE, Jali::OWNED);
  }

  //! Number of ghost cells in the mesh
  int num_ghost_cells() const {
    return jali_mesh_.num_entities(Jali::CELL, Jali::GHOST);
  }

  //! Number of ghost nodes in the mesh
  int num_ghost_nodes() const {
    return jali_mesh_.num_entities(Jali::NODE, Jali::GHOST);
  }

  //! Number of items of given entity
  int num_entities(int const entity) const {
    return jali_mesh_.num_entities((Jali::Entity_kind)entity, Jali::ALL);
  }

  //! Iterators on mesh entity - begin
  Portage::counting_iterator begin(int const entity) const {
    return Portage::make_counting_iterator(0);
  }

  //! Iterator on mesh entity - end
  Portage::counting_iterator end(int const entity) const {
    return (Portage::make_counting_iterator(0) + num_entities(entity));
  }

  //! Get list of nodes for a cell
  void cell_get_nodes(int cellid, std::vector<int> *nodes) const {
    jali_mesh_.cell_get_nodes(cellid, nodes);
  }


  //! Get node connected neighbors of cell
  void cell_get_node_adj_cells(int const cellid, int const ptype,
                               std::vector<int> *adjcells) const {
    jali_mesh_.cell_get_node_adj_cells(cellid, (Jali::Parallel_type) ptype,
                                       adjcells);
  }

  //! \brief Get "adjacent" nodes of given node
  //!
  //! Get "adjacent" nodes of given node - nodes that share a common
  //! cell with given node
  void node_get_cell_adj_nodes(int const nodeid, int const ptype,
                               std::vector<int> *adjnodes) const {
    adjnodes->clear();

    Jali::Entity_ID_List nodecells;
    jali_mesh_.node_get_cells(nodeid, (Jali::Parallel_type) ptype, &nodecells);

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

  //! \brief Get adjacent "dual cells" of a given "dual cell"
  void dual_cell_get_node_adj_cells(int const nodeid, int const ptype,
                                    std::vector<int> *adjnodes) const {
    node_get_cell_adj_nodes(nodeid,ptype,adjnodes);
  }
    

  //! 1D version of coords of a node
  void node_get_coordinates(int const nodeid, double *x) const {
    JaliGeometry::Point p;
    jali_mesh_.node_get_coordinates(nodeid, &p);
    assert(p.dim() == 1);
    *x = p[0];
  }

  //! 2D version of coords of a node
  void node_get_coordinates(int const nodeid, 
                            std::pair<double,double> *xy) const {
    JaliGeometry::Point p;
    jali_mesh_.node_get_coordinates(nodeid, &p);
    assert(p.dim() == 2);
    xy->first = p[0];
    xy->second = p[1];
  }

  //! 3D version of coords of a node
  void node_get_coordinates(int const nodeid, 
                            std::tuple<double,double,double> *xyz) const {
    JaliGeometry::Point p;
    jali_mesh_.node_get_coordinates(nodeid, &p);
    assert(p.dim() == 3);
    std::get<0>(*xyz) = p[0];
    std::get<1>(*xyz) = p[1];
    std::get<2>(*xyz) = p[2];
  }

  //! 1D version of coords of nodes of a cell

  void cell_get_coordinates(int const cellid, std::vector<double> *xlist) const {
    assert(jali_mesh_.space_dimension() == 1);

    std::vector<JaliGeometry::Point> plist;
    jali_mesh_.cell_get_coordinates(cellid, &plist);

    // should convert to a std::for_each or std::transform
    xlist->resize(plist.size());
    std::vector<JaliGeometry::Point>::iterator itp = plist.begin();
    std::vector<double>::iterator itx = xlist->begin();
    while (itp != plist.end()) {
      JaliGeometry::Point p = *itp;
      *itx = p[0];
      ++itp;
      ++itx;
    }
  }

  //! 2D version of coords of nodes of a cell

  void cell_get_coordinates(int const cellid, 
                            std::vector<std::pair<double,double> > *xylist) const {
    assert(jali_mesh_.space_dimension() == 2);

    std::vector<JaliGeometry::Point> plist;
    jali_mesh_.cell_get_coordinates(cellid, &plist);

    // should convert to a std::for_each or std::transform
    xylist->resize(plist.size());
    std::vector<JaliGeometry::Point>::iterator itp = plist.begin();
    std::vector<std::pair<double,double> >::iterator itx = xylist->begin();
    while (itp != plist.end()) {
      JaliGeometry::Point p = *itp;
      *itx = std::pair<double,double>(p[0],p[1]);
      ++itp;
      ++itx;
    }
  }

  //! 3D version of coords of nodes of a cell

  void cell_get_coordinates(int const cellid, 
                            std::vector<std::tuple<double,double,double> > *xyzlist) const {
    assert(jali_mesh_.space_dimension() == 3);

    std::vector<JaliGeometry::Point> plist;
    jali_mesh_.cell_get_coordinates(cellid, &plist);

    // should convert to a std::for_each or std::transform

    xyzlist->resize(plist.size());
    std::vector<JaliGeometry::Point>::iterator itp = plist.begin();
    std::vector<std::tuple<double,double,double> >::iterator itx = xyzlist->begin();
    while (itp != plist.end()) {
      JaliGeometry::Point p = *itp;
      *itx = std::tuple<double,double,double>(p[0],p[1],p[2]);
      ++itp;
      ++itx;
    }

  }

  //! 2D version of coords of nodes of a dual cell
  // Input is the node ID 'nodeid', and it returns the vertex coordinates of
  // the dual cell around this node in `xylist`. The vertices are ordered CCW.
  // For boundary node 'nodeid', the first vertex is the node itself, this
  // uniquely determines the 'xylist' vector. For node 'nodeid' not on a
  // boundary, the vector 'xylist' starts with a random vertex, but it is still
  // ordered CCW. Use the 'dual_cell_coordinates_canonical_rotation' to rotate
  // the 'xylist' into a canonical (unique) form.

  void dual_cell_get_coordinates(int const nodeid,
                    std::vector<std::pair<double,double> > *xylist) const {
    assert(jali_mesh_.space_dimension() == 2);

    Jali::Entity_ID cornerid, wedgeid, wedgeid0;
    Jali::Entity_ID_List cornerids, wedgeids;
    std::vector<JaliGeometry::Point> wcoords; // (node, edge midpoint, centroid)

    // Start with an arbitrary corner
    jali_mesh_.node_get_corners(nodeid, Jali::ALL, &cornerids);
    cornerid = cornerids[0];

    // Process this corner
    jali_mesh_.corner_get_wedges(cornerid, &wedgeids);
    order_wedges_ccw(&wedgeids);
    jali_mesh_.wedge_get_coordinates(wedgeids[0], &wcoords);
    xylist->push_back({wcoords[2].x(), wcoords[2].y()}); // centroid
    jali_mesh_.wedge_get_coordinates(wedgeids[1], &wcoords);
    xylist->push_back({wcoords[1].x(), wcoords[1].y()}); // edge midpoint

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
        xylist->push_back({wcoords[2].x(), wcoords[2].y()}); // centroid
        jali_mesh_.wedge_get_coordinates(wedgeids[1], &wcoords);
        xylist->push_back({wcoords[1].x(), wcoords[1].y()}); // edge midpoint
        wedgeid = jali_mesh_.wedge_get_opposite_wedge(wedgeids[1]);
    }

    if (wedgeid == -1) {
        // This is a boundary node, go the other way in a CW manner to get all
        // the coordinates and include the node (nodid) itself
        jali_mesh_.wedge_get_coordinates(wedgeid0, &wcoords);
        xylist->insert(xylist->begin(), {wcoords[1].x(), wcoords[1].y()}); // edge midpoint

        wedgeid = jali_mesh_.wedge_get_opposite_wedge(wedgeid0);
        // We must encounter the other boundary, so we only test for wedgeid ==
        // -1
        while (wedgeid != -1) {
            cornerid = jali_mesh_.wedge_get_corner(wedgeid);
            jali_mesh_.corner_get_wedges(cornerid, &wedgeids);
            order_wedges_ccw(&wedgeids);
            assert(wedgeids[1] == wedgeid);
            jali_mesh_.wedge_get_coordinates(wedgeids[1], &wcoords);
            xylist->insert(xylist->begin(), {wcoords[2].x(), wcoords[2].y()}); // centroid
            jali_mesh_.wedge_get_coordinates(wedgeids[0], &wcoords);
            xylist->insert(xylist->begin(), {wcoords[1].x(), wcoords[1].y()}); // edge midpoint
            wedgeid = jali_mesh_.wedge_get_opposite_wedge(wedgeids[0]);
        }

        // Include the node itself
        xylist->insert(xylist->begin(), {wcoords[0].x(), wcoords[0].y()}); // node
    }
  }

  void order_wedges_ccw(Jali::Entity_ID_List *wedgeids) const {
    assert(wedgeids->size() == 2);
    std::vector<JaliGeometry::Point> wcoords;
    jali_mesh_.wedge_get_coordinates((*wedgeids)[0], &wcoords);

    // Ensure (*wedgeids)[0] is the first wedge in a CCW direction
    if (not ccw(
            {wcoords[0].x(), wcoords[0].y()},
            {wcoords[1].x(), wcoords[1].y()},
            {wcoords[2].x(), wcoords[2].y()})) {
        std::swap((*wedgeids)[0], (*wedgeids)[1]);
    }
  }

  // Returns true if the three 2D points (p1, p2, p3) are a counter-clockwise
  // turn, otherwise returns false (corresponding to clockwise or collinear)
  bool ccw(const std::pair<double, double> p1,
          const std::pair<double, double> p2,
          const std::pair<double, double> p3) const {
      return (std::get<0>(p2) - std::get<0>(p1)) *
          (std::get<1>(p3) - std::get<1>(p1)) -
          (std::get<1>(p2) - std::get<1>(p1)) *
          (std::get<0>(p3) - std::get<0>(p1)) > 0;
  }

  std::vector<std::pair<double, double> > cellToXY(Jali::Entity_ID cellID) const{
    std::vector<std::pair<double, double> > cellPoints;
    cell_get_coordinates(cellID, &cellPoints);
    return cellPoints;
  }



 private:
  Jali::Mesh const & jali_mesh_;

}; // class Jali_Mesh_Wrapper


struct pointsToXY
{
  pointsToXY() { }
  std::vector<std::pair<double,double> > operator()(const std::vector<JaliGeometry::Point> ptList){    
    std::vector<std::pair<double, double> > xyList;
    std::for_each(ptList.begin(), ptList.end(), [&xyList](JaliGeometry::Point pt){xyList.emplace_back(pt.x(), pt.y());});								     
    return xyList;
  }
};


#endif // JALI_MESH_WRAPPER_H_
