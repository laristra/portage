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
#include <boost/iterator/counting_iterator.hpp>

#include "Mesh.hh"                      // Jali mesh header


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

  //! Number of ghost cells in the mesh
  int num_ghost_cells() const {
    return jali_mesh_.num_entities(Jali::CELL, Jali::GHOST);
  }

  //! Number of items of given entity
  int num_entities(int const entity) const {
    return jali_mesh_.num_entities((Jali::Entity_kind)entity, Jali::ALL);
  }

  //! Iterators on mesh entity - begin
  boost::counting_iterator<int> begin(int const entity) const {
    return boost::make_counting_iterator<int>(0);
  }

  //! Iterator on mesh entity - end
  boost::counting_iterator<int> end(int const entity) const {
    return (boost::make_counting_iterator<int>(0) + num_entities(entity));
  }

  //! Get list of nodes for a cell
  void cell_get_nodes(int cellid, std::vector<int> *nodes) const {
    jali_mesh_.cell_get_nodes(cellid, nodes);
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
  // the dual cell around this node in `xylist`.

  void dual_cell_get_coordinates(int const nodeid,
                    std::vector<std::pair<double,double> > *xylist) const {
    assert(jali_mesh_.space_dimension() == 2);

    Jali::Entity_ID_List cornerids;
    jali_mesh_.node_get_corners(nodeid, Jali::ALL, &cornerids);
    for (const auto cornerid : cornerids) {
        std::vector<JaliGeometry::Point> cncoords;
        jali_mesh_.corner_get_coordinates(cornerid, &cncoords);
        xylist->push_back({cncoords[1].x(), cncoords[1].y()}); // Edge midpoint
        xylist->push_back({cncoords[2].x(), cncoords[2].y()}); // Centroid
    }
    std::pair<double, double> center_node;
    node_get_coordinates(nodeid, &center_node);
    order_points(xylist, center_node);
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

  // Orders points in xylist in a CCW manner wrt. center_node
  void order_points(std::vector<std::pair<double, double> > *xylist,
          const std::pair<double, double> &center_node) const {
      std::sort(xylist->begin(), xylist->end(), [&](std::pair<double, double> a,
                  std::pair<double, double> b) {
              return ccw(b, center_node, a);
              });
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

struct cellToXY
{
  Jali_Mesh_Wrapper const & mesh;
  cellToXY(const Jali_Mesh_Wrapper & mesh): mesh(mesh){}
  std::vector<std::pair<double, double> > operator()(const Jali::Entity_ID cellID){
    std::vector<std::pair<double, double> > cellPoints;
    mesh.cell_get_coordinates(cellID, &cellPoints);
    return cellPoints;
  }
};

#endif // JALI_MESH_WRAPPER_H_
