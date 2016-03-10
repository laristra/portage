/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef FLECSI_MESH_WRAPPER_H_
#define FLECSI_MESH_WRAPPER_H_

namespace std
{
    typedef decltype(nullptr) nullptr_t;
}


#include <cassert>
#include <algorithm>
#include <array>

#include "flecsi/specializations/burton/burton.h"

#include "portage/support/portage.h"
#include "portage/support/Point.h"

using mesh_t = flecsi::burton_mesh_t;

using real_t = flecsi::mesh_t::real_t;
using vertex_t = flecsi::mesh_t::vertex_t;

namespace Portage {

  /*!
    @brief Helper routine to make a cartesian 2d grid in FleCSI
    @param[in] xmin,xmax,ymin,ymax The extents of the mesh
    @param[in] ncellsx,ncellsy The number of cells in each direction
    @param[in,out] mesh The flecsi::burton_mesh_t mesh object
  */
  void make_mesh_cart2d(const real_t xmin, const real_t xmax,
                        const real_t ymin, const real_t ymax,
                        const int ncellsx, const int ncellsy,
                        mesh_t & mesh) {
    // grid spacing
    auto dx = (xmax - xmin) / real_t(ncellsx);
    auto dy = (ymax - ymin) / real_t(ncellsy);

    //  mesh_t mesh;
    auto num_verts = (ncellsx+1)*(ncellsy+1);
    // this initializes storage for the number of vertices
    mesh.init_parameters(num_verts);

    // create the verts
    std::vector<vertex_t*> verts;
    for (auto j = 0; j < ncellsy+1; ++j) {
      for (auto i = 0; i < ncellsx+1; ++i) {
        auto vert = mesh.create_vertex(
            {xmin + dx*real_t(i), ymin + dy*real_t(j)});
        verts.push_back(vert);
      }
    }

    // create the cells
    auto ncellsx1 = ncellsx+1;
    for (auto j = 0; j < ncellsy; ++j) {
      for (auto i = 0; i < ncellsx; ++i) {
        // go over verts in counter-clockwise fashion
        auto c = mesh.create_cell({verts[i + j*ncellsx1],
                verts[i + 1 + j*ncellsx1],
                verts[i + 1 + (j + 1)*ncellsx1],
                verts[i + (j + 1)*ncellsx1]});
      }
    }

    // this setups up the faces, edges, wedges, etc and connectivity info
    mesh.init();

  }


/*!
  @class Flecsi_Mesh_Wrapper flecsi_mesh_wrapper.h
  @brief Flecsi_Mesh_Wrapper implements mesh methods for Flecsi

  Flecsi_Mesh_Wrapper implements methods required for Portage mesh
  queries for the Flecsi mesh infrastructure
*/

class Flecsi_Mesh_Wrapper {
 public:

  //! Constructor
  Flecsi_Mesh_Wrapper(mesh_t & mesh) :
      flecsi_mesh_(mesh)
  {}

  //! Copy constructor
  Flecsi_Mesh_Wrapper(Flecsi_Mesh_Wrapper const & inmesh) :
      flecsi_mesh_(inmesh.flecsi_mesh_)
  {}

  //! Assignment operator (disabled)
  Flecsi_Mesh_Wrapper & operator=(Flecsi_Mesh_Wrapper const &) = delete;

  //! Empty destructor
  ~Flecsi_Mesh_Wrapper() {};


  //! Dimension of space or mesh points
  int space_dimension() const {
    return flecsi_mesh_.dimension();
  }

  //! Number of owned cells in the mesh
  int num_owned_cells() const {
    return flecsi_mesh_.num_cells();
  }

  //! Number of owned nodes in the mesh
  int num_owned_nodes() const {
    return flecsi_mesh_.num_vertices();
  }

  // NOTE: I don't know where to get the ghosts yet, so setting this to 0.
  //! Number of ghost cells in the mesh
  int num_ghost_cells() const {
    return 0;
  }

  // NOTE: I don't know where to get the ghosts yet, so setting this to 0.
  //! Number of ghost nodes in the mesh
  int num_ghost_nodes() const {
    return 0;
  }

  //! Number of items of given entity
  int num_entities(Entity_kind const entity) const {
    switch(entity) {
      case NODE :
        return flecsi_mesh_.num_vertices();
        break;
      case EDGE :
        return flecsi_mesh_.num_edges();
        break;
        // This case does not work in FleCSI because only 2d implemented and
        // a 2d face is an edge
      /* case FACE : */
      /*   return flecsi_mesh_.num_faces(); */
      /*   break; */
      case CELL :
        return flecsi_mesh_.num_cells();
        break;
      case WEDGE :
        return flecsi_mesh_.num_wedges();
        break;
      case CORNER :
        return flecsi_mesh_.num_corners();
        break;
      default :
        assert(false && "Error: invalid entity request in num_entities");
    }
  }

  //! Iterators on mesh entity - begin
  counting_iterator begin(Entity_kind const entity) const {
    int start_index = 0;
    return make_counting_iterator(start_index);
  }

  //! Iterator on mesh entity - end
  counting_iterator end(Entity_kind const entity) const {
    int start_index = 0;
    return (make_counting_iterator(start_index) + num_entities(entity));
  }

  //! Get list of nodes for a cell
  void cell_get_nodes(int cellid, std::vector<int> *nodes) const {
    auto thisCell = flecsi_mesh_.cells()[cellid];
    nodes->clear();
    for (auto v : flecsi_mesh_.vertices(thisCell))
      nodes->emplace_back(v.global_id());
  }

  //! Get node connected neighbors of cell
  void cell_get_node_adj_cells(int const cellid,
                               Parallel_type const ptype,
                               std::vector<int> *adjcells) const {
    // I think this is supposed to work, but doesn't:
    // flecsi_mesh_.cells(mycell);
    auto thisCell = flecsi_mesh_.cells()[cellid];
    adjcells->clear();
    // Loop over all nodes of this cell
    for (auto node : flecsi_mesh_.vertices(thisCell)) {
      // Loop over all cells associated with this node
      for (auto cell : flecsi_mesh_.cells(node)) {
        if (cell != thisCell)
          adjcells->emplace_back(cell.global_id());
      }
    }
  }

  //! @brief Get "adjacent" nodes of given node
  //!
  //! Get "adjacent" nodes of given node - nodes that share a common
  //! cell with given node
  void node_get_cell_adj_nodes(int const nodeid,
                               Parallel_type const ptype,
                               std::vector<int> *adjnodes) const {
    auto thisNode = flecsi_mesh_.vertices()[nodeid];
    adjnodes->clear();
    // Loop over cells associated with this node
    for (auto cell : flecsi_mesh_.cells(thisNode)) {
      // Loop over nodes of this cell
      for (auto node : flecsi_mesh_.vertices(cell)) {
        if (thisNode != node)
          adjnodes->emplace_back(node.global_id());
      }
    }
  }

  //! @brief Get adjacent "dual cells" of a given "dual cell"
  void dual_cell_get_node_adj_cells(int const nodeid,
                                    Parallel_type const ptype,
                                    std::vector<int> *adjnodes) const {
    auto thisNode = flecsi_mesh_.vertices()[nodeid];
    adjnodes->clear();
    // Loop over cells associated with this node
    for (auto cell : flecsi_mesh_.cells(thisNode)) {
      // Loop over the nodes associated with this cell
      for (auto node : flecsi_mesh_.vertices(cell)) {
        if (thisNode != node)
          adjnodes->emplace_back(node.global_id());
      }
    }
  }

  //! 1D version of coords of a node
  // NOTE: I don't think FleCSI handles 1d!!!
  void node_get_coordinates(int const nodeid, double *x) const {
    auto thisVert = flecsi_mesh_.vertices()[nodeid];
    *x = thisVert->coordinates()[0];
  }

  //! 2D version of coords of a node
  void node_get_coordinates(int const nodeid,
                            std::pair<double,double> *xy) const {
    auto thisVert = flecsi_mesh_.vertices()[nodeid];
    auto coords = thisVert->coordinates();
    xy->first = coords[0];
    xy->second = coords[1];
  }

  //! 3D version of coords of a node
  // NOTE: FleCSI doesn't have 3d yet!!!
  void node_get_coordinates(int const nodeid,
                            std::tuple<double,double,double> *xyz) const {
    auto thisVert = flecsi_mesh_.vertices()[nodeid];
    auto coords = thisVert->coordinates();
    std::get<0>(*xyz) = coords[0];
    std::get<1>(*xyz) = coords[1];
    std::get<2>(*xyz) = coords[2];
  }

  //! 1D version of coords of nodes of a cell
  // NOTE: I don't think FleCSI handles 1d!!!
  void cell_get_coordinates(int const cellid, std::vector<double> *xlist)
      const {
    assert(space_dimension() == 1);

    // Get this cell object
    auto thisCell = flecsi_mesh_.cells()[cellid];

    // Loop over the vertices of this cell to get their coordinates
    xlist->clear();
    auto theseVerts = flecsi_mesh_.vertices(thisCell);
    for (auto v : theseVerts) {
      auto coords = v->coordinates();
      xlist->emplace_back(coords[0]);
    }
  }

  //! 2D version of coords of nodes of a cell

  void cell_get_coordinates(int const cellid,
                            std::vector<std::pair<double,double> > *xylist)
      const {
    assert(space_dimension() == 2);

    // Get this cell object
    auto thisCell = flecsi_mesh_.cells()[cellid];

    // Loop over the vertices of this cell to get their coordinates
    xylist->clear();
    auto theseVerts = flecsi_mesh_.vertices(thisCell);
    for (auto v : theseVerts) {
      auto coords = v->coordinates();
      xylist->emplace_back(coords[0], coords[1]);
    }
  }

  //! 3D version of coords of nodes of a cell
  // NOTE: FleCSI doesn't have 3d yet!!!
  void cell_get_coordinates(int const cellid,
                            std::vector<std::tuple<double,double,double> >
                            *xyzlist) const {
    assert(space_dimension() == 3);

    // Get this cell object
    auto thisCell = flecsi_mesh_.cells()[cellid];

    // Loop over the vertices of this cell to get their coordinates
    xyzlist->clear();
    auto theseVerts = flecsi_mesh_.vertices(thisCell);
    for (auto v : theseVerts) {
      auto coords = v->coordinates();
      xyzlist->emplace_back(coords[0], coords[1], coords[2]);
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
    assert(space_dimension() == 2);
    auto thisNode = flecsi_mesh_.vertices()[nodeid];
    xylist->clear();
    // Loop over the corners associated with this node
    for (auto corner : flecsi_mesh_.corners(thisNode)) {
      std::vector<int> wedgeIds;
      // The two wedges in this corner
      for (auto wedge : flecsi_mesh_.wedges(corner))
        wedgeIds.emplace_back(wedge.global_id());
      order_wedges_ccw(&wedgeIds);
      // Get the coordinates of each wedgeId, given the new ordering
      for (auto wId : wedgeIds) {
        auto w = flecsi_mesh_.wedges()[wId];
        // push back centroid and then edge midpoint
        auto cc = w->cell()->centroid();
        xylist->emplace_back(cc[0], cc[1]);
        auto nc = w->edge()->midpoint();
        xylist->emplace_back(nc[0], nc[1]);
      }
    }
    // @TODO worry about boundary cases

  }

  void order_wedges_ccw(std::vector<int> *wedgeids) const {
    assert(wedgeids->size() == 2);
    auto firstWedge = flecsi_mesh_.wedges()[(*wedgeids)[0]];
    // NOTE: This mimics the order of Jali in 2d: node, face/edge, cell
    std::vector<Point2> wcoords;
    auto nc = firstWedge->vertex()->coordinates();
    wcoords.emplace_back(nc[0], nc[1]);
    auto ec = firstWedge->edge()->midpoint();
    wcoords.emplace_back(ec[0], ec[1]);
    auto cc = firstWedge->cell()->centroid();
    wcoords.emplace_back(cc[0], cc[1]);

    // Ensure (*wedgeids)[0] is the first wedge in a CCW direction
    if (not ccw(
            {wcoords[0][0], wcoords[0][1]},
            {wcoords[1][0], wcoords[1][1]},
            {wcoords[2][0], wcoords[2][1]})) {
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

  std::vector<std::pair<double, double>>
      cellToXY(int const cellID) const {
    std::vector<std::pair<double, double> > cellPoints;
    cell_get_coordinates(cellID, &cellPoints);
    return cellPoints;
  }

  // NOTE: FleCSI doesn't have 3D yet!!!
  std::vector<std::tuple<double, double, double>>
      cellToXYZ(int const cellID) const {
    std::vector<std::tuple<double, double, double>> cellPoints;
    cell_get_coordinates(cellID, &cellPoints);
    return cellPoints;
  }

  // NOTE: This ASSUMES 3D - the "3" is for the 3 coordinates of a spatial point
  // and the "4" is for the four points that make up a wedge (tet) in 3d
  // NOTE: FleCSI doesn't have 3D yet!!!
  void wedges_get_coordinates(int const cellID,
      std::vector<std::array<std::array<double, 3>, 4>> *wcoords) const {
    assert(space_dimension() == 3);

    assert(false && "FleCSI 3D not implemented");
  }


  // NOTE: FleCSI doesn't have 3D yet!!!
  //! 3D version of coords of nodes of a dual cell
  // Input is the node ID 'nodeid', and it returns the vertex coordinates of
  // the dual cell around this node in `xyzlist`.  The vertices are NOT ordered
  // in any particular way

  void dual_cell_get_coordinates(int const nodeid,
        			 std::vector<std::tuple<double,double,double> > *xyzlist) const {
    assert(space_dimension() == 3);

    assert(false && "FleCSI 3D not implemented");
  }

  // NOTE: This ASSUMES 3D - the "3" is for the 3 coordinates of a spatial point
  // and the "4" is for the four points that make up a wedge (tet) in 3d
  // NOTE: FleCSI doesn't have 3D yet!!!
  // Get the coordinates of the wedges of the dual mesh
  void dual_wedges_get_coordinates(int const nodeID,
      std::vector<std::array<std::array<double, 3>, 4>> *wcoords) const {
    assert(space_dimension() == 3);

    assert(false && "FleCSI 3D not implemented yet");
  }


  /// @brief Centroid of a cell
  //
  // Return the centroid of a cell - THIS ROUTINE IS VIOLATING THE
  // CONVENTION THAT NODE_GET_COORDINATES AND CELL_GET_COORDINATES
  // USES FOR THE VARIABLE TYPE OF THE RETURN COORDINATES BECAUSE
  // BUILDING A GRADIENT OPERATOR WITH DIFFERENT TYPES FOR 2D
  // COORDINATES AND 3D COORDINATES IS VERY CONVOLUTED

  void cell_centroid(int const cellid,
                     std::vector<double> *centroid) const {
    auto thisCell = flecsi_mesh_.cells()[cellid];
    auto cntr = thisCell->centroid();
    centroid->clear();
    for (int i = 0; i < cntr.size(); ++i)
      centroid->emplace_back(cntr[i]);
  }

  /// @brief Centroid of a dual cell
  //
  // Centroid of a dual cell.

  //! \todo NOTE: THIS IS ASSUMED TO BE THE NODE COORDINATE BECAUSE
  //! THE NODAL VARIABLES LIVE THERE, BUT FOR DISTORTED GRIDS, THE
  //! NODE COORDINATED MAY NOT BE THE CENTROID OF THE DUAL CELL

  void dual_cell_centroid(int const nodeid,
                          std::vector<double> *centroid) const {
    auto thisNode = flecsi_mesh_.vertices()[nodeid];
    auto cc = thisNode->coordinates();
    centroid->clear();
    for (auto j = 0; j < cc.size(); ++j)
      centroid->emplace_back(cc[j]);
  }

 private:
  mesh_t & flecsi_mesh_;

}; // class Flecsi_Mesh_Wrapper


/* struct pointsToXY */
/* { */
/*   pointsToXY() { } */
/*   std::vector<std::pair<double,double> > operator()(const std::vector<FlecsiGeometry::Point> ptList){     */
/*     std::vector<std::pair<double, double> > xyList; */
/*     std::for_each(ptList.begin(), ptList.end(), [&xyList](FlecsiGeometry::Point pt){xyList.emplace_back(pt.x(), pt.y());});								      */
/*     return xyList; */
/*   } */
/* }; */

} // end namespace Portage

#endif // FLECSI_MESH_WRAPPER_H_
