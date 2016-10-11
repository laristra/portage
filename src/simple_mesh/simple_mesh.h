/*----------------------------------------------------------------------------*
 * Copyright (c) 2016 Los Alamos National Security, LLC
 * All rights reserved.
 *----------------------------------------------------------------------------*/

#ifndef SRC_SIMPLE_MESH_SIMPLE_MESH_H_
#define SRC_SIMPLE_MESH_SIMPLE_MESH_H_

#include <vector>
#include <memory>

#include "src/support/portage.h"

typedef std::vector<size_t> IDs;

namespace Portage {

class Simple_Mesh {
 public:
  // FIXME: account for non-3d meshes
Simple_Mesh(double x0, double y0, double z0,
            double x1, double y1, double z1,
            size_t nx, size_t ny, size_t nz) :
  nx_(nx), ny_(ny), nz_(nz),
      x0_(x0), y0_(y0), z0_(z0),
      x1_(x1), y1_(y1), z1_(z1) {
    num_cells_ = nx*ny*nz;
    num_nodes_ = (nx+1)*(ny+1)*(nz+1);
    num_faces_ = (nx_+1)*ny_*nz_ + nx_*(ny_+1)*nz_ + nx_*ny_*(nz_+1);

    // Construct the nodal coordinates from extents and number of nodes
    build_node_coords();

    // Build cell <--> node, cell --> face, face <--> node adjacencies
    build_cfn_adjacencies();

    // Build ownership information
    nodeids_owned_.resize(num_nodes_);
    for (size_t i(0); i < num_nodes_; ++i)
      nodeids_owned_[i] = i;
    nodeids_ghost_.resize(0);
    nodeids_all_ = nodeids_owned_;

    cellids_owned_.resize(num_cells_);
    for (size_t i(0); i < num_cells_; ++i)
      cellids_owned_[i] = i;
    cellids_ghost.resize(0);
    cellids_all_ = cellids_owned_;

    faceids_owned_.resize(num_faces_);
    for (size_t i(0); i < num_faces_; ++i)
      faceids_owned_[i] = i;
    faceids_ghost_.resize(0);
    faceids_all_ = faceids_owned_;

    // Build side, wedge, and corner adjacencies
    build_swc_adjacencies();
  }

  ~Simple_Mesh();

  //! Spatial dimension of posize_ts in the mesh
  inline size_t space_dimension() const {
    return spacedim;
  }

  size_t num_entities(const Entity_kind kind,
                      const Entity_type type) const {
    switch (kind) {
      case Entity_kind::NODE:
        switch (type) {
          case Entity_type::PARALLEL_OWNED:
            return nodeids_owned_.size();
          case Entity_type::PARALLEL_GHOST:
            return nodeids_ghost_.size();
          case Entity_type::ALL:
            return nodeids_all_.size();
          default:
            return 0;
        }
      case Entity_kind::CELL:
        switch (type) {
          case Entity_type::PARALLEL_OWNED:
            return cellids_owned_.size();
          case Entity_type::PARALLEL_GHOST:
            return cellids_ghost_.size();
          case Entity_type::ALL:
            return cellids_all_.size();
          default:
            return 0;
        }
      default:
        return 0;
    }
  }

  template<Entity_type type = Entity_type::ALL>
      const std::vector<size_t> & cells() const;

  double cell_centroid(const size_t cellid) const;
  double cell_volume(const size_t celliD) const;
  double corner_volume(const size_t cornerid) const;


  void cell_get_faces(const size_t cellid,
                      std::vector<size_t> *faces) const {
    auto offset = faces_per_cell_*cellid;
    faces->clear();
    for (size_t i(0); i < faces_per_cell_; ++i)
      faces->push_back(cell_to_face_[i+offset]);
  }
  void cell_get_nodes(const size_t cellid,
                      std::vector<size_t> *nodes) const;
  void cell_get_wedges(const size_t cellid,
                       std::vector<size_t> *wedges) const;
  void cell_get_sides(const size_t cellid,
                      std::vector<size_t> *sides) const;
  void node_get_cells(const size_t nodeid,
                      const Entity_type ptype,
                      std::vector<size_t> *cells) const;
  void node_get_corners(const size_t nodeid,
                        const Entity_type ptype,
                        std::vector<size_t> *corners) const;
  void node_get_wedges(const size_t nodeid,
                       const Entity_type ptype,
                       std::vector<size_t> *wedges) const;
  void corner_get_wedges(const size_t cornerid,
                         std::vector<size_t> *wedges) const;
  size_t wedge_get_edge(const size_t wedgeid) const;
  size_t wedge_get_face(const size_t wedgeid) const;
  size_t wedge_get_cell(const size_t wedgeid) const;
  size_t wedge_get_corner(const size_t wedgeid) const;
  size_t wedge_get_opposite_wedge(const size_t wedgeid) const;


  void cell_get_node_adj_cells(const size_t cellid,
                               const Entity_type ptype,
                               std::vector<size_t> *adjcells) const;
  void node_get_cell_adj_nodes(size_t const nodeid,
                               Entity_type const ptype,
                               std::vector<size_t> *adjnodes) const;


  void cell_get_coordinates(const size_t cellid,
                            std::vector<Posize_t3> *pplist) const;
  void wedge_get_coordinates(const size_t wedgeid,
                             std::vector<Posize_t3> *pplist) const;
  void side_get_coordinates(const size_t sideid,
                            std::vector<Posize_t3> *pplist) const;
  void node_get_coordinates(const size_t nodeid,
                            Posize_t3 *pp) const;

  // Wrapper needs to implement these
  size_t num_owned_cells() const;
  size_t num_owned_nodes() const;
  size_t num_ghost_cells() const;
  size_t num_ghost_nodes() const;

  // counting_iterator begin(Entity_kind const entity) const;
  // counting_iterator end(Entity_kind const entity) const;
  double dual_cell_volume(size_t const nodeids) const;
  double dual_cell_get_node_adj_cells(size_t const nodeid,
                                      Entity_type const ptype,
                                      std::vector<size_t> *adjnodes) const;

 private:
  void build_node_coords() {
    coordinates_.clear();
    nodeids_owned_.clear();

    double hx = (x1_ - x0_)/nx_;
    double hy = (y1_ - y0_)/ny_;
    double hz = (z1_ - z0_)/nz_;

    for (size_t iz(0); iz <= nz_; ++iz)
      for (size_t iy(0); iy <= ny_; ++iy)
        for (size_t ix(0); ix <= nx_; ++ix) {
          auto thisNode = node_index(ix, iy, iz);
          coordinates_.emplace_back(x0+ix*hx,
                                    y0+iy*hy,
                                    z0+iz*hz);
          nodeids_owned_.emplace_back(thisNode);
        }
    nodeids_ghost_.resize(0);
    nodeids_all_ = nodids_owned_;
  }

  void build_cfn_adjacencies() {
    // downward adjacencies
    cell_to_node.resize(nodes_per_cell_*num_cells_);
    cell_to_face.resize(faces_per_cell_*num_cells_);
    face_to_node.resize(nodes_per_face_*num_faces_);
    // upward adjacencies
    node_to_face.resize(faces_per_node*num_nodes_);
    node_to_cell.resize(cells_per_node_*num_nodes_);
    face_to_cell.resize(2*num_faces_);

    // this keeps track of how many cells we have so far per node
    std::vector<size_t> cells_at_node(num_nodes_);
    // likewise for faces we have per node
    std::vector<size_t> faces_at_node(num_nodes_);

    // cell adjacencies
    for (size_t iz(0); iz < nz; ++iz)
      for (size_t iy(0); iy < ny; ++iy)
        for (size_t ix(0); ix < nx; ++ix) {
          auto thisCell = cell_index_(ix, iy, iz);
          auto cstart = nodes_per_cell_ * thisCell;
          auto fstart = faces_per_cell_ * thisCell;

          /*
            Node ordering in hex:
                                  z
               7-------6          ^   y
              /|      /|          |  /
             / |     / |          | /
            4--|----5  |          |/
            |  3----|--2          +-------> x
            | /     | /
            0-------1
           */
          cell_to_node[cstart  ] = node_index(ix_, iy_, iz_);
          cell_to_node[cstart+1] = node_index(ix_+1, iy_, iz_);
          cell_to_node[cstart+2] = node_index(ix_+1, iy_+1, iz_);
          cell_to_node[cstart+3] = node_index(ix_, iy_+1, iz_);
          cell_to_node[cstart+4] = node_index(ix_, iy_, iz_+1);
          cell_to_node[cstart+5] = node_index(ix_+1, iy_, iz_+1);
          cell_to_node[cstart+6] = node_index(ix_+1, iy_+1, iz_+1);
          cell_to_node[cstart+7] = node_index(ix_, iy_+1, iz_+1);

          // loop over the 8 nodes attached to this cell, and assign
          // this cell to its cell connectivity
          // order shouldn't matter here
          for (size_t iiz(iz); iiz <= iz+1; ++iiz)
            for (size_t iiy(iy); iiy <= iy+1; ++iiy)
              for (size_t iix(ix); iix <= ix+1; ++iix) {
                auto thisNode = node_index(iix, iiy, iiz);
                auto cnstart = cells_per_node_ * thisNode;
                auto & c_at_n = cells_at_node[thisNode];
                node_to_cell[thisNode+c_at_n] = thisCell;
                c_at_n++;
              }

          /*
            Face ordering in hex:
                                 z
               +-------+         ^   y
              /|      /|         |  /
             / | 5  2/ |         | /
            +--|----+ 1|         |/
            |3 +0---|--+         +-------> x
            | /  4  | /
            +-------+
           */
          cell_to_face[fstart  ] = xzface_index(ix, iy, iz);    // front
          cell_to_face[fstart+1] = yzface_index(ix+1, iy, iz);  // right
          cell_to_face[fstart+2] = xzface_index(ix, iy+1, iz);  // back
          cell_to_face[fstart+3] = yzface_index(ix, iy, iz);    // left
          cell_to_face[fstart+4] = xyface_index(ix, iy, iz);    // bottom
          cell_to_face[fstart+5] = xyface_index(ix, iy, iz+1);  // top

          // go over the 6 faces attached to this cell, and assign
          // this cell to its cell connectivity
          auto cfstart = 2*xzface_index_(ix, iy, iz) + 1;  // +y side of front
          face_to_cell[cfstart] = thisCell;
          cfstart = 2*yzface_index(ix+1, iy, iz);          // -x side of right
          face_to_cell[cfstart] = thisCell;
          cfstart = 2*xzface_index(ix, iy+1, iz);          // -y side of back
          face_to_cell[cfstart] = thisCell;
          cfstart = 2*yzface_index(ix, iy, iz) + 1;        // +x side of left
          face_to_cell[cfstart] = thisCell;
          cfstart = 2*xyface_index(ix, iy, iz) + 1;        // +z side of bottom
          face_to_cell[cfstart] = thisCell;
          cfstart = 2*xyface_index(ix, iy, iz+1);          // -z side of top
          face_to_cell[cfstart] = thisCell;
        }

    // face adjacencies
    /* xy faces
         3-------2
        /       /
       /       /
      0-------1
    */
    for (size_t iz(0); iz <= nz_; ++iz)
      for (size_t iy(0); iy < ny_; ++iy)
        for (size_t ix(0); ix < nx_; ++ix) {
          auto thisFace = xyface_index(ix, iy, iz);
          auto nstart = nodes_per_face_ * thisFace;
          face_to_node[nstart  ] = node_index(ix, iy, iz);
          face_to_node[nstart+1] = node_index(ix+1, iy, iz);
          face_to_node[nstart+2] = node_index(ix+1, iy+1, iz);
          face_to_node[nstart+3] = node_index(ix, iy+1, iz);
          // loop over the 4 nodes attached to this face
          // and assign this face to their face connectivity
          // order shouldn't matter here
          for (size_t iiy(iy); iiy <= iy+1; ++iiy)
            for (size_t iix(ix); iix <= ix+1; ++iix) {
              auto thisNode = node_index(iix, iiy, iz);
              auto fnstart = faces_per_node_ * thisNode;
              auto & f_at_n = faces_at_node[thisNode];
              node_to_face[thisNode+f_at_n] = thisFace;
              f_at_n++;
              }
        }
    /* xz faces
       3-------2
       |       |
       |       |
       |       |
       0-------1
     */
    for (size_t iz(0); iz < nz_; ++iz)
      for (size_t iy(0); iy <= ny_; ++iy)
        for (size_t ix(0); ix < nx_; ++ix) {
          auto nstart = nodes_per_face_ * xzface_index(ix, iy, iz);
          face_to_node[nstart  ] = node_index(ix, iy, iz);
          face_to_node[nstart+1] = node_index(ix+1, iy, iz);
          face_to_node[nstart+2] = node_index(ix+1, iy, iz+1);
          face_to_node[nstart+3] = node_index(ix, iy, iz+1);
          // loop over the 4 nodes attached to this face
          // and assign this face to their face connectivity
          // order shouldn't matter here
          for (size_t iiz(iz); iiz <= iz+1; ++iiz)
            for (size_t iix(ix); iix <= ix+1; ++iix) {
              auto thisNode = node_index(iix, iy, iiz);
              auto fnstart = faces_per_node_ * thisNode;
              auto & f_at_n = faces_at_node[thisNode];
              node_to_face[thisNode+f_at_n] = thisFace;
              f_at_n++;
              }
        }
    /* yz faces
          2
         /|
        / |
       3  |
       |  1
       | /
       |/
       0
     */
    for (size_t iz(0); iz < nz_; ++iz)
      for (size_t iy(0); iy < ny_; ++iy)
        for (size_t ix(0); ix <= nx_; ++ix) {
          auto nstart = nodes_per_face_ * yzface_index(ix, iy, iz);
          face_to_node[nstart  ] = node_index(ix, iy, iz);
          face_to_node[nstart+1] = node_index(ix, iy+1, iz);
          face_to_node[nstart+2] = node_index(ix, iy+1, iz+1);
          face_to_node[nstart+3] = node_index(ix, iy, iz+1);
          // loop over the 4 nodes attached to this face
          // and assign this face to their face connectivity
          // order shouldn't matter here
          for (size_t iiz(iz); iiz <= iz+1; ++iiz)
            for (size_t iiy(iy); iiy <= iy+1; ++iiy) {
              auto thisNode = node_index(ix, iiy, iiz);
              auto fnstart = faces_per_node_ * thisNode;
              auto & f_at_n = faces_at_node[thisNode];
              node_to_face[thisNode+f_at_n] = thisFace;
              f_at_n++;
              }
        }
  }

  void build_swc_adjacencies() {
    cell_side_ids.resize(num_cells);
    // each side is associated with 1 face and 1 edge
    auto sides_per_cell = faces_per_cell_*edges_per_face;
    sideids_owned_.resize(num_cells_*sides_per_cell);
    sideids_ghost_.resize(0);
    sizeids_all.resize(num_cells_*sides_per_cell);
    for (auto const & c : cells()) {
      cell_side_ids[c].reserve(sides_per_cell);
      /* std::vector<size_t> cfaces; */
      /* cell_get_faces(c, &cfaces); */
    }
  }

  /***********************************************************************
   * DATA - FIXME: removem 3d assumptions
   **********************************************************************/

  size_t spacedim = 3;

  // number of cells in the three coordinate directions
  size_t nx_, ny_, nz;
  // coordinates of lower left front and upper right back of brick
  double x0_, x1_, y0_, y1_, z0_, z1;

  // node positions
  std::vector<Point<3>> coordinates;

  inline size_t node_index_(size_t i, size_t j, size_t k) const;
  inline size_t xyface_index_(size_t i, size_t j, size_t k) const;
  inline size_t yzface_index_(size_t i, size_t j, size_t k) const;
  inline size_t xzface_index_(size_t i, size_t j, size_t k) const;
  inline size_t cell_index_(size_t i, size_t j, size_t k) const;

  // hard coded to 3d hexes for now
  size_t nodes_per_face = 4;
  size_t nodes_per_cell = 8;
  size_t edges_per_face = 4;  // needed for sides
  size_t faces_per_cell = 6;
  size_t cells_per_node = 8;
  size_t faces_per_node = 12;

  size_t num_cells;
  size_t num_nodes;
  size_t num_faces;

  // One-to-many lookups
  std::vector<IDs> node_to_cells;
  std::vector<IDs> node_to_corners;
  std::vector<IDs> node_to_wedges;
  std::vector<IDs> corner_to_wedges;
  std::vector<IDs> cell_to_nodes;
  std::vector<IDs> cell_to_faces;
  std::vector<IDs> cell_to_corners;
  std::vector<IDs> cell_to_sides;

  /* std::vector<size_t> cell_to_face; */
  /* std::vector<size_t> cell_to_node; */
  /* std::vector<size_t> face_to_node; */
  /* std::vector<size_t> face_to_cell; */
  /* std::vector<size_t> node_to_face; */
  /* std::vector<size_t> node_to_cell; */

  // Some geometric quantities

  std::vector<double> cell_volumes, face_areas, edge_lengths,
      side_volumes, corner_volumes;
  std::vector<Posize_t3> cell_centroids, face_centroids, face_normal0,
      face_normal1, edge_vectors, edge_centroids;

  // outward facing normal from side to side in adjacent cell
  std::vector<Posize_t3> side_outward_facet_normal;
  // Normal of the common facet of the two wedges - normal posize_ts out
  // of wedge 0 of side size_to wedge 1
  std::vector<Posize_t3> side_mid_facet_normal;

  // Entity lists

  std::vector<size_t> nodeids_owned_, nodeids_ghost_, nodeids_all;
  std::vector<size_t> edgeids_owned_, edgeids_ghost_, edgeids_all;
  std::vector<size_t> faceids_owned_, faceids_ghost_, faceids_all;
  std::vector<size_t> sideids_owned_, sideids_ghost_,
    sideids_boundary_ghost_, sideids_all;
  std::vector<size_t> wedgeids_owned_, wedgeids_ghost_,
    wedgeids_boundary_ghost_, wedgeids_all;
  std::vector<size_t> cornerids_owned_, cornerids_ghost_,
    cornerids_boundary_ghost_, cornerids_all;
  std::vector<size_t> cellids_owned_, cellids_ghost_,
    cellids_boundary_ghost_, cellids_all;
  std::vector<size_t> dummy_list;  // for unspecialized cases

  // Type info for essential entities - sides, wedges and corners will
  // get their type from their owning cell

  std::vector<Entity_type> cell_type;
  std::vector<Entity_type> face_type;  // if faces requested
  std::vector<Entity_type> edge_type;  // if edges requested
  std::vector<Entity_type> node_type;

  // Some standard topological relationships that are cached. The rest
  // are computed on the fly or obtained from the derived class

  std::vector<size_t_List> cell_face_ids;
  std::vector<size_t_List> face_cell_ids;
  std::vector<size_t_List> cell_edge_ids;
  std::vector<size_t_List> face_edge_ids;
  std::vector<std::array<size_t, 2>> edge_node_ids;

  // Topological relationships involving standard and non-standard
  // entities (sides, corners and wedges). The non-standard entities
  // may be required for polyhedral elements and more accurate
  // discretizations.
  //
  // 1D:
  // A side is a line segment from a node to the cell. Wedges and
  // corners are the same as sides.
  //
  // 2D:
  // A side is a triangle formed by the two nodes of an edge/face and
  // the cell center. A wedge is half of a side formed by one node of
  // the edge, the edge center and the cell center. A corner is a
  // quadrilateral formed by the two wedges in a cell at a node
  //
  // 3D:
  // A side is a tet formed by the two nodes of an edge, a face center
  // and a cell center. A wedge is half a side, formed by a node of
  // the edge, the edge center, the face center and the cell center. A
  // corner is formed by all the wedges of a cell at a node.

  // Sides
  std::vector<size_t> side_cell_id;
  std::vector<size_t> side_face_id;
  std::vector<size_t> side_edge_id;
  std::vector<bool> side_edge_use;  // true: side, edge - p0, p1 match
  std::vector<std::array<size_t, 2>> side_node_ids;
  std::vector<size_t> side_opp_side_id;

  // Wedges - most wedge info is derived from sides
  std::vector<size_t> wedge_corner_id;

  // some other one-many adjacencies
  std::vector<std::vector<size_t>> cell_side_ids;
  std::vector<std::vector<size_t>> cell_corner_ids;
  //  std::vector<std::vector<size_t>> edge_side_ids;
  std::vector<std::vector<size_t>> node_corner_ids;
  std::vector<std::vector<size_t>> corner_wedge_ids;
}  // Simple_Mesh

template<> inline
    const std::vector<size_t> &
    Simple_Mesh::cells<Entity_type::ALL>() const {
  return cellids_all_;
}

size_t Simple_Mesh::node_index_(size_t i, size_t j, size_t k) const {
  return i + j*(nx_+1) + k*(nx_+1)*(ny_+1);
}
size_t Simple_Mesh::cell_index_(size_t i, size_t j, size_t k) const {
  return i + j*nx_ + k*nx_*ny_;
}
size_t Simple_Mesh::xyface_index_(size_t i, size_t j, size_t k) const {
  return i + j*nx_ + k*nx_*ny_;
}
size_t Simple_Mesh::xzface_index_(size_t i, size_t j, size_t k) const {
  return i + j*nx_ + k*nx_*(ny_+1) + xyface_index(0, 0, nz_+1);
}
size_t Simple_Mesh::cell_index_(size_t i, size_t j, size_t k) const {
  return i + j*(nx_+1) + k*(nx_+1)*ny_ + xzface_index(0, 0, nz_);
}

}  // namespace Portage

#endif  // SRC_SIMPLE_MESH_SIMPLE_MESH_H_
