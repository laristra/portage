/*----------------------------------------------------------------------------*
 * Copyright (c) 2016 Los Alamos National Security, LLC
 * All rights reserved.
 *----------------------------------------------------------------------------*/

#ifndef SRC_SIMPLE_MESH_SIMPLE_MESH_H_
#define SRC_SIMPLE_MESH_SIMPLE_MESH_H_

#include <vector>
#include <memory>

#include "portage/support/portage.h"
#include "portage/support/Point.h"

typedef int ID;

namespace Portage {

class Simple_Mesh {
 public:
  // FIXME: account for non-3d meshes
Simple_Mesh(double x0, double y0, double z0,
            double x1, double y1, double z1,
            ID nx, ID ny, ID nz) :
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

    // Build ownership information - no ghosts in simple_mesh
    nodeids_owned_.resize(num_nodes_);
    for (ID i(0); i < num_nodes_; ++i)
      nodeids_owned_[i] = i;
    nodeids_ghost_.resize(0);
    nodeids_all_ = nodeids_owned_;

    cellids_owned_.resize(num_cells_);
    for (ID i(0); i < num_cells_; ++i)
      cellids_owned_[i] = i;
    cellids_ghost_.resize(0);
    cellids_all_ = cellids_owned_;

    faceids_owned_.resize(num_faces_);
    for (ID i(0); i < num_faces_; ++i)
      faceids_owned_[i] = i;
    faceids_ghost_.resize(0);
    faceids_all_ = faceids_owned_;
  }

  ~Simple_Mesh() {}

  //! Spatial dimension of poIDs in the mesh
  inline ID space_dimension() const {
    return spacedim;
  }

  ID num_entities(const Entity_kind kind,
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
      case Entity_kind::FACE:
        switch (type) {
          case Entity_type::PARALLEL_OWNED:
            return faceids_owned_.size();
          case Entity_type::PARALLEL_GHOST:
            return faceids_ghost_.size();
          case Entity_type::ALL:
            return faceids_all_.size();
          default:
            return 0;
        }
      default:
        return 0;
    }
  }

  // @TODO: replace with std::copy?
  void cell_get_faces_and_dirs(const ID cellid,
                               std::vector<ID> *faces,
                               std::vector<int> *fdirs) const {
    auto offset = faces_per_cell_*cellid;
    faces->clear();
    fdirs->clear();
    for (ID i(0); i < faces_per_cell_; ++i) {
      faces->push_back(cell_to_face_[i+offset]);
      fdirs->push_back(cell_face_dirs_[i+offset]);
    }
  }
  // @TODO: replace with std::copy?
  void cell_get_nodes(const ID cellid,
                      std::vector<ID> *nodes) const {
    auto offset = nodes_per_cell_*cellid;
    nodes->clear();
    for (ID i(0); i < nodes_per_cell_; ++i)
      nodes->push_back(cell_to_node_[i+offset]);
  }
  // @TODO: replace with std::copy?
  void face_get_nodes(const ID faceid,
                      std::vector<ID> *nodes) const {
    auto offset = nodes_per_face_*faceid;
    nodes->clear();
    for (ID i(0); i < nodes_per_face_; ++i)
      nodes->push_back(face_to_node_[i+offset]);
  }
  // @TODO: replace with std::copy?
  void node_get_cells(const ID nodeid,
                      std::vector<ID> *cells) const {
    auto offset = cells_per_node_aug_*nodeid;
    cells->clear();
    std::cout << "   " << offset << " " << node_to_cell_[offset] << std::endl;
    for (ID i(0); i < node_to_cell_[offset]; ++i) {
      std::cout << "  pushing " << node_to_cell_[offset+i+1] << std::endl;
      cells->push_back(node_to_cell_[i+offset+1]);
    }
  }

  // General specification - specialization follows at bottom of file
  // @TODO throw error/exception
  template<long D>
  void node_get_coordinates(const ID nodeid,
                            Point<D> *pp) const {
    assert(D == space_dimension());
  }

 private:
  void build_node_coords() {
    coordinates_.clear();
    //    nodeids_owned_.clear();

    double hx = (x1_ - x0_)/nx_;
    double hy = (y1_ - y0_)/ny_;
    double hz = (z1_ - z0_)/nz_;

    for (ID iz(0); iz <= nz_; ++iz)
      for (ID iy(0); iy <= ny_; ++iy)
        for (ID ix(0); ix <= nx_; ++ix) {
          coordinates_.emplace_back(x0_+ix*hx,
                                    y0_+iy*hy,
                                    z0_+iz*hz);
        }
  }

  void build_cfn_adjacencies() {
    // downward adjacencies
    cell_to_node_.resize(nodes_per_cell_*num_cells_);
    cell_to_face_.resize(faces_per_cell_*num_cells_);
    cell_face_dirs_.resize(faces_per_cell_*num_cells_);
    face_to_node_.resize(nodes_per_face_*num_faces_);
    // upward adjacencies
    node_to_face_.resize(faces_per_node_aug_*num_nodes_);
    node_to_cell_.resize(cells_per_node_aug_*num_nodes_);
    face_to_cell_.resize(2*num_faces_);

    // likewise for faces we have per node
    std::vector<ID> faces_at_node(num_nodes_);

    // cell adjacencies
    for (ID iz(0); iz < nz_; ++iz)
      for (ID iy(0); iy < ny_; ++iy)
        for (ID ix(0); ix < nx_; ++ix) {
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
          cell_to_node_[cstart  ] = node_index_(ix, iy, iz);
          cell_to_node_[cstart+1] = node_index_(ix+1, iy, iz);
          cell_to_node_[cstart+2] = node_index_(ix+1, iy+1, iz);
          cell_to_node_[cstart+3] = node_index_(ix, iy+1, iz);
          cell_to_node_[cstart+4] = node_index_(ix, iy, iz+1);
          cell_to_node_[cstart+5] = node_index_(ix+1, iy, iz+1);
          cell_to_node_[cstart+6] = node_index_(ix+1, iy+1, iz+1);
          cell_to_node_[cstart+7] = node_index_(ix, iy+1, iz+1);

          // loop over the 8 nodes attached to this cell, and assign
          // this cell to its cell connectivity
          // order shouldn't matter here
          for (ID iiz(iz); iiz <= iz+1; ++iiz)
            for (ID iiy(iy); iiy <= iy+1; ++iiy)
              for (ID iix(ix); iix <= ix+1; ++iix) {
                auto thisNode = node_index_(iix, iiy, iiz);
                auto cnstart = cells_per_node_aug_ * thisNode;
                auto & c_at_n = node_to_cell_[cnstart];
                node_to_cell_[cnstart+c_at_n+1] = thisCell;
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
          cell_to_face_[fstart  ] = xzface_index_(ix, iy, iz);    // front
          cell_face_dirs_[fstart  ] = -1;
          cell_to_face_[fstart+1] = yzface_index_(ix+1, iy, iz);  // right
          cell_face_dirs_[fstart+1] = 1;
          cell_to_face_[fstart+2] = xzface_index_(ix, iy+1, iz);  // back
          cell_face_dirs_[fstart+2] = 1;
          cell_to_face_[fstart+3] = yzface_index_(ix, iy, iz);    // left
          cell_face_dirs_[fstart+3] = -1;
          cell_to_face_[fstart+4] = xyface_index_(ix, iy, iz);    // bottom
          cell_face_dirs_[fstart+4] = -1;
          cell_to_face_[fstart+5] = xyface_index_(ix, iy, iz+1);  // top
          cell_face_dirs_[fstart+5] = 1;

          // go over the 6 faces attached to this cell, and assign
          // this cell to its cell connectivity
          auto cfstart = 2*xzface_index_(ix, iy, iz) + 1;  // +y side of front
          face_to_cell_[cfstart] = thisCell;
          cfstart = 2*yzface_index_(ix+1, iy, iz);          // -x side of right
          face_to_cell_[cfstart] = thisCell;
          cfstart = 2*xzface_index_(ix, iy+1, iz);          // -y side of back
          face_to_cell_[cfstart] = thisCell;
          cfstart = 2*yzface_index_(ix, iy, iz) + 1;        // +x side of left
          face_to_cell_[cfstart] = thisCell;
          cfstart = 2*xyface_index_(ix, iy, iz) + 1;        // +z side of bottom
          face_to_cell_[cfstart] = thisCell;
          cfstart = 2*xyface_index_(ix, iy, iz+1);          // -z side of top
          face_to_cell_[cfstart] = thisCell;
        }

    // face adjacencies
    /* xy faces
         3-------2
        /       /
       /       /
      0-------1
    */
    for (ID iz(0); iz <= nz_; ++iz)
      for (ID iy(0); iy < ny_; ++iy)
        for (ID ix(0); ix < nx_; ++ix) {
          auto thisFace = xyface_index_(ix, iy, iz);
          auto nstart = nodes_per_face_ * thisFace;
          face_to_node_[nstart  ] = node_index_(ix, iy, iz);
          face_to_node_[nstart+1] = node_index_(ix+1, iy, iz);
          face_to_node_[nstart+2] = node_index_(ix+1, iy+1, iz);
          face_to_node_[nstart+3] = node_index_(ix, iy+1, iz);
          // loop over the 4 nodes attached to this face
          // and assign this face to their face connectivity
          // order shouldn't matter here
          for (ID iiy(iy); iiy <= iy+1; ++iiy)
            for (ID iix(ix); iix <= ix+1; ++iix) {
              auto thisNode = node_index_(iix, iiy, iz);
              auto fnstart = faces_per_node_aug_ * thisNode;
              auto & f_at_n = node_to_face_[fnstart];
              node_to_face_[fnstart+f_at_n+1] = thisFace;
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
    for (ID iz(0); iz < nz_; ++iz)
      for (ID iy(0); iy <= ny_; ++iy)
        for (ID ix(0); ix < nx_; ++ix) {
          auto thisFace = xzface_index_(ix, iy, iz);
          auto nstart = nodes_per_face_ * thisFace;
          face_to_node_[nstart  ] = node_index_(ix, iy, iz);
          face_to_node_[nstart+1] = node_index_(ix+1, iy, iz);
          face_to_node_[nstart+2] = node_index_(ix+1, iy, iz+1);
          face_to_node_[nstart+3] = node_index_(ix, iy, iz+1);
          // loop over the 4 nodes attached to this face
          // and assign this face to their face connectivity
          // order shouldn't matter here
          for (ID iiz(iz); iiz <= iz+1; ++iiz)
            for (ID iix(ix); iix <= ix+1; ++iix) {
              auto thisNode = node_index_(iix, iy, iiz);
              auto fnstart = faces_per_node_aug_ * thisNode;
              auto & f_at_n = node_to_face_[fnstart];
              node_to_face_[fnstart+f_at_n+1] = thisFace;
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
    for (ID iz(0); iz < nz_; ++iz)
      for (ID iy(0); iy < ny_; ++iy)
        for (ID ix(0); ix <= nx_; ++ix) {
          auto thisFace = yzface_index_(ix, iy, iz);
          auto nstart = nodes_per_face_ * thisFace;
          face_to_node_[nstart  ] = node_index_(ix, iy, iz);
          face_to_node_[nstart+1] = node_index_(ix, iy+1, iz);
          face_to_node_[nstart+2] = node_index_(ix, iy+1, iz+1);
          face_to_node_[nstart+3] = node_index_(ix, iy, iz+1);
          // loop over the 4 nodes attached to this face
          // and assign this face to their face connectivity
          // order shouldn't matter here
          for (ID iiz(iz); iiz <= iz+1; ++iiz)
            for (ID iiy(iy); iiy <= iy+1; ++iiy) {
              auto thisNode = node_index_(ix, iiy, iiz);
              auto fnstart = faces_per_node_aug_ * thisNode;
              auto & f_at_n = node_to_face_[fnstart];
              node_to_face_[fnstart+f_at_n+1] = thisFace;
              f_at_n++;
              }
        }
  }

  /***********************************************************************
   * DATA - FIXME: removem 3d assumptions
   **********************************************************************/

  ID spacedim = 3;

  // number of cells in the three coordinate directions
  ID nx_, ny_, nz_;
  // coordinates of lower left front and upper right back of brick
  double x0_, x1_, y0_, y1_, z0_, z1_;

  // node positions
  std::vector<Point<3>> coordinates_;

  // hard coded to 3d hexes for now
  ID nodes_per_face_ = 4;
  ID nodes_per_cell_ = 8;
  ID edges_per_face_ = 4;  // needed for sides
  ID faces_per_cell_ = 6;
  ID cells_per_node_aug_ = 9;   // 1 entry for the num cells actually attached
  ID faces_per_node_aug_ = 13;  // 1 entry for the num faces actually attached

  ID num_cells_;
  ID num_nodes_;
  ID num_faces_;

  std::vector<ID> cell_to_face_;
  std::vector<int> cell_face_dirs_;
  std::vector<ID> cell_to_node_;
  std::vector<ID> face_to_node_;
  std::vector<ID> face_to_cell_;
  std::vector<ID> node_to_face_;
  std::vector<ID> node_to_cell_;

  // Entity lists

  std::vector<ID> nodeids_owned_, nodeids_ghost_, nodeids_all_;
  std::vector<ID> faceids_owned_, faceids_ghost_, faceids_all_;
  std::vector<ID> cellids_owned_, cellids_ghost_,
    cellids_boundary_ghost_, cellids_all_;

  // helper functions for looking up indices
  ID node_index_(ID i, ID j, ID k) const {
    return i + j*(nx_+1) + k*(nx_+1)*(ny_+1);
  }
  ID cell_index_(ID i, ID j, ID k) const {
    return i + j*nx_ + k*nx_*ny_;
  }
  ID xyface_index_(ID i, ID j, ID k) const {
    return i + j*nx_ + k*nx_*ny_;
  }
  ID xzface_index_(ID i, ID j, ID k) const {
    return i + j*nx_ + k*nx_*(ny_+1) + xyface_index_(0, 0, nz_+1);
  }
  ID yzface_index_(ID i, ID j, ID k) const {
    return i + j*(nx_+1) + k*(nx_+1)*ny_ + xzface_index_(0, 0, nz_);
  }
};  // class Simple_Mesh

// Specializations
template<>
void Simple_Mesh::node_get_coordinates<3>(const ID nodeid,
                                          Point<3> *pp) const {
  *pp = coordinates_[nodeid];
}



}  // namespace Portage

#endif  // SRC_SIMPLE_MESH_SIMPLE_MESH_H_
